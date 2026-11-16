"""Phase 1: does flow history carry information that local variables do not?

Fits a ladder of models to the wing boundary-layer statistics and compares
them under leave-one-Reynolds-number-out validation:

    M0  constant           -- the floor any model must beat
    M1  ZPG correlation    -- canonical behaviour, no pressure gradient at all
    M2  linear, local      -- beta, Re_tau, up/Ue, dbeta/dx
    M3  trees, local       -- same inputs, nonlinear
    M4  trees, + history   -- adds the memory-kernel average of upstream beta

M3 -> M4 is the question. If history carries nothing, the paper's practice of
characterising a station by its local beta and Re_tau is sufficient and the
project has no premise.

    python scripts/phase1_baselines.py
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wingtbl.features import (  # noqa: E402
    HISTORY_FEATURES, LOCAL_FEATURES, TARGETS, build_table,
)
from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.validation import (  # noqa: E402
    blocked_streamwise, leave_one_case_out, summarise_errors,
)

SUCTION = ["top1n", "top2n", "top4n", "top10n"]

# W1 is not a well-behaved TBL (paper section 3.1) and its last station has
# separated. It is kept in the table for diagnosis but never trained on and
# never scored, so it cannot flatter or distort the comparison.
EXCLUDE_CASES = ["W1"]

KAPPA, C_CF = 0.384, 4.127


def zpg_cf(Re_theta: np.ndarray) -> np.ndarray:
    """Coles-Fernholz skin friction for a canonical ZPG boundary layer."""
    return 2.0 * ((1.0 / KAPPA) * np.log(np.maximum(Re_theta, 1.01)) + C_CF) ** -2


def zpg_reference(df: pd.DataFrame, target: str) -> np.ndarray | None:
    """What a zero-pressure-gradient layer would do at the same Re_theta.

    Only defined for the two targets with a well-established ZPG
    correlation. ``H`` for a ZPG layer sits near 1.4 over this Reynolds-number
    range and ``Ue_plus`` follows from Cf; the fluctuation quantities have no
    correlation of comparable standing, so M1 is not scored on them.
    """
    if target == "Cf":
        return zpg_cf(df["Re_theta"].to_numpy())
    if target == "Ue_plus":
        return np.sqrt(2.0 / zpg_cf(df["Re_theta"].to_numpy()))
    return None


def fit_predict(model, Xtr, ytr, Xte, log_target: bool):
    """Fit, optionally in log space.

    Every target here is a positive quantity whose dependence on beta and
    Reynolds number is expected to be closer to multiplicative than additive,
    and several span more than a decade across the wing. Fitting log(y) makes
    a power law a *linear* hypothesis, which is both the right inductive bias
    and the form symbolic regression will later be asked to recover.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if log_target:
            model.fit(Xtr, np.log(ytr))
            return np.exp(model.predict(Xte))
        model.fit(Xtr, ytr)
        return model.predict(Xte)


def ridge():
    return make_pipeline(StandardScaler(),
                         RidgeCV(alphas=np.logspace(-4, 4, 40)))


def trees():
    return HistGradientBoostingRegressor(
        max_iter=400, learning_rate=0.06, max_leaf_nodes=7,
        min_samples_leaf=5, l2_regularization=1.0, random_state=0)


# Placebo: as many columns as HISTORY_FEATURES, carrying no information.
# Any model gains flexibility from extra inputs, so the history result is
# only meaningful if it beats this. Added to the table in main().
PLACEBO_FEATURES = [f"placebo{i}" for i in range(len(HISTORY_FEATURES))]


def make_models():
    """The ladder, paired so that local vs +history is compared within a
    model family rather than across families."""
    return {
        "M2 linear, local": (ridge, LOCAL_FEATURES),
        "M3 linear, +history": (ridge, LOCAL_FEATURES + HISTORY_FEATURES),
        "Mp linear, +placebo": (ridge, LOCAL_FEATURES + PLACEBO_FEATURES),
        "M4 trees, local": (trees, LOCAL_FEATURES),
        "M5 trees, +history": (trees, LOCAL_FEATURES + HISTORY_FEATURES),
    }


# All targets are made positive before the log fit; uv_outer is negative
# throughout the wing, so its magnitude is modelled and the sign restored.
NEGATIVE_TARGETS = {"uv_outer"}


def run_scheme(df: pd.DataFrame, splitter, scheme_name: str) -> pd.DataFrame:
    models = make_models()
    rows = []
    for target in TARGETS:
        sign = -1.0 if target in NEGATIVE_TARGETS else 1.0
        y_all = sign * df[target].to_numpy()
        assert (y_all > 0).all(), f"{target} is not sign-definite"

        for fold, tr, te in splitter(df):
            ytr, yte = y_all[tr], y_all[te]

            # A constant in log space is the geometric mean, the right floor
            # for a multiplicative target.
            rows.append(dict(scheme=scheme_name, target=target, fold=fold,
                             model="M0 constant",
                             **summarise_errors(yte, np.full(yte.size,
                                                             np.exp(np.log(ytr).mean())))))

            ref = zpg_reference(df.iloc[te], target)
            if ref is not None:
                rows.append(dict(scheme=scheme_name, target=target, fold=fold,
                                 model="M1 ZPG correlation",
                                 **summarise_errors(yte, sign * ref)))

            for name, (make, feats) in models.items():
                X = df[feats].to_numpy()
                pred = fit_predict(make(), X[tr], ytr, X[te], log_target=True)
                rows.append(dict(scheme=scheme_name, target=target, fold=fold,
                                 model=name, **summarise_errors(yte, pred)))
    return pd.DataFrame(rows)


def pivot(res: pd.DataFrame, scheme: str) -> pd.DataFrame:
    sub = res[res.scheme == scheme]
    # Median across folds of the per-fold median relative error.
    return (sub.groupby(["target", "model"]).median_rel.median()
            .unstack("model") * 100).round(1)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--memory", type=float, default=20.0,
                   help="memory length as a multiple of local delta99")
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    raw = load_wing(args.data)
    df = build_table([raw[v] for v in SUCTION], memory_in_delta99=args.memory)
    df = df[~df.case.isin(EXCLUDE_CASES)].reset_index(drop=True)

    rng = np.random.default_rng(0)
    for c in PLACEBO_FEATURES:
        df[c] = rng.standard_normal(len(df))

    print(f"training table: {len(df)} stations, cases {sorted(set(df.case))}, "
          f"memory = {args.memory:g} * delta99\n")

    res = pd.concat([
        run_scheme(df, lambda d: leave_one_case_out(d), "leave-one-Re-out"),
        run_scheme(df, lambda d: blocked_streamwise(d, n_blocks=4),
                   "blocked streamwise"),
    ], ignore_index=True)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "phase1_results.csv", index=False)

    lines = ["# Phase 1 baselines\n",
             f"Median absolute relative error (%), median across folds. "
             f"Memory length {args.memory:g} x delta99. "
             f"W1 excluded (not a well-behaved TBL).\n"]
    for scheme in ("leave-one-Re-out", "blocked streamwise"):
        tbl = pivot(res, scheme)
        print(f"=== {scheme} ===")
        print(tbl.to_string(), "\n")
        lines += [f"\n## {scheme}\n", tbl.to_markdown(), ""]

    # The headline comparison, within each model family.
    loo = res[res.scheme == "leave-one-Re-out"]
    gains = {}
    for fam, (a, b) in {
        "linear": ("M2 linear, local", "M3 linear, +history"),
        "trees": ("M4 trees, local", "M5 trees, +history"),
        "placebo": ("M2 linear, local", "Mp linear, +placebo"),
    }.items():
        ea = loo[loo.model == a].groupby("target").median_rel.median()
        eb = loo[loo.model == b].groupby("target").median_rel.median()
        gains[fam] = ((ea - eb) / ea * 100).round(1)
    gain = pd.DataFrame(gains)
    print("=== error reduction from history, leave-one-Re-out (%) ===")
    print(gain.to_string(), "\n")
    lines += ["\n## Error reduction from adding history (leave-one-Re-out)\n",
              "Percent reduction in median relative error when the memory-kernel\n"
              "features are added, within each model family. Positive = history helps.\n",
              gain.to_markdown(), ""]

    # Which held-out Reynolds number is hardest.
    best = loo[loo.model == "M3 linear, +history"]
    perfold = (best.pivot_table(index="target", columns="fold",
                                values="median_rel") * 100).round(1)
    print("=== M3 linear, +history: per-fold median rel. error (%) ===")
    print(perfold.to_string())
    lines += ["\n## Per-fold error, M3 linear + history (%)\n",
              "Column = the Reynolds-number case held out.\n",
              perfold.to_markdown(), ""]

    (out / "phase1_results.md").write_text("\n".join(lines))
    print(f"\nwrote {out/'phase1_results.md'} and {out/'phase1_results.csv'}")


if __name__ == "__main__":
    main()
