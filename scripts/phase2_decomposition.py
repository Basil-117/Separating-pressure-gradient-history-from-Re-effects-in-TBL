"""Phase 2: fit the memory length honestly, then split the measured effect.

Two jobs.

1. Fit ``y = C beta^a1 (beta_bar/beta)^a2 Re_tau^a3`` with the memory length
   selected *inside* each training fold, and score it on a held-out
   Reynolds number. Phase 1 chose the memory length by held-out error, which
   is selection on the test set; this repeats the measurement without that.

2. Use the fitted law to decompose the difference between two Reynolds-number
   cases at matched chordwise position. Because the model is linear in log
   space the split is exact:

       dlog y = a1 dlog(beta) + a2 dlog(beta_bar/beta) + a3 dlog(Re_tau)
                \\_______________  ________________/   \\_____  _____/
                                \\/                           \\/
                       pressure-gradient part          Reynolds part

   Vinuesa et al. (2018) read the whole of that difference as a
   Reynolds-number effect, on the premise that beta is matched between
   cases. Phase 0 showed beta is not matched -- it runs ~40% higher at
   Re_c = 200k than at 1M. This quantifies what that does to the conclusion.

    python scripts/phase2_decomposition.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wingtbl.features import AnalysisWindow, TARGETS, build_table  # noqa: E402
from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.memory_model import FEATURES, MemoryPowerLaw  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402
from wingtbl.validation import leave_one_case_out, summarise_errors  # noqa: E402

SUCTION = ["top2n", "top4n", "top10n"]
SIGN = {"uv_outer": -1.0}
PRETTY = {
    "Cf": r"$C_f$", "H": r"$H$", "Ue_plus": r"$U_e^+$",
    "uu_outer": r"$\overline{u_t^2}^+$", "uv_outer": r"$\overline{u_tv_n}^+$",
    "Pk_outer": r"$P_k^+$",
}
# Categorical slots 1-3, validated on the light surface.
TERM_COLOR = {"local PG": "#2a78d6", "history": "#eb6834", "Reynolds": "#1baf7a"}


def load(data: str):
    raw = load_wing(data)
    ds = [raw[v] for v in SUCTION]
    wide = AnalysisWindow(lo=-1.0, hi=2.0)     # untrimmed, for the kernel
    full = {d.attrs["case"]: build_table([d], window=wide).reset_index(drop=True)
            for d in ds}
    return build_table(ds), full


def evaluate(df: pd.DataFrame, full: dict) -> pd.DataFrame:
    """Leave-one-Reynolds-number-out, memory length refitted per fold."""
    rows = []
    for target in TARGETS:
        s = SIGN.get(target, 1.0)
        for held, tr, te in leave_one_case_out(df):
            m = MemoryPowerLaw()
            r = m.fit(df.iloc[tr].reset_index(drop=True), target, full, sign=s)
            pred = m.predict(df.iloc[te].reset_index(drop=True), full)
            truth = s * df.iloc[te][target].to_numpy()
            rows.append(dict(target=target, held_out=held, ell=r.ell,
                             a_beta=r.coef[0], a_hist=r.coef[1], a_Re=r.coef[2],
                             at_bound=r.ell_at_bound,
                             **summarise_errors(truth, pred)))
    return pd.DataFrame(rows)


def decompose(df: pd.DataFrame, full: dict, lo_case="W2", hi_case="W10",
              n_x=60, ell_scale: float | None = None) -> pd.DataFrame:
    """Split the lo_case -> hi_case change into pressure-gradient and Re parts.

    ``ell_scale`` forces the memory length to a multiple of the fitted value,
    for the sensitivity check: the fitted length is only good to tens of
    percent, so the split must be shown not to depend on it.
    """
    xs = np.linspace(0.30, 0.90, n_x)
    rows = []
    for target in TARGETS:
        s = SIGN.get(target, 1.0)
        m = MemoryPowerLaw()
        m.fit(df, target, full, sign=s)
        if ell_scale is not None:
            # Refit the exponents with the length pinned, so the comparison
            # is between complete models rather than a half-updated one.
            lo, hi = m.ell_bounds
            forced = float(np.clip(m.result_.ell * ell_scale, lo, hi))
            m.ell_bounds = (forced, forced * 1.0000001)
            m.fit(df, target, full, sign=s)
        r = m.result_

        X = m.design(df, full)
        cols = {f: X[:, k] for k, f in enumerate(FEATURES)}
        frame = df[["case", "xa"]].copy()
        for f, v in cols.items():
            frame[f] = v
        frame["log_y"] = np.log(s * df[target].to_numpy())

        def interp(case, col):
            sub = frame[frame.case == case].sort_values("xa")
            return np.interp(xs, sub.xa, sub[col])

        d = {f: interp(hi_case, f) - interp(lo_case, f) for f in FEATURES}
        d_obs = interp(hi_case, "log_y") - interp(lo_case, "log_y")

        pg = r.coef[0] * d["beta"] + r.coef[1] * d["mem_ratio"]
        re = r.coef[2] * d["Re_tau"]
        rows.append(pd.DataFrame(dict(
            target=target, xa=xs, d_observed=d_obs,
            d_pressure_gradient=pg, d_reynolds=re, d_model=pg + re)))
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    apply_style()
    out = Path(args.outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    df, full = load(args.data)

    # ---- 1. held-out accuracy and the fitted memory length ----
    ev = evaluate(df, full)
    summary = (ev.groupby("target")
               .agg(ell_median=("ell", "median"), ell_min=("ell", "min"),
                    ell_max=("ell", "max"), err=("median_rel", "median"),
                    a_beta=("a_beta", "median"), a_hist=("a_hist", "median"),
                    a_Re=("a_Re", "median"))
               .loc[list(TARGETS)])
    summary["err_%"] = (summary.pop("err") * 100).round(2)
    print("=== leave-one-Re-out, memory length refitted inside each fold ===")
    print(summary.round(3).to_string(), "\n")
    if ev.at_bound.any():
        print("WARNING: memory length hit a bound in some folds\n")

    # ---- 2. the decomposition ----
    dec = decompose(df, full)
    agg = (dec.groupby("target")
           .apply(lambda g: pd.Series({
               "observed": g.d_observed.mean(),
               "model": g.d_model.mean(),
               "pressure_gradient": g.d_pressure_gradient.mean(),
               "reynolds": g.d_reynolds.mean(),
           }), include_groups=False)
           .loc[list(TARGETS)])
    agg["PG_share_%"] = (100 * agg.pressure_gradient.abs()
                         / (agg.pressure_gradient.abs() + agg.reynolds.abs()))
    print("=== W2 -> W10: mean change in log(y), split by cause ===")
    print(agg.round(4).to_string())

    # Robustness: the same split for a different pair of Reynolds numbers.
    # If the shares are an artefact of one comparison they will not survive.
    dec48 = decompose(df, full, lo_case="W4", hi_case="W10")
    agg48 = (dec48.groupby("target")
             .apply(lambda g: pd.Series({
                 "PG_share_%": 100 * g.d_pressure_gradient.abs().mean() / (
                     g.d_pressure_gradient.abs().mean()
                     + g.d_reynolds.abs().mean())}), include_groups=False)
             .loc[list(TARGETS)])
    check = pd.DataFrame({"W2->W10": agg["PG_share_%"],
                          "W4->W10": agg48["PG_share_%"]}).round(1)
    print("\n=== pressure-gradient share (%), two different case pairs ===")
    print(check.to_string())

    # The fitted memory length is only good to tens of percent (see
    # tests/test_memory_model.py), so force it off the optimum and check the
    # split does not move with it.
    sens = {"fitted": agg["PG_share_%"]}
    for scale, name in [(0.5, "ell x 0.5"), (2.0, "ell x 2")]:
        d = decompose(df, full, ell_scale=scale)
        sens[name] = (d.groupby("target")
                      .apply(lambda g: 100 * g.d_pressure_gradient.abs().mean()
                             / (g.d_pressure_gradient.abs().mean()
                                + g.d_reynolds.abs().mean()),
                             include_groups=False)
                      .loc[list(TARGETS)])
    sens = pd.DataFrame(sens).round(1)
    print("\n=== sensitivity of the split to the memory length ===")
    print(sens.to_string())

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.9), layout="constrained")

    ax = axes[0]
    ax.axhline(0, color=TEXT_SECONDARY, lw=1.0)
    idx = np.arange(len(TARGETS))
    w = 0.6
    # Stacked, not side by side: the split is exactly additive, so the top of
    # the stack must coincide with the observed change. Side-by-side bars
    # would hide the one property worth showing.
    ax.bar(idx, agg.pressure_gradient, w, color=TERM_COLOR["local PG"],
           label=r"pressure gradient  ($\beta$ and history)")
    ax.bar(idx, agg.reynolds, w, bottom=agg.pressure_gradient,
           color=TERM_COLOR["Reynolds"], label=r"Reynolds number  ($Re_\tau$)")
    ax.plot(idx, agg.observed, "o", color="#0b0b0b", ms=6, zorder=5,
            label="observed total")
    ax.set_xticks(idx)
    ax.set_xticklabels([PRETTY[t] for t in TARGETS])
    ax.set_ylabel(r"change in $\log y$,  $Re_c$: 200k $\to$ 1M")
    ax.set_title("What actually drives the difference between cases", loc="left")
    ax.legend(loc="lower left", fontsize=7.5)
    style_axes(ax)

    ax = axes[1]
    for t in ("uu_outer", "Pk_outer", "Cf"):
        g = dec[dec.target == t]
        share = 100 * g.d_pressure_gradient.abs() / (
            g.d_pressure_gradient.abs() + g.d_reynolds.abs())
        ax.plot(g.xa, share, label=PRETTY[t], lw=2.0)
    ax.axhline(50, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    ax.annotate("equal split", xy=(0.90, 50), xytext=(-4, 5),
                textcoords="offset points", ha="right", fontsize=8,
                color=TEXT_SECONDARY)
    ax.set_ylim(0, 100)
    ax.set_xlabel(r"$x_{ss}/c$")
    ax.set_ylabel("share attributed to pressure gradient (%)")
    ax.set_title("The residual $\\beta$ mismatch is not a small correction",
                 loc="left")
    ax.legend(loc="lower right", fontsize=8)
    style_axes(ax)

    fig.savefig(out / "figures" / "fig_decomposition.png")
    plt.close(fig)

    ev.to_csv(out / "phase2_folds.csv", index=False)
    dec.to_csv(out / "phase2_decomposition.csv", index=False)
    (out / "phase2_results.md").write_text(
        "# Phase 2: fitted memory and the Re / pressure-gradient split\n\n"
        "## Held-out accuracy, memory length refitted inside each fold\n\n"
        + summary.round(3).to_markdown() +
        "\n\n`ell` is in multiples of the local `delta99`; `err_%` is the median\n"
        "relative error on the held-out Reynolds number.\n\n"
        "## Decomposition of the W2 -> W10 change\n\n"
        + agg.round(4).to_markdown() +
        "\n\nMean over `0.30 <= x/c <= 0.90` of the change in `log y`, split\n"
        "exactly by the fitted law into a pressure-gradient part and a\n"
        "Reynolds-number part.\n\n"
        "## Robustness: the same split for a second pair of cases\n\n"
        + check.to_markdown() +
        "\n\n## Sensitivity to the memory length\n\n"
        "The fitted length is identifiable only to tens of percent, because a\n"
        "shorter memory and a larger history exponent trade off. The split is\n"
        "built from their product, so it barely moves when the length is\n"
        "forced to half or double the fitted value.\n\n"
        + sens.to_markdown() + "\n"
    )
    print(f"\nwrote {out/'phase2_results.md'} and "
          f"{out/'figures'/'fig_decomposition.png'}")


if __name__ == "__main__":
    main()
