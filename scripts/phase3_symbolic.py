"""Phase 3b: is the three-term power law the right complexity?

Searches exhaustively over expressions of the form

    log y = c0 + sum_i c_i t_i(beta, beta_bar/beta, Re_tau)

drawing the terms t_i from a small library, and builds a Pareto front of
complexity against error. Exhaustive rather than genetic: with a handful of
variables and a complexity cap the whole space is enumerable, which removes
the random seed, the population size and the stopping rule that a genetic
search would otherwise need defending. Within the space searched the result
is optimal rather than merely good.

The library is restricted to terms with a physical reading -- powers of the
three variables, a Reynolds-dependent exponent, curvature in log beta -- so
that whatever wins can be written down as a formula rather than admired as
a fit.

**Selection and validation are separated.** The front is chosen by
leave-one-Reynolds-number-out on the NACA4412, which is model selection and
not an unbiased error estimate. The chosen expressions are then evaluated
once on the NACA0012, which took no part in selection.

    python scripts/phase3_symbolic.py
"""

from __future__ import annotations

import argparse
from itertools import combinations
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
from wingtbl.memory_model import kernel_average  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402
from wingtbl.validation import leave_one_case_out, summarise_errors  # noqa: E402

SIGN = {"uv_outer": -1.0}
WIDE = AnalysisWindow(lo=-1.0, hi=2.0)
MAX_TERMS = 4
ELL_GRID = (0.2, 0.45, 1.0, 3.0)

# Term library. Each entry maps (Lb, Lm, LR, m) -> a column, where
# Lb = log(beta), Lm = log(beta_bar/beta), LR = log(Re_tau), m = beta_bar/beta.
LIBRARY: dict[str, callable] = {
    "log b":        lambda Lb, Lm, LR, m: Lb,
    "log m":        lambda Lb, Lm, LR, m: Lm,
    "log R":        lambda Lb, Lm, LR, m: LR,
    "log^2 b":      lambda Lb, Lm, LR, m: Lb ** 2,
    "log^2 R":      lambda Lb, Lm, LR, m: LR ** 2,
    "log b log R":  lambda Lb, Lm, LR, m: Lb * LR,
    "log m log R":  lambda Lb, Lm, LR, m: Lm * LR,
    "log b log m":  lambda Lb, Lm, LR, m: Lb * Lm,
    "m - 1":        lambda Lb, Lm, LR, m: m - 1.0,
}
NAMES = list(LIBRARY)

# What the three-term power law of Phase 2 corresponds to, for reference.
PHASE2_TERMS = ("log b", "log m", "log R")


def columns(df: pd.DataFrame, full: dict, ell: float) -> dict[str, np.ndarray]:
    beta_bar = np.concatenate([
        kernel_average(full[c]["xt"].to_numpy(), full[c]["beta"].to_numpy(),
                       full[c]["delta99"].to_numpy(), ell)[
            np.isin(full[c]["station"].to_numpy(),
                    df[df.case == c]["station"].to_numpy())]
        for c in dict.fromkeys(df["case"])
    ])
    beta = df["beta"].to_numpy()
    Lb, LR = np.log(beta), np.log(df["Re_tau"].to_numpy())
    m = beta_bar / beta
    Lm = np.log(m)
    return {k: f(Lb, Lm, LR, m) for k, f in LIBRARY.items()}


def lstsq_fit(X: np.ndarray, y: np.ndarray, alpha: float = 1e-8):
    Xc, yc = X - X.mean(0), y - y.mean()
    coef = np.linalg.solve(Xc.T @ Xc + alpha * np.eye(X.shape[1]), Xc.T @ yc)
    return coef, float(y.mean() - X.mean(0) @ coef)


def score_subset(df, full, target, terms, ell, sign):
    """Leave-one-Reynolds-number-out error for one expression."""
    cols = columns(df, full, ell)
    X = np.column_stack([cols[t] for t in terms])
    y = np.log(sign * df[target].to_numpy())
    errs = []
    for _, tr, te in leave_one_case_out(df):
        coef, b = lstsq_fit(X[tr], y[tr])
        pred = np.exp(X[te] @ coef + b)
        errs.append(summarise_errors(np.exp(y[te]), pred)["median_rel"])
    return float(np.median(errs))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--data0012", default="data/raw/naca0012.mat")
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    apply_style()
    out = Path(args.outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    raw = load_wing(args.data)
    suction = [raw[v] for v in ("top2n", "top4n", "top10n")]
    n0012 = load_wing(args.data0012)["top4n12"]
    full = {d.attrs["case"]: build_table([d], window=WIDE).reset_index(drop=True)
            for d in suction}
    g = build_table([n0012], window=WIDE).reset_index(drop=True)
    g["case"] = "N0012"
    full["N0012"] = g
    df = build_table(suction)
    df12 = build_table([n0012])
    df12["case"] = "N0012"

    subsets = [c for k in range(1, MAX_TERMS + 1)
               for c in combinations(NAMES, k)]
    print(f"searching {len(subsets)} expressions x {len(ELL_GRID)} memory "
          f"lengths x {len(TARGETS)} targets\n")

    front_rows, best_rows = [], []
    for target in TARGETS:
        sign = SIGN.get(target, 1.0)
        results = []
        for ell in ELL_GRID:
            cols = columns(df, full, ell)
            y = np.log(sign * df[target].to_numpy())
            for terms in subsets:
                X = np.column_stack([cols[t] for t in terms])
                errs = []
                for _, tr, te in leave_one_case_out(df):
                    coef, b = lstsq_fit(X[tr], y[tr])
                    errs.append(summarise_errors(
                        np.exp(y[te]), np.exp(X[te] @ coef + b))["median_rel"])
                results.append((len(terms), float(np.median(errs)), terms, ell))

        res = pd.DataFrame(results, columns=["k", "err", "terms", "ell"])
        # Pareto front: for each complexity, the best expression.
        front = res.loc[res.groupby("k").err.idxmin()].sort_values("k")
        # Keep only genuinely non-dominated points.
        keep, best_so_far = [], np.inf
        for _, r in front.iterrows():
            if r.err < best_so_far:
                keep.append(r)
                best_so_far = r.err
        front = pd.DataFrame(keep)

        for _, r in front.iterrows():
            front_rows.append(dict(target=target, k=int(r.k), err=r.err,
                                   ell=r.ell, terms=" + ".join(r.terms)))

        # Phase 2's own expression, for comparison at k = 3.
        p2 = res[(res.terms == PHASE2_TERMS)].err.min()

        def on_naca0012(terms, ell):
            """Fit on the NACA4412, score once on the unseen geometry."""
            cols_tr = columns(df, full, ell)
            cols_te = columns(df12, full, ell)
            X_tr = np.column_stack([cols_tr[t] for t in terms])
            X_te = np.column_stack([cols_te[t] for t in terms])
            coef, b = lstsq_fit(X_tr, np.log(sign * df[target].to_numpy()))
            e = summarise_errors(sign * df12[target].to_numpy(),
                                 np.exp(X_te @ coef + b))["median_rel"]
            return e, coef

        # Validate the selected expression once on the unseen geometry, and
        # put the three-term power law on exactly the same footing -- the
        # only comparison that is not contaminated by the search.
        chosen = front.iloc[-1]
        err12, coef = on_naca0012(chosen.terms, chosen.ell)
        ell_p2 = res[res.terms == PHASE2_TERMS].sort_values("err").iloc[0].ell
        p2_err12, _ = on_naca0012(PHASE2_TERMS, ell_p2)

        best_rows.append(dict(
            target=target, k=int(chosen.k), ell=chosen.ell,
            terms=" + ".join(chosen.terms), loo_err=chosen.err,
            phase2_loo=p2, naca0012_err=err12, phase2_naca0012=p2_err12,
            formula=" ".join(f"{c:+.3f}*[{t}]" for c, t in zip(coef, chosen.terms))))

    front_df = pd.DataFrame(front_rows)
    best = pd.DataFrame(best_rows).set_index("target").loc[list(TARGETS)]

    show = best[["k", "ell"]].copy()
    show["searched_loo_%"] = (best.loo_err * 100).round(2)
    show["searched_0012_%"] = (best.naca0012_err * 100).round(2)
    show["powerlaw_loo_%"] = (best.phase2_loo * 100).round(2)
    show["powerlaw_0012_%"] = (best.phase2_naca0012 * 100).round(2)
    show["optimism_x"] = (best.naca0012_err / best.loo_err).round(1)
    print("=== best searched expression vs the three-term power law ===")
    print(show.to_string(), "\n")
    print("terms chosen:")
    for t in TARGETS:
        print(f"  {t:10s} {best.loc[t,'terms']}")
    # Which library terms the search keeps choosing. A term picked for every
    # target, out of nine available and six independent selections, is far
    # more likely to be structure than noise.
    counts = pd.Series(
        [t for row in best.terms for t in row.split(" + ")]
    ).value_counts()
    print("\nterm frequency across the six targets (out of 6):")
    print(counts.to_string())

    wins = int((best.naca0012_err < best.phase2_naca0012).sum())
    print(f"\nOn the held-out geometry the searched expression beats the "
          f"three-term power law on {wins}/{len(TARGETS)} targets.")
    print(f"Selection optimism (held-out / selection error): median "
          f"{np.median(best.naca0012_err / best.loo_err):.1f}x")

    # ---- Pareto figure ----
    fig, ax = plt.subplots(figsize=(6.6, 4.1), layout="constrained")
    pretty = {"Cf": r"$C_f$", "H": r"$H$", "Ue_plus": r"$U_e^+$",
              "uu_outer": r"$\overline{u_t^2}^+$",
              "uv_outer": r"$\overline{u_tv_n}^+$", "Pk_outer": r"$P_k^+$"}
    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
    for (t, c) in zip(TARGETS, palette):
        s = front_df[front_df.target == t].sort_values("k")
        ax.plot(s.k, 100 * s.err, "o-", color=c, ms=5, label=pretty[t],
                markeredgecolor="#fcfcfb", markeredgewidth=1.0)
    ax.axvline(3, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    ax.annotate("three-term power law", xy=(3, ax.get_ylim()[0]),
                xytext=(-6, 12), textcoords="offset points", fontsize=8,
                color=TEXT_SECONDARY, va="bottom", ha="right")
    ax.set_yscale("log")
    ax.set_xticks(range(1, MAX_TERMS + 1))
    ax.set_xlabel("number of terms")
    ax.set_ylabel("leave-one-$Re$-out median error (%)")
    ax.set_title("Accuracy against complexity", loc="left")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    style_axes(ax)
    fig.savefig(out / "figures" / "fig_pareto.png")
    plt.close(fig)

    front_df.to_csv(out / "phase3_pareto.csv", index=False)
    (out / "phase3_symbolic.md").write_text(
        "# Phase 3b: symbolic search\n\n"
        f"Exhaustive over all {len(subsets)} expressions of up to "
        f"{MAX_TERMS} terms from a {len(NAMES)}-term library, at each of "
        f"{len(ELL_GRID)} memory lengths.\n\n"
        "## Best expression per target\n\n" + show.to_markdown() +
        "\n\n`searched_loo_%` was used for selection and is therefore\n"
        "optimistic. `searched_0012_%` and `powerlaw_0012_%` are both fitted\n"
        "on the NACA4412 and scored once on the NACA0012, which took no part\n"
        "in selection -- the only fair comparison of the two forms.\n"
        "`optimism_x` is how much worse the honest number is than the one\n"
        "selection reported.\n\n"
        "Terms chosen:\n\n"
        + "\n".join(f"- `{t}`: {best.loc[t,'terms']}" for t in TARGETS) + "\n\n"
        "## Term frequency across the six targets\n\n"
        + counts.to_frame("times chosen (of 6)").to_markdown() + "\n\n"
        "## Pareto front\n\n" + front_df.to_markdown(index=False) + "\n"
    )
    print(f"\nwrote {out/'phase3_symbolic.md'}, "
          f"{out/'figures'/'fig_pareto.png'}")


if __name__ == "__main__":
    main()
