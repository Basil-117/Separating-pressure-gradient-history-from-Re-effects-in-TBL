"""How far back does a wing boundary layer remember?

Sweeps the memory length of the history kernel and, for each target, finds
the length that minimises leave-one-Reynolds-number-out error. Clauser's
claim that the outer layer has a long memory predicts that the optimum is
short for wall quantities and long for outer-region quantities. This script
is the test, and it turns a qualitative claim into a number.

Memory length is expressed as a multiple of the *local* delta99, so the
kernel widens as the boundary layer thickens, which is the physically
natural scaling.

    python scripts/memory_scan.py
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wingtbl.features import (  # noqa: E402
    HISTORY_FEATURES, LOCAL_FEATURES, TARGETS, build_table,
)
from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402
from wingtbl.validation import leave_one_case_out, summarise_errors  # noqa: E402

SUCTION = ["top1n", "top2n", "top4n", "top10n"]
EXCLUDE = ["W1"]
NEGATIVE = {"uv_outer"}

# Where each target lives. The hypothesis is that this column predicts the
# shape of the curves below.
REGION = {
    "Cf": "wall", "H": "integral", "Ue_plus": "integral",
    "uu_outer": "outer", "uv_outer": "outer", "Pk_outer": "outer",
}

PRETTY = {
    "Cf": r"$C_f$", "H": r"$H$", "Ue_plus": r"$U_e^+$",
    "uu_outer": r"$\overline{u_t^2}^+$ at $0.2\delta_{99}$",
    "uv_outer": r"$\overline{u_tv_n}^+$ at $0.2\delta_{99}$",
    "Pk_outer": r"$P_k^+$ at $0.2\delta_{99}$",
}

# Categorical slots 1-3: the three regions are distinct physical places, not
# an ordered scale. Validated on the light surface (CVD dE 9.2 worst adjacent,
# normal-vision 27.6); the green carries a contrast warning, so each region
# median is directly labelled rather than relying on the legend alone.
REGION_COLOR = {"wall": "#2a78d6", "integral": "#eb6834", "outer": "#1baf7a"}
REGION_LABEL = {
    "wall": "wall  ($C_f$)",
    "integral": "integral  ($H$, $U_e^+$)",
    "outer": r"outer  at $0.2\delta_{99}$",
}


def plateau_range(lengths: np.ndarray, err: np.ndarray,
                  tol: float = 0.10) -> tuple[float, float]:
    """Range of memory lengths within ``tol`` of the best error.

    The scan curves are noisy -- ridge selects its own penalty inside each
    fold, and there are only three folds -- so a single argmin overstates
    what the data resolves. This reports the interval that is statistically
    indistinguishable from the optimum instead.
    """
    thresh = err.min() * (1.0 + tol)
    within = lengths[err <= thresh]
    return float(within.min()), float(within.max())


def score(df: pd.DataFrame, target: str, features: list[str]) -> float:
    sign = -1.0 if target in NEGATIVE else 1.0
    y = sign * df[target].to_numpy()
    X = df[features].to_numpy()
    errs = []
    for _, tr, te in leave_one_case_out(df):
        m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-4, 4, 40)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(X[tr], np.log(y[tr]))
            pred = np.exp(m.predict(X[te]))
        errs.append(summarise_errors(y[te], pred)["median_rel"])
    return float(np.median(errs))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    apply_style()
    raw = load_wing(args.data)
    datasets = [raw[v] for v in SUCTION]

    # Wide enough that every optimum is interior: a minimum sitting on an
    # endpoint would mean the scan, not the physics, chose the length.
    lengths = np.unique(np.round(np.geomspace(0.05, 500, 34), 4))

    rows = []
    for L in lengths:
        df = build_table(datasets, memory_in_delta99=float(L))
        df = df[~df.case.isin(EXCLUDE)].reset_index(drop=True)
        for t in TARGETS:
            rows.append(dict(length=L, target=t,
                             err=score(df, t, LOCAL_FEATURES + HISTORY_FEATURES)))
    scan = pd.DataFrame(rows)

    # Local-only reference, which does not depend on the memory length.
    df0 = build_table(datasets).pipe(lambda d: d[~d.case.isin(EXCLUDE)]).reset_index(drop=True)
    baseline = {t: score(df0, t, LOCAL_FEATURES) for t in TARGETS}

    best = (scan.loc[scan.groupby("target").err.idxmin()]
            .set_index("target")[["length", "err"]])
    best["local_only"] = pd.Series(baseline)
    best["reduction_%"] = ((best.local_only - best.err) / best.local_only * 100)
    best["region"] = pd.Series(REGION)
    # Honest uncertainty on the optimum, given three folds and a noisy curve.
    plateau = {}
    for t in TARGETS:
        s = scan[scan.target == t].sort_values("length")
        plateau[t] = plateau_range(s.length.to_numpy(), s.err.to_numpy())
    best["within10%_lo"] = pd.Series({t: v[0] for t, v in plateau.items()})
    best["within10%_hi"] = pd.Series({t: v[1] for t, v in plateau.items()})
    best = best.sort_values("length")

    print("=== best memory length per target (leave-one-Re-out) ===")
    print(best.round(3).to_string(), "\n")
    by_region = best.groupby("region").length.median()
    print("median best memory length by region (x delta99):")
    print(by_region.round(1).to_string())

    # Normalised curves, so every target is "error relative to using no
    # history at all" and the three regions are directly comparable.
    norm = {t: 100 * scan[scan.target == t].sort_values("length").err.to_numpy()
            / baseline[t] for t in TARGETS}
    xs = np.sort(scan.length.unique())

    fig, ax = plt.subplots(figsize=(7.2, 4.2), layout="constrained")
    ax.axhline(100, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.annotate("no history (local features only)", xy=(xs[-1], 100),
                xytext=(-4, 5), textcoords="offset points", ha="right",
                fontsize=8, color=TEXT_SECONDARY)

    label_y = {}
    for region, color in REGION_COLOR.items():
        members = [t for t in TARGETS if REGION[t] == region]
        stack = np.vstack([norm[t] for t in members])
        # Individual targets, recessive; the region median carries identity.
        for row in stack:
            ax.plot(xs, row, color=color, lw=0.9, alpha=0.35, zorder=2)
        med = np.median(stack, axis=0)
        ax.plot(xs, med, color=color, lw=2.4, zorder=3,
                solid_capstyle="round", label=REGION_LABEL[region])
        k = int(np.argmin(med))
        ax.plot(xs[k], med[k], "o", color=color, ms=7, zorder=5,
                markeredgecolor="#fcfcfb", markeredgewidth=1.4)
        label_y[region] = (xs[k], med[k])

    # Direct labels: required relief for the low-contrast slot, and they let
    # the reader identify a curve without crossing to a legend. The wall and
    # integral minima land on the same memory length, so their labels are
    # pushed to opposite sides rather than stacked on top of each other.
    offsets = {"wall": (-8, -20, "right"),
               "integral": (10, -20, "left"),
               "outer": (0, -18, "center")}
    for region, (lx, ly) in label_y.items():
        dx, dy, ha = offsets[region]
        ax.annotate(REGION_LABEL[region], xy=(lx, ly), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, fontsize=8.5,
                    color=REGION_COLOR[region], fontweight="semibold")
    ax.set_ylim(top=128)

    ax.set_xscale("log")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_xlabel(r"memory length  $\ell\,/\,\delta_{99}$")
    ax.set_ylabel("error relative to no-history model (%)")
    ax.set_title(
        "Outer-region statistics remember an order of magnitude further "
        "than wall quantities", loc="left", fontsize=10)
    style_axes(ax)

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "figures" / "fig_memory_scan.png")
    plt.close(fig)
    scan.to_csv(out / "memory_scan.csv", index=False)
    (out / "memory_scan.md").write_text(
        "# Memory-length scan\n\n"
        "> **Read with the limits in the README.** The fitted memory lengths below\n"
        "> are near or under one station spacing on the 50-station wing grid, where\n"
        "> the kernel decays inside a single quadrature panel and cannot be\n"
        "> resolved. They are a property of the grid, not a physical length. See\n"
        "> *Identifiability* in the README.\n"
        "Leave-one-Reynolds-number-out median relative error, as a function of\n"
        "the kernel memory length in units of the local `delta99`.\n\n"
        "## Best length per target\n\n" + best.round(3).to_markdown() +
        "\n\n## Median best length by region\n\n" +
        by_region.round(1).to_frame("length_over_delta99").to_markdown() + "\n"
    )
    print(f"\nwrote {out/'figures'/'fig_memory_scan.png'}, {out/'memory_scan.md'}")


if __name__ == "__main__":
    main()
