"""The residual pressure-gradient variation between the wing cases.

Vinuesa et al. (2018) compare the four wings at "approximately the same"
beta(x) and read the remaining differences as Reynolds-number effects. This
script shows that beta is in fact matched only to tens of percent, and that
the mismatch is systematic in Re -- so it runs in the same direction as the
effect being measured.

    python scripts/beta_confound.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.style import (  # noqa: E402
    CASE_COLOR, CASE_LABEL, TEXT_SECONDARY, apply_style, style_axes,
)

SUCTION = {"W2": "top2n", "W4": "top4n", "W10": "top10n"}
XLO, XHI = 0.25, 0.95


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--outdir", default="paper/figures")
    args = p.parse_args()

    apply_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = load_wing(args.data)
    x = np.linspace(XLO, XHI, 200)
    beta = {c: np.interp(x, raw[v].xa.values, raw[v].beta_db.values)
            for c, v in SUCTION.items()}

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), layout="constrained")

    # Left: beta on a log axis -- how the paper plots it.
    ax = axes[0]
    for c in ("W2", "W4", "W10"):
        ax.plot(x, beta[c], color=CASE_COLOR[c], label=CASE_LABEL[c])
    ax.set_yscale("log")
    ax.set_xlabel(r"$x_{ss}/c$")
    ax.set_ylabel(r"$\beta$")
    ax.set_title("As plotted in the paper (log axis)", loc="left")
    ax.legend(loc="upper left")
    style_axes(ax)

    # Right: the same three curves as a ratio to the highest-Re case.
    ax = axes[1]
    for c in ("W2", "W4"):
        ax.plot(x, beta[c] / beta["W10"], color=CASE_COLOR[c],
                label=f"{CASE_LABEL[c]} / {CASE_LABEL['W10']}")
    ax.axhline(1.0, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    ax.annotate("perfectly matched", xy=(XHI, 1.0), xytext=(-4, 5),
                textcoords="offset points", ha="right", va="bottom",
                fontsize=8, color=TEXT_SECONDARY)
    ax.set_xlabel(r"$x_{ss}/c$")
    ax.set_ylabel(r"$\beta / \beta_{Re_c=1\mathrm{M}}$")
    ax.set_ylim(0.95, 1.82)
    ax.set_title("The same data as a ratio", loc="left")
    ax.legend(loc="upper left")
    style_axes(ax)

    fig.suptitle(
        r"The matched-$\beta$ premise holds only to tens of percent",
        fontsize=10.5, fontweight="semibold",
    )
    out = outdir / "fig_beta_confound.png"
    fig.savefig(out)
    plt.close(fig)

    r2 = beta["W2"] / beta["W10"]
    r4 = beta["W4"] / beta["W10"]
    print(f"over {XLO} <= x/c <= {XHI}:")
    print(f"  W2 / W10 : mean {r2.mean():.3f}  range {r2.min():.3f}-{r2.max():.3f}")
    print(f"  W4 / W10 : mean {r4.mean():.3f}  range {r4.min():.3f}-{r4.max():.3f}")
    print(f"  beta decreases monotonically with Re at every station: "
          f"{bool(np.all((beta['W2'] > beta['W4']) & (beta['W4'] > beta['W10'])))}")
    print("wrote", out)


if __name__ == "__main__":
    main()
