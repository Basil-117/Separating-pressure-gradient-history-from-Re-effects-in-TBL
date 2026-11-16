"""Reproduce the figures and tables of Vinuesa et al. (2018) from the WING data.

This is the credibility anchor for the project: before modelling anything, show
that the published results come back out of the published files.

    python scripts/reproduce_paper.py --data data/raw/naca4412.mat \
        --outdir paper/figures

Produces:
    fig05_bl_development.png   beta, Re_theta, Re_tau   (paper Fig. 5)
    fig06_delta99.png          delta99/c                (paper Fig. 6)
    fig08_cf_H.png             Cf, H                    (paper Fig. 8)
    reproduction_report.md     tables 3 and 4, with errors
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
    CASE_COLOR, CASE_DASH, CASE_LABEL, CASE_ORDER, apply_style, style_axes,
)

SUCTION_VARS = {"W1": "top1n", "W2": "top2n", "W4": "top4n", "W10": "top10n"}

# Published values, Tables 3 (x/c = 0.4) and 4 (x/c = 0.7).
PAPER_TABLES = {
    0.4: {
        "W1": dict(beta=1.6, Re_tau=76, Re_theta=276, Cf=2.8e-3, H=2.04),
        "W2": dict(beta=0.67, Re_tau=160, Re_theta=452, Cf=4.4e-3, H=1.68),
        "W4": dict(beta=0.67, Re_tau=241, Re_theta=759, Cf=3.9e-3, H=1.61),
        "W10": dict(beta=0.58, Re_tau=449, Re_theta=1465, Cf=3.4e-3, H=1.50),
    },
    0.7: {
        "W1": dict(beta=3.7, Re_tau=136, Re_theta=566, Cf=3.0e-3, H=1.92),
        "W2": dict(beta=2.6, Re_tau=228, Re_theta=901, Cf=3.0e-3, H=1.73),
        "W4": dict(beta=2.4, Re_tau=360, Re_theta=1506, Cf=2.8e-3, H=1.68),
        "W10": dict(beta=2.0, Re_tau=671, Re_theta=2877, Cf=2.5e-3, H=1.58),
    },
}

# Lift, drag and efficiency, Table 2. Reported for completeness: these come
# from integrating forces over the whole wing, which the profile database does
# not contain, so they are quoted rather than recomputed.
PAPER_TABLE2 = {
    "W1": dict(Cl=0.8315, Cd=0.0242),
    "W2": dict(Cl=0.8854, Cd=0.0185),
    "W4": dict(Cl=0.9041, Cd=0.0162),
    "W10": dict(Cl=0.9112, Cd=0.0143),
}


def series(ds, name):
    """One per-station quantity, with H derived from the shipped thicknesses."""
    if name == "H":
        return (ds["delta_star"] / ds["theta"]).values
    return ds[name].values


def plot_panel(ax, cases, quantity, ylabel, logy=False):
    for cn in CASE_ORDER:
        ds = cases[cn]
        x = ds.xa.values
        y = np.array(series(ds, quantity), dtype=float)
        if logy:
            # beta is negative over the first few stations (the suction side
            # still accelerates there) and at W1's last station, where the
            # mean flow has separated so u_tau -> 0 and beta changes sign.
            # A log axis cannot show those; mask rather than let matplotlib
            # draw a spurious vertical excursion.
            y = np.where(y > 0, y, np.nan)
        ax.plot(
            x, y,
            color=CASE_COLOR[cn],
            dashes=CASE_DASH[cn] or (None, None),
            label=CASE_LABEL[cn],
            solid_capstyle="round",
        )
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(r"$x_{ss}/c$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0.2, 1.0)
    style_axes(ax)


def fig_bl_development(cases, outdir: Path) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.5), layout="constrained")
    plot_panel(axes[0], cases, "beta_db", r"$\beta$", logy=True)
    plot_panel(axes[1], cases, "Re_theta", r"$Re_\theta$")
    plot_panel(axes[2], cases, "Re_tau", r"$Re_\tau$")
    axes[0].set_ylim(1e-2, 3e2)
    axes[0].set_title("Clauser pressure-gradient parameter", loc="left")
    axes[1].set_title("Momentum-thickness Reynolds number", loc="left")
    axes[2].set_title("Friction Reynolds number", loc="left")
    axes[0].legend(loc="upper left", ncol=1)
    fig.suptitle(
        "Boundary-layer development on the suction side  "
        "(reproduction of Vinuesa et al. 2018, Fig. 5)",
        fontsize=10.5, fontweight="semibold",
    )
    out = outdir / "fig05_bl_development.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def fig_delta99(cases, outdir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    plot_panel(ax, cases, "delta99", r"$\delta_{99}/c$")
    ax.set_title("Boundary-layer thickness (paper Fig. 6)", loc="left")
    ax.legend(loc="upper left")
    out = outdir / "fig06_delta99.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def fig_cf_H(cases, outdir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), layout="constrained")
    plot_panel(axes[0], cases, "Cf", r"$C_f$")
    plot_panel(axes[1], cases, "H", r"$H$")
    axes[0].set_title("Skin-friction coefficient", loc="left")
    axes[1].set_title("Shape factor", loc="left")
    axes[0].set_ylim(0, 8e-3)
    axes[1].set_ylim(1.2, 3.0)
    axes[1].legend(loc="upper left")
    fig.suptitle(
        "Reproduction of Vinuesa et al. 2018, Fig. 8",
        fontsize=10.5, fontweight="semibold",
    )
    out = outdir / "fig08_cf_H.png"
    fig.savefig(out)
    plt.close(fig)
    return out


# Maxima quoted in the text of section 3.1. The paper restricts its analysis
# to x_ss/c < 0.95, and the Re_theta figures match that cutoff rather than the
# last station in the database (x/c = 0.983).
PAPER_MAXIMA = {
    "Re_theta_at_0.95": {"W1": 1050, "W2": 1730, "W4": 2930, "W10": 6000},
    "Re_tau_max": {"W1": 139, "W2": 231, "W4": 366, "W10": 707},
}


def maxima_report(cases) -> tuple[str, np.ndarray]:
    lines, errs = [], []
    lines.append("\n### Maxima quoted in section 3.1\n")
    lines.append("| quantity | case | paper | reproduced | diff |")
    lines.append("|---|---|---|---|---|")
    for cn in CASE_ORDER:
        ds = cases[cn]
        x = ds.xa.values
        got = float(np.interp(0.95, x, ds["Re_theta"].values))
        ref = PAPER_MAXIMA["Re_theta_at_0.95"][cn]
        d = 100 * (got - ref) / ref
        errs.append(abs(d))
        lines.append(f"| Re_theta at x/c=0.95 | {cn} | {ref} | {got:.0f} | {d:+.1f}% |")
    for cn in CASE_ORDER:
        ds = cases[cn]
        got = float(ds["Re_tau"].values.max())
        ref = PAPER_MAXIMA["Re_tau_max"][cn]
        d = 100 * (got - ref) / ref
        errs.append(abs(d))
        lines.append(f"| max Re_tau | {cn} | {ref} | {got:.0f} | {d:+.1f}% |")
    return "\n".join(lines), np.array(errs)


def table_report(cases) -> tuple[str, np.ndarray]:
    lines, errs = [], []
    for xt, table in PAPER_TABLES.items():
        lines.append(f"\n### x_ss/c = {xt} "
                     f"(paper Table {3 if xt == 0.4 else 4})\n")
        lines.append("| case | quantity | paper | reproduced | diff |")
        lines.append("|---|---|---|---|---|")
        for cn in CASE_ORDER:
            ds = cases[cn]
            x = ds.xa.values
            for q, fmt in [("beta", "{:.3g}"), ("Re_tau", "{:.0f}"),
                           ("Re_theta", "{:.0f}"), ("Cf", "{:.3e}"),
                           ("H", "{:.3f}")]:
                key = "beta_db" if q == "beta" else q
                got = float(np.interp(xt, x, series(ds, key)))
                ref = table[cn][q]
                d = 100 * (got - ref) / ref
                errs.append(abs(d))
                lines.append(
                    f"| {cn} | {q} | {fmt.format(ref)} | {fmt.format(got)} "
                    f"| {d:+.1f}% |"
                )
    return "\n".join(lines), np.array(errs)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--outdir", default="paper/figures")
    args = p.parse_args()

    apply_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = load_wing(args.data)
    cases = {cn: raw[v] for cn, v in SUCTION_VARS.items()}

    made = [fig_bl_development(cases, outdir),
            fig_delta99(cases, outdir),
            fig_cf_H(cases, outdir)]

    body, errs_t = table_report(cases)
    mbody, errs_m = maxima_report(cases)
    body = body + "\n" + mbody
    errs = np.concatenate([errs_t, errs_m])
    summary = (
        f"- quantities compared: **{errs.size}**\n"
        f"- median absolute error: **{np.median(errs):.2f}%**\n"
        f"- mean absolute error: **{errs.mean():.2f}%**\n"
        f"- within 2%: **{(errs < 2).sum()}/{errs.size}**, "
        f"within 5%: **{(errs < 5).sum()}/{errs.size}**\n"
    )
    report = (
        "# Reproduction of Vinuesa et al. (2018)\n\n"
        "Values are interpolated to the exact chordwise station quoted in the\n"
        "paper. Residual differences in `beta` are dominated by the paper\n"
        "quoting it to two significant figures.\n\n"
        "## Summary\n\n" + summary + "\n## Tables\n" + body + "\n"
    )
    (outdir.parent / "reproduction_report.md").write_text(report)

    print(summary)
    for m in made:
        print("wrote", m)
    print("wrote", outdir.parent / "reproduction_report.md")


if __name__ == "__main__":
    main()
