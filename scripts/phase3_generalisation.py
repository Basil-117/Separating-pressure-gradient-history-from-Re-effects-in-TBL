"""Phase 3a: does the fitted law transfer off the wing it was fitted to?

Every conclusion so far rests on three Reynolds numbers of one airfoil at one
angle of attack. Three tests, in increasing order of how much they ask:

    A  interpolate in Re     hold out Re_c = 400k, train on 200k and 1M
    B  extrapolate in beta   train on x/c <= 0.70 (beta <~ 2.5),
                             predict x/c > 0.70 (beta up to 46)
    C  change the geometry   train on NACA4412 (5 deg), predict NACA0012 (0 deg)

Test C is the real one. The NACA0012 is a different section at a different
angle of attack with a milder pressure gradient (beta reaches 3.3 rather
than 46), so nothing about it is in the training distribution except the
physics the law is supposed to have captured.

A fourth test is not possible and the reason is worth stating: the pressure
side of the wing is a *favourable* pressure gradient, beta < 0 at 88-95% of
its stations, and a law containing beta^a has no real value there. The form
is an adverse-pressure-gradient law and does not continue through beta = 0.

    python scripts/phase3_generalisation.py
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
from wingtbl.memory_model import MemoryPowerLaw  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402
from wingtbl.validation import summarise_errors  # noqa: E402

SIGN = {"uv_outer": -1.0}
PRETTY = {
    "Cf": r"$C_f$", "H": r"$H$", "Ue_plus": r"$U_e^+$",
    "uu_outer": r"$\overline{u_t^2}^+$", "uv_outer": r"$\overline{u_tv_n}^+$",
    "Pk_outer": r"$P_k^+$",
}
WIDE = AnalysisWindow(lo=-1.0, hi=2.0)
# Categorical slots 1-3, validated on the light surface.
TEST_COLOR = {"A interpolate in Re": "#2a78d6",
              "B extrapolate in beta": "#eb6834",
              "C change the geometry": "#1baf7a"}


def tables(data_4412: str, data_0012: str):
    raw = load_wing(data_4412)
    suction = [raw[v] for v in ("top2n", "top4n", "top10n")]
    n0012 = load_wing(data_0012)["top4n12"]

    full = {d.attrs["case"]: build_table([d], window=WIDE).reset_index(drop=True)
            for d in suction}
    # The NACA0012 ships one case whose label collides with nothing, but give
    # it an explicit name so the two airfoils cannot be confused downstream.
    g = build_table([n0012], window=WIDE).reset_index(drop=True)
    g["case"] = "N0012"
    full["N0012"] = g

    df4412 = build_table(suction)
    df0012 = build_table([n0012])
    df0012["case"] = "N0012"
    return df4412, df0012, full


def run(train: pd.DataFrame, test: pd.DataFrame, full: dict, label: str):
    rows = []
    for target in TARGETS:
        s = SIGN.get(target, 1.0)
        if test[target].isna().all():
            continue
        m = MemoryPowerLaw()
        r = m.fit(train.reset_index(drop=True), target, full, sign=s)
        pred = m.predict(test.reset_index(drop=True), full)
        truth = s * test[target].to_numpy()
        rows.append(dict(test=label, target=target, ell=r.ell,
                         a_beta=r.coef[0], a_hist=r.coef[1], a_Re=r.coef[2],
                         **summarise_errors(truth, pred)))
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="data/raw/naca4412.mat")
    p.add_argument("--data0012", default="data/raw/naca0012.mat")
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    apply_style()
    out = Path(args.outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    df4412, df0012, full = tables(args.data, args.data0012)

    rows = []
    # A: hold out the middle Reynolds number.
    rows += run(df4412[df4412.case != "W4"], df4412[df4412.case == "W4"],
                full, "A interpolate in Re")
    # B: extrapolate in pressure-gradient strength.
    rows += run(df4412[df4412.xa <= 0.70], df4412[df4412.xa > 0.70],
                full, "B extrapolate in beta")
    # C: a different airfoil at a different angle of attack.
    rows += run(df4412, df0012, full, "C change the geometry")

    res = pd.DataFrame(rows)
    tbl = (res.pivot(index="target", columns="test", values="median_rel") * 100
           ).loc[list(TARGETS)].round(2)
    print("=== median relative error (%) ===")
    print(tbl.to_string(), "\n")

    # Context for test B: how far outside the training range it reaches.
    b_tr = df4412[df4412.xa <= 0.70].beta
    b_te = df4412[df4412.xa > 0.70].beta
    print(f"test B: trained on beta in [{b_tr.min():.2f}, {b_tr.max():.2f}], "
          f"tested on [{b_te.min():.2f}, {b_te.max():.2f}]")
    print(f"test C: NACA0012 beta in [{df0012.beta.min():.2f}, "
          f"{df0012.beta.max():.2f}], Re_tau in [{df0012.Re_tau.min():.0f}, "
          f"{df0012.Re_tau.max():.0f}]")

    # ---- figure: predicted vs true on the hardest test ----
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), layout="constrained")

    ax = axes[0]
    idx = np.arange(len(TARGETS))
    w = 0.26
    for k, test in enumerate(tbl.columns):
        ax.bar(idx + (k - 1) * w, tbl[test].to_numpy(), w,
               color=TEST_COLOR[test], label=test)
    ax.set_xticks(idx)
    ax.set_xticklabels([PRETTY[t] for t in TARGETS])
    ax.set_yscale("log")
    ax.set_ylabel("median relative error (%)")
    ax.set_title("Generalisation, hardest test last", loc="left")
    ax.legend(fontsize=7.5, loc="upper left")
    style_axes(ax)

    # Parity on the cross-geometry test, the one that matters.
    ax = axes[1]
    shown = ["Cf", "H", "uv_outer"]
    for t, c in zip(shown, ("#2a78d6", "#eb6834", "#1baf7a")):
        s = SIGN.get(t, 1.0)
        m = MemoryPowerLaw()
        m.fit(df4412, t, full, sign=s)
        pred = m.predict(df0012.reset_index(drop=True), full)
        truth = s * df0012[t].to_numpy()
        # Normalised by the training geometric mean so three quantities with
        # different units share one pair of axes.
        ref = np.exp(np.log(s * df4412[t].to_numpy()).mean())
        ax.plot(truth / ref, pred / ref, "o", ms=4.5, color=c, alpha=0.8,
                label=PRETTY[t], markeredgecolor="#fcfcfb", markeredgewidth=0.6)
    lim = np.array([0.35, 2.6])
    ax.plot(lim, lim, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel("NACA0012 truth  (/ NACA4412 mean)")
    ax.set_ylabel("predicted from NACA4412")
    ax.set_title("Trained on one airfoil, tested on another", loc="left")
    ax.legend(fontsize=8, loc="upper left")
    style_axes(ax)

    fig.savefig(out / "figures" / "fig_generalisation.png")
    plt.close(fig)

    res.to_csv(out / "phase3_generalisation.csv", index=False)
    (out / "phase3_generalisation.md").write_text(
        "# Phase 3a: generalisation\n\n"
        "Median relative error (%), by test.\n\n" + tbl.to_markdown() +
        "\n\n- **A** holds out Re_c = 400k and trains on 200k and 1M.\n"
        f"- **B** trains on `x/c <= 0.70` (beta <= {b_tr.max():.1f}) and predicts\n"
        f"  `x/c > 0.70` (beta up to {b_te.max():.1f}).\n"
        "- **C** trains on the NACA4412 at 5 deg and predicts the NACA0012 at\n"
        "  0 deg, a different section with a far milder pressure gradient.\n\n"
        "The pressure side cannot be tested: it is a favourable pressure\n"
        "gradient with beta < 0 at 88-95% of stations, and a law containing\n"
        "beta^a has no real value there. The fitted form is an\n"
        "adverse-pressure-gradient law and does not continue through beta = 0.\n"
    )
    print(f"\nwrote {out/'phase3_generalisation.md'}")


if __name__ == "__main__":
    main()
