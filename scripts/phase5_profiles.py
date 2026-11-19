"""Phase 5: show the flow, not just the regression.

Two figures the paper needs.

**The non-equilibrium parameter itself.** Lambda is the central quantity of
the model and appeared in no figure. Plotted against position it separates
the two flow classes immediately: the wing's Lambda is positive throughout
and reaches 1.4, while all five flat plates sit within 0.17 of zero.

**A controlled profile comparison.** Plotting profiles at arbitrary stations
proves nothing. Instead we search for a wing station and a flat-plate station
at the *same* beta and the *same* Re_tau, and compare the profiles there. If
they still differ, the difference cannot be the local pressure gradient and
cannot be Reynolds number -- the remaining difference is Lambda. This is the
model's claim tested directly in the flow rather than inferred from a fit.

    python scripts/phase5_profiles.py
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.phase4_apg import (  # noqa: E402
    APG_COLOR, APG_WINDOW, apg_tables, third_term, wing_tables,
)
from wingtbl.apg_loader import load_apg  # noqa: E402
from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402

WING_COLOR = "#0d366b"
PLATE_COLOR = "#D85A30"


def find_pair(dw, da, tol=0.08):
    """Wing and flat-plate stations at matched beta and Re_tau.

    Among all pairs inside the tolerance, take the one whose Lambda differs
    most: that is the cleanest isolation of the one variable left free.
    """
    best = None
    for _, w in dw.iterrows():
        cand = da[(np.abs(da.beta - w.beta) / w.beta < tol)
                  & (np.abs(da.Re_tau - w.Re_tau) / w.Re_tau < tol)]
        for _, a in cand.iterrows():
            # Prefer a large Lambda separation, but break ties towards the
            # closest match in beta and Re_tau.
            quality = abs(w.Lam - a.Lam) - 0.5 * (
                abs(w.beta - a.beta) / w.beta + abs(w.Re_tau - a.Re_tau) / w.Re_tau)
            if best is None or quality > best[0]:
                best = (quality, w, a)
    return best[1], best[2]


def profile(ds, station):
    s = ds.isel(station=station)
    yn = s.yn.values
    ok = np.isfinite(yn) & (yn >= 0)
    return dict(
        yn=yn[ok], U=s.U.values[ok], uu=s.uu.values[ok], uv=s.uv.values[ok],
        u_tau=float(s.u_tau), nu=ds.attrs["nu"], d99=float(s.delta99),
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--wing", default="data/raw/naca4412.mat")
    p.add_argument("--apg", default="data/raw/apg.mat")
    p.add_argument("--stride", type=int, default=2)
    p.add_argument("--outdir", default="paper/figures")
    args = p.parse_args()

    warnings.filterwarnings("ignore")
    apply_style()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    dw, fw = wing_tables(args.wing)
    da, fa, ds_a = apg_tables(args.apg, args.stride)
    dw["Lam"] = third_term(dw, fw, "dbeta_dx")
    da["Lam"] = third_term(da, fa, "dbeta_dx")

    # ---------------- figure 1: Lambda ----------------
    fig, ax = plt.subplots(figsize=(6.6, 3.9), layout="constrained")
    ax.axhline(0, color=TEXT_SECONDARY, lw=1.0, ls=(0, (4, 3)))
    for d in ds_a:
        sub = da[da.case == d.attrs["case"]].sort_values("xa")
        xn = (sub.xa - APG_WINDOW.lo) / (APG_WINDOW.hi - APG_WINDOW.lo)
        ax.plot(xn, sub.Lam, color=APG_COLOR[d.attrs["case"]], lw=1.8,
                label=f"flat plate {d.attrs['case']}")
    for case, dash in [("W4", (None, None)), ("W10", (4, 2))]:
        sub = dw[dw.case == case].sort_values("xa")
        xn = (sub.xa - 0.25) / 0.70
        ax.plot(xn, sub.Lam, color=WING_COLOR, lw=2.4, dashes=dash,
                label=f"wing {case}")
    ax.annotate("equilibrium", xy=(1.0, 0), xytext=(-4, -14),
                textcoords="offset points", ha="right", fontsize=8,
                color=TEXT_SECONDARY)
    ax.set_xlabel("position through the analysis window")
    ax.set_ylabel(r"$\Lambda = (\delta_{99}/\beta)\,\mathrm{d}\beta/\mathrm{d}x_t$")
    ax.set_title("How far each layer sits from equilibrium", loc="left")
    ax.legend(fontsize=7.5, ncol=2, loc="upper left")
    style_axes(ax)
    fig.savefig(out / "fig_lambda.png")
    plt.close(fig)

    # ---------------- figure 2: matched-condition profiles ----------------
    w, a = find_pair(dw, da)
    raw_w = load_wing(args.wing)
    wvar = {"W2": "top2n", "W4": "top4n", "W10": "top10n"}[w.case]
    apg = load_apg(args.apg, stride=args.stride)
    avar = {d.attrs["case"]: k for k, d in apg.items()}[a.case]

    pw = profile(raw_w[wvar], int(w.station))
    pa = profile(apg[avar], int(a.station))

    print("matched pair")
    print(f"  wing {w.case} at x/c = {w.xa:.3f}: beta = {w.beta:.2f}, "
          f"Re_tau = {w.Re_tau:.0f}, Lambda = {w.Lam:+.3f}")
    print(f"  flat plate {a.case} at x = {a.xa:.0f}: beta = {a.beta:.2f}, "
          f"Re_tau = {a.Re_tau:.0f}, Lambda = {a.Lam:+.3f}")
    print(f"  mismatch: beta {100*abs(w.beta-a.beta)/w.beta:.1f}%, "
          f"Re_tau {100*abs(w.Re_tau-a.Re_tau)/w.Re_tau:.1f}%")

    fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.7), layout="constrained")

    ax = axes[0]
    for pr, c, lab in [(pw, WING_COLOR, f"wing {w.case}"),
                       (pa, PLATE_COLOR, f"flat plate {a.case}")]:
        yp = pr["yn"] * pr["u_tau"] / pr["nu"]
        m = yp > 0.8
        ax.semilogx(yp[m], pr["U"][m] / pr["u_tau"], color=c, lw=2.0, label=lab)
    ax.set_xlabel(r"$y^+$")
    ax.set_ylabel(r"$U^+$")
    ax.set_xlim(1, 1200)
    ax.set_title("Mean velocity, inner units", loc="left")
    ax.legend(fontsize=8, loc="upper left")
    style_axes(ax)

    for ax, key, lab, sign in [
        (axes[1], "uu", r"$\overline{u^2}^+$", 1.0),
        (axes[2], "uv", r"$-\overline{uv}^+$", -1.0),
    ]:
        for pr, c in [(pw, WING_COLOR), (pa, PLATE_COLOR)]:
            ax.plot(pr["yn"] / pr["d99"], sign * pr[key] / pr["u_tau"] ** 2,
                    color=c, lw=2.0)
        ax.set_xlabel(r"$y_n/\delta_{99}$")
        ax.set_ylabel(lab)
        ax.set_xlim(0, 1.6)
        style_axes(ax)
    axes[1].set_title("Streamwise stress", loc="left")
    axes[2].set_title("Shear stress", loc="left")

    fig.suptitle(
        rf"Matched $\beta \approx {w.beta:.1f}$ and $Re_\tau \approx {w.Re_tau:.0f}$, "
        rf"but $\Lambda = {w.Lam:+.2f}$ against ${a.Lam:+.2f}$",
        fontsize=10.5, fontweight="semibold")
    fig.savefig(out / "fig_matched_profiles.png")
    plt.close(fig)

    # Quantify the difference, so the caption can state it.
    yg = np.linspace(0.05, 1.0, 40)
    uu_w = np.interp(yg, pw["yn"] / pw["d99"], pw["uu"] / pw["u_tau"] ** 2)
    uu_a = np.interp(yg, pa["yn"] / pa["d99"], pa["uu"] / pa["u_tau"] ** 2)
    uv_w = np.interp(yg, pw["yn"] / pw["d99"], -pw["uv"] / pw["u_tau"] ** 2)
    uv_a = np.interp(yg, pa["yn"] / pa["d99"], -pa["uv"] / pa["u_tau"] ** 2)
    print(f"\n  outer-region difference at matched conditions:")
    print(f"    uu+  wing/plate at y/d99=0.2: "
          f"{np.interp(0.2,yg,uu_w):.2f} vs {np.interp(0.2,yg,uu_a):.2f} "
          f"({100*(np.interp(0.2,yg,uu_w)/np.interp(0.2,yg,uu_a)-1):+.0f}%)")
    print(f"    -uv+ wing/plate at y/d99=0.2: "
          f"{np.interp(0.2,yg,uv_w):.2f} vs {np.interp(0.2,yg,uv_a):.2f} "
          f"({100*(np.interp(0.2,yg,uv_w)/np.interp(0.2,yg,uv_a)-1):+.0f}%)")
    print(f"\nwrote {out/'fig_lambda.png'}, {out/'fig_matched_profiles.png'}")


if __name__ == "__main__":
    main()
