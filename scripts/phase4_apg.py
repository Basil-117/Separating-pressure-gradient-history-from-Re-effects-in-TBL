"""Phase 4: the flat-plate test, and what it says about the history term.

Everything before this came from wings, where each geometry supplies exactly
one beta(x) history. The flat-plate database of Bobke et al. (2017) supplies
five deliberately different ones -- two at nominally constant beta, three
following a power-law freestream -- at 1570 stations per case instead of 50.

It settles three questions, and the answers revise the earlier story.

1. **Was the memory length ever resolved?** No. On the wing grid the fitted
   ell was 0.37-0.56 times the spacing *between stations*, so the kernel
   decayed inside a single quadrature panel. The flat plates put ell/dx above
   20, and are the first data here on which a memory length is measurable.

2. **What is the third term actually measuring?** The local rate of change of
   beta, not a memory. On the wing, replacing the kernel with a plain
   (d beta/dx) term gives interchangeable accuracy. The kernel at sub-grid
   ell is a finite-difference gradient wearing a convolution's clothes.

3. **Does that matter for the headline result?** The split stays
   qualitatively the same -- outer-region Reynolds stresses are
   pressure-gradient dominated, wall and integral quantities are
   Reynolds-number dominated -- but the shares move, and one of them crosses
   50%. Reported both ways.

The physical reading is coherent: the term measures *departure from
equilibrium*. The wing's beta climbs from 0.3 to 46 across the window, so it
is far from equilibrium and the term is essential. The flat plates were
built to sit near equilibrium, their beta is nearly flat, and the term
correctly contributes almost nothing.

    python scripts/phase4_apg.py
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

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wingtbl.apg_loader import load_apg  # noqa: E402
from wingtbl.features import AnalysisWindow, TARGETS, build_table  # noqa: E402
from wingtbl.loader import load_wing  # noqa: E402
from wingtbl.memory_model import kernel_average  # noqa: E402
from wingtbl.style import apply_style, style_axes, TEXT_SECONDARY  # noqa: E402
from wingtbl.validation import leave_one_case_out, summarise_errors  # noqa: E402

SIGN = {"uv_outer": -1.0}
PRETTY = {
    "Cf": r"$C_f$", "H": r"$H$", "Ue_plus": r"$U_e^+$",
    "uu_outer": r"$\overline{u^2}^+$", "uv_outer": r"$\overline{uv}^+$",
    "Pk_outer": r"$P_k^+$",
}

# Below x ~ 700 the pressure gradient is still being established; above
# x ~ 1900 beta collapses towards zero as the outflow boundary is
# approached. Between them each case sits on its intended history.
APG_WINDOW = AnalysisWindow(lo=700.0, hi=1900.0)
WING_WINDOW = AnalysisWindow(lo=0.25, hi=0.95)
WIDE = AnalysisWindow(lo=-1e9, hi=1e9)

APG_COLOR = {"b1": "#85B7EB", "b2": "#378ADD", "m13": "#F0997B",
             "m16": "#D85A30", "m18": "#993C1D"}


def apg_tables(path: str, stride: int):
    apg = load_apg(path, stride=stride)
    datasets = list(apg.values())
    full = {d.attrs["case"]: build_table([d], window=WIDE).reset_index(drop=True)
            for d in datasets}
    return build_table(datasets, window=APG_WINDOW), full, datasets


def wing_tables(path: str):
    raw = load_wing(path)
    datasets = [raw[v] for v in ("top2n", "top4n", "top10n")]
    full = {d.attrs["case"]: build_table([d], window=WIDE).reset_index(drop=True)
            for d in datasets}
    return build_table(datasets, window=WING_WINDOW), full


def resolution_ratio(full: dict, df: pd.DataFrame, ell_over_d99: float) -> float:
    """``ell / dx`` at a given memory length: below 1 the kernel is unresolved."""
    dxs, d99s = [], []
    for c in dict.fromkeys(df.case):
        sub = df[df.case == c]
        dxs.append(np.median(np.diff(np.sort(sub.xt.to_numpy()))))
        d99s.append(np.median(sub.delta99.to_numpy()))
    return float(ell_over_d99 * np.median(d99s) / np.median(dxs))


def third_term(df: pd.DataFrame, full: dict, kind: str,
               ell: float = 0.5) -> np.ndarray:
    """The non-local feature, either as a kernel or as a plain gradient."""
    out = []
    for c in dict.fromkeys(df.case):
        g = full[c]
        keep = np.isin(g.station.to_numpy(), df[df.case == c].station.to_numpy())
        beta = g.beta.to_numpy()
        if kind == "memory":
            bb = kernel_average(g.xt.to_numpy(), beta, g.delta99.to_numpy(), ell)
            out.append(np.log(bb / beta)[keep])
        elif kind == "dbeta_dx":
            # (d beta / dx) * delta99 / beta -- dimensionless, and the
            # first-order Taylor limit of the kernel as ell -> 0.
            v = np.gradient(beta, g.xt.to_numpy()) * g.delta99.to_numpy() / beta
            out.append(v[keep])
        else:
            raise ValueError(kind)
    return np.concatenate(out)


def design(df, full, kind, ell=0.5):
    cols = [np.log(df.beta.to_numpy()), np.log(df.Re_tau.to_numpy())]
    if kind != "none":
        cols.insert(1, third_term(df, full, kind, ell))
    return np.column_stack(cols)


def _solve(X, y, alpha=1e-8):
    Xc = X - X.mean(0)
    coef = np.linalg.solve(Xc.T @ Xc + alpha * np.eye(X.shape[1]),
                           Xc.T @ (y - y.mean()))
    return coef, float(y.mean() - X.mean(0) @ coef)


def score(df, full, target, kind, ell=0.5) -> float:
    s = SIGN.get(target, 1.0)
    y = np.log(s * df[target].to_numpy())
    X = design(df, full, kind, ell)
    errs = []
    for _, tr, te in leave_one_case_out(df):
        coef, b = _solve(X[tr], y[tr])
        errs.append(summarise_errors(np.exp(y[te]),
                                     np.exp(X[te] @ coef + b))["median_rel"])
    return 100 * float(np.median(errs))


def pg_share(df, full, target, kind, lo="W2", hi="W10") -> float:
    """Share of the lo->hi change attributed to the pressure gradient."""
    s = SIGN.get(target, 1.0)
    X = design(df, full, kind)
    y = np.log(s * df[target].to_numpy())
    coef, _ = _solve(X, y)
    xs = np.linspace(0.30, 0.90, 60)

    def interp(case, col):
        m = (df.case == case).to_numpy()
        o = np.argsort(df.xa.to_numpy()[m])
        return np.interp(xs, df.xa.to_numpy()[m][o], X[m, col][o])

    d = [interp(hi, k) - interp(lo, k) for k in range(3)]
    pg = np.abs(coef[0] * d[0] + coef[1] * d[1]).mean()
    re = np.abs(coef[2] * d[2]).mean()
    return 100 * pg / (pg + re)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--apg", default="data/raw/apg.mat")
    p.add_argument("--wing", default="data/raw/naca4412.mat")
    p.add_argument("--stride", type=int, default=4)
    p.add_argument("--outdir", default="paper")
    args = p.parse_args()

    warnings.filterwarnings("ignore")
    apply_style()
    out = Path(args.outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    df_a, full_a, ds_a = apg_tables(args.apg, args.stride)
    df_w, full_w = wing_tables(args.wing)
    print(f"flat plates: {len(df_a)} stations, {df_a.case.nunique()} histories")
    print(f"wing:        {len(df_w)} stations, {df_w.case.nunique()} cases\n")

    # ---- 1. was the memory length ever resolved? ----
    rw = resolution_ratio(full_w, df_w, 0.5)
    ra = resolution_ratio(full_a, df_a, 0.5)
    print("=== resolution of a memory length of 0.5 delta99 ===")
    print(f"  wing        ell/dx = {rw:5.2f}   {'UNRESOLVED' if rw < 1 else 'ok'}")
    print(f"  flat plates ell/dx = {ra:5.2f}   {'UNRESOLVED' if ra < 1 else 'ok'}\n")

    # ---- 2. does the third term earn its place, and as what? ----
    rows = []
    for label, (df, full) in [("wing", (df_w, full_w)),
                              ("flat plate", (df_a, full_a))]:
        for t in TARGETS:
            if df[t].isna().all():
                continue
            rows.append(dict(
                dataset=label, target=t,
                local_only=score(df, full, t, "none"),
                memory=score(df, full, t, "memory"),
                dbeta_dx=score(df, full, t, "dbeta_dx"),
            ))
    comp = pd.DataFrame(rows)
    for label in ("wing", "flat plate"):
        sub = comp[comp.dataset == label].set_index("target")[
            ["local_only", "memory", "dbeta_dx"]]
        print(f"=== {label}: median rel. error (%), leave-one-case-out ===")
        print(sub.round(2).to_string(), "\n")

    # ---- 3. does the headline split survive the reinterpretation? ----
    shares = pd.DataFrame({
        "memory": {t: pg_share(df_w, full_w, t, "memory") for t in TARGETS},
        "dbeta_dx": {t: pg_share(df_w, full_w, t, "dbeta_dx") for t in TARGETS},
    }).loc[list(TARGETS)]
    print("=== pressure-gradient share of the W2 -> W10 change (%) ===")
    print(shares.round(1).to_string())

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.9), layout="constrained")

    ax = axes[0]
    for d in ds_a:
        x = d.xa.values
        m = (x >= APG_WINDOW.lo) & (x <= APG_WINDOW.hi)
        ax.plot((x[m] - APG_WINDOW.lo) / (APG_WINDOW.hi - APG_WINDOW.lo),
                d.beta_db.values[m], color=APG_COLOR[d.attrs["case"]], lw=2.0,
                label=d.attrs["case"])
    dw = load_wing(args.wing)["top4n"]
    x = dw.xa.values
    m = (x >= 0.25) & (x <= 0.95)
    ax.plot((x[m] - 0.25) / 0.70, dw.beta_db.values[m], color="#0b0b0b",
            lw=2.2, dashes=(4, 2), label="wing")
    ax.set_yscale("log")
    ax.set_xlabel("position through the analysis window")
    ax.set_ylabel(r"$\beta$")
    ax.set_title("Near equilibrium, and far from it", loc="left")
    ax.legend(fontsize=7.5, ncol=2, loc="upper left")
    style_axes(ax)

    ax = axes[1]
    idx = np.arange(len(TARGETS))
    w = 0.38
    for lab, (df, full), off, col in [
        ("wing", (df_w, full_w), -w / 2, "#0d366b"),
        ("flat plate", (df_a, full_a), w / 2, "#86b6ef"),
    ]:
        gain = []
        for t in TARGETS:
            base = score(df, full, t, "none")
            best = min(score(df, full, t, "memory"), score(df, full, t, "dbeta_dx"))
            gain.append(100 * (base - best) / base)
        ax.bar(idx + off, gain, w, color=col, label=lab)
    ax.axhline(0, color=TEXT_SECONDARY, lw=1.0)
    ax.set_xticks(idx)
    ax.set_xticklabels([PRETTY[t] for t in TARGETS])
    ax.set_ylabel("error removed by the third term (%)")
    ax.set_title("The term earns its place only off equilibrium", loc="left")
    ax.set_ylim(top=112)
    ax.legend(fontsize=8, loc="upper right")
    style_axes(ax)

    fig.savefig(out / "figures" / "fig_apg.png")
    plt.close(fig)

    comp.to_csv(out / "phase4_comparison.csv", index=False)
    (out / "phase4_results.md").write_text(
        "# Phase 4: the flat-plate test\n\n"
        "## Was the memory length ever resolved?\n\n"
        f"Ratio of a 0.5 `delta99` memory length to the station spacing:\n"
        f"wing **{rw:.2f}**, flat plates **{ra:.2f}**. Below 1 the kernel decays\n"
        "inside a single quadrature panel, so the wing fit reported a length\n"
        "its own grid could not represent.\n\n"
        "## Does the third term earn its place?\n\n"
        + comp.round(2).to_markdown(index=False) +
        "\n\nOn the wing it removes most of the error; on the near-equilibrium\n"
        "flat plates it does almost nothing, and the kernel and the plain\n"
        "`d beta/dx` gradient are interchangeable throughout.\n\n"
        "## Does the headline split survive?\n\n"
        + shares.round(1).to_markdown() +
        "\n\nThe ordering is unchanged, but `uu_outer` crosses 50%, so that\n"
        "share is parameterisation dependent and must be quoted as a range.\n"
    )
    print(f"\nwrote {out/'phase4_results.md'}, {out/'figures'/'fig_apg.png'}")


if __name__ == "__main__":
    main()
