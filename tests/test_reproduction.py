"""Regression tests: the published results must come back out of the files.

These are the tests that stop a refactor from quietly changing the physics.
Every tolerance below is justified in docs/DATA_NOTES.md.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from wingtbl import blquantities as blq
from wingtbl.loader import load_wing, parse_case_name

DATA = Path(__file__).resolve().parents[1] / "data" / "raw"
NACA4412 = DATA / "naca4412.mat"
NACA0012 = DATA / "naca0012.mat"

pytestmark = pytest.mark.skipif(
    not NACA4412.exists(),
    reason="WING data not present; see docs/DATA_NOTES.md section 7",
)

SUCTION = {"W1": "top1n", "W2": "top2n", "W4": "top4n", "W10": "top10n"}

# Vinuesa et al. (2018), Tables 3 and 4.
TABLES = {
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


@pytest.fixture(scope="module")
def cases():
    raw = load_wing(NACA4412)
    return {cn: raw[v] for cn, v in SUCTION.items()}


# -- loader ---------------------------------------------------------------

def test_parse_case_name():
    assert parse_case_name("top10n") == ("top", 1_000_000)
    assert parse_case_name("bot1n") == ("bot", 100_000)
    # NACA0012 appends the airfoil to the variable name.
    assert parse_case_name("top4n12") == ("top", 400_000)
    assert parse_case_name("nonsense") is None
    assert parse_case_name("top3n") is None


def test_all_cases_present():
    raw = load_wing(NACA4412)
    assert len(raw) == 8
    assert set(raw) == {f"{s}{r}n" for s in ("top", "bot") for r in (1, 2, 4, 10)}
    assert set(load_wing(NACA0012)) == {"top4n12"}


def test_nu_matches_reynolds_number(cases):
    """nu = 1/Re_c for unit chord and freestream -- checks the Re code."""
    for ds in cases.values():
        assert np.isclose(ds.attrs["nu"] * ds.attrs["Re_c"], 1.0, rtol=1e-6)


def test_pressure_side_budgets_are_nan_not_crash():
    """The distributed file has empty TKE budgets on the pressure side."""
    ds = load_wing(NACA4412)["bot4n"]
    assert np.isnan(ds["Pk"].values).all()
    assert np.isfinite(ds["U"].values).all()


# -- the velocity-scale finding (DATA_NOTES section 1) --------------------

@pytest.mark.parametrize("var", list(SUCTION.values()) + ["bot4n"])
def test_U_ref_is_the_shipped_velocity_scale(var):
    """Cf, Re_delta* and Re_theta are all built on max(U), not on Ue."""
    ds = load_wing(NACA4412)[var]
    nu = ds.attrs["nu"]
    U_ref = ds["U_ref"].values
    ok = ds["Cf"].values > 0

    # W1 separates at its last station, where Cf <= 0 and the inversion is
    # undefined; `ok` masks it out, so suppress the warning from the masked
    # elements rather than computing sqrt of a negative.
    with np.errstate(invalid="ignore", divide="ignore"):
        from_Cf = ds["u_tau"].values * np.sqrt(2.0 / ds["Cf"].values)
    from_Reds = ds["Re_delta_star"].values * nu / ds["delta_star"].values
    from_Reth = ds["Re_theta"].values * nu / ds["theta"].values

    for implied in (from_Cf, from_Reds, from_Reth):
        assert np.allclose(implied[ok], U_ref[ok], rtol=1e-6)


def test_Ue_field_would_be_wrong(cases):
    """Guard the finding itself: Ue is materially different from U_ref.

    If a future release of the database fixed this, the test fails loudly
    rather than letting the loader keep applying an obsolete correction.
    """
    ds = cases["W4"]
    ratio = (ds["Ue"] / ds["U_ref"]).values
    assert ratio.max() < 0.995
    assert ratio.min() > 0.85
    assert 0.96 < np.median(ratio) < 0.99


def test_thicknesses_recomputed_with_U_ref(cases):
    """delta* and theta follow from U_ref, integrated to delta99."""
    ds = cases["W4"]
    err_ds, err_th = [], []
    for i in range(ds.sizes["station"]):
        s = ds.isel(station=i)
        yn, U = s.yn.values, s.U.values
        Uref, d99 = float(s.U_ref), float(s.delta99)
        err_ds.append(abs(blq.displacement_thickness(yn, U, Uref, d99)
                          / float(s.delta_star) - 1))
        err_th.append(abs(blq.momentum_thickness(yn, U, Uref, d99)
                          / float(s.theta) - 1))
    assert np.median(err_ds) < 0.01
    assert np.median(err_th) < 0.01
    assert max(err_ds) < 0.03
    assert max(err_th) < 0.05


def test_wall_shear_recovers_u_tau(cases):
    """u_tau from the near-wall slope matches the shipped value."""
    ds = cases["W4"]
    nu = ds.attrs["nu"]
    errs = []
    for i in range(ds.sizes["station"]):
        s = ds.isel(station=i)
        tau = blq.wall_shear_stress(s.yn.values, s.U.values, nu=nu)
        errs.append(abs(blq.friction_velocity(tau) / float(s.u_tau) - 1))
    assert np.median(errs) < 0.005
    assert max(errs) < 0.02


# -- the paper ------------------------------------------------------------

@pytest.mark.parametrize("x_target", [0.4, 0.7])
@pytest.mark.parametrize("case", list(SUCTION))
def test_paper_tables(cases, x_target, case):
    ds = cases[case]
    x = ds.xa.values
    ref = TABLES[x_target][case]

    got = {
        "beta": np.interp(x_target, x, ds["beta_db"].values),
        "Re_tau": np.interp(x_target, x, ds["Re_tau"].values),
        "Re_theta": np.interp(x_target, x, ds["Re_theta"].values),
        "Cf": np.interp(x_target, x, ds["Cf"].values),
        "H": np.interp(x_target, x, (ds["delta_star"] / ds["theta"]).values),
    }
    for q, want in ref.items():
        # beta is quoted to two significant figures in the paper, so it
        # carries a looser tolerance than the rest.
        tol = 0.08 if q == "beta" else 0.05
        assert abs(got[q] - want) / want < tol, (
            f"{case} {q} at x/c={x_target}: paper {want}, got {got[q]:.4g}"
        )


@pytest.mark.parametrize(
    "case,want", [("W1", 139), ("W2", 231), ("W4", 366), ("W10", 707)]
)
def test_max_Re_tau(cases, case, want):
    got = cases[case]["Re_tau"].values.max()
    assert abs(got - want) / want < 0.02


@pytest.mark.parametrize(
    "case,want", [("W1", 1050), ("W2", 1730), ("W4", 2930), ("W10", 6000)]
)
def test_Re_theta_at_analysis_cutoff(cases, case, want):
    """The paper's quoted maxima correspond to x/c = 0.95, its analysis
    cutoff -- not to the last stored station."""
    ds = cases[case]
    got = np.interp(0.95, ds.xa.values, ds["Re_theta"].values)
    assert abs(got - want) / want < 0.05


def test_residual_beta_variation_is_substantial(cases):
    """The matched-beta premise holds only to tens of percent.

    The paper describes the Re_c >= 200k cases as subjected to "approximately
    the same" beta(x), with differences "on the order of 10%" only beyond
    x/c > 0.9. In the distributed data beta decreases systematically with
    Reynolds number across the whole analysis window: W2 runs ~37% above W10
    on average. The paper's own Table 4 agrees (beta = 2.6 vs 2.0 at
    x/c = 0.7, a ratio of 1.30); the apparent collapse in its Fig. 5 is a
    consequence of the four-decade log axis.

    This is the confound the project exists to quantify, so it is pinned
    here: if a data revision removes it, the motivation changes.
    """
    x = np.linspace(0.4, 0.9, 60)
    beta = {c: np.interp(x, cases[c].xa.values, cases[c]["beta_db"].values)
            for c in ("W2", "W4", "W10")}

    # Monotone in Reynolds number at every station.
    assert np.all(beta["W2"] > beta["W4"])
    assert np.all(beta["W4"] > beta["W10"])

    ratio = beta["W2"] / beta["W10"]
    assert 1.3 < ratio.mean() < 1.45
    assert ratio.min() > 1.15          # never close to matched
    assert ratio.max() < 1.6


def test_W1_last_station_is_separated(cases):
    """W1 separates at the final station; it must not reach a model."""
    assert cases["W1"]["beta_db"].values[-1] < -100


# -- integration helper ---------------------------------------------------

def test_integrate_to_is_continuous_in_the_upper_limit():
    """The partial final panel matters: snapping to grid points would make
    delta* jump as delta99 crosses a node."""
    yn = np.linspace(0, 1, 11)
    U = yn.copy()
    # Straddle the node at 0.55. The integral must change by about
    # f(0.55) * dL, not by a whole panel.
    dL = 2e-4
    a = blq._integrate_to(yn, U, 0.55 - dL / 2)
    b = blq._integrate_to(yn, U, 0.55 + dL / 2)
    assert abs(b - a) == pytest.approx(0.55 * dL, rel=0.02)
    # Exact for a linear integrand.
    assert np.isclose(blq._integrate_to(yn, U, 1.0), 0.5, rtol=1e-12)
    assert np.isclose(blq._integrate_to(yn, U, 0.5), 0.125, rtol=1e-12)
