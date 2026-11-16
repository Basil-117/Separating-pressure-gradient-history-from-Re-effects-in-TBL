"""Tests for the fitted-memory power law.

The test that matters is :func:`test_recovers_a_known_law`: synthetic data is
generated from a power law with a chosen memory length, and the fitter has to
find both the exponents and the length. Everything else in Phase 2 -- the
held-out accuracy, the Reynolds-vs-pressure-gradient split -- is only worth
reading if that works.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from wingtbl.features import memory_average
from wingtbl.memory_model import (
    FEATURES, MemoryPowerLaw, kernel_average, trapezoid_matrix,
)


# -- quadrature and kernel ------------------------------------------------

def test_trapezoid_matrix_integrates_exactly():
    xt = np.linspace(0, 2, 41)
    T = trapezoid_matrix(xt)
    f = 3 * xt + 1
    # Row i integrates from xt[0] to xt[i]; exact for a linear integrand.
    for i in (1, 10, 40):
        want = 1.5 * xt[i] ** 2 + xt[i]
        assert T[i] @ f == pytest.approx(want, rel=1e-12)


def test_trapezoid_matrix_is_lower_triangular():
    T = trapezoid_matrix(np.linspace(0, 1, 12))
    assert np.allclose(np.triu(T, k=1), 0.0)
    assert T[0].sum() == 0.0          # nothing upstream of the first station


def test_kernel_matches_the_reference_implementation():
    """The model's vectorised kernel and the loop in features.py must agree."""
    xt = np.linspace(0, 1, 60)
    beta = np.sin(3 * xt) + 1.5
    d99 = np.full(60, 0.02)
    for L in (0.5, 5.0, 20.0):
        a = memory_average(xt, beta, L * d99)
        b = kernel_average(xt, beta, d99, L)
        assert np.allclose(a, b, atol=1e-12)


def test_constant_beta_is_a_fixed_point():
    xt = np.linspace(0, 1, 30)
    beta = np.full(30, 1.7)
    for L in (0.05, 1.0, 25.0):
        assert np.allclose(kernel_average(xt, beta, np.full(30, 0.03), L), 1.7)


def test_short_memory_collapses_onto_local_beta():
    """At short memory the history feature must switch itself off, which is
    what makes the ratio parameterisation well behaved."""
    xt = np.linspace(0, 1, 120)
    beta = 0.5 + 2 * xt
    d99 = np.full(120, 0.01)
    bb = kernel_average(xt, beta, d99, 0.02)
    assert np.allclose(bb[2:], beta[2:], rtol=2e-3)


# -- the fitter -----------------------------------------------------------

def _synthetic(ell_true: float, coef_true: tuple[float, float, float],
               n_case: int = 3, n: int = 50, seed: int = 0):
    """Three cases with wing-like beta(x) and Re_tau, and a target built
    from a known power law."""
    rng = np.random.default_rng(seed)
    frames, full = [], {}
    for c in range(n_case):
        case = f"S{c}"
        xt = np.linspace(0.0, 0.9, n)
        # beta rises steeply towards the trailing edge, as on the wing, and
        # its magnitude differs between cases, as it does between Reynolds
        # numbers in the real data.
        beta = (0.1 + 3.0 * xt ** 2) * (1.0 + 0.3 * c)
        d99 = 0.004 * (1 + 3 * xt) * (1.0 - 0.15 * c)
        Re_tau = 150.0 * (2.0 ** c) * (1 + 2 * xt)

        bb = kernel_average(xt, beta, d99, ell_true)
        logy = (0.4 + coef_true[0] * np.log(beta)
                + coef_true[1] * np.log(bb / beta)
                + coef_true[2] * np.log(Re_tau)
                + 0.002 * rng.standard_normal(n))

        d = pd.DataFrame(dict(case=case, station=np.arange(n), xa=xt, xt=xt,
                              beta=beta, delta99=d99, Re_tau=Re_tau,
                              y=np.exp(logy)))
        full[case] = d
        frames.append(d[(d.xa >= 0.2) & (d.xa <= 0.85)])
    return pd.concat(frames, ignore_index=True), full


def _noiseless(ell_true: float, coef_true: tuple[float, float, float]):
    df, full = _synthetic(ell_true, coef_true)
    for d in full.values():
        bb = kernel_average(d.xt.to_numpy(), d.beta.to_numpy(),
                            d.delta99.to_numpy(), ell_true)
        d["y"] = np.exp(0.4
                        + coef_true[0] * np.log(d.beta.to_numpy())
                        + coef_true[1] * np.log(bb / d.beta.to_numpy())
                        + coef_true[2] * np.log(d.Re_tau.to_numpy()))
    df = pd.concat([d[(d.xa >= 0.2) & (d.xa <= 0.85)] for d in full.values()],
                   ignore_index=True)
    return df, full


@pytest.mark.parametrize("ell_true", [0.5, 3.0, 12.0])
def test_recovers_a_known_law_exactly_without_noise(ell_true):
    """With clean data the fitter must be essentially exact. This separates
    'the optimiser works' from 'the parameters are identifiable', which are
    different questions and have different answers here."""
    coef_true = (0.6, -2.5, -0.3)
    df, full = _noiseless(ell_true, coef_true)
    r = MemoryPowerLaw(alpha=1e-10).fit(df, "y", full)

    assert r.ell == pytest.approx(ell_true, rel=2e-3)
    assert not r.ell_at_bound
    for got, want in zip(r.coef, coef_true):
        assert got == pytest.approx(want, abs=2e-3)


@pytest.mark.parametrize("ell_true", [0.5, 3.0, 12.0])
def test_memory_length_is_only_weakly_identified_under_noise(ell_true):
    """Documents a real limitation, so it cannot be quietly lost.

    0.2% noise in log space -- far below the scatter in the LES statistics --
    moves the fitted memory length by tens of percent, because a shorter
    memory shrinks log(beta_bar/beta) and a larger history exponent
    compensates almost exactly. The length is good to a factor of order
    1.3, not to two digits, and the project must not quote it as if it were.
    """
    df, full = _synthetic(ell_true, (0.6, -2.5, -0.3))
    r = MemoryPowerLaw(alpha=1e-8).fit(df, "y", full)
    assert r.ell == pytest.approx(ell_true, rel=0.45)
    # The terms that do not trade off are still recovered well.
    assert r.coef[0] == pytest.approx(0.6, abs=0.05)
    assert r.coef[2] == pytest.approx(-0.3, abs=0.05)


@pytest.mark.parametrize("ell_true", [0.5, 3.0, 12.0])
def test_history_contribution_survives_the_tradeoff(ell_true):
    """The quantity the decomposition actually uses is the *product* of the
    history exponent and the history feature, not either alone. The
    exponent-length trade-off cancels in that product, which is why the
    Reynolds-vs-pressure-gradient split is trustworthy even though the
    memory length on its own is not."""
    coef_true = (0.6, -2.5, -0.3)
    df, full = _synthetic(ell_true, coef_true)

    df_c, full_c = _noiseless(ell_true, coef_true)
    truth = MemoryPowerLaw(alpha=1e-10)
    truth.fit(df_c, "y", full_c)

    fitted = MemoryPowerLaw(alpha=1e-8)
    fitted.fit(df, "y", full)

    k = FEATURES.index("mem_ratio")
    a = truth.attribute(df_c, full_c)[f"phi_{FEATURES[k]}"].to_numpy()
    b = fitted.attribute(df, full)[f"phi_{FEATURES[k]}"].to_numpy()
    # Same stations in both frames, so compare pointwise.
    assert np.max(np.abs(a - b)) < 0.06
    assert np.corrcoef(a, b)[0, 1] > 0.99


def test_recovered_model_predicts_a_held_out_case():
    df, full = _synthetic(3.0, (0.6, -2.5, -0.3))
    m = MemoryPowerLaw(alpha=1e-8)
    tr = df[df.case != "S2"].reset_index(drop=True)
    te = df[df.case == "S2"].reset_index(drop=True)
    m.fit(tr, "y", full)
    pred = m.predict(te, full)
    rel = np.abs(pred / te.y.to_numpy() - 1)
    assert np.median(rel) < 0.02


def test_no_history_dependence_is_detected_as_such():
    """If the target does not depend on history, the fitted history exponent
    must be near zero -- otherwise the decomposition would invent an effect."""
    df, full = _synthetic(3.0, (0.6, 0.0, -0.3))
    r = MemoryPowerLaw(alpha=1e-8).fit(df, "y", full)
    assert abs(r.coef[1]) < 0.1


def test_attribution_is_exact():
    """Shapley values plus baseline must reconstruct the prediction exactly;
    for a linear model there is no sampling error to hide behind."""
    df, full = _synthetic(3.0, (0.6, -2.5, -0.3))
    m = MemoryPowerLaw()
    m.fit(df, "y", full)
    phi = m.attribute(df, full)
    total = phi[[f"phi_{f}" for f in FEATURES]].sum(axis=1) + phi["baseline"]
    assert np.allclose(np.exp(total.to_numpy()), m.predict(df, full), rtol=1e-10)


def test_attribution_is_centred():
    """Each feature's contribution averages to zero over the fitting set, so
    shares are read against a common baseline rather than an arbitrary one."""
    df, full = _synthetic(3.0, (0.6, -2.5, -0.3))
    m = MemoryPowerLaw()
    m.fit(df, "y", full)
    phi = m.attribute(df, full)
    for f in FEATURES:
        assert abs(phi[f"phi_{f}"].mean()) < 1e-10


def test_fit_before_predict_is_an_error():
    df, full = _synthetic(1.0, (0.5, -1.0, -0.2))
    with pytest.raises(RuntimeError):
        MemoryPowerLaw().predict(df, full)
