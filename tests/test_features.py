"""Tests for the history features and the validation splitters.

The memory kernel is the one genuinely new piece of machinery in the
project, so its limiting behaviour is pinned here: if it stops reducing to
the local value at short memory, or to the upstream mean at long memory, the
"memory length" it reports means nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from wingtbl.features import accumulated_beta, memory_average, outer_value
from wingtbl.validation import (
    blocked_streamwise, leave_one_case_out, relative_error, summarise_errors,
)


# -- memory kernel --------------------------------------------------------

@pytest.fixture
def ramp():
    xt = np.linspace(0.0, 1.0, 101)
    return xt, xt.copy()          # beta rises linearly from 0 to 1


def test_short_memory_recovers_the_local_value(ramp):
    xt, beta = ramp
    out = memory_average(xt, beta, length=1e-4)
    assert np.allclose(out[1:], beta[1:], atol=5e-3)


def test_long_memory_approaches_the_upstream_mean(ramp):
    xt, beta = ramp
    out = memory_average(xt, beta, length=1e6)
    # With a kernel far wider than the domain, every upstream point is
    # weighted equally, so the result is the running mean of beta.
    running = np.concatenate([[beta[0]], [beta[: i + 1].mean()
                                          for i in range(1, beta.size)]])
    assert np.allclose(out[1:], running[1:], atol=1e-3)


def test_memory_average_is_bounded_by_the_signal(ramp):
    """A weighted mean can never leave the range of what it averages."""
    xt, beta = ramp
    for L in (0.01, 0.1, 1.0, 10.0):
        out = memory_average(xt, beta, length=L)
        for i in range(beta.size):
            assert beta[: i + 1].min() - 1e-9 <= out[i] <= beta[: i + 1].max() + 1e-9


def test_constant_beta_is_a_fixed_point(ramp):
    """Averaging a constant must return that constant at every length."""
    xt, _ = ramp
    beta = np.full_like(xt, 2.5)
    for L in (0.001, 0.5, 100.0):
        assert np.allclose(memory_average(xt, beta, L), 2.5)


def test_memory_lags_a_rising_signal(ramp):
    """The point of the feature: with a real memory the weighted average sits
    below a rising beta, and further below as the memory lengthens."""
    xt, beta = ramp
    short = memory_average(xt, beta, 0.02)
    long_ = memory_average(xt, beta, 0.5)
    assert np.all(short[5:] <= beta[5:] + 1e-12)
    assert np.all(long_[5:] < short[5:])


def test_per_station_memory_length_is_accepted(ramp):
    xt, beta = ramp
    ell = np.linspace(0.01, 0.5, xt.size)
    out = memory_average(xt, beta, ell)
    assert out.shape == beta.shape
    assert np.isfinite(out).all()


def test_box_kernel_converges_to_a_plain_upstream_mean():
    """The box kernel is a step, so quadrature on a fixed grid mis-weights
    the panel straddling the window edge by about half a panel. That is a
    discretisation effect, not a bug, so the test is that the error shrinks
    with the grid rather than that it is zero.

    The exponential kernel used by the model is smooth and has no such edge,
    which is one reason it is the default.
    """
    L, x_eval = 0.25, 0.75

    def err(n):
        xt = np.linspace(0, 1, n)
        beta = np.sin(4 * xt) + 2
        i = int(np.argmin(np.abs(xt - x_eval)))
        got = memory_average(xt, beta, L, kernel="box")[i]
        m = (xt >= xt[i] - L) & (xt <= xt[i])
        want = np.trapezoid(beta[m], xt[m]) / (xt[m][-1] - xt[m][0])
        return abs(got - want)

    coarse, fine = err(201), err(1601)
    assert coarse < 0.01                     # small to begin with
    assert fine < coarse / 4                 # and first-order convergent


def test_rejects_nonpositive_memory(ramp):
    xt, beta = ramp
    with pytest.raises(ValueError):
        memory_average(xt, beta, 0.0)


def test_accumulated_beta_is_monotone_for_positive_beta(ramp):
    xt, beta = ramp
    acc = accumulated_beta(xt, beta, delta99=np.ones_like(xt))
    assert np.all(np.diff(acc) >= -1e-12)
    assert acc[0] == 0.0
    # integral of x from 0 to 1 is 1/2
    assert acc[-1] == pytest.approx(0.5, rel=1e-6)


def test_outer_value_interpolates():
    yn = np.linspace(0, 1, 11)
    f = 3.0 * yn
    assert outer_value(yn, f, delta99=1.0, frac=0.2) == pytest.approx(0.6, rel=1e-9)


# -- splitters ------------------------------------------------------------

@pytest.fixture
def toy():
    return pd.DataFrame({
        "case": ["A"] * 10 + ["B"] * 10 + ["C"] * 10,
        "xa": list(np.linspace(0.25, 0.95, 10)) * 3,
    })


def test_leave_one_case_out_partitions(toy):
    folds = list(leave_one_case_out(toy))
    assert [f[0] for f in folds] == ["A", "B", "C"]
    for _, tr, te in folds:
        assert set(tr).isdisjoint(te)
        assert len(tr) + len(te) == len(toy)
        assert toy.iloc[te].case.nunique() == 1


def test_blocked_streamwise_holds_out_a_contiguous_stretch(toy):
    for _, tr, te in blocked_streamwise(toy, n_blocks=4):
        assert set(tr).isdisjoint(te)
        held = toy.iloc[te].xa
        kept = toy.iloc[tr].xa
        # the held-out x-range must not be interleaved with the training one
        assert not ((kept > held.min()) & (kept < held.max())).any()
        # and it must span every case, unlike leave-one-case-out
        assert toy.iloc[te].case.nunique() == 3


def test_relative_error_handles_negative_targets():
    e = relative_error(np.array([-2.0]), np.array([-1.0]))
    assert e[0] == pytest.approx(0.5)


def test_summarise_errors_perfect_prediction():
    y = np.array([1.0, 2.0, 3.0])
    s = summarise_errors(y, y)
    assert s["median_rel"] == 0.0
    assert s["r2"] == pytest.approx(1.0)
    assert s["n"] == 3
