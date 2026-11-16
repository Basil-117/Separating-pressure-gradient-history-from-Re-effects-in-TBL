"""Features and targets for the history-aware model.

The question the features are built to answer: how much of a boundary layer's
departure from canonical behaviour is set by the *local* pressure gradient,
how much by the pressure gradient it has already been through, and how much
by Reynolds number?

That requires a history variable that is more than "beta somewhere upstream".
:func:`memory_average` convolves beta(x) with a decaying kernel whose length
is a free parameter, so the model can learn how far back the boundary layer
remembers instead of being told. Clauser's observation that the outer layer
has a long memory then becomes a number, and a number that can differ between
inner-layer and outer-layer quantities.

All streamwise integration is along surface arc length ``xt``, never ``x/c``:
on the forward part of the wing the surface slope is far from 1, and using
``x/c`` would silently rescale the pressure gradient.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Sequence

import numpy as np
import pandas as pd
import xarray as xr

__all__ = [
    "AnalysisWindow",
    "DEFAULT_WINDOW",
    "memory_average",
    "accumulated_beta",
    "outer_value",
    "build_table",
    "LOCAL_FEATURES",
    "HISTORY_FEATURES",
    "TARGETS",
]


@dataclass(frozen=True)
class AnalysisWindow:
    """Chordwise range the model is fitted on.

    Lower bound: past the tripping at x/c = 0.1 and past the stations where
    the suction side is still accelerating. Upper bound: beyond this the
    structures detach from the wall and the boundary-layer description stops
    meaning much, so the paper restricts its own analysis the same way.
    """

    lo: float = 0.25
    hi: float = 0.95


DEFAULT_WINDOW = AnalysisWindow()


# --------------------------------------------------------------------------
# history
# --------------------------------------------------------------------------

def memory_average(
    xt: np.ndarray,
    beta: np.ndarray,
    length: np.ndarray | float,
    kernel: Literal["exp", "box"] = "exp",
) -> np.ndarray:
    """Kernel-weighted average of the upstream pressure gradient.

    .. math::
        \\bar\\beta(x;\\ell) =
        \\frac{\\int_{x_0}^{x} \\beta(\\xi)\\,K\\!\\left(\\frac{x-\\xi}{\\ell}\\right)d\\xi}
             {\\int_{x_0}^{x} K\\!\\left(\\frac{x-\\xi}{\\ell}\\right)d\\xi}

    Normalising by the kernel's own integral over the *available* upstream
    distance matters: without it, stations near the start of the data would
    be assigned a small history simply because less of the wing lies upstream
    of them, and the model would learn chordwise position rather than memory.
    The normalisation makes ``beta_bar`` a weighted mean of beta, so it is
    directly comparable with the local beta at the same station.

    Parameters
    ----------
    xt
        Surface arc length, increasing.
    beta
        Clauser parameter at each station.
    length
        Memory length in the same units as ``xt``. A scalar gives a fixed
        length; an array (one per station) lets the length scale with a local
        quantity such as ``delta99``, which is the physically natural choice.
    kernel
        ``exp`` for exponential decay, ``box`` for a plain upstream mean over
        one memory length (kept as the simple comparator).

    Returns
    -------
    Weighted upstream mean of beta at each station. The first station has no
    upstream data, so it is given its own local value.
    """
    xt = np.asarray(xt, dtype=float)
    beta = np.asarray(beta, dtype=float)
    n = xt.size
    ell = np.full(n, float(length)) if np.isscalar(length) else np.asarray(length, float)
    if np.any(ell <= 0):
        raise ValueError("memory length must be positive")

    out = np.empty(n)
    for i in range(n):
        d = xt[i] - xt[: i + 1]          # >= 0, distance upstream
        if kernel == "exp":
            w = np.exp(-d / ell[i])
        elif kernel == "box":
            w = (d <= ell[i]).astype(float)
        else:
            raise ValueError(f"unknown kernel {kernel!r}")
        if i == 0:
            out[i] = beta[0]
            continue
        num = np.trapezoid(w * beta[: i + 1], xt[: i + 1])
        den = np.trapezoid(w, xt[: i + 1])
        # den -> 0 when the memory is far shorter than the first grid spacing;
        # the weighted mean is then just the local value.
        out[i] = num / den if den > 1e-14 else beta[i]
    return out


def accumulated_beta(xt: np.ndarray, beta: np.ndarray,
                     delta99: np.ndarray) -> np.ndarray:
    """Pressure gradient accumulated upstream, in units of local thickness.

    A memory-length-free history measure in the spirit of the integral
    history parameter used in recent experimental work: the total upstream
    beta, scaled by how thick the layer is now. Included as the baseline that
    the learned kernel has to beat.
    """
    xt = np.asarray(xt, float)
    beta = np.asarray(beta, float)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (beta[1:] + beta[:-1]) * np.diff(xt))])
    return cum / np.asarray(delta99, float)


# --------------------------------------------------------------------------
# outer-region sampling
# --------------------------------------------------------------------------

def outer_value(yn: np.ndarray, f: np.ndarray, delta99: float,
                frac: float = 0.2) -> float:
    """Interpolate a profile at ``frac`` * delta99.

    0.2 is the paper's sampling height for the wing-to-ZPG ratios: the end of
    the overlap region for Re_theta < 6000.
    """
    return float(np.interp(frac * delta99, np.asarray(yn), np.asarray(f)))


# --------------------------------------------------------------------------
# table assembly
# --------------------------------------------------------------------------

LOCAL_FEATURES = ["beta", "log_Re_tau", "up_over_Ue", "dbeta_dx_scaled"]
HISTORY_FEATURES = ["beta_bar", "beta_bar_over_beta", "beta_acc"]

TARGETS = [
    "H",            # shape factor
    "Cf",           # skin friction
    "Ue_plus",      # inner-scaled edge velocity
    "uu_outer",     # tangential fluctuation at 0.2 delta99, inner-scaled
    "uv_outer",     # Reynolds shear stress there, inner-scaled
    "Pk_outer",     # TKE production there, inner-scaled
]


def _case_frame(ds: xr.Dataset, memory_in_delta99: float,
                window: AnalysisWindow) -> pd.DataFrame:
    n = ds.sizes["station"]
    xa = ds.xa.values
    xt = ds.xt.values
    nu = ds.attrs["nu"]

    beta = ds.beta_db.values
    d99 = ds.delta99.values
    u_tau = ds.u_tau.values
    U_ref = ds.U_ref.values
    dstar = ds.delta_star.values

    # History is integrated over the whole available chord, then the window
    # is applied. Trimming first would throw away the upstream record that
    # the history variable exists to capture.
    ell = memory_in_delta99 * d99
    beta_bar = memory_average(xt, beta, ell)
    beta_acc = accumulated_beta(xt, beta, d99)

    # Local rate of change, made dimensionless with the local thickness.
    dbeta_dx_scaled = np.gradient(beta, xt) * d99

    # up/Ue: outer-scaled pressure-gradient magnitude. beta = (up/u_tau)^2,
    # so up = u_tau sqrt(beta), defined only where beta > 0.
    with np.errstate(invalid="ignore"):
        up_over_Ue = np.where(beta > 0, u_tau * np.sqrt(np.abs(beta)) / U_ref, np.nan)

    rows = []
    for i in range(n):
        yn = ds.yn.values[i]
        ut2 = u_tau[i] ** 2
        rec = {
            "case": ds.attrs["case"],
            "airfoil": ds.attrs["airfoil"],
            "side": ds.attrs["side"],
            "Re_c": ds.attrs["Re_c"],
            "station": i,
            "xa": xa[i],
            "xt": xt[i],
            # features
            "beta": beta[i],
            "Re_tau": ds.Re_tau.values[i],
            "log_Re_tau": np.log(ds.Re_tau.values[i]),
            "Re_theta": ds.Re_theta.values[i],
            "up_over_Ue": up_over_Ue[i],
            "dbeta_dx_scaled": dbeta_dx_scaled[i],
            "beta_bar": beta_bar[i],
            "beta_bar_over_beta": beta_bar[i] / beta[i] if beta[i] > 0 else np.nan,
            "beta_acc": beta_acc[i],
            # targets
            "H": dstar[i] / ds.theta.values[i],
            "Cf": ds.Cf.values[i],
            "Ue_plus": U_ref[i] / u_tau[i],
            "uu_outer": outer_value(yn, ds.uu.values[i], d99[i]) / ut2,
            "uv_outer": outer_value(yn, ds.uv.values[i], d99[i]) / ut2,
            "Pk_outer": outer_value(yn, ds.Pk.values[i], d99[i]) * nu / u_tau[i] ** 4,
            # context
            "delta99": d99[i],
            "u_tau": u_tau[i],
        }
        rows.append(rec)

    df = pd.DataFrame(rows)
    return df[(df.xa >= window.lo) & (df.xa <= window.hi)].reset_index(drop=True)


def build_table(
    datasets: Iterable[xr.Dataset],
    memory_in_delta99: float = 20.0,
    window: AnalysisWindow = DEFAULT_WINDOW,
) -> pd.DataFrame:
    """Assemble the modelling table across cases.

    Parameters
    ----------
    datasets
        Per-case datasets from :func:`wingtbl.loader.load_wing`.
    memory_in_delta99
        Memory length as a multiple of the local ``delta99``. The default is
        a starting guess only -- Phase 2 fits it, and the point of the project
        is that its best value differs between inner and outer quantities.
    window
        Chordwise range to keep.
    """
    frames = [_case_frame(ds, memory_in_delta99, window) for ds in datasets]
    return pd.concat(frames, ignore_index=True)
