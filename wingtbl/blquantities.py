"""Boundary-layer integral quantities for pressure-gradient TBLs.

Everything here works on a single wall-normal profile expressed in the local
tangential / wall-normal frame, which is how the WING database stores its
statistics.

On the boundary-layer edge
--------------------------
Determining the edge of a pressure-gradient TBL is genuinely ambiguous: the
tangential velocity is not constant for yn > delta99, so "where the velocity
stops changing" is not well defined. The WING database ships ``Ue`` and
``delta99`` computed with the diagnostic-scaling method of

    R. Vinuesa, A. Bobke, R. Orlu and P. Schlatter, "On determining
    characteristic length scales in pressure-gradient turbulent boundary
    layers", Phys. Fluids 28, 055101 (2016).

Those shipped values are canonical and are what :func:`bl_parameters` uses by
default. The alternative estimators in this module exist for the sensitivity
study, *not* to replace the published values -- they are independent and
deliberately simple, so that disagreement between them is informative.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Callable

import numpy as np

__all__ = [
    "BLParameters",
    "cumulative_trapezoid",
    "displacement_thickness",
    "momentum_thickness",
    "wall_shear_stress",
    "friction_velocity",
    "interp_at_edge",
    "edge_by_velocity_fraction",
    "edge_by_turbulence_intensity",
    "clauser_beta",
    "pressure_velocity",
    "bl_parameters",
]


# --------------------------------------------------------------------------
# integration helpers
# --------------------------------------------------------------------------

def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral of ``y`` over ``x``, starting at 0."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def _integrate_to(yn: np.ndarray, f: np.ndarray, upper: float) -> float:
    """Trapezoidal integral of ``f`` from yn[0] to ``upper``, interpolating
    the final partial panel rather than snapping to the nearest grid point.

    Snapping would make the result jump as ``upper`` crosses a grid point,
    which matters because delta99 lands between points and we compare
    integrals across cases with different grids.
    """
    yn = np.asarray(yn, dtype=float)
    f = np.asarray(f, dtype=float)
    if upper <= yn[0]:
        return 0.0
    if upper >= yn[-1]:
        return float(np.trapezoid(f, yn))

    k = int(np.searchsorted(yn, upper))  # yn[k-1] < upper <= yn[k]
    total = float(np.trapezoid(f[:k], yn[:k])) if k > 1 else 0.0
    # partial panel from yn[k-1] to upper
    t = (upper - yn[k - 1]) / (yn[k] - yn[k - 1])
    f_up = f[k - 1] + t * (f[k] - f[k - 1])
    total += 0.5 * (f[k - 1] + f_up) * (upper - yn[k - 1])
    return total


# --------------------------------------------------------------------------
# integral thicknesses
# --------------------------------------------------------------------------

def displacement_thickness(yn: np.ndarray, U: np.ndarray, Ue: float,
                           delta99: float) -> float:
    """delta* = int_0^delta99 (1 - U/Ue) dyn.

    The upper limit is delta99, not infinity. In a PG TBL the integrand does
    not vanish above the edge (U keeps changing), so integrating to the top of
    the domain would not converge to anything meaningful. This matches the
    paper, which states that the momentum thickness is "obtained through
    integration up to yn = delta99".
    """
    return _integrate_to(yn, 1.0 - np.asarray(U) / Ue, delta99)


def momentum_thickness(yn: np.ndarray, U: np.ndarray, Ue: float,
                       delta99: float) -> float:
    """theta = int_0^delta99 (U/Ue)(1 - U/Ue) dyn."""
    u = np.asarray(U) / Ue
    return _integrate_to(yn, u * (1.0 - u), delta99)


# --------------------------------------------------------------------------
# wall quantities
# --------------------------------------------------------------------------

def wall_shear_stress(yn: np.ndarray, U: np.ndarray, nu: float,
                      rho: float = 1.0, n_points: int = 4) -> float:
    """tau_w = rho * nu * dU/dyn at the wall.

    Fits a polynomial with no constant term (U(0) = 0 by no-slip) through the
    first ``n_points`` off-wall points and takes its slope at the wall. A
    plain two-point difference is noticeably biased when the first point sits
    at yn+ ~ 1, which it does here.
    """
    yn = np.asarray(yn, dtype=float)
    U = np.asarray(U, dtype=float)

    # Drop a yn == 0 point if present; no-slip is imposed through the fit.
    mask = yn > 0
    y = yn[mask][:n_points]
    u = U[mask][:n_points]
    if y.size < 2:
        raise ValueError("need at least two off-wall points for dU/dy")

    # U = a1*y + a2*y^2 + ... , no constant term.
    powers = np.arange(1, min(3, y.size) + 1)
    A = y[:, None] ** powers[None, :]
    coef, *_ = np.linalg.lstsq(A, u, rcond=None)
    dUdy_wall = coef[0]
    return float(rho * nu * dUdy_wall)


def friction_velocity(tau_w: float, rho: float = 1.0) -> float:
    return float(np.sqrt(abs(tau_w) / rho)) * (1.0 if tau_w >= 0 else -1.0)


# --------------------------------------------------------------------------
# edge estimators (for the sensitivity study -- see module docstring)
# --------------------------------------------------------------------------

def interp_at_edge(yn: np.ndarray, f: np.ndarray, delta: float) -> float:
    """Linearly interpolate ``f`` at wall distance ``delta``."""
    return float(np.interp(delta, np.asarray(yn), np.asarray(f)))


def edge_by_velocity_fraction(yn: np.ndarray, U: np.ndarray, Ue: float,
                              fraction: float = 0.99) -> float:
    """First wall distance where U reaches ``fraction`` * Ue.

    Takes the *first* upward crossing, scanning outward. In an APG TBL the
    profile can wobble around the threshold near the edge, and taking the last
    crossing would place the edge far out in the freestream.
    """
    yn = np.asarray(yn, dtype=float)
    U = np.asarray(U, dtype=float)
    target = fraction * Ue
    above = np.flatnonzero(U >= target)
    if above.size == 0:
        return float("nan")
    k = int(above[0])
    if k == 0:
        return float(yn[0])
    u0, u1 = U[k - 1], U[k]
    if u1 == u0:
        return float(yn[k])
    t = (target - u0) / (u1 - u0)
    return float(yn[k - 1] + t * (yn[k] - yn[k - 1]))


def edge_by_turbulence_intensity(yn: np.ndarray, u_rms: np.ndarray,
                                 Ue: float, threshold: float = 0.01) -> float:
    """Outermost wall distance where u_rms/Ue still exceeds ``threshold``.

    A crude but fully independent estimator: it asks where the turbulence
    dies out rather than where the mean velocity flattens, so it does not
    inherit the mean-profile ambiguity that makes the edge hard to define
    under a pressure gradient. Used only as a cross-check.
    """
    yn = np.asarray(yn, dtype=float)
    ti = np.asarray(u_rms, dtype=float) / Ue
    above = np.flatnonzero(ti > threshold)
    if above.size == 0:
        return float("nan")
    k = int(above[-1])
    if k >= yn.size - 1:
        return float(yn[-1])
    t0, t1 = ti[k], ti[k + 1]
    if t1 == t0:
        return float(yn[k])
    t = (threshold - t0) / (t1 - t0)
    return float(yn[k] + t * (yn[k + 1] - yn[k]))


# --------------------------------------------------------------------------
# pressure-gradient parameters
# --------------------------------------------------------------------------

def pressure_velocity(delta_star: float, dPe_dxt: float,
                      rho: float = 1.0) -> float:
    """Mellor & Gibson pressure velocity, up = sqrt((delta*/rho) dPe/dxt).

    Returns NaN under a favourable pressure gradient (dPe/dxt < 0), where the
    square root is imaginary -- up is only defined for decelerating flow.
    """
    arg = delta_star * dPe_dxt / rho
    return float(np.sqrt(arg)) if arg >= 0 else float("nan")


def clauser_beta(delta_star: float, tau_w: float, dPe_dxt: float) -> float:
    """beta = (delta* / tau_w) * dPe/dxt."""
    return float(delta_star * dPe_dxt / tau_w)


def dPe_dxt_from_edge_velocity(Ue: np.ndarray, xt: np.ndarray,
                               rho: float = 1.0) -> np.ndarray:
    """dPe/dxt = -rho * Ue * dUe/dxt, from Bernoulli along the edge streamline.

    ``xt`` must be arc length along the surface, not x/c, or the gradient is
    wrong by the local surface slope.
    """
    Ue = np.asarray(Ue, dtype=float)
    xt = np.asarray(xt, dtype=float)
    return -rho * Ue * np.gradient(Ue, xt)


# --------------------------------------------------------------------------
# the lot, for one profile
# --------------------------------------------------------------------------

@dataclass
class BLParameters:
    """Integral parameters of one wall-normal profile."""

    delta99: float
    delta_star: float
    theta: float
    H: float
    u_tau: float
    tau_w: float
    Ue: float
    Cf: float
    Re_tau: float
    Re_theta: float
    Re_delta_star: float
    beta: float
    up_over_Ue: float
    Ue_plus: float

    def as_dict(self) -> dict:
        return asdict(self)


def bl_parameters(yn: np.ndarray, U: np.ndarray, *, nu: float, Ue: float,
                  delta99: float, dPe_dxt: float, rho: float = 1.0,
                  tau_w: float | None = None) -> BLParameters:
    """Assemble the integral parameters of a single profile.

    Parameters
    ----------
    yn, U
        Wall-normal coordinate (starting at the wall) and tangential mean
        velocity, in the local frame.
    nu
        Kinematic viscosity, in the same units as the data.
    Ue, delta99
        Edge velocity and boundary-layer thickness. Pass the values shipped
        with the database; see the module docstring.
    dPe_dxt
        Streamwise pressure gradient at the edge, along the surface.
    tau_w
        Wall-shear stress. Computed from the profile if not given.
    """
    if tau_w is None:
        tau_w = wall_shear_stress(yn, U, nu=nu, rho=rho)
    u_tau = friction_velocity(tau_w, rho=rho)

    delta_star = displacement_thickness(yn, U, Ue, delta99)
    theta = momentum_thickness(yn, U, Ue, delta99)
    up = pressure_velocity(delta_star, dPe_dxt, rho=rho)

    return BLParameters(
        delta99=float(delta99),
        delta_star=delta_star,
        theta=theta,
        H=delta_star / theta if theta else float("nan"),
        u_tau=u_tau,
        tau_w=float(tau_w),
        Ue=float(Ue),
        Cf=2.0 * (u_tau / Ue) ** 2,
        Re_tau=delta99 * u_tau / nu,
        Re_theta=Ue * theta / nu,
        Re_delta_star=Ue * delta_star / nu,
        beta=clauser_beta(delta_star, tau_w, dPe_dxt),
        up_over_Ue=up / Ue,
        Ue_plus=Ue / u_tau,
    )
