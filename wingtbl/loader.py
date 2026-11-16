"""Load the VinuesaLab WING database into labelled xarray datasets.

The published files are MATLAB structs whose variable names encode the case:
``top``/``bot`` for the suction/pressure side, then a digit for the Reynolds
number (1 -> 100k, 2 -> 200k, 4 -> 400k, 10 -> 1M), then a trailing ``n``.
Each struct is a 1 x N array with one entry per streamwise station.

This module turns that into one :class:`xarray.Dataset` per case, with
dimensions ``(station, j)``: ``station`` indexes chordwise position and ``j``
indexes the wall-normal grid. ``yn`` has to stay a 2-D coordinate because
each station has its own wall-normal grid.

Source
------
R. Vinuesa, P. S. Negi, M. Atzori, A. Hanifi, D. S. Henningson and
P. Schlatter, "Turbulent boundary layers around wing sections up to
Re_c = 1,000,000", Int. J. Heat Fluid Flow 72, 86-99 (2018).

A. Tanarro, R. Vinuesa and P. Schlatter, "Effect of adverse pressure
gradients on turbulent wing boundary layers", J. Fluid Mech. 883, A8 (2020).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr
from scipy.io import loadmat

__all__ = ["CASES", "parse_case_name", "load_case", "load_wing", "concat_cases"]


# Reynolds-number code used in the MATLAB variable names.
_RE_CODE = {"1": 100_000, "2": 200_000, "4": 400_000, "10": 1_000_000}

# Case labels as used in the paper (Table 2).
CASES = {
    ("top", 100_000): "W1",
    ("top", 200_000): "W2",
    ("top", 400_000): "W4",
    ("top", 1_000_000): "W10",
}

_SIDE = {"top": "suction", "bot": "pressure"}

# Profile fields: (station, j)
_PROFILE_FIELDS = {
    "U": ("tangential mean velocity", "U_inf"),
    "V": ("wall-normal mean velocity", "U_inf"),
    "uu": ("tangential velocity variance <u_t u_t>", "U_inf^2"),
    "vv": ("wall-normal velocity variance <v_n v_n>", "U_inf^2"),
    "ww": ("spanwise velocity variance <w w>", "U_inf^2"),
    "uv": ("Reynolds shear stress <u_t v_n>", "U_inf^2"),
    "P": ("mean pressure", "rho U_inf^2"),
    "prms": ("rms pressure fluctuation", "rho U_inf^2"),
    "Pk": ("TKE production", "U_inf^3 / c"),
    "Dk": ("TKE dissipation", "U_inf^3 / c"),
    "Tk": ("TKE turbulent transport", "U_inf^3 / c"),
    "VDk": ("TKE viscous diffusion", "U_inf^3 / c"),
    "VPGk": ("TKE velocity-pressure-gradient term", "U_inf^3 / c"),
}

# Per-station scalars shipped with the database.
_STATION_FIELDS = {
    "ut": ("friction velocity u_tau", "U_inf"),
    "delta99": ("99% boundary-layer thickness", "c"),
    "deltas": ("displacement thickness delta*", "c"),
    "theta": ("momentum thickness", "c"),
    "Ue": ("boundary-layer edge velocity", "U_inf"),
    "Ret": ("friction Reynolds number delta99 u_tau / nu", "1"),
    "Reds": ("Reynolds number based on delta*", "1"),
    "Reth": ("Reynolds number based on momentum thickness", "1"),
    "Cf": ("skin-friction coefficient 2 (u_tau/Ue)^2", "1"),
    "beta": ("Clauser pressure-gradient parameter", "1"),
}

# Renamed on load so the names say what they are.
_RENAME = {
    "ut": "u_tau",
    "deltas": "delta_star",
    "Ret": "Re_tau",
    "Reds": "Re_delta_star",
    "Reth": "Re_theta",
    "beta": "beta_db",
}


def parse_case_name(name: str) -> tuple[str, int] | None:
    """``'top10n'`` -> ``('top', 1_000_000)``. Returns None if it isn't a case.

    The NACA0012 file names its single case ``top4n12``, where the trailing
    ``12`` is the airfoil rather than part of the Reynolds-number code, so the
    digits after ``n`` are matched and discarded.
    """
    m = re.fullmatch(r"(top|bot)(\d+)n(\d*)", name)
    if not m:
        return None
    side, code = m.group(1), m.group(2)
    if code not in _RE_CODE:
        return None
    return side, _RE_CODE[code]


def _as_profile(entry, field: str, n_j: int) -> np.ndarray:
    """Pull one profile field, padding absent fields with NaN.

    The TKE budget terms are empty on the pressure side of the NACA4412 file,
    so those fields have to survive as all-NaN rather than blowing up.
    """
    a = np.atleast_1d(np.asarray(entry[field], dtype=float).squeeze())
    if a.size == 0:
        return np.full(n_j, np.nan)
    if a.size != n_j:
        out = np.full(n_j, np.nan)
        out[: a.size] = a[:n_j]
        return out
    return a


def load_case(mat: dict, varname: str, *, airfoil: str) -> xr.Dataset:
    """Build the dataset for one MATLAB variable, e.g. ``top4n``."""
    parsed = parse_case_name(varname)
    if parsed is None:
        raise ValueError(f"{varname!r} is not a WING case variable")
    side, Re_c = parsed

    struct = np.atleast_1d(mat[varname])
    n_station = struct.size
    n_j = max(np.atleast_1d(e["yn"]).size for e in struct)

    nu = float(np.asarray(struct[0]["nu"]).squeeze())

    coords = {
        "station": np.arange(n_station),
        "xa": ("station", np.array([float(e["xa"]) for e in struct])),
        "ya": ("station", np.array([float(e["ya"]) for e in struct])),
        "yn": (
            ("station", "j"),
            np.stack([_as_profile(e, "yn", n_j) for e in struct]),
        ),
    }

    data_vars = {}
    for f, (desc, units) in _PROFILE_FIELDS.items():
        arr = np.stack([_as_profile(e, f, n_j) for e in struct])
        data_vars[f] = (("station", "j"), arr, {"long_name": desc, "units": units})

    for f, (desc, units) in _STATION_FIELDS.items():
        arr = np.array([float(np.asarray(e[f]).squeeze()) for e in struct])
        name = _RENAME.get(f, f)
        data_vars[name] = (("station",), arr, {"long_name": desc, "units": units})

    ds = xr.Dataset(data_vars, coords=coords)

    # Global-frame coordinates of the profile points, kept for plotting the
    # profiles on the airfoil.
    for f in ("x", "y"):
        ds[f"{f}_global"] = (
            ("station", "j"),
            np.stack([_as_profile(e, f, n_j) for e in struct]),
        )

    ds["tke"] = 0.5 * (ds["uu"] + ds["vv"] + ds["ww"])
    ds["tke"].attrs = {"long_name": "turbulent kinetic energy", "units": "U_inf^2"}

    # The velocity scale behind every shipped integral quantity is the peak
    # tangential velocity of the profile, NOT the `Ue` field. See
    # docs/DATA_NOTES.md section 1 -- using `Ue` here costs 9-13% on the
    # thicknesses and ~6% on Cf.
    ds["U_ref"] = ds["U"].max(dim="j")
    ds["U_ref"].attrs = {
        "long_name": "velocity scale used by the shipped integral quantities",
        "units": "U_inf",
        "definition": "max over the stored wall-normal range of U",
    }

    # Surface arc length, measured from the first stored station. Needed for
    # any streamwise derivative: differentiating along x/c instead would be
    # wrong by the local surface slope, which is far from 1 near the nose.
    xa_v, ya_v = np.asarray(coords["xa"][1]), np.asarray(coords["ya"][1])
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xa_v), np.diff(ya_v)))])
    ds = ds.assign_coords(xt=("station", s))
    ds["xt"].attrs = {
        "long_name": "arc length along the surface from the first station",
        "units": "c",
    }

    ds.attrs.update(
        airfoil=airfoil,
        side=_SIDE[side],
        side_code=side,
        Re_c=Re_c,
        nu=nu,
        case=CASES.get((side, Re_c), f"{side}{Re_c}"),
        matlab_var=varname,
        aoa_deg=5.0 if airfoil == "naca4412" else 0.0,
        note=(
            "Statistics in the local tangential/wall-normal frame. "
            "Ue and delta99 from Vinuesa et al., Phys. Fluids 28, 055101 (2016). "
            "Lengths normalised by chord, velocities by freestream."
        ),
    )
    # Cross-check: the file's nu should equal 1/Re_c for unit chord and
    # freestream. A mismatch would mean our Re-code mapping is wrong.
    if not np.isclose(nu * Re_c, 1.0, rtol=1e-6):
        ds.attrs["warning"] = f"nu*Re_c = {nu * Re_c:.6g}, expected 1"

    return ds


def load_wing(path: str | Path, airfoil: str | None = None) -> dict[str, xr.Dataset]:
    """Load every case in a WING ``.mat`` file.

    Returns a dict keyed by MATLAB variable name, e.g. ``{'top4n': ...}``.
    """
    path = Path(path)
    if airfoil is None:
        airfoil = path.stem.lower()

    mat = loadmat(path, struct_as_record=True, squeeze_me=True)
    out = {}
    for key in mat:
        if key.startswith("__"):
            continue
        if parse_case_name(key) is None:
            continue
        out[key] = load_case(mat, key, airfoil=airfoil)
    if not out:
        raise ValueError(f"no WING case variables found in {path}")
    return dict(sorted(out.items()))


def concat_cases(datasets: Iterable[xr.Dataset], dim: str = "case") -> xr.Dataset:
    """Stack per-case datasets along a new ``case`` dimension.

    Per-case scalars that live in ``attrs`` (Re_c, nu, side) are promoted to
    coordinates, because once the cases sit in one object the attrs of the
    first would otherwise silently stand for all of them.
    """
    dss = list(datasets)
    labels, Re_c, nus, sides = [], [], [], []
    for ds in dss:
        labels.append(ds.attrs["case"])
        Re_c.append(ds.attrs["Re_c"])
        nus.append(ds.attrs["nu"])
        sides.append(ds.attrs["side"])

    out = xr.concat(dss, dim=dim, combine_attrs="drop_conflicts")
    out = out.assign_coords(
        {
            dim: labels,
            "Re_c": (dim, np.array(Re_c)),
            "nu": (dim, np.array(nus)),
            "side": (dim, np.array(sides)),
        }
    )
    return out
