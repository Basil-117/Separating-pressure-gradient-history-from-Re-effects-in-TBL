"""Load the VinuesaLab APG flat-plate database.

Well-resolved LES of zero-to-adverse-pressure-gradient turbulent boundary
layers on a flat plate, from

    A. Bobke, R. Vinuesa, R. Orlu and P. Schlatter, "History effects and near
    equilibrium in adverse-pressure-gradient turbulent boundary layers",
    J. Fluid Mech. 820, 667-692 (2017).

Why this database matters here: each wing geometry supplies exactly one
``beta(x)`` history, so a history term fitted on wings alone could be fitting
the shape of that single curve. These cases were constructed to have
*different* histories at comparable Reynolds numbers -- two at nominally
constant beta, three following a power-law freestream -- which is what makes
them a real test rather than another sample of the same thing.

Case naming
-----------
``b1n``, ``b2n``
    Nominally constant ``beta`` of 1 and 2.
``m13n``, ``m16n``, ``m18n``
    Power-law freestream ``U_e ~ x^m`` with ``m`` = -0.13, -0.16, -0.18, so
    ``beta`` grows along the plate.

.. warning::
   **The velocity scale differs from the WING database.** Here every shipped
   integral quantity is built on the ``Ue`` field; in WING they are built on
   ``max(U)`` and the ``Ue`` field is something else entirely (see
   docs/DATA_NOTES.md section 1). Verified by inverting ``Cf``, ``Reds`` and
   ``Reth`` for the velocity each implies: all three agree with ``Ue`` to six
   decimal places across all five cases. Reusing the wing convention here
   would be silently wrong.

Units are inlet-based rather than chord-based: ``x`` runs from about 200 to
2500, ``nu`` is 1/450, and lengths are in units of the inlet displacement
thickness. Nothing downstream depends on the choice, because every feature
the model uses is dimensionless.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.io import loadmat

__all__ = ["APG_CASES", "parse_apg_case", "load_apg"]

# Pretty labels and the nominal pressure-gradient history of each case.
APG_CASES = {
    "b1n": ("b1", "constant beta = 1"),
    "b2n": ("b2", "constant beta = 2"),
    "m13n": ("m13", "power law, m = -0.13"),
    "m16n": ("m16", "power law, m = -0.16"),
    "m18n": ("m18", "power law, m = -0.18"),
}

_PROFILE_FIELDS = {
    "U": ("streamwise mean velocity", "U_inf"),
    "V": ("wall-normal mean velocity", "U_inf"),
    "uu": ("streamwise velocity variance", "U_inf^2"),
    "vv": ("wall-normal velocity variance", "U_inf^2"),
    "ww": ("spanwise velocity variance", "U_inf^2"),
    "uv": ("Reynolds shear stress", "U_inf^2"),
    "P": ("mean pressure", "rho U_inf^2"),
    "prms": ("rms pressure fluctuation", "rho U_inf^2"),
    "Pk": ("TKE production", "U_inf^3 / delta*_in"),
    "Dk": ("TKE dissipation", "U_inf^3 / delta*_in"),
    "Tk": ("TKE turbulent transport", "U_inf^3 / delta*_in"),
    "VDk": ("TKE viscous diffusion", "U_inf^3 / delta*_in"),
    "VPGk": ("TKE velocity-pressure-gradient term", "U_inf^3 / delta*_in"),
    "Ck": ("TKE convection", "U_inf^3 / delta*_in"),
}

_STATION_FIELDS = {
    "ut": "u_tau",
    "delta99": "delta99",
    "deltas": "delta_star",
    "theta": "theta",
    "Ue": "Ue",
    "Ret": "Re_tau",
    "Reds": "Re_delta_star",
    "Reth": "Re_theta",
    "Cf": "Cf",
    "beta": "beta_db",
}


def parse_apg_case(name: str) -> tuple[str, str] | None:
    """``'m13n'`` -> ``('m13', 'power law, m = -0.13')``."""
    return APG_CASES.get(name)


def _stack(struct, field: str, n_j: int) -> np.ndarray:
    out = np.full((struct.size, n_j), np.nan)
    for i, e in enumerate(struct):
        a = np.atleast_1d(np.asarray(e[field], dtype=float).squeeze())
        out[i, : min(a.size, n_j)] = a[:n_j]
    return out


def _load_case(mat: dict, var: str, stride: int) -> xr.Dataset:
    parsed = parse_apg_case(var)
    if parsed is None:
        raise ValueError(f"{var!r} is not an APG case variable")
    label, history = parsed

    struct = np.atleast_1d(mat[var])[::stride]
    n_j = max(np.atleast_1d(e["y"]).size for e in struct)
    nu = float(np.asarray(struct[0]["nu"]).squeeze())

    coords = {
        "station": np.arange(struct.size),
        # Flat plate, so streamwise distance and surface arc length coincide.
        # Both names are kept so the wing feature builder works unchanged.
        "xa": ("station", np.array([float(e["x"]) for e in struct])),
        "xt": ("station", np.array([float(e["x"]) for e in struct])),
        "yn": (("station", "j"), _stack(struct, "y", n_j)),
    }

    data_vars = {
        f: (("station", "j"), _stack(struct, f, n_j),
            {"long_name": d, "units": u})
        for f, (d, u) in _PROFILE_FIELDS.items()
    }
    for src, dst in _STATION_FIELDS.items():
        data_vars[dst] = (
            ("station",),
            np.array([float(np.asarray(e[src]).squeeze()) for e in struct]),
        )

    ds = xr.Dataset(data_vars, coords=coords)
    ds["tke"] = 0.5 * (ds["uu"] + ds["vv"] + ds["ww"])

    # The velocity scale. Unlike WING this IS the Ue field; see the module
    # warning. Stored under the same name so downstream code is shared.
    ds["U_ref"] = ds["Ue"].copy()
    ds["U_ref"].attrs = {
        "long_name": "velocity scale used by the shipped integral quantities",
        "units": "U_inf",
        "definition": "the Ue field (NOT max(U), unlike the WING database)",
    }

    ds.attrs.update(
        airfoil="flat_plate",
        side="flat_plate",
        side_code="fp",
        case=label,
        history=history,
        nu=nu,
        Re_c=int(round(1.0 / nu)),
        matlab_var=var,
        velocity_scale="Ue_field",
        aoa_deg=0.0,
        note=(
            "Bobke et al., J. Fluid Mech. 820, 667-692 (2017). "
            "Lengths in units of the inlet displacement thickness."
        ),
    )
    return ds


def load_apg(path: str | Path, stride: int = 1) -> dict[str, xr.Dataset]:
    """Load every case in the APG ``.mat`` file.

    Parameters
    ----------
    stride
        Keep every ``stride``-th streamwise station. The file ships 1570
        stations per case; the full set is useful for resolving the memory
        kernel but makes the n-by-n convolution matrices large, so a stride
        is offered for exploratory work. Default keeps everything.
    """
    mat = loadmat(Path(path), struct_as_record=True, squeeze_me=True)
    out = {}
    for key in mat:
        if not key.startswith("__") and parse_apg_case(key):
            out[key] = _load_case(mat, key, stride)
    if not out:
        raise ValueError(f"no APG case variables found in {path}")
    return dict(sorted(out.items()))
