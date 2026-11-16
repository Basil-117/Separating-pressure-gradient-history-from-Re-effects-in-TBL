"""Introspect a MATLAB .mat file and print its structure.

Handles both the classic (<= v7) format via scipy.io and the HDF5-based
v7.3 format via h5py, so it works whatever the WING files turn out to be.

Usage:
    python -m wingtbl.matinspect path/to/file.mat [--depth 4] [--max-fields 40]
"""

from __future__ import annotations

import argparse
import os
from typing import Any

import numpy as np


def _fmt_array(a: np.ndarray, name: str = "") -> str:
    kind = a.dtype.kind
    base = f"ndarray{a.shape} {a.dtype}"
    if kind in "fiu" and a.size:
        finite = a[np.isfinite(a)] if kind == "f" else a
        if finite.size:
            base += f"  min={finite.min():.6g} max={finite.max():.6g}"
    return base


def _is_struct(a: Any) -> bool:
    return isinstance(a, np.ndarray) and a.dtype.names is not None


def _walk_scipy(obj: Any, name: str, depth: int, max_depth: int,
                max_fields: int, indent: int = 0) -> None:
    pad = "  " * indent
    if depth > max_depth:
        print(f"{pad}{name}: ... (max depth)")
        return

    if _is_struct(obj):
        fields = list(obj.dtype.names)
        shown = fields[:max_fields]
        print(f"{pad}{name}: struct{obj.shape} with {len(fields)} fields")
        # A 1xN struct array = one entry per streamwise station. Show entry 0
        # as the template, and say how many entries there are.
        flat = obj.ravel()
        if flat.size == 0:
            return
        entry = flat[0]
        if flat.size > 1:
            print(f"{pad}  ({flat.size} entries; showing entry 0)")
        for f in shown:
            _walk_scipy(entry[f], f, depth + 1, max_depth, max_fields, indent + 1)
        if len(fields) > max_fields:
            print(f"{pad}  ... and {len(fields) - max_fields} more fields")
        return

    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            flat = obj.ravel()
            print(f"{pad}{name}: object array{obj.shape}")
            if flat.size:
                _walk_scipy(flat[0], f"{name}[0]", depth + 1, max_depth,
                            max_fields, indent + 1)
            return
        # Squeeze the 1x1 wrappers MATLAB leaves everywhere.
        sq = np.squeeze(obj)
        if sq.ndim == 0:
            print(f"{pad}{name}: scalar {sq.dtype} = {sq.item():.6g}"
                  if sq.dtype.kind in "fiu" else f"{pad}{name}: {sq!r}")
        else:
            print(f"{pad}{name}: {_fmt_array(sq, name)}")
        return

    print(f"{pad}{name}: {type(obj).__name__} = {obj!r}"[:200])


def inspect_scipy(path: str, max_depth: int, max_fields: int) -> bool:
    from scipy.io import loadmat
    try:
        mat = loadmat(path, struct_as_record=True, squeeze_me=False)
    except NotImplementedError:
        return False  # v7.3 / HDF5
    except Exception as exc:  # noqa: BLE001
        print(f"scipy.io.loadmat failed: {exc}")
        return False

    keys = [k for k in mat if not k.startswith("__")]
    print(f"format: MATLAB <= v7 (scipy.io)")
    print(f"top-level variables ({len(keys)}): {', '.join(sorted(keys))}\n")
    for k in sorted(keys):
        _walk_scipy(mat[k], k, 0, max_depth, max_fields)
        print()
    return True


def inspect_h5(path: str, max_depth: int, max_fields: int) -> None:
    import h5py

    print("format: MATLAB v7.3 (HDF5)")

    def visit(name: str, obj: Any) -> None:
        depth = name.count("/")
        if depth > max_depth:
            return
        pad = "  " * depth
        if isinstance(obj, h5py.Dataset):
            print(f"{pad}{name}: dataset{obj.shape} {obj.dtype}")
        else:
            print(f"{pad}{name}/  ({len(obj)} members)")

    with h5py.File(path, "r") as f:
        f.visititems(visit)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path")
    p.add_argument("--depth", type=int, default=4)
    p.add_argument("--max-fields", type=int, default=60)
    args = p.parse_args()

    size_mb = os.path.getsize(args.path) / 1e6
    print(f"file: {args.path}  ({size_mb:.1f} MB)\n")

    if not inspect_scipy(args.path, args.depth, args.max_fields):
        inspect_h5(args.path, args.depth, args.max_fields)


if __name__ == "__main__":
    main()
