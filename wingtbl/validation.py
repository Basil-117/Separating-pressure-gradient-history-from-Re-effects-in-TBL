"""Cross-validation schemes that do not leak.

The obvious mistake with this dataset is a random train/test split. Adjacent
stations on the same wing are nearly the same boundary layer: a split that
puts x/c = 0.60 in train and x/c = 0.61 in test reports an error close to the
interpolation error of the grid, which says nothing about whether a model
generalises. Both splitters here keep correlated samples together.

:func:`leave_one_case_out` is the honest test and the one to quote. A model
that is only asked to interpolate between Reynolds numbers it has already
seen has not shown that it separates Reynolds-number effects from
pressure-gradient effects -- it may simply have memorised each case.
"""

from __future__ import annotations

from typing import Iterator, Sequence

import numpy as np
import pandas as pd

__all__ = [
    "leave_one_case_out",
    "blocked_streamwise",
    "relative_error",
    "summarise_errors",
]


def leave_one_case_out(
    df: pd.DataFrame, case_col: str = "case",
    cases: Sequence[str] | None = None,
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """Yield ``(held_out_case, train_idx, test_idx)``.

    Parameters
    ----------
    cases
        Restrict which cases are held out. Useful for excluding a case from
        being a *test* fold while still allowing it in training, or the
        reverse.
    """
    all_cases = list(dict.fromkeys(df[case_col]))
    for c in (cases if cases is not None else all_cases):
        test = np.flatnonzero((df[case_col] == c).to_numpy())
        train = np.flatnonzero((df[case_col] != c).to_numpy())
        if test.size and train.size:
            yield c, train, test


def blocked_streamwise(
    df: pd.DataFrame, n_blocks: int = 5, case_col: str = "case",
    x_col: str = "xa",
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """Contiguous chordwise blocks, held out from every case at once.

    Tests extrapolation in pressure-gradient strength rather than across
    Reynolds number: the model must predict a stretch of the wing it has
    never seen, at any Reynolds number.
    """
    x = df[x_col].to_numpy()
    edges = np.quantile(x, np.linspace(0, 1, n_blocks + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    for b in range(n_blocks):
        mask = (x >= edges[b]) & (x < edges[b + 1])
        test = np.flatnonzero(mask)
        train = np.flatnonzero(~mask)
        if test.size and train.size:
            lo, hi = x[test].min(), x[test].max()
            yield f"x/c {lo:.2f}-{hi:.2f}", train, test


def relative_error(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    """Signed relative error, guarded against targets that cross zero.

    ``uv_outer`` is negative throughout and ``Cf`` approaches zero near
    separation, so a bare ``(pred - true) / true`` is unstable. Scaling by
    ``|true|`` keeps the sign meaningful and the magnitude finite.
    """
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    denom = np.maximum(np.abs(y_true), 1e-12)
    return (y_pred - y_true) / denom


def summarise_errors(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Median and 90th-percentile absolute relative error, plus R^2.

    The median is quoted rather than the mean because the trailing-edge
    stations, where inner scaling diverges as u_tau -> 0, otherwise dominate
    every average.
    """
    e = np.abs(relative_error(y_true, y_pred))
    y_true = np.asarray(y_true, float)
    ss_res = float(np.sum((y_true - np.asarray(y_pred, float)) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    return {
        "median_rel": float(np.median(e)),
        "p90_rel": float(np.percentile(e, 90)),
        "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        "n": int(y_true.size),
    }
