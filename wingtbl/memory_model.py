"""A power law with a fitted memory.

The model is deliberately the smallest thing that can answer the question:

.. math::
    y \\;=\\; C\\,\\beta^{a_1}\\;\\bar\\beta(\\ell)^{a_2}\\;Re_\\tau^{a_3}

so each target is a product of a *local* pressure-gradient term, a *history*
term and a *Reynolds-number* term, and the memory length ``ell`` of the
history kernel is fitted rather than chosen. Three things follow from that
choice:

* it is linear in log space, so Shapley attribution is exact and analytic
  rather than sampled -- see :meth:`MemoryPowerLaw.attribute`;
* the fitted exponents are already a symbolic expression, which is what
  Phase 3 refines;
* with roughly 40 usable stations per Reynolds number, three coefficients
  and one length is close to the most this data can support. Phase 1 found
  that ridge in log space beat gradient boosting on every target, which is
  the same lesson.

Fitting uses variable projection: for a fixed ``ell`` the model is linear and
solved exactly by regularised least squares, leaving a well-behaved
one-dimensional search over ``ell``. That is global over the searched
interval and has no learning rate, no initialisation and no stopping rule --
all of which would otherwise need defending.

**Memory length is selected inside the training fold.** A scan that picks
``ell`` by held-out error has already used the test set; :meth:`fit` runs its
own leave-one-case-out on the training cases instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

__all__ = [
    "trapezoid_matrix",
    "kernel_average",
    "MemoryPowerLaw",
    "FitResult",
    "ELL_BOUNDS",
]

# Memory length in multiples of the local delta99.
#
# Upper bound: past roughly 30 delta99 the kernel spans the whole measured
# chord, so beta_bar stops being a memory and becomes a constant; it also
# reaches back into the favourable-pressure-gradient region near the nose,
# where beta < 0 drives the weighted mean negative and the power law is
# undefined. Lower bound: below the station spacing the kernel cannot be
# resolved and beta_bar collapses onto the local beta.
ELL_BOUNDS = (0.02, 30.0)

# The history feature is the *ratio* beta_bar / beta, not beta_bar itself.
#
# Parameterised with beta and beta_bar as separate powers, the two are nearly
# collinear -- at short memory beta_bar -> beta -- and the fit returns huge
# equal-and-opposite exponents (beta^-6.9 beta_bar^+6.9) that are unstable
# and say nothing. Since beta^a1 beta_bar^a2 = beta^(a1+a2) (beta_bar/beta)^a2,
# the ratio form is an exact reparameterisation of the same model class, but
# the two terms now mean different things: `beta` carries the magnitude of
# the local pressure gradient and `mem_ratio` carries how far the upstream
# history departs from it. At short memory mem_ratio -> 1 and its log -> 0,
# so the history term switches itself off rather than fighting the local one.
FEATURES = ("beta", "mem_ratio", "Re_tau")
PRETTY_FEATURE = {
    "beta": r"local $\beta$",
    "mem_ratio": r"history $\bar\beta/\beta$",
    "Re_tau": r"$Re_\tau$",
}
GROUPS = {"beta": "local PG", "mem_ratio": "history", "Re_tau": "Reynolds"}


def trapezoid_matrix(xt: np.ndarray) -> np.ndarray:
    """Lower-triangular quadrature weights.

    ``T[i, j]`` is the trapezoidal weight of node ``j`` when integrating from
    ``xt[0]`` to ``xt[i]``. Precomputed because it does not depend on the
    memory length, which is what makes the one-dimensional search cheap.
    """
    xt = np.asarray(xt, dtype=float)
    n = xt.size
    T = np.zeros((n, n))
    for i in range(1, n):
        d = np.diff(xt[: i + 1])
        w = np.zeros(i + 1)
        w[:-1] += d / 2.0
        w[1:] += d / 2.0
        T[i, : i + 1] = w
    return T


def kernel_average(xt: np.ndarray, beta: np.ndarray, delta99: np.ndarray,
                   ell_over_delta: float, T: np.ndarray | None = None) -> np.ndarray:
    """Exponentially weighted mean of the upstream pressure gradient.

    The memory length at each station is ``ell_over_delta * delta99``, so the
    kernel widens as the layer thickens -- the physically natural scaling,
    and the one that makes the fitted number comparable between cases at
    different Reynolds numbers.
    """
    xt = np.asarray(xt, float)
    beta = np.asarray(beta, float)
    ell = float(ell_over_delta) * np.asarray(delta99, float)
    if T is None:
        T = trapezoid_matrix(xt)

    d = xt[:, None] - xt[None, :]          # distance upstream, negative above
    W = np.where(d >= 0, np.exp(-np.clip(d, 0, None) / ell[:, None]), 0.0)
    num = (T * W * beta[None, :]).sum(axis=1)
    den = (T * W).sum(axis=1)
    out = np.where(den > 1e-14, num / np.where(den > 1e-14, den, 1.0), beta)
    return out


@dataclass
class FitResult:
    """What a fit learned."""

    ell: float
    coef: np.ndarray                     # exponents, in FEATURES order
    intercept: float
    feature_means: np.ndarray            # of log-features, for attribution
    ell_at_bound: bool = False
    condition: float = float("nan")      # of the log design matrix
    inner_scores: dict = field(default_factory=dict)

    def as_formula(self, target: str) -> str:
        terms = " ".join(
            f"{name}^{c:+.3f}" for name, c in zip(FEATURES, self.coef)
        )
        return (f"{target} = {np.exp(self.intercept):.4g} * {terms}"
                f"   [ell = {self.ell:.2f} d99]")


class MemoryPowerLaw:
    """Fit ``y = C beta^a1 beta_bar(ell)^a2 Re_tau^a3`` to grouped data.

    Parameters
    ----------
    alpha
        Ridge penalty applied in log space to the exponents (not the
        intercept). Small by default: with three physically chosen features
        there is little to shrink, and heavy shrinkage would bias the
        exponents that the attribution then reports.
    """

    def __init__(self, alpha: float = 1e-3, ell_bounds: tuple = ELL_BOUNDS):
        self.alpha = alpha
        self.ell_bounds = ell_bounds
        self.result_: FitResult | None = None

    # -- internals --------------------------------------------------------

    @staticmethod
    def _design(df: pd.DataFrame, beta_bar: np.ndarray) -> np.ndarray:
        beta = df["beta"].to_numpy()
        return np.column_stack([
            np.log(beta),
            np.log(beta_bar / beta),      # memory ratio; -> 0 at short memory
            np.log(df["Re_tau"].to_numpy()),
        ])

    def _beta_bar(self, groups: Sequence[dict], ell: float) -> np.ndarray:
        return np.concatenate([
            kernel_average(g["xt"], g["beta"], g["delta99"], ell, g["T"])[g["keep"]]
            for g in groups
        ])

    def _solve(self, X: np.ndarray, logy: np.ndarray) -> tuple[np.ndarray, float]:
        """Ridge in log space, intercept left unpenalised."""
        Xc = X - X.mean(axis=0)
        yc = logy - logy.mean()
        n_f = X.shape[1]
        A = Xc.T @ Xc + self.alpha * np.eye(n_f)
        coef = np.linalg.solve(A, Xc.T @ yc)
        intercept = logy.mean() - X.mean(axis=0) @ coef
        return coef, float(intercept)

    def _inner_score(self, groups_tr: Sequence[dict], logy_tr_by_case: Sequence[np.ndarray],
                     dfs_tr: Sequence[pd.DataFrame], ell: float) -> float:
        """Leave-one-case-out error *within* the training set.

        This is what selects ``ell``. Using the training residual instead
        would favour ever-longer memories, because a wider kernel smooths
        beta_bar towards a constant the model can always absorb.
        """
        n_case = len(groups_tr)
        if n_case < 2:
            # Nothing to hold out; fall back to the training residual.
            bb = self._beta_bar(groups_tr, ell)
            X = self._design(pd.concat(dfs_tr, ignore_index=True), bb)
            y = np.concatenate(logy_tr_by_case)
            coef, b = self._solve(X, y)
            return float(np.mean((X @ coef + b - y) ** 2))

        errs = []
        for h in range(n_case):
            tr_idx = [k for k in range(n_case) if k != h]
            bb_tr = self._beta_bar([groups_tr[k] for k in tr_idx], ell)
            X_tr = self._design(pd.concat([dfs_tr[k] for k in tr_idx],
                                          ignore_index=True), bb_tr)
            y_tr = np.concatenate([logy_tr_by_case[k] for k in tr_idx])
            coef, b = self._solve(X_tr, y_tr)

            bb_h = self._beta_bar([groups_tr[h]], ell)
            X_h = self._design(dfs_tr[h], bb_h)
            errs.append(np.mean((X_h @ coef + b - logy_tr_by_case[h]) ** 2))
        return float(np.mean(errs))

    @staticmethod
    def _groups(df: pd.DataFrame, full: dict[str, pd.DataFrame]) -> tuple:
        """Per-case arrays.

        The kernel integrates over the *whole* chord available for a case,
        then keeps only the in-window stations. Integrating over the window
        alone would discard exactly the upstream record the history term is
        supposed to carry.
        """
        groups, dfs, cases = [], [], list(dict.fromkeys(df["case"]))
        for c in cases:
            whole = full[c]
            sub = df[df["case"] == c]
            keep = np.isin(whole["station"].to_numpy(), sub["station"].to_numpy())
            groups.append(dict(
                xt=whole["xt"].to_numpy(),
                beta=whole["beta"].to_numpy(),
                delta99=whole["delta99"].to_numpy(),
                T=trapezoid_matrix(whole["xt"].to_numpy()),
                keep=keep,
            ))
            dfs.append(sub.reset_index(drop=True))
        return groups, dfs, cases

    # -- API --------------------------------------------------------------

    def fit(self, df: pd.DataFrame, target: str,
            full: dict[str, pd.DataFrame], sign: float = 1.0) -> FitResult:
        """Fit to ``df``; ``full`` maps case -> its untrimmed station table."""
        groups, dfs, _ = self._groups(df, full)
        logy_by_case = [np.log(sign * d[target].to_numpy()) for d in dfs]

        lo, hi = self.ell_bounds
        obj = lambda t: self._inner_score(groups, logy_by_case, dfs, float(np.exp(t)))

        # Coarse grid first: the inner score is not always unimodal, and a
        # bare local search would land in whichever basin it started in.
        grid = np.linspace(np.log(lo), np.log(hi), 40)
        vals = [obj(t) for t in grid]
        k = int(np.argmin(vals))
        a = grid[max(k - 1, 0)]
        b = grid[min(k + 1, grid.size - 1)]
        best_t = (minimize_scalar(obj, bounds=(a, b), method="bounded",
                                  options={"xatol": 1e-3}).x
                  if b > a else grid[k])
        ell = float(np.exp(best_t))

        bb = self._beta_bar(groups, ell)
        X = self._design(pd.concat(dfs, ignore_index=True), bb)
        y = np.concatenate(logy_by_case)
        coef, intercept = self._solve(X, y)

        self.result_ = FitResult(
            ell=ell, coef=coef, intercept=intercept,
            feature_means=X.mean(axis=0),
            ell_at_bound=bool(ell < lo * 1.05 or ell > hi * 0.95),
            # Reported so that collinearity between the terms, which is what
            # made the un-ratioed parameterisation useless, stays visible.
            condition=float(np.linalg.cond(X - X.mean(axis=0))),
            inner_scores=dict(zip(np.exp(grid).round(3), np.round(vals, 6))),
        )
        return self.result_

    def design(self, df: pd.DataFrame, full: dict[str, pd.DataFrame]) -> np.ndarray:
        if self.result_ is None:
            raise RuntimeError("fit first")
        groups, dfs, _ = self._groups(df, full)
        bb = self._beta_bar(groups, self.result_.ell)
        return self._design(pd.concat(dfs, ignore_index=True), bb)

    def predict(self, df: pd.DataFrame, full: dict[str, pd.DataFrame]) -> np.ndarray:
        r = self.result_
        return np.exp(self.design(df, full) @ r.coef + r.intercept)

    def attribute(self, df: pd.DataFrame,
                  full: dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Exact Shapley values, in log space.

        For a linear model the Shapley value of feature *k* is
        ``a_k (x_k - E[x_k])``, with the expectation over the training
        distribution. No sampling, no approximation, no background-set
        choice to defend -- which is worth more here than a tree model's
        extra flexibility.

        Returned values are additive contributions to ``log y``, so a value
        of 0.1 means that feature multiplies the prediction by ``exp(0.1)``.
        """
        r = self.result_
        X = self.design(df, full)
        phi = (X - r.feature_means) * r.coef
        out = pd.DataFrame(phi, columns=[f"phi_{f}" for f in FEATURES])
        out["baseline"] = r.intercept + r.feature_means @ r.coef
        return out
