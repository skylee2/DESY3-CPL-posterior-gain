"""Weighted posterior summaries with explicit normalization conventions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _validated(samples: ArrayLike, weights: ArrayLike | None = None) -> tuple[NDArray, NDArray]:
    x = np.asarray(samples, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    if x.ndim != 2 or x.shape[0] < 2:
        raise ValueError("samples must have shape (n_samples, n_parameters), with n_samples >= 2")
    w = np.ones(x.shape[0], dtype=float) if weights is None else np.asarray(weights, dtype=float)
    if w.ndim != 1 or len(w) != len(x):
        raise ValueError("weights must be one-dimensional and match the number of samples")
    mask = np.all(np.isfinite(x), axis=1) & np.isfinite(w) & (w >= 0)
    x, w = x[mask], w[mask]
    if len(x) < 2 or not np.any(w > 0):
        raise ValueError("fewer than two finite, positive-weight samples remain")
    return x, w


def normalized_weights(weights: ArrayLike) -> NDArray:
    w = np.asarray(weights, dtype=float)
    total = np.sum(w)
    if w.ndim != 1 or np.any(~np.isfinite(w)) or np.any(w < 0) or total <= 0:
        raise ValueError("weights must be finite, non-negative, one-dimensional, and have positive sum")
    return w / total


def effective_sample_size(weights: ArrayLike) -> float:
    """Kish effective sample size, invariant under an overall weight rescaling."""
    w = normalized_weights(weights)
    return float(1.0 / np.sum(w * w))


def weighted_mean(samples: ArrayLike, weights: ArrayLike | None = None) -> NDArray:
    x, w = _validated(samples, weights)
    return np.average(x, axis=0, weights=w)


def weighted_covariance(
    samples: ArrayLike,
    weights: ArrayLike | None = None,
    *,
    normalization: str = "population",
) -> NDArray:
    """Return a symmetric weighted covariance.

    ``population`` is the posterior-moment convention (divide by sum of
    weights). ``unbiased`` applies the reliability-weight correction
    1/(1-sum(normalized_weight**2)).
    """
    x, w = _validated(samples, weights)
    wn = w / np.sum(w)
    mean = np.sum(x * wn[:, None], axis=0)
    dx = x - mean
    cov = (dx * wn[:, None]).T @ dx
    if normalization == "unbiased":
        denom = 1.0 - np.sum(wn * wn)
        if denom <= 0:
            raise ValueError("unbiased covariance is undefined for effective sample size <= 1")
        cov /= denom
    elif normalization != "population":
        raise ValueError("normalization must be 'population' or 'unbiased'")
    return 0.5 * (cov + cov.T)


def weighted_quantile(values: ArrayLike, quantiles: ArrayLike, weights: ArrayLike | None = None) -> NDArray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("values must be one-dimensional")
    x, w = _validated(values[:, None], weights)
    v = x[:, 0]
    order = np.argsort(v)
    v, w = v[order], w[order]
    cdf = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    q = np.asarray(quantiles, dtype=float)
    if np.any((q < 0) | (q > 1)):
        raise ValueError("quantiles must lie in [0, 1]")
    return np.interp(q, cdf, v, left=v[0], right=v[-1])
