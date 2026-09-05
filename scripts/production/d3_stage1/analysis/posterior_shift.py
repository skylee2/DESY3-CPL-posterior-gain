"""Posterior-location diagnostics, kept separate from covariance contraction."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .information_modes import validate_symmetric_positive_definite


def mean_shift(mu1: ArrayLike, mu2: ArrayLike, covariance: ArrayLike | None = None) -> dict[str, object]:
    first, second = np.asarray(mu1, dtype=float), np.asarray(mu2, dtype=float)
    if first.shape != second.shape or first.ndim != 1:
        raise ValueError("means must be one-dimensional vectors with matching shape")
    delta = second - first
    result: dict[str, object] = {"delta": delta, "euclidean_norm": float(np.linalg.norm(delta))}
    if covariance is not None:
        cov = validate_symmetric_positive_definite(covariance, name="shift covariance")
        if cov.shape[0] != len(delta):
            raise ValueError("covariance dimension does not match means")
        d2 = float(delta @ np.linalg.solve(cov, delta))
        result.update({"mahalanobis_squared": d2, "mahalanobis_distance": float(np.sqrt(d2))})
    return result


def lcdm_displacement(mean: ArrayLike, covariance: ArrayLike, lcdm: ArrayLike = (-1.0, 0.0)) -> float:
    mean, lcdm = np.asarray(mean, dtype=float), np.asarray(lcdm, dtype=float)
    if mean.shape != lcdm.shape:
        raise ValueError("mean and LCDM reference point must have matching shape")
    cov = validate_symmetric_positive_definite(covariance, name="covariance")
    delta = mean - lcdm
    return float(np.sqrt(delta @ np.linalg.solve(cov, delta)))
