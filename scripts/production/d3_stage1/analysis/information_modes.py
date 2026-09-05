"""Ordinary covariance/information eigensystems and ellipse summaries."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh
from scipy.stats import chi2


def validate_symmetric_positive_definite(matrix: ArrayLike, *, name: str = "matrix") -> NDArray:
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"{name} must be square")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains non-finite values")
    scale = max(1.0, float(np.max(np.abs(a))))
    if not np.allclose(a, a.T, rtol=1e-10, atol=1e-12 * scale):
        raise ValueError(f"{name} is not symmetric")
    a = 0.5 * (a + a.T)
    if np.min(eigh(a, eigvals_only=True)) <= 0:
        raise ValueError(f"{name} is not positive definite")
    return a


def matrix_diagnostics(matrix: ArrayLike, *, name: str = "matrix") -> dict[str, float | bool]:
    a = np.asarray(matrix, dtype=float)
    symmetric = bool(a.ndim == 2 and a.shape[0] == a.shape[1] and np.allclose(a, a.T))
    eig = eigh(0.5 * (a + a.T), eigvals_only=True) if symmetric else np.array([np.nan])
    return {
        "symmetric": symmetric,
        "positive_definite": bool(np.all(eig > 0)),
        "condition_number": float(np.linalg.cond(a)),
        "minimum_eigenvalue": float(np.min(eig)),
    }


def fisher_from_covariance(covariance: ArrayLike) -> NDArray:
    c = validate_symmetric_positive_definite(covariance, name="covariance")
    # solve is more stable than forming an explicit element-by-element inverse.
    return np.linalg.solve(c, np.eye(c.shape[0]))


def ordinary_eigensystem(information: ArrayLike) -> tuple[NDArray, NDArray]:
    f = validate_symmetric_positive_definite(information, name="information matrix")
    values, vectors = eigh(f)
    order = np.argsort(values)[::-1]
    return values[order], vectors[:, order]


def effective_dimension(eigenvalues: ArrayLike) -> float:
    values = np.asarray(eigenvalues, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)) or np.any(values < 0) or not np.any(values > 0):
        raise ValueError("eigenvalues must be a finite, non-negative vector with positive sum")
    return float(np.sum(values) ** 2 / np.sum(values**2))


def gaussian_fom(covariance: ArrayLike) -> float:
    """Coordinate-dependent Gaussian FoM proportional to inverse ellipse area."""
    c = validate_symmetric_positive_definite(covariance, name="covariance")
    return float(1.0 / np.sqrt(np.linalg.det(c)))


def ellipse_properties(covariance: ArrayLike, *, confidence: float = 0.68) -> dict[str, object]:
    c = validate_symmetric_positive_definite(covariance, name="covariance")
    if c.shape != (2, 2):
        raise ValueError("ellipse properties require a 2x2 covariance")
    values, vectors = eigh(c)
    order = np.argsort(values)[::-1]  # major, minor covariance axes
    values, vectors = values[order], vectors[:, order]
    scale = float(np.sqrt(chi2.ppf(confidence, df=2)))
    semi_axes = scale * np.sqrt(values)
    angle = float(np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0])))
    angle = ((angle + 90.0) % 180.0) - 90.0
    return {
        "confidence": confidence,
        "semi_axes": semi_axes,
        "axis_ratio_minor_to_major": float(semi_axes[1] / semi_axes[0]),
        "orientation_degrees": angle,
        "area": float(np.pi * semi_axes[0] * semi_axes[1]),
        "major_vector": vectors[:, 0],
        "minor_vector": vectors[:, 1],
    }


def ordinary_mode_metrics(covariance: ArrayLike) -> dict[str, object]:
    """Collect ordinary eigenspectrum and 2D ellipse diagnostics."""
    c = validate_symmetric_positive_definite(covariance, name="covariance")
    f = fisher_from_covariance(c)
    values, vectors = ordinary_eigensystem(f)
    result: dict[str, object] = {
        "information": f,
        "eigenvalues": values,
        "eigenvectors": vectors,
        "eigenvalue_ratio_weak_to_strong": float(values[-1] / values[0]),
        "effective_dimension": effective_dimension(values),
        "condition_number": float(np.linalg.cond(c)),
    }
    if c.shape == (2, 2):
        result["ellipse"] = ellipse_properties(c)
    return result
