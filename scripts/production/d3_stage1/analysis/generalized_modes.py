"""Stable generalized information-gain eigenmodes."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh

from .information_modes import validate_symmetric_positive_definite


def generalized_eigensystem(f_new: ArrayLike, f_ref: ArrayLike) -> tuple[NDArray, NDArray]:
    """Solve F_new v = Lambda F_ref v, returning descending gains.

    SciPy's symmetric-definite solver normalizes vectors so V.T F_ref V = I.
    Eigenvectors are directions in the supplied parameter coordinates; their
    component values therefore inherit the chosen parameter units.
    """
    new = validate_symmetric_positive_definite(f_new, name="new information matrix")
    ref = validate_symmetric_positive_definite(f_ref, name="reference information matrix")
    if new.shape != ref.shape:
        raise ValueError("new and reference matrices must have the same shape")
    values, vectors = eigh(new, ref, check_finite=True)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    if np.any(values <= 0):
        raise ValueError("generalized eigenvalues must be positive")
    return values, vectors


def determinant_identity(f_new: ArrayLike, f_ref: ArrayLike, gains: ArrayLike) -> dict[str, float | bool]:
    new = validate_symmetric_positive_definite(f_new, name="new information matrix")
    ref = validate_symmetric_positive_definite(f_ref, name="reference information matrix")
    gains = np.asarray(gains, dtype=float)
    product = float(np.prod(gains))
    _, logdet_new = np.linalg.slogdet(new)
    _, logdet_ref = np.linalg.slogdet(ref)
    ratio = float(np.exp(logdet_new - logdet_ref))
    return {
        "generalized_product": product,
        "information_determinant_ratio": ratio,
        "relative_error": float(abs(product - ratio) / ratio),
        "passes": bool(np.isclose(product, ratio, rtol=1e-8, atol=1e-10)),
    }


def fom_identity(gains: ArrayLike, fom_new: float, fom_ref: float) -> dict[str, float | bool]:
    if fom_new <= 0 or fom_ref <= 0:
        raise ValueError("FoMs must be positive")
    product = float(np.prod(np.asarray(gains, dtype=float)))
    ratio_squared = float((fom_new / fom_ref) ** 2)
    return {
        "generalized_product": product,
        "fom_ratio_squared": ratio_squared,
        "relative_error": float(abs(product - ratio_squared) / ratio_squared),
        "passes": bool(np.isclose(product, ratio_squared, rtol=5e-2)),
    }


def incremental_gain_diagnostics(gains: ArrayLike) -> dict[str, object]:
    """Summarize how the information increments Lambda_i - 1 are distributed."""
    values = np.asarray(gains, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("generalized eigenvalues must be a finite positive vector")
    delta = values - 1.0
    denominator = float(np.sum(delta**2))
    participation = float(np.sum(delta) ** 2 / denominator) if denominator > 0 else float("nan")
    return {"delta": delta, "r_gain": participation}


def mode_angles(vectors: ArrayLike, reference_vectors: ArrayLike) -> NDArray:
    """Acute Euclidean angles (degrees) between corresponding displayed directions."""
    v = np.asarray(vectors, dtype=float)
    r = np.asarray(reference_vectors, dtype=float)
    if v.shape != r.shape:
        raise ValueError("vector matrices must have the same shape")
    vn = v / np.linalg.norm(v, axis=0)
    rn = r / np.linalg.norm(r, axis=0)
    return np.degrees(np.arccos(np.clip(np.abs(np.sum(vn * rn, axis=0)), 0.0, 1.0)))


def pairwise_mode_angles(vectors: ArrayLike, reference_vectors: ArrayLike) -> NDArray:
    """All acute Euclidean angles (degrees) between two sets of displayed directions."""
    v = np.asarray(vectors, dtype=float)
    r = np.asarray(reference_vectors, dtype=float)
    if v.ndim != 2 or r.ndim != 2 or v.shape[0] != r.shape[0]:
        raise ValueError("vector matrices must have the same row dimension")
    vn = v / np.linalg.norm(v, axis=0)
    rn = r / np.linalg.norm(r, axis=0)
    return np.degrees(np.arccos(np.clip(np.abs(vn.T @ rn), 0.0, 1.0)))
