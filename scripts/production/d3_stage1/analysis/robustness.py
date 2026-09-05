"""Prior-boundary and non-Gaussianity robustness estimators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import gaussian_kde

from .generalized_modes import determinant_identity, generalized_eigensystem
from .information_modes import fisher_from_covariance
from .weighted_stats import effective_sample_size, weighted_covariance, weighted_mean, weighted_quantile


@dataclass
class CovarianceEstimate:
    label: str
    covariance: NDArray
    mean: NDArray
    retained_samples: int
    retained_weight: float
    effective_sample_size: float


def prior_boundary_mask(samples: ArrayLike, *, margin: float = 0.0, w0_index: int = 0, wa_index: int = 1) -> NDArray:
    """Keep CPL samples interior to w0 + wa < -margin."""
    x = np.asarray(samples, dtype=float)
    if margin < 0:
        raise ValueError("margin must be non-negative")
    return x[:, w0_index] + x[:, wa_index] < -margin


def covariance_with_mask(samples: ArrayLike, weights: ArrayLike, mask: ArrayLike, label: str) -> CovarianceEstimate:
    x, w, mask = np.asarray(samples, float), np.asarray(weights, float), np.asarray(mask, bool)
    if len(mask) != len(x) or len(w) != len(x) or np.count_nonzero(mask) < 2:
        raise ValueError("mask must retain at least two samples")
    retained_weight = float(np.sum(w[mask]) / np.sum(w))
    subset_weights = np.asarray(w[mask], dtype=float).copy()
    subset_total = subset_weights.sum()
    if not np.isfinite(subset_total) or subset_total <= 0:
        raise ValueError("retained weights must have a finite positive sum")
    subset_weights /= subset_total
    return CovarianceEstimate(
        label=label,
        covariance=weighted_covariance(x[mask], w[mask]),
        mean=weighted_mean(x[mask], w[mask]),
        retained_samples=int(np.count_nonzero(mask)),
        retained_weight=retained_weight,
        effective_sample_size=effective_sample_size(subset_weights),
    )


def full_covariance(samples: ArrayLike, weights: ArrayLike) -> CovarianceEstimate:
    x = np.asarray(samples, float)
    return covariance_with_mask(x, weights, np.ones(len(x), bool), "full weighted posterior")


def boundary_excluded_covariance(samples: ArrayLike, weights: ArrayLike, *, margin: float) -> CovarianceEstimate:
    return covariance_with_mask(samples, weights, prior_boundary_mask(samples, margin=margin), f"w0+wa < -{margin:g}")


def hpd_interior_covariance(
    samples: ArrayLike,
    weights: ArrayLike,
    log_posterior: ArrayLike,
    *,
    retained_mass: float = 0.68,
) -> CovarianceEstimate:
    """Local covariance of the highest-log-posterior samples containing a target mass."""
    x, w, logp = np.asarray(samples, float), np.asarray(weights, float), np.asarray(log_posterior, float)
    if not 0 < retained_mass < 1 or len(x) != len(w) or len(x) != len(logp):
        raise ValueError("invalid retained_mass or mismatched arrays")
    order = np.argsort(logp)[::-1]
    cumulative = np.cumsum(w[order]) / np.sum(w)
    count = max(2, int(np.searchsorted(cumulative, retained_mass, side="left") + 1))
    mask = np.zeros(len(x), dtype=bool)
    mask[order[:count]] = True
    return covariance_with_mask(x, w, mask, f"{retained_mass:.0%} HPD interior")


def local_mode_covariance(
    samples: ArrayLike,
    weights: ArrayLike,
    log_posterior: ArrayLike,
    *,
    log_posterior_drop: float = 2.30 / 2.0,
) -> CovarianceEstimate:
    logp = np.asarray(log_posterior, float)
    mask = logp >= np.max(logp) - log_posterior_drop
    return covariance_with_mask(samples, weights, mask, f"local mode: delta logP <= {log_posterior_drop:g}")


def kde_contour_axes(
    samples: ArrayLike,
    weights: ArrayLike,
    *,
    credible_mass: float = 0.68,
    grid_size: int = 160,
) -> dict[str, object]:
    """Estimate axes from points near a weighted KDE HPD contour in 2D."""
    x, w = np.asarray(samples, float), np.asarray(weights, float)
    if x.ndim != 2 or x.shape[1] != 2 or len(w) != len(x):
        raise ValueError("KDE contour estimator requires 2D samples and matching weights")
    qx = weighted_quantile(x[:, 0], [0.005, 0.995], w)
    qy = weighted_quantile(x[:, 1], [0.005, 0.995], w)
    gx, gy = np.meshgrid(np.linspace(*qx, grid_size), np.linspace(*qy, grid_size))
    points = np.vstack([gx.ravel(), gy.ravel()])
    kde = gaussian_kde(x.T, weights=w)
    density = kde(points).reshape(gx.shape)
    flat = density.ravel()
    order = np.argsort(flat)[::-1]
    cumulative = np.cumsum(flat[order]) / np.sum(flat)
    threshold = flat[order[np.searchsorted(cumulative, credible_mass)]]
    tolerance = max(0.02, 3.0 / grid_size)
    contour = points[:, np.abs(density.ravel() / threshold - 1.0) < tolerance].T
    if len(contour) < 8:
        raise RuntimeError("too few grid points near the KDE contour; increase grid_size")
    center = np.mean(contour, axis=0)
    values, vectors = np.linalg.eigh(np.cov((contour - center).T))
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    return {
        "credible_mass": credible_mass,
        "threshold": float(threshold),
        "center": center,
        "orientation_degrees": float(((angle + 90) % 180) - 90),
        "axis_ratio_minor_to_major": float(np.sqrt(values[1] / values[0])),
        "contour_equivalent_covariance": vectors @ np.diag(values) @ vectors.T,
        "major_vector": vectors[:, 0],
        "minor_vector": vectors[:, 1],
        "contour_points": contour,
    }


def compare_covariance_estimators(
    reference: dict[str, ArrayLike], combined: dict[str, ArrayLike]
) -> dict[str, dict[str, object]]:
    """Recompute generalized gains for covariance estimators with matching labels."""
    labels = sorted(set(reference) & set(combined))
    if not labels:
        raise ValueError("reference and combined estimators have no matching labels")
    output: dict[str, dict[str, object]] = {}
    for label in labels:
        f_ref = fisher_from_covariance(reference[label])
        f_new = fisher_from_covariance(combined[label])
        gains, vectors = generalized_eigensystem(f_new, f_ref)
        output[label] = {
            "generalized_eigenvalues": gains,
            "generalized_eigenvectors": vectors,
            "determinant_check": determinant_identity(f_new, f_ref, gains),
        }
    return output
