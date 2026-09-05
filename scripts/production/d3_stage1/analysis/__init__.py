"""DES Y6 CPL information-mode analysis package."""

from .generalized_modes import generalized_eigensystem
from .information_modes import effective_dimension, ellipse_properties, fisher_from_covariance, ordinary_eigensystem
from .posterior_shift import mean_shift
from .weighted_stats import weighted_covariance, weighted_mean

__all__ = [
    "weighted_mean",
    "weighted_covariance",
    "fisher_from_covariance",
    "ordinary_eigensystem",
    "generalized_eigensystem",
    "effective_dimension",
    "ellipse_properties",
    "mean_shift",
]
