"""First-pass marginalized (w0, wa) validation pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .config import ChainConfig
from .generalized_modes import (
    determinant_identity,
    fom_identity,
    generalized_eigensystem,
    incremental_gain_diagnostics,
    pairwise_mode_angles,
)
from .information_modes import fisher_from_covariance, gaussian_fom, ordinary_mode_metrics
from .io_chains import read_chain
from .posterior_shift import mean_shift
from .weighted_stats import weighted_covariance, weighted_mean, weighted_quantile


def analyze_pair(ref_config: ChainConfig, new_config: ChainConfig, *, published_foms: tuple[float, float] | None = None) -> dict[str, object]:
    ref, new = read_chain(ref_config), read_chain(new_config)
    if ref.names != ("w0", "wa") or new.names != ("w0", "wa"):
        raise ValueError("first validation is intentionally restricted to canonical parameters ('w0', 'wa')")
    c_ref = weighted_covariance(ref.samples, ref.weights)
    c_new = weighted_covariance(new.samples, new.weights)
    mu_ref = weighted_mean(ref.samples, ref.weights)
    mu_new = weighted_mean(new.samples, new.weights)
    f_ref, f_new = fisher_from_covariance(c_ref), fisher_from_covariance(c_new)
    gains, gain_vectors = generalized_eigensystem(f_new, f_ref)
    ordinary_ref, ordinary_new = ordinary_mode_metrics(c_ref), ordinary_mode_metrics(c_new)
    geometry_axes = np.column_stack([
        ordinary_ref["ellipse"]["major_vector"],
        ordinary_ref["ellipse"]["minor_vector"],
    ])
    result: dict[str, object] = {
        "ref_metadata": ref.metadata,
        "new_metadata": new.metadata,
        "ref_mean": mu_ref,
        "new_mean": mu_new,
        "ref_intervals_16_50_84": np.vstack([weighted_quantile(ref.samples[:, i], [.16, .5, .84], ref.weights) for i in range(2)]),
        "new_intervals_16_50_84": np.vstack([weighted_quantile(new.samples[:, i], [.16, .5, .84], new.weights) for i in range(2)]),
        "ref_covariance": c_ref,
        "new_covariance": c_new,
        "ref_gaussian_fom": gaussian_fom(c_ref),
        "new_gaussian_fom": gaussian_fom(c_new),
        "ordinary_ref": ordinary_ref,
        "ordinary_new": ordinary_new,
        "generalized_eigenvalues": gains,
        "generalized_eigenvectors": gain_vectors,
        "incremental_gain": incremental_gain_diagnostics(gains),
        "pairwise_display_angles_to_geometry_axes_degrees": pairwise_mode_angles(gain_vectors, geometry_axes),
        "determinant_check": determinant_identity(f_new, f_ref, gains),
        "mean_shift": mean_shift(mu_ref, mu_new, c_ref),
        "hard_prior_w0_plus_wa": {
            "reference_maximum": float(np.max(np.sum(ref.samples, axis=1))),
            "new_maximum": float(np.max(np.sum(new.samples, axis=1))),
            "reference_weight_within_0p1": float(np.sum(ref.weights[np.sum(ref.samples, axis=1) > -0.1])),
            "new_weight_within_0p1": float(np.sum(new.weights[np.sum(new.samples, axis=1) > -0.1])),
            "warning": "The hard prior w0+wa<0 can truncate broad CPL posteriors; use the dedicated robustness estimators.",
        },
    }
    if published_foms is not None:
        result["published_fom_check"] = fom_identity(gains, published_foms[1], published_foms[0])
    return result


def _json_default(value: object) -> object:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    raise TypeError(type(value).__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("combined", type=Path)
    parser.add_argument("--w0-column", default="w0")
    parser.add_argument("--wa-column", default="wa")
    parser.add_argument("--weight-column", default="weight")
    parser.add_argument("--log-weight-column")
    parser.add_argument("--published-fom-reference", type=float)
    parser.add_argument("--published-fom-combined", type=float)
    parser.add_argument("--output", type=Path, default=Path("results/cpl_validation.json"))
    args = parser.parse_args()
    weight = None if args.log_weight_column else args.weight_column
    common = dict(parameter_map={"w0": args.w0_column, "wa": args.wa_column},
                  weight_column=weight, log_weight_column=args.log_weight_column)
    supplied_foms = (args.published_fom_reference, args.published_fom_combined)
    if (supplied_foms[0] is None) != (supplied_foms[1] is None):
        parser.error("supply both published FoMs or neither")
    published_foms = None if supplied_foms[0] is None else supplied_foms
    result = analyze_pair(
        ChainConfig(path=args.reference, **common),
        ChainConfig(path=args.combined, **common),
        published_foms=published_foms,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, default=_json_default, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
