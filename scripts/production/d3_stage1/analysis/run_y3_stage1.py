"""DES Y3 Stage-1 real-chain validation (not the primary science result)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .config import ChainConfig
from .generalized_modes import generalized_eigensystem, incremental_gain_diagnostics, pairwise_mode_angles
from .information_modes import fisher_from_covariance, ordinary_mode_metrics
from .io_chains import ChainData, read_chain
from .robustness import (
    boundary_excluded_covariance,
    full_covariance,
    hpd_interior_covariance,
    kde_contour_axes,
    local_mode_covariance,
)
from .validate_cpl import _json_default, analyze_pair
from .weighted_stats import weighted_covariance, weighted_quantile


LABEL = "STAGE 1 VALIDATION — NOT PRIMARY SCIENTIFIC RESULT"
PARAMETER_MAP = {
    "w0": "cosmological_parameters--w",
    "wa": "cosmological_parameters--wa",
}


def _config(path: Path) -> ChainConfig:
    return ChainConfig(
        path=path,
        parameter_map=PARAMETER_MAP,
        weight_column="weight",
        log_posterior_column="post",
    )


def _covariance_estimators(chain: ChainData) -> dict[str, object]:
    if chain.log_posterior is None:
        raise ValueError("Stage-1 robustness requires the stored post column")
    estimates = [full_covariance(chain.samples, chain.weights)]
    for margin in (0.05, 0.10, 0.20):
        estimates.append(boundary_excluded_covariance(chain.samples, chain.weights, margin=margin))
    estimates.append(hpd_interior_covariance(chain.samples, chain.weights, chain.log_posterior, retained_mass=0.68))
    try:
        estimates.append(local_mode_covariance(chain.samples, chain.weights, chain.log_posterior))
    except ValueError:
        pass
    return {item.label: item for item in estimates}


def _acute_angle_degrees(first: np.ndarray, second: np.ndarray) -> float:
    a, b = first / np.linalg.norm(first), second / np.linalg.norm(second)
    return float(np.degrees(np.arccos(np.clip(abs(a @ b), 0.0, 1.0))))


def _pivot_diagnostics(chain: ChainData) -> dict[str, object]:
    covariance = weighted_covariance(chain.samples, chain.weights)
    a_pivot = 1.0 + covariance[0, 1] / covariance[1, 1]
    w_pivot = chain.samples[:, 0] + (1.0 - a_pivot) * chain.samples[:, 1]
    return {
        "a_pivot": float(a_pivot),
        "z_pivot": float(1.0 / a_pivot - 1.0),
        "w_pivot_quantiles_16_50_84": weighted_quantile(
            w_pivot, [0.16, 0.50, 0.84], chain.weights
        ),
        "w0_quantiles_2p5_5": weighted_quantile(
            chain.samples[:, 0], [0.025, 0.05], chain.weights
        ),
    }


def _bootstrap(
    reference: ChainData,
    combined: ChainData,
    *,
    repetitions: int,
    seed: int,
) -> dict[str, object]:
    """Multiplier bootstrap of weighted posterior moments for numerical stability."""
    rng = np.random.default_rng(seed)
    full_ref = weighted_covariance(reference.samples, reference.weights)
    full_new = weighted_covariance(combined.samples, combined.weights)
    full_gains, full_vectors = generalized_eigensystem(
        fisher_from_covariance(full_new), fisher_from_covariance(full_ref)
    )
    records = []
    for _ in range(repetitions):
        # Exponential multipliers implement a Bayesian/multiplier bootstrap while
        # retaining the released nested-sampling weights.
        wr = reference.weights * rng.exponential(size=len(reference.weights))
        wn = combined.weights * rng.exponential(size=len(combined.weights))
        cr = weighted_covariance(reference.samples, wr)
        cn = weighted_covariance(combined.samples, wn)
        gains, vectors = generalized_eigensystem(fisher_from_covariance(cn), fisher_from_covariance(cr))
        ref_axes = np.column_stack([
            ordinary_mode_metrics(cr)["ellipse"]["major_vector"],
            ordinary_mode_metrics(cr)["ellipse"]["minor_vector"],
        ])
        records.append(
            np.r_[
                gains,
                [_acute_angle_degrees(vectors[:, i], full_vectors[:, i]) for i in range(2)],
                pairwise_mode_angles(vectors, ref_axes).ravel(),
                incremental_gain_diagnostics(gains)["r_gain"],
            ]
        )
    values = np.asarray(records)
    names = [
        "Lambda1", "Lambda2",
        "gain1_angle_to_full", "gain2_angle_to_full",
        "gain1_to_ref_major", "gain1_to_ref_minor",
        "gain2_to_ref_major", "gain2_to_ref_minor",
        "r_gain",
    ]
    return {
        "method": "independent exponential-multiplier bootstrap of released weighted samples",
        "repetitions": repetitions,
        "seed": seed,
        "full_sample_generalized_eigenvalues": full_gains,
        "quantiles_16_50_84": {
            name: np.quantile(values[:, index], [0.16, 0.50, 0.84])
            for index, name in enumerate(names)
        },
        "fraction_Lambda2_above_one": float(np.mean(values[:, 1] > 1.0)),
    }


def run(reference_path: Path, combined_path: Path, output_dir: Path, *, bootstrap: int, seed: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    ref, new = read_chain(_config(reference_path)), read_chain(_config(combined_path))
    core = analyze_pair(_config(reference_path), _config(combined_path), published_foms=None)
    estimates_ref, estimates_new = _covariance_estimators(ref), _covariance_estimators(new)
    robustness: dict[str, object] = {}
    for label in sorted(set(estimates_ref) & set(estimates_new)):
        er, en = estimates_ref[label], estimates_new[label]
        gains, vectors = generalized_eigensystem(
            fisher_from_covariance(en.covariance), fisher_from_covariance(er.covariance)
        )
        robustness[label] = {
            "reference": er.__dict__,
            "combined": en.__dict__,
            "generalized_eigenvalues": gains,
            "generalized_eigenvectors": vectors,
            "incremental_gain": incremental_gain_diagnostics(gains),
        }

    # KDE is used only as a non-Gaussian orientation diagnostic. A deterministic
    # weighted resample caps its cost and is not used for the Fisher calculation.
    rng = np.random.default_rng(seed)
    kde = {}
    for name, chain in (("reference", ref), ("combined", new)):
        count = min(5000, len(chain.samples))
        indices = rng.choice(len(chain.samples), size=count, replace=True, p=chain.weights)
        kde[name] = kde_contour_axes(
            chain.samples[indices], np.ones(count), credible_mass=0.68, grid_size=120
        )
        kde[name].pop("contour_points")

    payload = {
        "label": LABEL,
        "scientific_scope_warning": (
            "Adds BAO+RSD+Pantheon SN to DES Y3 3x2pt and reverses the conditioning of the "
            "planned primary geometry-to-LSS test. Do not use as the paper's primary result."
        ),
        "core": core,
        "published_quantity_validation": {
            "reference": _pivot_diagnostics(ref),
            "combined": _pivot_diagnostics(new),
            "targets": {
                "DES_Y3_3x2pt_z_pivot": 0.24,
                "DES_Y3_3x2pt_plus_BAO_RSD_SN_z_pivot": 0.21,
                "DES_Y3_3x2pt_w0_95_percent_lower_bound_approximately": -1.4,
            },
        },
        "robustness": robustness,
        "kde_orientation_diagnostic": kde,
        "bootstrap": _bootstrap(ref, new, repetitions=bootstrap, seed=seed),
    }
    (output_dir / "stage1_validation.json").write_text(
        json.dumps(payload, default=_json_default, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "README.md").write_text(
        f"# {LABEL}\n\n"
        "This directory contains a real-public-chain pipeline validation for DES Y3 3x2pt "
        "to DES Y3 3x2pt+BAO+RSD+Pantheon. It is not the primary geometry-to-growth result.\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("combined", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=400)
    parser.add_argument("--seed", type=int, default=220705766)
    args = parser.parse_args()
    run(args.reference, args.combined, args.output_dir, bootstrap=args.bootstrap, seed=args.seed)


if __name__ == "__main__":
    main()
