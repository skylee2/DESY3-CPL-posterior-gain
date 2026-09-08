#!/usr/bin/env python3
"""Reproducible weighted-Gaussian validation of the production gain pipeline.

Run from any directory with

    python3 scripts/validate_synthetic_gaussian.py

The default machine-readable result is written to
``results/validation/synthetic_gaussian_validation.json``. This is the current
standalone validation for manuscript v3. Earlier undocumented validation
attempts are superseded and are not reconstructed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import scipy


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_ROOT = PACKAGE_ROOT / "scripts" / "production" / "d3_stage1"
sys.path.insert(0, str(PRODUCTION_ROOT))

from analysis.generalized_modes import (  # noqa: E402
    determinant_identity,
    generalized_eigensystem,
    incremental_gain_diagnostics,
)
from analysis.information_modes import fisher_from_covariance  # noqa: E402
from analysis.weighted_stats import (  # noqa: E402
    effective_sample_size,
    normalized_weights,
    weighted_covariance,
)


SEED = 20260907
N_REF = 1_000_000
N_NEW = 1_000_000
MEAN_REF = np.array([0.0, 0.0])
MEAN_NEW = np.array([0.0, 0.0])
F_REF = np.diag([100.0, 1.0])
F_NEW = np.diag([1000.0, 5.0])
C_REF = np.linalg.inv(F_REF)
C_NEW = np.linalg.inv(F_NEW)
LOGNORMAL_SIGMA = 0.35

ANALYTIC_LAMBDA = np.array([10.0, 5.0])
ANALYTIC_DETERMINANT_RATIO = 50.0
ANALYTIC_R_GAIN = 169.0 / 97.0

# Monte Carlo tolerances are deliberately fixed in the archived script.
COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE = 5.0e-3
EIGENVALUE_RELATIVE_TOLERANCE = 7.5e-3
DETERMINANT_RELATIVE_TOLERANCE = 1.0e-2
R_GAIN_ABSOLUTE_TOLERANCE = 5.0e-3


def draw_gaussian(rng: np.random.Generator, n: int, mean: np.ndarray, covariance: np.ndarray) -> np.ndarray:
    """Draw a Gaussian sample using an explicit standard-normal/Cholesky map."""
    z = rng.standard_normal((n, len(mean)))
    return mean + z @ np.linalg.cholesky(covariance).T


def draw_weights(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return raw and normalized nonuniform weights independent of samples.

    The prescription is raw_w[i] = exp(sigma*z[i] - sigma**2/2), where
    z[i] are independent standard-normal variates and sigma=0.35.  The
    subtraction makes E[raw_w]=1 and has no effect after normalization.
    """
    z = rng.standard_normal(n)
    raw = np.exp(LOGNORMAL_SIGMA * z - 0.5 * LOGNORMAL_SIGMA**2)
    return raw, normalized_weights(raw)


def relative_frobenius_error(measured: np.ndarray, target: np.ndarray) -> float:
    return float(np.linalg.norm(measured - target) / np.linalg.norm(target))


def array_sha256(values: np.ndarray) -> str:
    """Hash a canonical little-endian float64 representation of an array."""
    canonical = np.asarray(values, dtype="<f8")
    return hashlib.sha256(canonical.tobytes(order="C")).hexdigest()


def test_record(name: str, passed: bool, **details: object) -> dict[str, object]:
    return {"name": name, "passed": bool(passed), **details}


def run_validation() -> dict[str, object]:
    # Spawn order is part of the reproducibility prescription: ref samples,
    # new samples, ref weights, then new weights.  Streams are independent.
    streams = [np.random.Generator(np.random.PCG64(s)) for s in np.random.SeedSequence(SEED).spawn(4)]
    ref_samples = draw_gaussian(streams[0], N_REF, MEAN_REF, C_REF)
    new_samples = draw_gaussian(streams[1], N_NEW, MEAN_NEW, C_NEW)
    ref_raw_weights, ref_weights = draw_weights(streams[2], N_REF)
    new_raw_weights, new_weights = draw_weights(streams[3], N_NEW)

    # This is the same production-function path used for DES chain pairs.
    c_ref_est = weighted_covariance(ref_samples, ref_weights, normalization="population")
    c_new_est = weighted_covariance(new_samples, new_weights, normalization="population")
    f_ref_est = fisher_from_covariance(c_ref_est)
    f_new_est = fisher_from_covariance(c_new_est)
    gains, _ = generalized_eigensystem(f_new_est, f_ref_est)
    determinant = determinant_identity(f_new_est, f_ref_est, gains)
    gain_diagnostics = incremental_gain_diagnostics(gains)

    neff_ref = effective_sample_size(ref_weights)
    neff_new = effective_sample_size(new_weights)
    covariance_error_ref = relative_frobenius_error(c_ref_est, C_REF)
    covariance_error_new = relative_frobenius_error(c_new_est, C_NEW)
    determinant_ratio = float(determinant["information_determinant_ratio"])
    r_gain = float(gain_diagnostics["r_gain"])

    exact_gains, _ = generalized_eigensystem(F_NEW, F_REF)
    tests = [
        test_record(
            "analytic_covariances_are_matrix_inverses",
            np.allclose(C_REF @ F_REF, np.eye(2), rtol=0.0, atol=1e-14)
            and np.allclose(C_NEW @ F_NEW, np.eye(2), rtol=0.0, atol=1e-14),
            absolute_tolerance=1e-14,
        ),
        test_record(
            "analytic_generalized_eigenvalues",
            np.allclose(exact_gains, ANALYTIC_LAMBDA, rtol=0.0, atol=1e-13),
            measured=exact_gains.tolist(),
            target=ANALYTIC_LAMBDA.tolist(),
            absolute_tolerance=1e-13,
        ),
        test_record(
            "reference_weights_normalized",
            np.isclose(np.sum(ref_weights), 1.0, rtol=0.0, atol=5e-15),
            measured_sum=float(np.sum(ref_weights)),
            target=1.0,
            absolute_tolerance=5e-15,
        ),
        test_record(
            "new_weights_normalized",
            np.isclose(np.sum(new_weights), 1.0, rtol=0.0, atol=5e-15),
            measured_sum=float(np.sum(new_weights)),
            target=1.0,
            absolute_tolerance=5e-15,
        ),
        test_record(
            "synthetic_weights_are_nonuniform",
            np.ptp(ref_raw_weights) > 0.0 and np.ptp(new_raw_weights) > 0.0,
            reference_raw_weight_range=[float(np.min(ref_raw_weights)), float(np.max(ref_raw_weights))],
            new_raw_weight_range=[float(np.min(new_raw_weights)), float(np.max(new_raw_weights))],
        ),
        test_record(
            "effective_sample_sizes_are_valid",
            1.0 < neff_ref < N_REF and 1.0 < neff_new < N_NEW,
            measured={"reference": neff_ref, "new": neff_new},
            required_ranges={"reference": [1.0, N_REF], "new": [1.0, N_NEW]},
        ),
        test_record(
            "reference_weighted_covariance",
            covariance_error_ref < COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE,
            measured_relative_frobenius_error=covariance_error_ref,
            tolerance=COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE,
        ),
        test_record(
            "new_weighted_covariance",
            covariance_error_new < COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE,
            measured_relative_frobenius_error=covariance_error_new,
            tolerance=COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE,
        ),
        test_record(
            "sample_generalized_eigenvalues",
            np.allclose(gains, ANALYTIC_LAMBDA, rtol=EIGENVALUE_RELATIVE_TOLERANCE, atol=0.0),
            measured=gains.tolist(),
            target=ANALYTIC_LAMBDA.tolist(),
            relative_tolerance=EIGENVALUE_RELATIVE_TOLERANCE,
        ),
        test_record(
            "production_determinant_closure",
            bool(determinant["passes"]),
            generalized_product=float(determinant["generalized_product"]),
            information_determinant_ratio=determinant_ratio,
            measured_relative_error=float(determinant["relative_error"]),
            production_relative_tolerance=1e-8,
            production_absolute_tolerance=1e-10,
        ),
        test_record(
            "sample_determinant_ratio",
            np.isclose(determinant_ratio, ANALYTIC_DETERMINANT_RATIO, rtol=DETERMINANT_RELATIVE_TOLERANCE, atol=0.0),
            measured=determinant_ratio,
            target=ANALYTIC_DETERMINANT_RATIO,
            relative_tolerance=DETERMINANT_RELATIVE_TOLERANCE,
        ),
        test_record(
            "sample_r_gain",
            np.isclose(r_gain, ANALYTIC_R_GAIN, rtol=0.0, atol=R_GAIN_ABSOLUTE_TOLERANCE),
            measured=r_gain,
            target=ANALYTIC_R_GAIN,
            absolute_tolerance=R_GAIN_ABSOLUTE_TOLERANCE,
        ),
    ]

    result = {
        "schema_version": 1,
        "description": "Sole standalone manuscript-v3 synthetic validation in release candidate v1.1.0; earlier undocumented validation attempts are superseded.",
        "software": {"numpy": np.__version__, "scipy": scipy.__version__},
        "rng": {
            "seed": SEED,
            "bit_generator": "PCG64",
            "stream_construction": "numpy.random.SeedSequence(seed).spawn(4)",
            "spawn_order": ["reference_samples", "new_samples", "reference_weights", "new_weights"],
        },
        "sample_sizes": {"reference": N_REF, "new": N_NEW},
        "means": {"reference": MEAN_REF.tolist(), "new": MEAN_NEW.tolist()},
        "analytic": {
            "F_ref": F_REF.tolist(),
            "F_new": F_NEW.tolist(),
            "C_ref": C_REF.tolist(),
            "C_new": C_NEW.tolist(),
            "generalized_eigenvalues": ANALYTIC_LAMBDA.tolist(),
            "determinant_ratio": ANALYTIC_DETERMINANT_RATIO,
            "r_gain": ANALYTIC_R_GAIN,
        },
        "weights": {
            "prescription": "raw_w[i] = exp(0.35*z[i] - 0.5*0.35^2), z[i] iid N(0,1), independently generated for each posterior and independently of samples; normalized_w = raw_w/sum(raw_w)",
            "lognormal_sigma": LOGNORMAL_SIGMA,
            "normalized_weight_sums": {
                "reference": float(np.sum(ref_weights)),
                "new": float(np.sum(new_weights)),
            },
            "normalized_weight_ranges": {
                "reference": [float(np.min(ref_weights)), float(np.max(ref_weights))],
                "new": [float(np.min(new_weights)), float(np.max(new_weights))],
            },
            "normalized_weight_float64_sha256": {
                "reference": array_sha256(ref_weights),
                "new": array_sha256(new_weights),
            },
            "effective_sample_size_definition": "1/sum(normalized_w^2)",
            "effective_sample_sizes": {"reference": neff_ref, "new": neff_new},
        },
        "production_pipeline": [
            "analysis.weighted_stats.weighted_covariance(normalization='population')",
            "analysis.information_modes.fisher_from_covariance",
            "analysis.generalized_modes.generalized_eigensystem",
            "analysis.generalized_modes.determinant_identity",
            "analysis.generalized_modes.incremental_gain_diagnostics",
        ],
        "estimated": {
            "C_ref": c_ref_est.tolist(),
            "C_new": c_new_est.tolist(),
            "F_ref": f_ref_est.tolist(),
            "F_new": f_new_est.tolist(),
            "Lambda1": float(gains[0]),
            "Lambda2": float(gains[1]),
            "generalized_eigenvalues": gains.tolist(),
            "generalized_eigenvalue_product": float(np.prod(gains)),
            "determinant_ratio": determinant_ratio,
            "covariance_determinant_ratio": float(np.linalg.det(c_ref_est) / np.linalg.det(c_new_est)),
            "r_gain": r_gain,
            "N_eff_ref": neff_ref,
            "N_eff_new": neff_new,
            "seed": SEED,
            "N_ref": N_REF,
            "N_new": N_NEW,
        },
        "tolerances": {
            "covariance_relative_frobenius": COVARIANCE_RELATIVE_FROBENIUS_TOLERANCE,
            "generalized_eigenvalue_relative": EIGENVALUE_RELATIVE_TOLERANCE,
            "determinant_ratio_relative": DETERMINANT_RELATIVE_TOLERANCE,
            "r_gain_absolute": R_GAIN_ABSOLUTE_TOLERANCE,
        },
        "tests": tests,
        "test_summary": {
            "passed": int(sum(test["passed"] for test in tests)),
            "total": len(tests),
            "all_passed": bool(all(test["passed"] for test in tests)),
        },
    }

    for test in tests:
        status = "PASS" if test["passed"] else "FAIL"
        print(f"[{status}] {test['name']}")
    print(json.dumps(result["estimated"], indent=2))
    print(f"Tests passed: {result['test_summary']['passed']}/{result['test_summary']['total']}")

    failed = [test["name"] for test in tests if not test["passed"]]
    assert not failed, "Synthetic validation failures: " + ", ".join(failed)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=PACKAGE_ROOT / "results" / "validation" / "synthetic_gaussian_validation.json",
        help="Machine-readable JSON output path.",
    )
    args = parser.parse_args()
    result = run_validation()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
