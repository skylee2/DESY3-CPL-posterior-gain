#!/usr/bin/env python3
"""Hard-prior directional robustness for a weighted CPL chain pair.

The chain reader is shared with the clean main-analysis reproducer.  All
covariance, inverse-covariance, generalized-eigensystem, incremental-gain, and
reference-axis calculations call the preserved production analysis modules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
PRODUCTION_ROOT = SCRIPT_DIR / "production" / "d3_stage1"
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(PRODUCTION_ROOT))

from analysis.generalized_modes import (  # noqa: E402
    determinant_identity,
    generalized_eigensystem,
    incremental_gain_diagnostics,
)
from analysis.information_modes import fisher_from_covariance, ordinary_mode_metrics  # noqa: E402
from analysis.robustness import boundary_excluded_covariance, full_covariance  # noqa: E402
from paper1_y3_reproduce import read_chain  # noqa: E402


MARGINS = (0.05, 0.10, 0.20)
CSV_FIELDS = [
    "selection",
    "margin",
    "reference_retained_weight_fraction",
    "new_retained_weight_fraction",
    "Lambda1",
    "Lambda2",
    "r_gain",
    "weak_axis_w0",
    "weak_axis_wa",
    "delta_phi_weak_deg",
    "delta_phi_strong_deg",
    "determinant_closure_relative_error",
]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def acute_angle_degrees(first: np.ndarray, second: np.ndarray) -> float:
    a = np.asarray(first, dtype=float)
    b = np.asarray(second, dtype=float)
    cosine = abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))


def analyze_estimates(reference: object, new: object, *, selection: str, margin: float | None) -> dict[str, object]:
    f_ref = fisher_from_covariance(reference.covariance)
    f_new = fisher_from_covariance(new.covariance)
    gains, vectors = generalized_eigensystem(f_new, f_ref)
    gain_summary = incremental_gain_diagnostics(gains)
    reference_modes = ordinary_mode_metrics(reference.covariance)
    weak_axis = np.asarray(reference_modes["ellipse"]["major_vector"], dtype=float)
    strong_axis = np.asarray(reference_modes["ellipse"]["minor_vector"], dtype=float)
    v1 = np.asarray(vectors[:, 0], dtype=float)
    determinant = determinant_identity(f_new, f_ref, gains)
    angle_weak = acute_angle_degrees(v1, weak_axis)
    angle_strong = acute_angle_degrees(v1, strong_axis)
    return {
        "selection": selection,
        "margin": margin,
        "reference_retained_weight_fraction": float(reference.retained_weight),
        "new_retained_weight_fraction": float(new.retained_weight),
        "reference_retained_samples": int(reference.retained_samples),
        "new_retained_samples": int(new.retained_samples),
        "reference_effective_sample_size": float(reference.effective_sample_size),
        "new_effective_sample_size": float(new.effective_sample_size),
        "reference_covariance": np.asarray(reference.covariance).tolist(),
        "new_covariance": np.asarray(new.covariance).tolist(),
        "Lambda1": float(gains[0]),
        "Lambda2": float(gains[1]),
        "r_gain": float(gain_summary["r_gain"]),
        "weak_axis": weak_axis.tolist(),
        "strong_axis": strong_axis.tolist(),
        "dominant_generalized_direction_v1": v1.tolist(),
        "delta_phi_weak_deg": angle_weak,
        "delta_phi_strong_deg": angle_strong,
        "determinant_closure": determinant,
    }


def run(reference_path: Path, new_path: Path, *, label: str) -> dict[str, object]:
    reference_chain, new_chain = read_chain(reference_path), read_chain(new_path)
    results = [
        analyze_estimates(
            full_covariance(reference_chain.theta, reference_chain.weight),
            full_covariance(new_chain.theta, new_chain.weight),
            selection="full",
            margin=None,
        )
    ]
    for margin in MARGINS:
        results.append(
            analyze_estimates(
                boundary_excluded_covariance(reference_chain.theta, reference_chain.weight, margin=margin),
                boundary_excluded_covariance(new_chain.theta, new_chain.weight, margin=margin),
                selection=f"w0+wa < -{margin:.2f}",
                margin=margin,
            )
        )

    checks = {
        "all_retained_weight_fractions_in_unit_interval": bool(
            all(
                0.0 < row["reference_retained_weight_fraction"] <= 1.0
                and 0.0 < row["new_retained_weight_fraction"] <= 1.0
                for row in results
            )
        ),
        "retained_weight_decreases_monotonically": bool(
            all(
                results[i]["reference_retained_weight_fraction"] >= results[i + 1]["reference_retained_weight_fraction"]
                and results[i]["new_retained_weight_fraction"] >= results[i + 1]["new_retained_weight_fraction"]
                for i in range(len(results) - 1)
            )
        ),
        "all_generalized_eigenvalues_positive": bool(
            all(row["Lambda1"] > 0.0 and row["Lambda2"] > 0.0 for row in results)
        ),
        "all_determinant_closures_pass": bool(
            all(row["determinant_closure"]["passes"] for row in results)
        ),
        "all_angles_acute": bool(
            all(
                0.0 <= row["delta_phi_weak_deg"] <= 90.0
                and 0.0 <= row["delta_phi_strong_deg"] <= 90.0
                for row in results
            )
        ),
        "weak_strong_angles_are_complementary": bool(
            all(abs(row["delta_phi_weak_deg"] + row["delta_phi_strong_deg"] - 90.0) < 1e-8 for row in results)
        ),
    }
    payload = {
        "schema_version": 1,
        "label": label,
        "inputs": {
            "reference_file": reference_path.name,
            "new_file": new_path.name,
            "reference_sha256": file_sha256(reference_path),
            "new_sha256": file_sha256(new_path),
        },
        "cuts": ["w0+wa < -0.05", "w0+wa < -0.10", "w0+wa < -0.20"],
        "production_pipeline": [
            "analysis.robustness.boundary_excluded_covariance",
            "analysis.weighted_stats.weighted_covariance (called by boundary_excluded_covariance)",
            "analysis.information_modes.fisher_from_covariance",
            "analysis.generalized_modes.generalized_eigensystem",
            "analysis.generalized_modes.determinant_identity",
            "analysis.generalized_modes.incremental_gain_diagnostics",
            "analysis.information_modes.ordinary_mode_metrics",
        ],
        "angle_definition": "acute Euclidean angle in native (w0,wa) coordinates between v1 and the major axis of the selected reference covariance",
        "results": results,
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
    }
    assert payload["all_checks_passed"], f"robustness checks failed for {label}: {checks}"
    return payload


def write_outputs(prefix: Path, payload: dict[str, object]) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    csv_path = prefix.with_suffix(".csv")
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for item in payload["results"]:
            writer.writerow(
                {
                    "selection": item["selection"],
                    "margin": item["margin"],
                    "reference_retained_weight_fraction": item["reference_retained_weight_fraction"],
                    "new_retained_weight_fraction": item["new_retained_weight_fraction"],
                    "Lambda1": item["Lambda1"],
                    "Lambda2": item["Lambda2"],
                    "r_gain": item["r_gain"],
                    "weak_axis_w0": item["weak_axis"][0],
                    "weak_axis_wa": item["weak_axis"][1],
                    "delta_phi_weak_deg": item["delta_phi_weak_deg"],
                    "delta_phi_strong_deg": item["delta_phi_strong_deg"],
                    "determinant_closure_relative_error": item["determinant_closure"]["relative_error"],
                }
            )
    print(f"{payload['label']} hard-prior directional robustness")
    for item in payload["results"]:
        print(
            f"{item['selection']:15s}  weights=({item['reference_retained_weight_fraction']:.9f},"
            f" {item['new_retained_weight_fraction']:.9f})  Lambda=({item['Lambda1']:.8f},"
            f" {item['Lambda2']:.8f})  r_gain={item['r_gain']:.8f}  "
            f"Delta_phi_weak={item['delta_phi_weak_deg']:.6f} deg"
        )
    print(f"All checks passed: {payload['all_checks_passed']}")
    print(f"Wrote {json_path} and {csv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    write_outputs(args.output_prefix, run(args.reference, args.new, label=args.label))


if __name__ == "__main__":
    main()
