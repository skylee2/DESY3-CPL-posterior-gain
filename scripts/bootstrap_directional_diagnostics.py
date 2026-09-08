#!/usr/bin/env python3
"""Directional weighted-multiplier bootstrap for a CPL chain pair.

This extends the archived BS -> BRS diagnostic prescription to other
transitions without modifying the historical BS -> BRS production script or
its outputs.  It imports the clean reproducer's chain reader, weighted
population covariance, generalized-gain, reference-axis, and acute-angle
functions so that the definitions are identical.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from paper1_y3_reproduce import (  # noqa: E402
    acute_angle_deg,
    generalized_gain,
    ordinary_axes,
    pair_diagnostics,
    percentile_summary,
    read_chain,
    weighted_mean_cov,
)


CSV_FIELDS = [
    "realization",
    "Lambda1",
    "Lambda2",
    "r_gain",
    "weak_axis_w0",
    "weak_axis_wa",
    "delta_phi_weak_deg",
    "delta_phi_strong_deg",
    "angle_weak_lt_5deg",
    "angle_weak_lt_10deg",
    "angle_weak_lt_45deg",
    "angle_weak_lt_angle_strong",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(values: np.ndarray) -> dict[str, float]:
    return percentile_summary(np.asarray(values, dtype=float))


def run_bootstrap(reference_path: Path, new_path: Path, *, label: str, nboot: int, seed: int) -> tuple[dict, list[dict]]:
    reference, new = read_chain(reference_path), read_chain(new_path)
    point = pair_diagnostics(reference, new)
    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    for realization in range(nboot):
        reference_weights = reference.weight * rng.exponential(size=len(reference.weight))
        new_weights = new.weight * rng.exponential(size=len(new.weight))
        _, c_ref = weighted_mean_cov(reference.theta, reference_weights)
        _, c_new = weighted_mean_cov(new.theta, new_weights)
        reference_axes = ordinary_axes(c_ref)
        gain = generalized_gain(c_ref, c_new)
        weak_axis = np.asarray(reference_axes["e_weak"], dtype=float)
        strong_axis = np.asarray(reference_axes["e_strong"], dtype=float)
        v1 = np.asarray(gain["v1"], dtype=float)
        angle_weak = acute_angle_deg(v1, weak_axis)
        angle_strong = acute_angle_deg(v1, strong_axis)
        rows.append(
            {
                "realization": realization,
                "Lambda1": float(gain["lambda1"]),
                "Lambda2": float(gain["lambda2"]),
                "r_gain": float(gain["r_gain_algebraic"]),
                "weak_axis_w0": float(weak_axis[0]),
                "weak_axis_wa": float(weak_axis[1]),
                "delta_phi_weak_deg": angle_weak,
                "delta_phi_strong_deg": angle_strong,
                "angle_weak_lt_5deg": int(angle_weak < 5.0),
                "angle_weak_lt_10deg": int(angle_weak < 10.0),
                "angle_weak_lt_45deg": int(angle_weak < 45.0),
                "angle_weak_lt_angle_strong": int(angle_weak < angle_strong),
            }
        )

    arrays = {name: np.asarray([row[name] for row in rows]) for name in CSV_FIELDS[1:]}
    lambda2_ge_one = arrays["Lambda2"] >= 1.0
    all_lambda2_ge_one = bool(np.all(lambda2_ge_one))
    reported_r_gain = arrays["r_gain"] if all_lambda2_ge_one else arrays["r_gain"][lambda2_ge_one]

    indicators = {
        "delta_phi_weak_lt_5deg": arrays["angle_weak_lt_5deg"].astype(bool),
        "delta_phi_weak_lt_10deg": arrays["angle_weak_lt_10deg"].astype(bool),
        "delta_phi_weak_lt_45deg": arrays["angle_weak_lt_45deg"].astype(bool),
        "delta_phi_weak_lt_delta_phi_strong": arrays["angle_weak_lt_angle_strong"].astype(bool),
    }
    angle_counts = {
        name: {
            "count": int(np.count_nonzero(indicator)),
            "total": int(nboot),
            "fraction": float(np.mean(indicator)),
            "percent": float(100.0 * np.mean(indicator)),
        }
        for name, indicator in indicators.items()
    }

    checks = {
        "positive_generalized_eigenvalues": bool(np.all(arrays["Lambda2"] > 0.0)),
        "acute_angles_in_range": bool(
            np.all((arrays["delta_phi_weak_deg"] >= 0.0) & (arrays["delta_phi_weak_deg"] <= 90.0))
            and np.all((arrays["delta_phi_strong_deg"] >= 0.0) & (arrays["delta_phi_strong_deg"] <= 90.0))
        ),
        "orthogonal_axis_angle_complementarity": bool(
            np.allclose(
                arrays["delta_phi_weak_deg"] + arrays["delta_phi_strong_deg"],
                90.0,
                rtol=0.0,
                atol=1e-8,
            )
        ),
        "indicator_columns_match_angles": bool(
            np.array_equal(arrays["angle_weak_lt_5deg"].astype(bool), arrays["delta_phi_weak_deg"] < 5.0)
            and np.array_equal(arrays["angle_weak_lt_10deg"].astype(bool), arrays["delta_phi_weak_deg"] < 10.0)
            and np.array_equal(arrays["angle_weak_lt_45deg"].astype(bool), arrays["delta_phi_weak_deg"] < 45.0)
            and np.array_equal(
                arrays["angle_weak_lt_angle_strong"].astype(bool),
                arrays["delta_phi_weak_deg"] < arrays["delta_phi_strong_deg"],
            )
        ),
    }

    payload = {
        "schema_version": 1,
        "label": label,
        "method": "independent exponential weighted-multiplier bootstrap",
        "definitions": {
            "reference": "first input chain; its resampled covariance defines the weak/strong axes",
            "weights": "stored normalized posterior weights multiplied by iid Exponential(scale=1) variates separately for reference and new chains",
            "weak_axis": "oriented eigenvector of the largest eigenvalue of each resampled reference covariance",
            "angles": "acute Euclidean angles in native (w0, wa) coordinates",
            "r_gain": "((Lambda1-1)+(Lambda2-1))^2 / ((Lambda1-1)^2+(Lambda2-1)^2)",
        },
        "inputs": {
            "reference_file": reference_path.name,
            "new_file": new_path.name,
            "reference_sha256": sha256(reference_path),
            "new_sha256": sha256(new_path),
            "reference_retained_samples": int(len(reference.weight)),
            "new_retained_samples": int(len(new.weight)),
        },
        "nboot": int(nboot),
        "seed": int(seed),
        "point": {
            "Lambda1": point["lambda1"],
            "Lambda2": point["lambda2"],
            "r_gain": point["r_gain_algebraic"],
            "weak_axis": point["e_weak_ref"],
            "delta_phi_weak_deg": point["delta_phi_weak_deg"],
            "delta_phi_strong_deg": point["delta_phi_strong_deg"],
        },
        "bootstrap_16_50_84": {
            "Lambda1": summarize(arrays["Lambda1"]),
            "Lambda2": summarize(arrays["Lambda2"]),
            "r_gain_unconditional": summarize(arrays["r_gain"]),
            "r_gain_reported": summarize(reported_r_gain),
            "delta_phi_weak_deg": summarize(arrays["delta_phi_weak_deg"]),
            "delta_phi_strong_deg": summarize(arrays["delta_phi_strong_deg"]),
            "weak_axis_w0": summarize(arrays["weak_axis_w0"]),
            "weak_axis_wa": summarize(arrays["weak_axis_wa"]),
        },
        "lambda2_ge_one": {
            "count": int(np.count_nonzero(lambda2_ge_one)),
            "total": int(nboot),
            "fraction": float(np.mean(lambda2_ge_one)),
            "all_realizations": all_lambda2_ge_one,
            "r_gain_reporting": "unconditional" if all_lambda2_ge_one else "conditional_on_Lambda2_ge_1",
            "r_gain_reported_count": int(len(reported_r_gain)),
        },
        "angle_indicators": angle_counts,
        "corr_Lambda1_delta_phi_weak": float(
            np.corrcoef(arrays["Lambda1"], arrays["delta_phi_weak_deg"])[0, 1]
        ),
        "checks": checks,
        "all_checks_passed": bool(all(checks.values())),
    }
    assert payload["all_checks_passed"], f"internal diagnostic check failed for {label}: {checks}"
    return payload, rows


def write_outputs(prefix: Path, payload: dict, rows: list[dict]) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_suffix(".csv")
    json_path = prefix.with_suffix(".json")
    text_path = prefix.with_suffix(".txt")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    def triplet(name: str) -> str:
        item = payload["bootstrap_16_50_84"][name]
        return f"{item['q16']:.10f} / {item['median']:.10f} / {item['q84']:.10f}"

    lines = [
        f"{payload['label']} directional weighted-multiplier bootstrap",
        "",
        f"nboot = {payload['nboot']}",
        f"seed = {payload['seed']}",
        f"reference samples = {payload['inputs']['reference_retained_samples']}",
        f"new samples = {payload['inputs']['new_retained_samples']}",
        "",
        "Bootstrap 16/50/84:",
        f"  Lambda1 = {triplet('Lambda1')}",
        f"  Lambda2 = {triplet('Lambda2')}",
        f"  r_gain ({payload['lambda2_ge_one']['r_gain_reporting']}) = {triplet('r_gain_reported')}",
        f"  Delta phi_weak = {triplet('delta_phi_weak_deg')} deg",
        f"  Delta phi_strong = {triplet('delta_phi_strong_deg')} deg",
        "",
        "Fractions and counts:",
        f"  Lambda2 >= 1 = {payload['lambda2_ge_one']['fraction']*100:.2f}% ({payload['lambda2_ge_one']['count']}/{payload['nboot']})",
    ]
    for name, item in payload["angle_indicators"].items():
        lines.append(f"  {name} = {item['percent']:.2f}% ({item['count']}/{item['total']})")
    lines.extend(
        [
            f"  Corr(Lambda1, Delta phi_weak) = {payload['corr_Lambda1_delta_phi_weak']:.8f}",
            "",
            f"All internal checks passed = {payload['all_checks_passed']}",
        ]
    )
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Wrote {csv_path}, {json_path}, and {text_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("new", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--nboot", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    payload, rows = run_bootstrap(
        args.reference,
        args.new,
        label=args.label,
        nboot=args.nboot,
        seed=args.seed,
    )
    write_outputs(args.output_prefix, payload, rows)


if __name__ == "__main__":
    main()
