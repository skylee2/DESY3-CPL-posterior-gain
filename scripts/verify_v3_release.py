#!/usr/bin/env python3
"""Verify every manuscript-v3-facing result archived in this release."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CHECKS: list[dict[str, object]] = []


def record(name: str, passed: bool, actual: object, expected: object, tolerance: float | None = None) -> None:
    CHECKS.append(
        {"name": name, "passed": bool(passed), "actual": actual, "expected": expected, "tolerance": tolerance}
    )


def scalar(name: str, actual: float, expected: float, tolerance: float) -> None:
    record(name, abs(float(actual) - expected) <= tolerance, float(actual), expected, tolerance)


def load(relative: str) -> dict[str, object]:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def angle(a: list[float], b: list[float]) -> float:
    aa, bb = np.asarray(a), np.asarray(b)
    cosine = abs(float(aa @ bb)) / (np.linalg.norm(aa) * np.linalg.norm(bb))
    return float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))


def check_point_results() -> None:
    expected = {
        "bs_brs_point.json": (3.1499, 1.0067, 3.1709, 1.0062, 4.632, 0.25),
        "br_brs_point.json": (6.7343, 1.4972, 10.0828, 1.1721, 3.345, 1.382),
        "d3_d3brs_point.json": (37.2203, 7.06185, 262.8443, 1.32560, 50.82, 0.699),
    }
    # Lambda2 is printed to four decimals for BS and BR and five for D3.
    # The common tolerance therefore follows the least precise manuscript entry.
    tolerances = (5e-5, 5e-5, 5e-5, 5e-5, 5e-3, 5e-3)
    keys = ("lambda1", "lambda2", "det_ratio", "r_gain_algebraic", "delta_phi_weak_deg", "D_ref")
    for filename, values in expected.items():
        point = load(f"results/points/{filename}")
        for key, target, tolerance in zip(keys, values, tolerances):
            scalar(f"{filename}:{key}", point[key], target, tolerance)
        scalar(f"{filename}:determinant_closure", point["lambda_product"] - point["det_ratio"], 0.0, 2e-12)
    d3 = load("results/points/d3_d3brs_point.json")
    scalar("D3:v2_weak_axis_angle_deg", angle(d3["v2"], d3["e_weak_ref"]), 2.23, 5e-3)


def check_bs_bootstrap() -> None:
    with (ROOT / "results/bootstrap/bs_brs_bootstrap_v072.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    arrays = {key: np.asarray([float(row[key]) for row in rows]) for key in (
        "Lambda1", "Lambda2", "r_gain_algebraic", "delta_phi_weak_deg", "delta_phi_strong_deg"
    )}
    q1 = np.percentile(arrays["Lambda1"], [16, 50, 84])
    q2 = np.percentile(arrays["Lambda2"], [16, 50, 84])
    mask = arrays["Lambda2"] >= 1.0
    qr = np.percentile(arrays["r_gain_algebraic"][mask], [16, 50, 84])
    qa = np.percentile(arrays["delta_phi_weak_deg"], [16, 50, 84])
    for name, measured, target, tolerance in (
        ("Lambda1_median", q1[1], 3.142, 5e-4), ("Lambda1_minus", q1[1]-q1[0], 0.101, 5e-4),
        ("Lambda1_plus", q1[2]-q1[1], 0.112, 5e-4), ("Lambda2_median", q2[1], 1.006, 5e-4),
        ("Lambda2_minus", q2[1]-q2[0], 0.033, 5e-4), ("Lambda2_plus", q2[2]-q2[1], 0.034, 5e-4),
        ("r_gain_median", qr[1], 1.023, 5e-4), ("r_gain_minus", qr[1]-qr[0], 0.015, 5e-4),
        ("r_gain_plus", qr[2]-qr[1], 0.025, 5e-4), ("angle_median", qa[1], 4.65, 5e-3),
        ("angle_minus", qa[1]-qa[0], 0.26, 5e-3), ("angle_plus", qa[2]-qa[1], 0.27, 5e-3),
    ):
        scalar(f"BS_bootstrap:{name}", measured, target, tolerance)
    scalar("BS_bootstrap:P_Lambda2_ge_1_percent", 100*np.mean(mask), 56.9, 0.05)
    scalar("BS_bootstrap:P_angle_lt_5_percent", 100*np.mean(arrays["delta_phi_weak_deg"] < 5), 91.6, 0.05)
    scalar("BS_bootstrap:P_angle_lt_10_percent", 100*np.mean(arrays["delta_phi_weak_deg"] < 10), 100.0, 0.0)
    scalar("BS_bootstrap:P_angle_lt_45_percent", 100*np.mean(arrays["delta_phi_weak_deg"] < 45), 100.0, 0.0)
    scalar("BS_bootstrap:corr", np.corrcoef(arrays["Lambda1"], arrays["delta_phi_weak_deg"])[0, 1], -0.205, 5e-4)


def check_json_bootstrap(prefix: str, expected: dict[str, tuple[float, float, float]]) -> None:
    data = load(f"results/bootstrap/{prefix}.json")
    for field, (median, minus, plus) in expected.items():
        item = data["bootstrap_16_50_84"][field]
        tolerance = 5e-4 if field != "delta_phi_weak_deg" else 5e-3
        scalar(f"{prefix}:{field}_median", item["median"], median, tolerance)
        scalar(f"{prefix}:{field}_minus", item["minus"], minus, tolerance)
        scalar(f"{prefix}:{field}_plus", item["plus"], plus, tolerance)
    record(f"{prefix}:internal_checks", data["all_checks_passed"] is True, data["all_checks_passed"], True)


def check_other_bootstraps() -> None:
    check_json_bootstrap("br_brs_directional_bootstrap", {
        "Lambda1": (6.718, 0.205, 0.209), "Lambda2": (1.500, 0.053, 0.054),
        "r_gain_reported": (1.174, 0.019, 0.020), "delta_phi_weak_deg": (3.35, 0.16, 0.16),
    })
    check_json_bootstrap("d3_d3brs_directional_bootstrap", {
        "Lambda1": (37.341, 1.128, 1.174), "Lambda2": (7.036, 0.200, 0.221),
        "r_gain_reported": (1.323, 0.013, 0.017), "delta_phi_weak_deg": (50.84, 2.97, 3.40),
    })
    br = load("results/bootstrap/br_brs_directional_bootstrap.json")
    d3 = load("results/bootstrap/d3_d3brs_directional_bootstrap.json")
    scalar("BR_bootstrap:P_angle_lt_45_percent", br["angle_indicators"]["delta_phi_weak_lt_45deg"]["percent"], 100.0, 0.0)
    scalar("D3_bootstrap:P_angle_lt_45_percent", d3["angle_indicators"]["delta_phi_weak_lt_45deg"]["percent"], 4.0, 0.0)
    for prefix, data, expected in (("BR", br, (100.0, 100.0, 100.0)), ("D3", d3, (0.0, 0.0, 100.0))):
        scalar(f"{prefix}_bootstrap:P_angle_lt_5_percent",
               data["angle_indicators"]["delta_phi_weak_lt_5deg"]["percent"], expected[0], 0.0)
        scalar(f"{prefix}_bootstrap:P_angle_lt_10_percent",
               data["angle_indicators"]["delta_phi_weak_lt_10deg"]["percent"], expected[1], 0.0)
        scalar(f"{prefix}_bootstrap:P_Lambda2_ge_1_percent",
               100.0 * data["lambda2_ge_one"]["fraction"], expected[2], 0.0)
    scalar("BR_bootstrap:corr", br["corr_Lambda1_delta_phi_weak"], -0.102, 5e-4)
    scalar("D3_bootstrap:corr", d3["corr_Lambda1_delta_phi_weak"], 0.048, 5e-4)


def check_robustness() -> None:
    for prefix, targets in {
        "bs_brs_hard_prior_directional": [
            (0.96896, 0.999994, 2.9656, 1.0070, 1.0071, 4.707),
            (0.93834, 0.999978, 2.8356, 1.0093, 1.0101, 4.909),
            (0.87058, 0.998611, 2.5882, 1.0173, 1.0217, 5.292),
        ],
        "br_brs_hard_prior_directional": [
            (0.999783, 0.999994, 6.7216, 1.4975, 1.1726, 3.351),
            (0.999535, 0.999978, 6.7097, 1.4976, 1.1730, 3.356),
            (0.998664, 0.998611, 6.7201, 1.4988, 1.1731, 3.360),
        ],
    }.items():
        data = load(f"results/robustness/{prefix}.json")
        record(f"{prefix}:internal_checks", data["all_checks_passed"] is True, data["all_checks_passed"], True)
        keys = ("reference_retained_weight_fraction", "new_retained_weight_fraction",
                "Lambda1", "Lambda2", "r_gain", "delta_phi_weak_deg")
        for row, target in zip(data["results"][1:], targets):
            for key, value in zip(keys, target):
                # Table 3 prints retained weights to five or six decimal places.
                tolerance = 5e-6 if "weight_fraction" in key else (5e-4 if key == "delta_phi_weak_deg" else 5e-5)
                scalar(f"{prefix}:{row['selection']}:{key}", row[key], value, tolerance)
    d3 = load("results/robustness/d3_d3brs_robustness.json")
    d3_targets = {
        "full": (37.220, 7.062, 1.326), "w0+wa<-0.05": (37.059, 7.011, 1.324),
        "w0+wa<-0.1": (36.790, 6.927, 1.322), "w0+wa<-0.2": (36.491, 6.730, 1.315),
        "hpd_0.68": (35.009, 8.143, 1.402),
    }
    for label, target in d3_targets.items():
        row = d3[label]
        for key, value in zip(("lambda1", "lambda2", "r_gain_algebraic"), target):
            scalar(f"D3_robustness:{label}:{key}", row[key], value, 5e-4)


def check_additional_manuscript_values() -> None:
    br = load("results/points/br_brs_point.json")
    bs = load("results/points/bs_brs_point.json")
    d3 = load("results/points/d3_d3brs_point.json")
    for label, point, expected in (
        ("BR", br, (8045, 2964.7, -0.7175, -0.5942, 18.74, 83.2)),
        ("BS", bs, (10110, 2907.3, -0.9814, 0.2819, 33.41, 67.2)),
        ("BRS", bs, (11174, 3265.7, -0.9726, 0.1599, 59.49, 59.2)),
    ):
        prefix = "ref" if label != "BRS" else "new"
        scalar(f"Table1:{label}:samples", point[f"n_{prefix}"], expected[0], 0.0)
        scalar(f"Table1:{label}:N_eff", point[f"neff_{prefix}"], expected[1], 0.05)
        scalar(f"Table1:{label}:mean_w0", point[f"mu_{prefix}"][0], expected[2], 5e-5)
        scalar(f"Table1:{label}:mean_wa", point[f"mu_{prefix}"][1], expected[3], 5e-5)
        scalar(f"Table1:{label}:FoM", point[f"fom_{prefix}"], expected[4], 5e-3)
        scalar(f"Table1:{label}:kappa", point[f"kappa_{prefix}"], expected[5], 0.05)
    for label, point, which, expected in (
        ("BR", br, "ref", [[0.035489, -0.115679], [-0.115679, 0.457333]]),
        ("BS", bs, "ref", [[0.005445, -0.020746], [-0.020746, 0.243606]]),
        ("BRS", bs, "new", [[0.005405, -0.019983], [-0.019983, 0.126154]]),
        ("D3", d3, "ref", [[0.101275, -0.219956], [-0.219956, 1.120123]]),
        ("D3+BRS", d3, "new", [[0.006214, -0.026916], [-0.026916, 0.156413]]),
    ):
        record(f"covariance:{label}", np.allclose(point[f"cov_{which}"], expected, atol=5e-7, rtol=0),
               point[f"cov_{which}"], expected, 5e-7)
    scalar("D3:N_eff_ref", d3["neff_ref"], 5881.15, 5e-3)
    scalar("D3:N_eff_new", d3["neff_new"], 3114.64, 5e-3)
    scalar("D3:kappa_ref", d3["kappa_ref"], 20.88, 5e-3)
    scalar("D3:kappa_new", d3["kappa_new"], 104.84, 5e-3)
    ratio_ref = np.linalg.eigvalsh(np.asarray(d3["cov_ref"]))[0] / np.linalg.eigvalsh(np.asarray(d3["cov_ref"]))[1]
    ratio_new = np.linalg.eigvalsh(np.asarray(d3["cov_new"]))[0] / np.linalg.eigvalsh(np.asarray(d3["cov_new"]))[1]
    scalar("D3:ordinary_axis_ratio_ref", ratio_ref, 0.04789, 5e-6)
    scalar("D3:ordinary_axis_ratio_new", ratio_new, 0.00954, 5e-6)
    for name, actual, expected in (
        ("v1", d3["v1"], [-0.26614, 0.13856]),
        ("v2", d3["v2"], [0.17449, -1.04925]),
        ("e_major", d3["e_weak_ref"], [-0.20239, 0.97931]),
        ("e_minor", d3["e_strong_ref"], [-0.97931, -0.20239]),
    ):
        actual_array, expected_array = np.asarray(actual), np.asarray(expected)
        difference = min(np.max(np.abs(actual_array - expected_array)), np.max(np.abs(actual_array + expected_array)))
        scalar(f"D3:{name}_up_to_sign", difference, 0.0, 5e-6)
    scalar("D3:v1_strong_axis_angle", angle(d3["v1"], d3["e_strong_ref"]), 39.18, 5e-3)
    scalar("D3:v2_strong_axis_angle", angle(d3["v2"], d3["e_strong_ref"]), 87.77, 5e-3)
    scalar("D3:delta_mu_w0", d3["delta_mu"][0], -0.16661, 5e-6)
    scalar("D3:delta_mu_wa", d3["delta_mu"][1], 0.73261, 5e-6)
    scalar("BS:delta1", bs["delta1"], 2.1499, 5e-5)
    scalar("BS:delta2", bs["delta2"], 0.0067, 5e-5)
    scalar("BS:radial_width_factor", bs["lambda1"] ** -0.5, 0.56, 5e-3)
    scalar("BS:weak_axis_w0_up_to_sign", abs(bs["e_weak_ref"][0]), 0.086, 5e-4)
    scalar("BS:weak_axis_wa_up_to_sign", abs(bs["e_weak_ref"][1]), 0.996, 5e-4)
    bs_robustness = load("results/robustness/bs_brs_hard_prior_directional.json")
    br_robustness = load("results/robustness/br_brs_hard_prior_directional.json")
    scalar("BS:strongest_cut_retained_weight_percent",
           100.0 * bs_robustness["results"][-1]["reference_retained_weight_fraction"], 87.1, 0.05)
    scalar("BS:strongest_cut_angle_change_deg",
           bs_robustness["results"][-1]["delta_phi_weak_deg"] - bs_robustness["results"][0]["delta_phi_weak_deg"],
           0.660, 5e-4)
    record("BR:strongest_cut_retained_weight_above_99p8_percent",
           100.0 * br_robustness["results"][-1]["reference_retained_weight_fraction"] > 99.8,
           100.0 * br_robustness["results"][-1]["reference_retained_weight_fraction"], ">99.8")
    local = load("results/production/d3_stage1/local_ess_bootstrap.json")
    stage1 = load("results/production/d3_stage1/stage1_validation.json")
    local_mode = stage1["robustness"]["local mode: delta logP <= 1.15"]
    scalar("D3:local_N_eff_ref", local_mode["reference"]["effective_sample_size"], 2.8229, 5e-5)
    scalar("D3:local_N_eff_new", local_mode["combined"]["effective_sample_size"], 9.3317, 5e-5)
    scalar("D3:local_bootstrap_Lambda2_q025", local["lambda2_q025_q16_q50_q84_q975"][0], 0.031, 5e-4)
    scalar("D3:local_bootstrap_Lambda2_q975", local["lambda2_q025_q16_q50_q84_q975"][4], 1.825, 5e-4)
    scalar("D3:local_bootstrap_repetitions", local["repetitions"], 20000, 0.0)
    reverse = load("results/points/brs_d3brs_point.json")
    scalar("BRS_to_D3BRS:Lambda1", reverse["lambda1"], 1.45, 5e-3)
    scalar("BRS_to_D3BRS:Lambda2", reverse["lambda2"], 0.79, 5e-3)
    scalar("BRS_to_D3BRS:determinant_gain", reverse["det_ratio"], 1.14, 5e-3)


def check_validation_and_provenance() -> None:
    synthetic = load("results/validation/synthetic_gaussian_validation.json")
    targets = {
        "Lambda1": 9.990717157020178, "Lambda2": 4.997205887416219,
        "determinant_ratio": 49.925670596571464, "r_gain": 1.7424337629753992,
        "N_eff_ref": 884973.3591494124, "N_eff_new": 884685.4314052031,
    }
    for key, target in targets.items():
        scalar(f"synthetic:{key}", synthetic["estimated"][key], target, 1e-10)
    record("synthetic:analytic_targets", synthetic["analytic"]["generalized_eigenvalues"] == [10.0, 5.0]
           and synthetic["analytic"]["determinant_ratio"] == 50.0,
           synthetic["analytic"], {"Lambda": [10.0, 5.0], "determinant_ratio": 50.0})
    record("synthetic:12_of_12", synthetic["test_summary"] == {"passed": 12, "total": 12, "all_passed": True},
           synthetic["test_summary"], {"passed": 12, "total": 12, "all_passed": True})
    core = load("results/production/d3_stage1/core_validation.json")
    public = load("results/production/d3_stage1/public_d3_closure.json")
    scalar("public_D3:z_p", public["z_pivot"], 0.24435, 5e-6)
    scalar("public_D3:Q_0.05_w0", public["w0_q025_q05_q16_q50_q84_q95_q975"][1], -1.39847, 5e-6)
    scalar("public_D3:determinant_closure_relative_error", core["determinant_check"]["relative_error"], 1.1e-15, 5e-17)
    scalar("public_D3:boundary_weight_within_0p1_percent",
           100.0 * core["hard_prior_w0_plus_wa"]["reference_weight_within_0p1"], 0.506, 5e-4)
    figure1 = load("results/validation/figure1_metadata.json")
    record("Figure1:replacement_status", figure1["generator_status"].startswith("reconstructed replacement"),
           figure1["generator_status"], "reconstructed replacement; original historical plotting script not retained")
    scalar("Figure1:Lambda1", figure1["Lambda1"], 3.1499, 5e-5)
    scalar("Figure1:Lambda2", figure1["Lambda2"], 1.0067, 5e-5)
    scalar("Figure1:v1_angle", figure1["v1_weak_axis_angle_deg"], 4.632, 5e-4)
    figure = ROOT / figure1["output_figure"]
    record("Figure1:hash", sha256(figure) == figure1["figure_sha256"], sha256(figure), figure1["figure_sha256"])
    figure2 = load("results/validation/figure2_metadata.json")
    record("Figure2:replacement_status", figure2["generator_status"].startswith("reconstructed replacement"),
           figure2["generator_status"], "reconstructed replacement; original historical plotting script not retained")
    scalar("Figure2:Lambda1", figure2["Lambda1"], 6.7343, 5e-5)
    scalar("Figure2:Lambda2", figure2["Lambda2"], 1.4972, 5e-5)
    scalar("Figure2:v1_angle", figure2["v1_weak_axis_angle_deg"], 3.345, 5e-4)
    figure = ROOT / figure2["output_figure"]
    record("Figure2:hash", sha256(figure) == figure2["figure_sha256"], sha256(figure), figure2["figure_sha256"])
    metadata = load("results/validation/figure3_metadata.json")
    scalar("Figure3:Lambda1", metadata["Lambda1"], 37.2203, 5e-5)
    scalar("Figure3:Lambda2", metadata["Lambda2"], 7.06185, 5e-6)
    scalar("Figure3:v1_angle", metadata["v1_weak_axis_angle_deg"], 50.82, 5e-3)
    scalar("Figure3:v2_angle", metadata["v2_weak_axis_angle_deg"], 2.23, 5e-3)
    figure = ROOT / metadata["output_figure"]
    record("Figure3:hash", sha256(figure) == metadata["figure_sha256"], sha256(figure), metadata["figure_sha256"])
    expected_hashes = {
        "br": "499c6902b3f08239f78f2201ef4f7e912c695f0c3517af8ed1056f9b8bb28f54",
        "bs": "d315b9813808cf46bb624138eda41d29e478623ba381867db435d3fdae69d127",
        "brs": "d7ae346ceacd3c313fe7f731c63fb62cd79eda62f24be5303dd6b1ce2c47aea6",
        "pbrs": "700d723a7017906690e02f9fc58021ea294b110f425876a3e4588a9c75e1aeed",
    }
    checksum_text = (ROOT / "input/external_chain_SHA256SUMS.txt").read_text(encoding="utf-8")
    record("provenance:external_hashes", all(value in checksum_text for value in expected_hashes.values()),
           "all present" if all(value in checksum_text for value in expected_hashes.values()) else "missing", "all present")
    raw_chains = [path.name for path in ROOT.rglob("chain_*.txt")]
    record("provenance:no_raw_historical_chains", not raw_chains, raw_chains, [])
    d3_point = load("results/points/d3_d3brs_point.json")
    reverse = load("results/points/brs_d3brs_point.json")
    d3_robustness = load("results/robustness/d3_d3brs_robustness.json")
    expected_pair = ("d3_w0wa_nla_realy3dat.txt", "d3_brs_w0wa_nla_realy3dat.txt")
    identifiers = [(d3_point["ref_file"], d3_point["new_file"]),
                   (reverse["new_file"],)]
    identifiers.extend((row["ref_file"], row["new_file"]) for row in d3_robustness.values())
    portable = identifiers[0] == expected_pair and identifiers[1] == (expected_pair[1],)
    portable = portable and all(pair == expected_pair for pair in identifiers[2:])
    portable = portable and all("/" not in item and "\\" not in item for group in identifiers for item in group)
    record("provenance:portable_D3_input_identifiers", portable, identifiers,
           "portable basenames without parent-repository paths")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/validation/v3_release_audit.json")
    args = parser.parse_args()
    check_point_results()
    check_bs_bootstrap()
    check_other_bootstraps()
    check_robustness()
    check_additional_manuscript_values()
    check_validation_and_provenance()
    payload = {
        "schema_version": 1,
        "authoritative_manuscript": "paper1_Y3_submission_finaljournal_v3.tex",
        "authoritative_manuscript_sha256": "56919dc92f8466f3c2b23ffb8c860e49b8d833f97c73ed913e1c85893d79e340",
        "release_candidate": "1.1.0",
        "summary": {
            "passed": sum(bool(item["passed"]) for item in CHECKS),
            "total": len(CHECKS),
            "all_passed": all(bool(item["passed"]) for item in CHECKS),
        },
        "checks": CHECKS,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for item in CHECKS:
        print(f"[{'PASS' if item['passed'] else 'FAIL'}] {item['name']}")
    print(f"Checks passed: {payload['summary']['passed']}/{payload['summary']['total']}")
    if not payload["summary"]["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
