#!/usr/bin/env python3
"""Reconstructed generator for manuscript Figure 1 (BS -> BRS geometry).

The original historical plotting script was not retained.  This replacement
uses only the validated, release-permitted point-result JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from scipy.stats import chi2


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def acute_angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    cosine = abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0))))


def covariance_ellipse(
    mean: np.ndarray,
    covariance: np.ndarray,
    probability: float,
    linestyle: str,
) -> Ellipse:
    values, vectors = np.linalg.eigh(covariance)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    scale = np.sqrt(chi2.ppf(probability, df=2))
    width, height = 2.0 * scale * np.sqrt(values)
    angle = float(np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0])))
    return Ellipse(
        mean,
        width,
        height,
        angle=angle,
        fill=False,
        color="black",
        linewidth=2.0,
        linestyle=linestyle,
    )


def direction_segment(center: np.ndarray, direction: np.ndarray, half_length: float) -> np.ndarray:
    unit = direction / np.linalg.norm(direction)
    return np.column_stack((center - half_length * unit, center + half_length * unit))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def generate(point_path: Path, output_path: Path, metadata_path: Path) -> dict[str, object]:
    point = json.loads(point_path.read_text(encoding="utf-8"))
    mu_ref = np.asarray(point["mu_ref"], dtype=float)
    mu_new = np.asarray(point["mu_new"], dtype=float)
    cov_ref = np.asarray(point["cov_ref"], dtype=float)
    cov_new = np.asarray(point["cov_new"], dtype=float)
    weak = np.asarray(point["e_weak_ref"], dtype=float)
    v1 = np.asarray(point["v1"], dtype=float)
    v2 = np.asarray(point["v2"], dtype=float)
    angle_v1 = acute_angle_deg(v1, weak)

    plt.rcParams.update({"font.size": 12, "axes.titlesize": 16, "axes.labelsize": 14})
    figure, axis = plt.subplots(figsize=(9.2, 7.6))
    for probability, ref_style, new_style in ((0.68, "-", "-."), (0.95, "--", ":")):
        axis.add_patch(covariance_ellipse(mu_ref, cov_ref, probability, ref_style))
        axis.add_patch(covariance_ellipse(mu_new, cov_new, probability, new_style))
    axis.scatter(*mu_ref, marker="o", s=55, color="tab:blue", zorder=5)
    axis.scatter(*mu_new, marker="s", s=55, color="tab:orange", zorder=5)

    weak_line = direction_segment(mu_ref, weak, 1.10)
    v1_line = direction_segment(mu_ref, v1, 1.10)
    v2_line = direction_segment(mu_ref, v2, 1.10)
    axis.plot(*weak_line, "--", color="tab:blue", linewidth=2.0)
    axis.plot(*v1_line, "-", color="tab:orange", linewidth=2.2)
    axis.plot(*v2_line, ":", color="tab:green", linewidth=2.0)

    handles = [
        Ellipse((0, 0), 1, 1, fill=False, color="black", linewidth=2.0, linestyle="-"),
        Ellipse((0, 0), 1, 1, fill=False, color="black", linewidth=2.0, linestyle="--"),
        Ellipse((0, 0), 1, 1, fill=False, color="black", linewidth=2.0, linestyle="-."),
        Ellipse((0, 0), 1, 1, fill=False, color="black", linewidth=2.0, linestyle=":"),
        plt.Line2D([], [], marker="o", linestyle="none", color="tab:blue", markersize=7),
        plt.Line2D([], [], marker="s", linestyle="none", color="tab:orange", markersize=7),
        plt.Line2D([], [], linestyle="--", color="tab:blue", linewidth=2.0),
        plt.Line2D([], [], linestyle="-", color="tab:orange", linewidth=2.2),
        plt.Line2D([], [], linestyle=":", color="tab:green", linewidth=2.0),
    ]
    labels = [
        "BS = BAO+Pantheon SN (68%)",
        "BS = BAO+Pantheon SN (95%)",
        "BRS = BAO+RSD+Pantheon SN (68%)",
        "BRS = BAO+RSD+Pantheon SN (95%)",
        "BS mean",
        "BRS mean",
        "BS weak covariance axis",
        rf"dominant $v_1$ ($\Lambda_1={point['lambda1']:.2f}$; $\Delta\phi_{{\rm weak}}={angle_v1:.2f}^\circ$)",
        rf"second $v_2$ ($\Lambda_2={point['lambda2']:.2f}$)",
    ]
    axis.legend(handles, labels, loc="lower left", fontsize=9.4, framealpha=0.90)
    axis.set_title(r"BS $\rightarrow$ BRS: RSD-induced generalized covariance gain")
    axis.set_xlabel(r"$w_0$")
    axis.set_ylabel(r"$w_a$")
    axis.set_xlim(-1.35, -0.55)
    axis.set_ylim(-1.35, 1.35)
    axis.grid(alpha=0.22)
    figure.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output_path,
        bbox_inches="tight",
        metadata={
            "Title": "BS to BRS generalized covariance geometry",
            "Author": "Seokcheon Lee",
            "Subject": (
                "Reconstructed replacement generator; original historical plotter not retained. "
                f"Lambda1={point['lambda1']:.8f}, Lambda2={point['lambda2']:.8f}, "
                f"dominant weak-axis angle={angle_v1:.8f} deg"
            ),
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(figure)

    metadata = {
        "schema_version": 1,
        "generator_status": "reconstructed replacement; original historical plotting script not retained",
        "source_point_result": str(point_path.relative_to(PACKAGE_ROOT)),
        "output_figure": str(output_path.relative_to(PACKAGE_ROOT)),
        "Lambda1": float(point["lambda1"]),
        "Lambda2": float(point["lambda2"]),
        "v1_weak_axis_angle_deg": angle_v1,
        "figure_sha256": sha256(output_path),
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--point", type=Path, default=PACKAGE_ROOT / "results" / "points" / "bs_brs_point.json"
    )
    parser.add_argument(
        "--output", type=Path, default=PACKAGE_ROOT / "figures" / "fig_BS_to_BRS_geometry_v07.pdf"
    )
    parser.add_argument(
        "--metadata", type=Path, default=PACKAGE_ROOT / "results" / "validation" / "figure1_metadata.json"
    )
    args = parser.parse_args()
    print(json.dumps(generate(args.point, args.output, args.metadata), indent=2))


if __name__ == "__main__":
    main()
