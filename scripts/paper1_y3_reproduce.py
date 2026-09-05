#!/usr/bin/env python3
"""
paper1_y3_reproduce.py

Clean reproducibility implementation for the CPL covariance-pair diagnostics
used in the DES Y3 manuscript.

It reads the actual CosmoSIS/PolyChord text chains, uses the stored *linear*
posterior weights, constructs the marginalized (w0, wa) weighted covariance,
and computes

    F_new v_i = Lambda_i F_ref v_i

with F = C^{-1}.

It also supports the weighted exponential-multiplier bootstrap used for the
BS -> BRS finite-chain stability test.

Important:
- Weighted covariance uses the population/second-moment convention in the
  manuscript: C = sum_i w_i (theta_i-mu)(theta_i-mu)^T, with sum_i w_i = 1.
  No Bessel correction is applied.
- scipy.linalg.eigh(F_new, F_ref) returns generalized eigenvectors normalized
  so v_i^T F_ref v_j = delta_ij.
- Euclidean alignment angles are descriptive quantities in native (w0, wa)
  coordinates.
- Exact bootstrap percentile values depend on the RNG seed/stream. If the
  historical production script used a different seed, point estimates will
  agree but bootstrap medians/percentiles can differ slightly.

Dependencies: numpy, scipy
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.linalg import eigh

W0_COL = "cosmological_parameters--w"
WA_COL = "cosmological_parameters--wa"
WEIGHT_COL = "weight"
POST_COL = "post"


@dataclass
class Chain2D:
    path: Path
    w0: np.ndarray
    wa: np.ndarray
    weight: np.ndarray
    post: Optional[np.ndarray] = None

    @property
    def theta(self) -> np.ndarray:
        return np.column_stack([self.w0, self.wa])


def _read_header_columns(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") and not line.startswith("##"):
                # First CosmoSIS output header line contains tab-separated names.
                text = line[1:].strip()
                if W0_COL in text and WA_COL in text and WEIGHT_COL in text:
                    return text.split()
    raise ValueError(f"Could not find chain column header in {path}")


def read_chain(
    path: str | Path,
    w0_col: str = W0_COL,
    wa_col: str = WA_COL,
    weight_col: str = WEIGHT_COL,
    post_col: str = POST_COL,
) -> Chain2D:
    path = Path(path)
    names = _read_header_columns(path)
    index = {name: i for i, name in enumerate(names)}

    for required in (w0_col, wa_col, weight_col):
        if required not in index:
            raise KeyError(f"{required!r} not found in {path.name}. Columns: {names}")

    data = np.loadtxt(path, comments="#", dtype=float, ndmin=2)

    w0 = data[:, index[w0_col]]
    wa = data[:, index[wa_col]]
    weight = data[:, index[weight_col]]
    post = data[:, index[post_col]] if post_col in index else None

    good = np.isfinite(w0) & np.isfinite(wa) & np.isfinite(weight) & (weight >= 0)
    if post is not None:
        post = post[good]
    w0, wa, weight = w0[good], wa[good], weight[good]

    if len(weight) == 0 or not np.isfinite(weight.sum()) or weight.sum() <= 0:
        raise ValueError(f"No valid positive total weight in {path}")

    return Chain2D(path=path, w0=w0, wa=wa, weight=weight, post=post)


def normalize_weight(weight: np.ndarray) -> np.ndarray:
    weight = np.asarray(weight, float)
    s = weight.sum()
    if s <= 0 or not np.isfinite(s):
        raise ValueError("Weights must have finite positive sum.")
    return weight / s


def weighted_mean_cov(theta: np.ndarray, weight: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Weighted population mean/covariance exactly as defined in the manuscript."""
    w = normalize_weight(weight)
    mu = np.sum(w[:, None] * theta, axis=0)
    dx = theta - mu
    cov = (dx * w[:, None]).T @ dx
    cov = 0.5 * (cov + cov.T)
    return mu, cov


def neff(weight: np.ndarray) -> float:
    w = normalize_weight(weight)
    return float(1.0 / np.sum(w * w))


def fom(cov: np.ndarray) -> float:
    return float(1.0 / np.sqrt(np.linalg.det(cov)))


def condition_number(cov: np.ndarray) -> float:
    return float(np.linalg.cond(cov))


def orient_axis(v: np.ndarray) -> np.ndarray:
    """Deterministic sign convention: positive wa component if possible."""
    v = np.asarray(v, float).copy()
    if v[1] < 0 or (v[1] == 0 and v[0] < 0):
        v *= -1
    return v


def ordinary_axes(cov: np.ndarray) -> Dict[str, np.ndarray | float]:
    vals, vecs = eigh(cov)  # ascending covariance eigenvalues
    strong = orient_axis(vecs[:, 0])  # small variance / short axis
    weak = orient_axis(vecs[:, -1])   # large variance / long axis
    return {
        "var_strong": float(vals[0]),
        "var_weak": float(vals[-1]),
        "e_strong": strong,
        "e_weak": weak,
    }


def generalized_gain(cov_ref: np.ndarray, cov_new: np.ndarray) -> Dict[str, np.ndarray | float]:
    Fref = np.linalg.inv(cov_ref)
    Fnew = np.linalg.inv(cov_new)

    # Symmetric-definite generalized eigenproblem:
    # Fnew v = lambda Fref v
    vals, vecs = eigh(Fnew, Fref)  # ascending, Fref-orthonormal
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    # Orient dominant vector toward reference weak axis for readable reporting.
    axes = ordinary_axes(cov_ref)
    weak = np.asarray(axes["e_weak"])
    if np.dot(vecs[:, 0], weak) < 0:
        vecs[:, 0] *= -1

    # The second sign is arbitrary; orient positive wa for reproducible text output.
    vecs[:, 1] = orient_axis(vecs[:, 1])

    lam1, lam2 = map(float, vals)
    d1, d2 = lam1 - 1.0, lam2 - 1.0
    rgain = float((d1 + d2) ** 2 / (d1 * d1 + d2 * d2)) if (d1*d1 + d2*d2) > 0 else np.nan

    return {
        "F_ref": Fref,
        "F_new": Fnew,
        "lambda1": lam1,
        "lambda2": lam2,
        "v1": vecs[:, 0],
        "v2": vecs[:, 1],
        "delta1": d1,
        "delta2": d2,
        "r_gain_algebraic": rgain,
        "det_ratio": float(np.linalg.det(cov_ref) / np.linalg.det(cov_new)),
    }


def acute_angle_deg(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    c = abs(float(np.dot(a, b))) / (np.linalg.norm(a) * np.linalg.norm(b))
    c = np.clip(c, -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def pair_diagnostics(ref: Chain2D, new: Chain2D) -> Dict:
    mu_ref, cov_ref = weighted_mean_cov(ref.theta, ref.weight)
    mu_new, cov_new = weighted_mean_cov(new.theta, new.weight)
    axes = ordinary_axes(cov_ref)
    gg = generalized_gain(cov_ref, cov_new)

    angle_weak = acute_angle_deg(gg["v1"], axes["e_weak"])
    angle_strong = acute_angle_deg(gg["v1"], axes["e_strong"])

    dmu = mu_new - mu_ref
    Dref = float(np.sqrt(dmu @ np.linalg.inv(cov_ref) @ dmu))

    return {
        "ref_file": str(ref.path),
        "new_file": str(new.path),
        "n_ref": int(len(ref.w0)),
        "n_new": int(len(new.w0)),
        "neff_ref": neff(ref.weight),
        "neff_new": neff(new.weight),
        "mu_ref": mu_ref.tolist(),
        "mu_new": mu_new.tolist(),
        "cov_ref": cov_ref.tolist(),
        "cov_new": cov_new.tolist(),
        "fom_ref": fom(cov_ref),
        "fom_new": fom(cov_new),
        "kappa_ref": condition_number(cov_ref),
        "kappa_new": condition_number(cov_new),
        "lambda1": gg["lambda1"],
        "lambda2": gg["lambda2"],
        "delta1": gg["delta1"],
        "delta2": gg["delta2"],
        "r_gain_algebraic": gg["r_gain_algebraic"],
        "det_ratio": gg["det_ratio"],
        "lambda_product": gg["lambda1"] * gg["lambda2"],
        "fom_ratio_squared": (fom(cov_new) / fom(cov_ref)) ** 2,
        "e_weak_ref": np.asarray(axes["e_weak"]).tolist(),
        "e_strong_ref": np.asarray(axes["e_strong"]).tolist(),
        "v1": np.asarray(gg["v1"]).tolist(),
        "v2": np.asarray(gg["v2"]).tolist(),
        "delta_phi_weak_deg": angle_weak,
        "delta_phi_strong_deg": angle_strong,
        "delta_mu": dmu.tolist(),
        "D_ref": Dref,
        "sigma_w0_ref": float(np.sqrt(cov_ref[0, 0])),
        "sigma_wa_ref": float(np.sqrt(cov_ref[1, 1])),
        "sigma_w0_new": float(np.sqrt(cov_new[0, 0])),
        "sigma_wa_new": float(np.sqrt(cov_new[1, 1])),
    }


def percentile_summary(x: np.ndarray) -> Dict[str, float]:
    q16, q50, q84 = np.percentile(np.asarray(x, float), [16, 50, 84])
    return {
        "median": float(q50),
        "minus": float(q50 - q16),
        "plus": float(q84 - q50),
        "q16": float(q16),
        "q84": float(q84),
    }


def multiplier_bootstrap(
    ref: Chain2D,
    new: Chain2D,
    nboot: int,
    seed: int,
    csv_path: Optional[Path] = None,
) -> Dict:
    rng = np.random.default_rng(seed)

    tr, tn = ref.theta, new.theta
    wr0, wn0 = np.asarray(ref.weight, float), np.asarray(new.weight, float)

    rows = []
    lam1s = np.empty(nboot)
    lam2s = np.empty(nboot)
    angles = np.empty(nboot)
    angle_strong = np.empty(nboot)
    rgains = np.empty(nboot)

    for b in range(nboot):
        wr = wr0 * rng.exponential(scale=1.0, size=len(wr0))
        wn = wn0 * rng.exponential(scale=1.0, size=len(wn0))

        _, Cr = weighted_mean_cov(tr, wr)
        _, Cn = weighted_mean_cov(tn, wn)

        axes = ordinary_axes(Cr)
        gg = generalized_gain(Cr, Cn)

        l1 = float(gg["lambda1"])
        l2 = float(gg["lambda2"])
        aw = acute_angle_deg(gg["v1"], axes["e_weak"])
        ast = acute_angle_deg(gg["v1"], axes["e_strong"])
        rg = float(gg["r_gain_algebraic"])

        lam1s[b], lam2s[b], angles[b], angle_strong[b], rgains[b] = l1, l2, aw, ast, rg
        rows.append((b, l1, l2, rg, aw, ast))

    if csv_path is not None:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["realization", "lambda1", "lambda2", "r_gain_algebraic",
                        "delta_phi_weak_deg", "delta_phi_strong_deg"])
            w.writerows(rows)

    keep = lam2s >= 1.0
    rg_cond = rgains[keep]

    return {
        "nboot": int(nboot),
        "seed": int(seed),
        "lambda1": percentile_summary(lam1s),
        "lambda2": percentile_summary(lam2s),
        "delta_phi_weak_deg": percentile_summary(angles),
        "p_lambda2_ge_1": float(np.mean(keep)),
        "p_angle_lt_5deg": float(np.mean(angles < 5.0)),
        "n_angle_ge_10deg": int(np.sum(angles >= 10.0)),
        "n_angle_ge_45deg": int(np.sum(angles >= 45.0)),
        "n_angle_lt_strong_angle": int(np.sum(angles < angle_strong)),
        "corr_lambda1_angle": float(np.corrcoef(lam1s, angles)[0, 1]),
        "r_gain_conditional_lambda2_ge_1": percentile_summary(rg_cond) if np.any(keep) else None,
        "n_r_gain_retained": int(np.sum(keep)),
    }


def subset_hard_boundary(chain: Chain2D, cut: float) -> Chain2D:
    mask = (chain.w0 + chain.wa) < cut
    post = chain.post[mask] if chain.post is not None else None
    return Chain2D(
        path=chain.path,
        w0=chain.w0[mask],
        wa=chain.wa[mask],
        weight=chain.weight[mask],
        post=post,
    )


def hpd_subset(chain: Chain2D, mass: float = 0.68) -> Chain2D:
    """
    Weighted HPD interior using the chain's 'post' column as log posterior density:
    sort by post descending and retain the smallest high-density set whose
    cumulative normalized posterior weight reaches 'mass'.
    """
    if chain.post is None:
        raise ValueError(f"{chain.path.name} has no '{POST_COL}' column; cannot construct HPD subset.")
    w = normalize_weight(chain.weight)
    order = np.argsort(chain.post)[::-1]
    cw = np.cumsum(w[order])
    k = int(np.searchsorted(cw, mass, side="left")) + 1
    idx = order[:k]
    return Chain2D(
        path=chain.path,
        w0=chain.w0[idx],
        wa=chain.wa[idx],
        weight=chain.weight[idx],
        post=chain.post[idx],
    )


def jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


def print_pair(label: str, d: Dict) -> None:
    print(f"\n=== {label} ===")
    print(f"N_eff(ref,new) = ({d['neff_ref']:.4f}, {d['neff_new']:.4f})")
    print(f"mu_ref = ({d['mu_ref'][0]:.6f}, {d['mu_ref'][1]:.6f})")
    print(f"mu_new = ({d['mu_new'][0]:.6f}, {d['mu_new'][1]:.6f})")
    print("C_ref =")
    print(np.array(d["cov_ref"]))
    print("C_new =")
    print(np.array(d["cov_new"]))
    print(f"FoM: {d['fom_ref']:.5f} -> {d['fom_new']:.5f}")
    print(f"kappa(C): {d['kappa_ref']:.4f} -> {d['kappa_new']:.4f}")
    print(f"(Lambda1,Lambda2) = ({d['lambda1']:.8f}, {d['lambda2']:.8f})")
    print(f"(delta1,delta2) = ({d['delta1']:.8f}, {d['delta2']:.8f})")
    print(f"r_gain(algebraic) = {d['r_gain_algebraic']:.8f}")
    print(f"Lambda1*Lambda2 = {d['lambda_product']:.10f}")
    print(f"det(Cref)/det(Cnew) = {d['det_ratio']:.10f}")
    print(f"(FoMnew/FoMref)^2 = {d['fom_ratio_squared']:.10f}")
    print(f"e_weak(ref) = {np.array(d['e_weak_ref'])}")
    print(f"v1 = {np.array(d['v1'])}")
    print(f"Delta phi_weak = {d['delta_phi_weak_deg']:.6f} deg")
    print(f"D_ref = {d['D_ref']:.6f}")


def cmd_point(args):
    ref = read_chain(args.ref)
    new = read_chain(args.new)
    d = pair_diagnostics(ref, new)
    print_pair(args.label, d)
    if args.output:
        Path(args.output).write_text(json.dumps(jsonable(d), indent=2), encoding="utf-8")


def cmd_bootstrap(args):
    ref = read_chain(args.ref)
    new = read_chain(args.new)
    point = pair_diagnostics(ref, new)
    print_pair(args.label + " point estimate", point)

    boot = multiplier_bootstrap(
        ref, new, nboot=args.nboot, seed=args.seed,
        csv_path=Path(args.csv) if args.csv else None
    )
    print("\n=== Bootstrap ===")
    print(json.dumps(boot, indent=2))
    if args.output:
        payload = {"point": point, "bootstrap": boot}
        Path(args.output).write_text(json.dumps(jsonable(payload), indent=2), encoding="utf-8")


def cmd_robustness(args):
    ref = read_chain(args.ref)
    new = read_chain(args.new)
    results = {"full": pair_diagnostics(ref, new)}

    for cut in (-0.05, -0.10, -0.20):
        rr = subset_hard_boundary(ref, cut)
        nn = subset_hard_boundary(new, cut)
        results[f"w0+wa<{cut}"] = pair_diagnostics(rr, nn)

    if args.hpd:
        rr = hpd_subset(ref, args.hpd_mass)
        nn = hpd_subset(new, args.hpd_mass)
        results[f"hpd_{args.hpd_mass:.2f}"] = pair_diagnostics(rr, nn)

    for key, d in results.items():
        print_pair(key, d)

    if args.output:
        Path(args.output).write_text(json.dumps(jsonable(results), indent=2), encoding="utf-8")


def build_parser():
    p = argparse.ArgumentParser(
        description="Reproduce generalized CPL covariance-pair diagnostics from DES Y3 chains."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("point", help="Point-estimate pair diagnostics.")
    q.add_argument("ref", help="Reference chain (e.g. BS).")
    q.add_argument("new", help="Updated chain (e.g. BRS).")
    q.add_argument("--label", default="pair")
    q.add_argument("--output", help="Write JSON.")
    q.set_defaults(func=cmd_point)

    q = sub.add_parser("bootstrap", help="Weighted exponential-multiplier bootstrap.")
    q.add_argument("ref")
    q.add_argument("new")
    q.add_argument("--label", default="BS -> BRS")
    q.add_argument("--nboot", type=int, default=1000)
    q.add_argument("--seed", type=int, default=20260904,
                   help="Reproducible RNG seed. Use the historical production seed if known.")
    q.add_argument("--csv", help="Write realization-level CSV.")
    q.add_argument("--output", help="Write point+bootstrap JSON.")
    q.set_defaults(func=cmd_bootstrap)

    q = sub.add_parser("robustness", help="Hard-boundary and optional weighted-HPD checks.")
    q.add_argument("ref")
    q.add_argument("new")
    q.add_argument("--hpd", action="store_true")
    q.add_argument("--hpd-mass", type=float, default=0.68)
    q.add_argument("--output", help="Write JSON.")
    q.set_defaults(func=cmd_robustness)

    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
