#!/usr/bin/env python3
"""Expanded weighted multiplier bootstrap for Paper 1 v07.1.

For every BS -> BRS bootstrap realization, recompute:
  Lambda_1, Lambda_2, r_gain, Delta phi_weak.

The weak-axis angle is recomputed against the *bootstrapped BS weak covariance axis*,
so the bootstrap tests the stability of the large-Lambda <-> weak-axis pairing itself.

Expected input files (CosmoSIS text chains):
  ...bs_w0wa_realy3dat.txt
  ...brs_w0wa_realy3dat.txt

Outputs:
  bs_brs_bootstrap_v071.csv
  bs_brs_bootstrap_summary_v071.txt
  bs_brs_bootstrap_summary_v071.tex

Place the generated .tex summary beside paper1_Y3_draft_v07_1.tex; the manuscript
will then automatically use the raw-chain bootstrap values instead of the fallback text.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

BS_DEFAULT = "chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.bs_w0wa_realy3dat.txt"
BRS_DEFAULT = "chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.brs_w0wa_realy3dat.txt"


def read_chain(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().lstrip("#")
    cols = header.split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    required = ["cosmological_parameters--w", "cosmological_parameters--wa", "weight"]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(f"Missing columns in {path.name}: {missing}")
    data = np.loadtxt(path, comments="#")
    x = np.column_stack([
        data[:, idx["cosmological_parameters--w"]],
        data[:, idx["cosmological_parameters--wa"]],
    ])
    w = data[:, idx["weight"]].astype(float)
    mask = np.isfinite(x).all(axis=1) & np.isfinite(w) & (w > 0)
    x, w = x[mask], w[mask]
    w /= w.sum()
    return x, w


def wcov(x: np.ndarray, w: np.ndarray):
    w = np.asarray(w, float)
    w = w / w.sum()
    mu = np.sum(x * w[:, None], axis=0)
    d = x - mu
    C = (d * w[:, None]).T @ d
    return mu, C


def gevp(C_ref: np.ndarray, C_new: np.ndarray):
    F_ref = np.linalg.inv(C_ref)
    F_new = np.linalg.inv(C_new)
    vals, vecs = eigh(F_new, F_ref)
    order = np.argsort(vals)[::-1]
    return vals[order], vecs[:, order]


def weak_axis(C: np.ndarray):
    vals, vecs = np.linalg.eigh(C)
    return vecs[:, np.argmax(vals)]


def acute_angle_deg(a: np.ndarray, b: np.ndarray):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    c = np.clip(abs(float(a @ b)), 0.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def r_gain(lam1: float, lam2: float):
    d1, d2 = lam1 - 1.0, lam2 - 1.0
    den = d1*d1 + d2*d2
    return np.nan if den == 0 else ((d1+d2)**2 / den)


def diagnostics(C_ref, C_new):
    lam, V = gevp(C_ref, C_new)
    wa = weak_axis(C_ref)
    angle = acute_angle_deg(V[:, 0], wa)
    rg = r_gain(lam[0], lam[1])
    return float(lam[0]), float(lam[1]), float(rg), float(angle)


def q16_50_84(a):
    return np.percentile(np.asarray(a, float), [16, 50, 84])


def fmt_interval(q, nd=3):
    lo, med, hi = q
    return f"{med:.{nd}f}^{{+{hi-med:.{nd}f}}}_{{-{med-lo:.{nd}f}}}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bs", type=Path, default=Path(BS_DEFAULT))
    ap.add_argument("--brs", type=Path, default=Path(BRS_DEFAULT))
    ap.add_argument("--nboot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260904)
    ap.add_argument("--outdir", type=Path, default=Path("."))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    x_bs, w_bs = read_chain(args.bs)
    x_brs, w_brs = read_chain(args.brs)
    mu_bs, C_bs = wcov(x_bs, w_bs)
    mu_brs, C_brs = wcov(x_brs, w_brs)
    point = diagnostics(C_bs, C_brs)

    rng = np.random.default_rng(args.seed)
    out = np.empty((args.nboot, 5), float)
    for b in range(args.nboot):
        m_bs = rng.exponential(scale=1.0, size=w_bs.size)
        m_brs = rng.exponential(scale=1.0, size=w_brs.size)
        wb_bs = w_bs * m_bs
        wb_brs = w_brs * m_brs
        wb_bs /= wb_bs.sum()
        wb_brs /= wb_brs.sum()
        _, Cb_bs = wcov(x_bs, wb_bs)
        _, Cb_brs = wcov(x_brs, wb_brs)
        l1, l2, rg, ang = diagnostics(Cb_bs, Cb_brs)
        out[b] = (b, l1, l2, rg, ang)

    names = ["rep", "Lambda1", "Lambda2", "r_gain", "delta_phi_weak_deg"]
    csv_path = args.outdir / "bs_brs_bootstrap_v071.csv"
    np.savetxt(csv_path, out, delimiter=",", header=",".join(names), comments="", fmt=["%d","%.12g","%.12g","%.12g","%.12g"])

    q_l1 = q16_50_84(out[:,1])
    q_l2 = q16_50_84(out[:,2])
    q_rg = q16_50_84(out[:,3])
    q_ang = q16_50_84(out[:,4])
    p_l2 = 100*np.mean(out[:,2] > 1.0)
    p_ang10 = 100*np.mean(out[:,4] < 10.0)
    p_ang5 = 100*np.mean(out[:,4] < 5.0)

    summary = f"""BS -> BRS expanded multiplier bootstrap (v07.1)\n
nboot = {args.nboot}\nseed = {args.seed}\nBS samples = {len(w_bs)}\nBRS samples = {len(w_brs)}\n
Full-chain point values:\n  Lambda1 = {point[0]:.8f}\n  Lambda2 = {point[1]:.8f}\n  r_gain  = {point[2]:.8f}\n  Delta phi_weak = {point[3]:.6f} deg\n
Bootstrap 16/50/84 percentiles:\n  Lambda1 = {q_l1[0]:.8f} / {q_l1[1]:.8f} / {q_l1[2]:.8f}\n  Lambda2 = {q_l2[0]:.8f} / {q_l2[1]:.8f} / {q_l2[2]:.8f}\n  r_gain  = {q_rg[0]:.8f} / {q_rg[1]:.8f} / {q_rg[2]:.8f}\n  Delta phi_weak [deg] = {q_ang[0]:.6f} / {q_ang[1]:.6f} / {q_ang[2]:.6f}\n
Fractions:\n  Lambda2 > 1 = {p_l2:.2f}%\n  Delta phi_weak < 5 deg = {p_ang5:.2f}%\n  Delta phi_weak < 10 deg = {p_ang10:.2f}%\n"""
    txt_path = args.outdir / "bs_brs_bootstrap_summary_v071.txt"
    txt_path.write_text(summary)

    # TeX macros: values are intended for math mode except the percentages.
    tex = f"""% Auto-generated by bs_brs_bootstrap_v071.py\n% nboot={args.nboot}, seed={args.seed}\n\\newcommand{{\\BSBootLone}}{{{fmt_interval(q_l1,3)}}}\n\\newcommand{{\\BSBootLtwo}}{{{fmt_interval(q_l2,3)}}}\n\\newcommand{{\\BSBootRgain}}{{{fmt_interval(q_rg,3)}}}\n\\newcommand{{\\BSBootAngle}}{{{fmt_interval(q_ang,2)}^\\circ}}\n\\newcommand{{\\BSBootPtwo}}{{{p_l2:.1f}\\%}}\n\\newcommand{{\\BSBootPangleTen}}{{{p_ang10:.1f}\\%}}\n\\newcommand{{\\BSPointAngle}}{{{point[3]:.3f}^\\circ}}\n"""
    tex_path = args.outdir / "bs_brs_bootstrap_summary_v071.tex"
    tex_path.write_text(tex)

    print(summary)
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {txt_path}")
    print(f"Wrote: {tex_path}")

if __name__ == "__main__":
    main()
