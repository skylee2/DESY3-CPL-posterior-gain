#!/usr/bin/env python3
"""Expanded weighted multiplier bootstrap for Paper 1 v07.2.

For every BS -> BRS realization recompute:
  Lambda_1, Lambda_2, algebraic r_gain,
  weak/strong-axis alignment angles, and the large-Lambda/weak-axis pairing.

The manuscript interprets r_gain as an effective rank only when Lambda_2 >= 1.
Accordingly, the reported r_gain interval is conditional on that subset, while
P(Lambda_2 >= 1) is reported separately.
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
    req = ["cosmological_parameters--w", "cosmological_parameters--wa", "weight"]
    miss = [c for c in req if c not in idx]
    if miss:
        raise ValueError(f"Missing columns in {path.name}: {miss}")
    data = np.loadtxt(path, comments="#")
    x = np.column_stack([data[:, idx[req[0]]], data[:, idx[req[1]]]])
    w = data[:, idx[req[2]]].astype(float)
    m = np.isfinite(x).all(axis=1) & np.isfinite(w) & (w > 0)
    x, w = x[m], w[m]
    w /= w.sum()
    return x, w

def wcov(x, w):
    w = np.asarray(w, float); w = w / w.sum()
    mu = np.sum(x*w[:,None], axis=0)
    d = x-mu
    return mu, (d*w[:,None]).T @ d

def gevp(Cref, Cnew):
    Fref, Fnew = np.linalg.inv(Cref), np.linalg.inv(Cnew)
    vals, vecs = eigh(Fnew, Fref)
    o = np.argsort(vals)[::-1]
    return vals[o], vecs[:,o]

def axes(C):
    vals, vecs = np.linalg.eigh(C)
    # largest covariance eigenvalue = weak/major axis; smallest = strong/minor
    return vecs[:, np.argmax(vals)], vecs[:, np.argmin(vals)]

def acute_deg(a,b):
    a=a/np.linalg.norm(a); b=b/np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(abs(float(a@b)),0,1))))

def r_gain(l1,l2):
    d1,d2=l1-1,l2-1
    den=d1*d1+d2*d2
    return np.nan if den == 0 else (d1+d2)**2/den

def diagnostics(Cref,Cnew):
    lam,V=gevp(Cref,Cnew)
    ew,es=axes(Cref)
    aw=acute_deg(V[:,0],ew); ast=acute_deg(V[:,0],es)
    return float(lam[0]),float(lam[1]),float(r_gain(lam[0],lam[1])),aw,ast

def q(a): return np.percentile(np.asarray(a,float),[16,50,84])
def fmt(qv,nd=3):
    lo,med,hi=qv
    return f"{med:.{nd}f}^{{+{hi-med:.{nd}f}}}_{{-{med-lo:.{nd}f}}}"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--bs',type=Path,default=Path(BS_DEFAULT))
    ap.add_argument('--brs',type=Path,default=Path(BRS_DEFAULT))
    ap.add_argument('--nboot',type=int,default=1000)
    ap.add_argument('--seed',type=int,default=20260904)
    ap.add_argument('--outdir',type=Path,default=Path('.'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    xbs,wbs=read_chain(args.bs); xbrs,wbrs=read_chain(args.brs)
    _,Cbs=wcov(xbs,wbs); _,Cbrs=wcov(xbrs,wbrs)
    point=diagnostics(Cbs,Cbrs)
    rng=np.random.default_rng(args.seed)
    out=np.empty((args.nboot,6),float)
    for b in range(args.nboot):
        wb=wbs*rng.exponential(size=wbs.size); wb/=wb.sum()
        wn=wbrs*rng.exponential(size=wbrs.size); wn/=wn.sum()
        _,Cb=wcov(xbs,wb); _,Cn=wcov(xbrs,wn)
        l1,l2,rg,aw,ast=diagnostics(Cb,Cn)
        out[b]=(b,l1,l2,rg,aw,ast)
    np.savetxt(args.outdir/'bs_brs_bootstrap_v072.csv',out,delimiter=',',
               header='rep,Lambda1,Lambda2,r_gain_algebraic,delta_phi_weak_deg,delta_phi_strong_deg',
               comments='',fmt=['%d','%.12g','%.12g','%.12g','%.12g','%.12g'])
    l1,l2,rg,aw,ast=out[:,1],out[:,2],out[:,3],out[:,4],out[:,5]
    valid=l2>=1.0
    ql1,ql2,qaw=q(l1),q(l2),q(aw)
    qrg=q(rg[valid]) if np.any(valid) else np.array([np.nan]*3)
    p2=100*np.mean(valid)
    p5=100*np.mean(aw<5); p10=100*np.mean(aw<10); p45=100*np.mean(aw<45)
    pcloser=100*np.mean(aw<ast)
    corr=float(np.corrcoef(l1,aw)[0,1])
    txt=f"""BS -> BRS expanded multiplier bootstrap (v07.2)\n
nboot = {args.nboot}\nseed = {args.seed}\nBS samples = {len(wbs)}\nBRS samples = {len(wbrs)}\n
Full-chain point values:\n  Lambda1 = {point[0]:.8f}\n  Lambda2 = {point[1]:.8f}\n  r_gain = {point[2]:.8f}\n  Delta phi_weak = {point[3]:.6f} deg\n  Delta phi_strong = {point[4]:.6f} deg\n
Bootstrap 16/50/84:\n  Lambda1 = {ql1[0]:.8f} / {ql1[1]:.8f} / {ql1[2]:.8f}\n  Lambda2 = {ql2[0]:.8f} / {ql2[1]:.8f} / {ql2[2]:.8f}\n  r_gain | Lambda2>=1 = {qrg[0]:.8f} / {qrg[1]:.8f} / {qrg[2]:.8f}\n  Delta phi_weak = {qaw[0]:.6f} / {qaw[1]:.6f} / {qaw[2]:.6f} deg\n
Fractions:\n  Lambda2 >= 1 = {p2:.2f}%\n  Delta phi_weak < 5 deg = {p5:.2f}%\n  Delta phi_weak < 10 deg = {p10:.2f}%\n  Delta phi_weak < 45 deg = {p45:.2f}%\n  Delta phi_weak < Delta phi_strong = {pcloser:.2f}%\n  Corr(Lambda1, Delta phi_weak) = {corr:.5f}\n"""
    (args.outdir/'bs_brs_bootstrap_summary_v072.txt').write_text(txt)
    tex=f"""% Auto-generated by bs_brs_bootstrap_v072.py\n% nboot={args.nboot}, seed={args.seed}\n\\newcommand{{\\BSBootLone}}{{{fmt(ql1,3)}}}\n\\newcommand{{\\BSBootLtwo}}{{{fmt(ql2,3)}}}\n\\newcommand{{\\BSBootRgainCond}}{{{fmt(qrg,3)}}}\n\\newcommand{{\\BSBootAngle}}{{{fmt(qaw,2)}^\\circ}}\n\\newcommand{{\\BSBootPtwo}}{{{p2:.1f}\\%}}\n\\newcommand{{\\BSBootPangleFive}}{{{p5:.1f}\\%}}\n\\newcommand{{\\BSBootPangleTen}}{{{p10:.1f}\\%}}\n\\newcommand{{\\BSBootPangleFortyFive}}{{{p45:.1f}\\%}}\n\\newcommand{{\\BSBootPweakCloser}}{{{pcloser:.1f}\\%}}\n\\newcommand{{\\BSBootCorrLoneAngle}}{{{corr:.3f}}}\n\\newcommand{{\\BSPointAngle}}{{{point[3]:.3f}^\\circ}}\n"""
    (args.outdir/'bs_brs_bootstrap_summary_v072.tex').write_text(tex)
    print(txt)

if __name__=='__main__': main()
