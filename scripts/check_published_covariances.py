#!/usr/bin/env python3
"""
check_published_covariances.py

Independent deterministic check using only the covariance matrices printed in
paper1_Y3_submission_finaljournal_v3.tex. Because those matrices are rounded to 6 decimals,
the result is expected to agree with, but not reproduce every last digit of,
the full-chain values.

Expected from rounded C_BS and C_BRS:
  Lambda1 ~ 3.1500
  Lambda2 ~ 1.00678
  e_weak^BS ~ (-0.086, 0.996) up to sign
"""

import numpy as np
from scipy.linalg import eigh

C_BS = np.array([
    [ 0.005445, -0.020746],
    [-0.020746,  0.243606],
], dtype=float)

C_BRS = np.array([
    [ 0.005405, -0.019983],
    [-0.019983,  0.126154],
], dtype=float)

F_BS = np.linalg.inv(C_BS)
F_BRS = np.linalg.inv(C_BRS)

lam, V = eigh(F_BRS, F_BS)
idx = np.argsort(lam)[::-1]
lam = lam[idx]
V = V[:, idx]

cvals, cvecs = eigh(C_BS)
eweak = cvecs[:, -1]
if eweak[1] < 0:
    eweak *= -1
if np.dot(V[:, 0], eweak) < 0:
    V[:, 0] *= -1

def acute_angle(a, b):
    c = abs(np.dot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))

d1, d2 = lam - 1.0
rgain = (d1 + d2)**2 / (d1*d1 + d2*d2)

print("C_BS =")
print(C_BS)
print("\nC_BRS =")
print(C_BRS)
print("\nGeneralized eigenvalues (descending):", lam)
print("Lambda1 =", lam[0])
print("Lambda2 =", lam[1])
print("Lambda1*Lambda2 =", np.prod(lam))
print("det(C_BS)/det(C_BRS) =", np.linalg.det(C_BS)/np.linalg.det(C_BRS))
print("r_gain =", rgain)
print("e_weak^BS =", eweak)
print("v1 =", V[:, 0])
print("Delta phi_weak [deg] =", acute_angle(V[:, 0], eweak))
print("FoM_BS =", 1/np.sqrt(np.linalg.det(C_BS)))
print("FoM_BRS =", 1/np.sqrt(np.linalg.det(C_BRS)))
