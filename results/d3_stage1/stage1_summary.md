# STAGE 1 VALIDATION — NOT PRIMARY SCIENTIFIC RESULT

This is a real-public-chain validation of

`DES Y3 3x2pt -> DES Y3 3x2pt + BAO + RSD + Pantheon SN`.

It is not the planned geometry-to-LSS comparison: the increment is the mixed
BAO+RSD+SN block and the conditioning is reversed relative to the primary test.

## Full weighted-posterior result

| Quantity | DES Y3 3x2pt | DES Y3 3x2pt+BRS+SN |
|---|---:|---:|
| samples | 20,671 | 12,468 |
| Kish effective sample size | 5,881.15 | 3,114.64 |
| mean w0 | -0.805919 | -0.972531 |
| mean wa | -0.920856 | -0.188243 |
| ordinary information eigenvalues | (17.9156, 0.857942) | (650.819, 6.20766) |
| weak/strong information ratio | 0.0478879 | 0.00953823 |
| effective information dimension | 1.09556 | 1.01907 |
| covariance condition number | 20.8821 | 104.841 |

Covariances in `(w0, wa)` order:

```text
C_ref = [[ 0.101274817, -0.219956411],
         [-0.219956411,  1.120122854]]

C_new = [[ 0.006214314, -0.026916196],
         [-0.026916196,  0.156413452]]
```

Generalized information gains, ordered largest first:

```text
Lambda = (37.2203426, 7.06184552)
delta  = (36.2203426, 6.06184552)
r_gain = 1.32560060
```

The `F_ref`-normalized generalized vectors returned by `scipy.linalg.eigh` are

```text
v1 = (-0.266136565,  0.138555374)
v2 = ( 0.174488239, -1.049249857)
```

The reference ellipse major/minor vectors are

```text
major = (-0.202389629,  0.979305079)
minor = (-0.979305079, -0.202389629)
```

Using the documented acute Euclidean display-angle convention, the matrix
`[gain mode] x [reference major, reference minor]` is

```text
[[50.82099, 39.17901],
 [ 2.23496, 87.76504]] degrees.
```

Determinant closure is satisfied at machine precision:

```text
Lambda1*Lambda2               = 262.84430948415763
det(F_new)/det(F_ref)         = 262.84430948415790
relative difference           = 1.08e-15
```

Posterior-center shift:

```text
Delta mean (w0, wa)           = (-0.1666120, 0.7326125)
reference-Mahalanobis distance = 0.698623
```

## Public-result closure

The covariance-derived pivot redshifts are `0.24435` for DES Y3 3x2pt and
`0.20785` for DES Y3 3x2pt+BRS+SN, reproducing the published values 0.24 and
0.21. The 5th percentile of the DES-only w0 chain is `-1.39847`, reproducing
the published approximate 95% lower limit `w0 >= -1.4`.

## Robustness and Lambda2

| Estimator | retained reference weight | retained combined weight | Lambda1 | Lambda2 | r_gain |
|---|---:|---:|---:|---:|---:|
| full | 1.000000 | 1.000000 | 37.2203 | 7.06185 | 1.32560 |
| w0+wa < -0.05 | 0.998149 | ~1 | 37.0590 | 7.01055 | 1.32436 |
| w0+wa < -0.10 | 0.994937 | ~1 | 36.7898 | 6.92659 | 1.32235 |
| w0+wa < -0.20 | 0.985943 | ~1 | 36.4911 | 6.73008 | 1.31470 |
| 68% log-posterior HPD interior | 0.680013 | 0.680239 | 35.0086 | 8.14342 | 1.40234 |
| local delta-logP <= 1.15 | 0.0000639 | 0.001050 | 13.4419 | 0.86236 | 0.97788 |

The local-mode subsets have actual Kish effective sample sizes, computed after
renormalizing the retained weights, of only `2.82286` (3 D3 samples) and
`9.33170` (16 D3+BRS+SN samples). They are not estimated by multiplying the
retained-weight fraction by the global ESS. A 20,000-replicate conditional
multiplier bootstrap gives Lambda2 16/50/84 percentiles
`0.184/0.595/1.166`, with only 23.57% above one. Thus the reported local
`Lambda2=0.86236` is not statistically resolved and is not probative. Boundary cuts,
the 68% HPD estimate, and a 400-replicate multiplier bootstrap are stable.
The bootstrap 16/50/84 percentiles are:

```text
Lambda1 = 36.214 / 37.341 / 38.515
Lambda2 =  6.837 /  7.036 /  7.257
r_gain  =  1.310 /  1.323 /  1.340
```

All 400 bootstrap realizations had `Lambda2 > 1`. Thus `Lambda2 > 1` is
numerically stable for this validation comparison, but it does not establish
the paper's geometry-to-LSS conclusion.

The weighted KDE 68% contour orientations are `-79.97 deg` (reference) and
`-79.51 deg` (combined), compared with covariance-ellipse orientations
`-78.32 deg` and `-80.14 deg`.

The hard boundary is sampled closely (`max(w0+wa)=-0.00176` and `-0.000308`),
but only 0.506% of the reference weight and effectively none of the combined
weight lies within 0.1 of it. This does not remove the need for the planned
prior-boundary tests in the geometry-dominated Stage-2 baseline.
