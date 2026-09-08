# Numerical provenance for release v1.1.1 and manuscript v3

All values below were recomputed from verified chain inputs or from the documented synthetic construction. JSON and CSV files retain more digits than the manuscript.

## Point estimates

| Transition | Lambda1 | Lambda2 | determinant gain | r_gain | dominant-angle to weak axis | D_ref |
|---|---:|---:|---:|---:|---:|---:|
| BS -> BRS | 3.1499 | 1.0067 | 3.1709 | 1.0062 | 4.632 deg | approximately 0.25 |
| BR -> BRS | 6.7343 | 1.4972 | 10.0828 | 1.1721 | 3.345 deg | 1.382 |
| D3 -> D3+BRS | 37.2203 | 7.06185 | 262.8443 | 1.32560 | 50.82 deg | 0.699 |

For D3 -> D3+BRS, the subdominant `v2` angle to the weak D3 axis is 2.23 deg. The small `v2` angle must not be attributed to dominant `v1`.

## Weighted multiplier bootstrap

Intervals are medians with 16th--84th percentiles. BS and BR use 1000 realizations; D3 uses 400.

| Diagnostic | BS -> BRS | BR -> BRS | D3 -> D3+BRS |
|---|---:|---:|---:|
| Lambda1 | 3.142 +0.112/-0.101 | 6.718 +0.209/-0.205 | 37.341 +1.174/-1.128 |
| Lambda2 | 1.006 +0.034/-0.033 | 1.500 +0.054/-0.053 | 7.036 +0.221/-0.200 |
| r_gain | 1.023 +0.025/-0.015 | 1.174 +0.020/-0.019 | 1.323 +0.017/-0.013 |
| dominant weak-axis angle | 4.65 +0.27/-0.26 deg | 3.35 +0.16/-0.16 deg | 50.84 +3.40/-2.97 deg |
| P(Lambda2 >= 1) | 56.9% | 100% | 100% |
| P(angle < 5 deg) | 91.6% | 100% | 0% |
| P(angle < 10 deg) | 100% | 100% | 0% |
| P(angle < 45 deg) | 100% | 100% | 4.0% |

The BS `r_gain` interval is conditional on Lambda2 >= 1; the other two are unconditional. Seeds and exact quantiles are retained in the JSON summaries.

## Synthetic validation

`scripts/validate_synthetic_gaussian.py` uses seed 20260907, NumPy PCG64 with four spawned streams, and lognormal importance weights with sigma 0.35. The analytic generalized eigenvalues are (10, 5), the determinant target is 50, and analytic `r_gain` is 169/97.

Recovered values are:

- Lambda1 = 9.990717157020178
- Lambda2 = 4.997205887416219
- determinant ratio = 49.925670596571464
- r_gain = 1.7424337629753992
- N_eff_ref = 884973.3591494124
- N_eff_new = 884685.4314052031

All 12 explicitly named checks pass. This is the sole standalone synthetic validation introduced in v1.1.0 and retained byte-identically in v1.1.1; earlier undocumented validation attempts are superseded.

## Public-chain checks

- pivot redshift `z_p = 0.2443508229`, reported as 0.24435;
- `Q_0.05(w0) = -1.3984674949`, reported as -1.39847;
- D3 -> D3+BRS determinant-closure relative error `1.0813134774e-15`, reported as approximately `1.1e-15`.

Table 3 hard-prior and HPD values are stored at full precision in `results/robustness/` and checked row by row by `scripts/verify_v3_release.py`.

## Other v3 numerical checks

The release-wide verifier also checks the manuscript's chain row counts, effective sample sizes, means, FoMs, condition numbers, and all five printed covariance matrices. For the D3 transition it checks the printed generalized vectors up to their arbitrary signs, the ordinary-axis ratios 0.04789 and 0.00954, the complete 2x2 angle matrix, and the posterior-center displacement vector.

The diagnostic local effective sample sizes are 2.8229 and 9.3317. The 20,000-realization local diagnostic gives the printed 95% Lambda2 range 0.031--1.825. The reverse BRS -> D3+BRS check gives generalized gains approximately (1.45, 0.79) and determinant gain approximately 1.14. These supporting values are stored in the D3 production products and `results/points/brs_d3brs_point.json`.

## Externally sourced DES Year 6 context

The DES Year 6 quantities quoted in manuscript Section 5 are literature-only contextual values from the DES Year 6 cosmology analysis (`DES:2026jmi`, arXiv:2605.27221): FoMs 11, 48, 61, 110, 202, and 222, and reported departure significances 3.2 sigma and 3.0 sigma. They are not DES Year 3 chain outputs and are not presented as locally reproduced measurements.

Two manuscript values are transparent algebraic summaries of those literature inputs: `(48/11)^2` is approximately 19.0, and the change from 202 to 222 is approximately 10%. Their provenance is therefore “derived locally from externally sourced published scalars,” not “reproduced from packaged chains.” The published inputs are recorded in `scripts/production/d3_stage1/analysis/config.py`; the significance values remain citation-backed literature statements rather than inputs to this package's numerical pipeline.

Every internally derived manuscript-facing DES Year 3 quantity is stored at full precision in the machine-readable files mapped by `docs/REPRODUCIBILITY_MAP.md`. The release-wide audit JSON records direct checks of the manuscript-rounded values; quantities described narratively, such as the 0.506% D3 boundary weight, the 0.660-degree BS angle change, and the approximately 0.56 radial-width factor, are traceable respectively to `results/production/d3_stage1/core_validation.json`, `results/robustness/bs_brs_hard_prior_directional.json`, and `results/points/bs_brs_point.json`.
