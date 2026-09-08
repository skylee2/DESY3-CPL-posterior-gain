# Validation protocol

The release candidate was validated in a fresh virtual environment with the dependencies in `requirements.txt`.

## Standalone synthetic validation

The validation constructs two diagonal two-dimensional Gaussian covariances from `F_ref = diag(100, 1)` and `F_new = diag(1000, 5)`, whose analytic generalized eigenvalues are exactly (10, 5). Independent Gaussian draws are converted to weighted samples using documented lognormal weights. The JSON output records the seed, bit generator, spawned stream identifiers, sample size, construction matrices, tolerances, effective sample sizes, and every pass/fail result.

The 12 executed checks are, exactly: (1) the two analytic covariance matrices are inverses of their precision matrices; (2) the analytic generalized eigenvalues are `(10, 5)`; (3) the reference weights sum to one; (4) the new weights sum to one; (5) both raw-weight samples are nonuniform; (6) both effective sample sizes lie strictly between one and the nominal sample size; (7) the reference weighted covariance meets its relative Frobenius-error tolerance; (8) the new weighted covariance meets the same tolerance; (9) both sample generalized eigenvalues meet their relative tolerance; (10) the production determinant identity closes; (11) the sample determinant ratio recovers the analytic target; and (12) the sample `r_gain` recovers the analytic target. This suite does not test generalized-eigenvector normalization, directional recovery, weak-axis recovery, or angle recovery. See `results/validation/synthetic_gaussian_validation.json` for the authoritative names, measurements, and tolerances.

## Release-wide audit

Run:

```bash
.venv/bin/python scripts/verify_v3_release.py
```

The resulting `results/validation/v3_release_audit.json` is machine-readable and must report `all_passed: true`. This is separate from each production script's internal checks.

## Figure validation

The original historical plotting scripts were not retained. `scripts/generate_figure1.py` and `scripts/generate_figure2.py` are reconstructed replacement generators, not original plotters. They consume only `results/points/bs_brs_point.json` and `results/points/br_brs_point.json`, respectively, and reproduce the validated scientific geometry of the final manuscript figures. Each records the packaged input, both generalized gains, the recomputed dominant weak-axis angle, its output path, and the PDF SHA-256 in `results/validation/figure1_metadata.json` or `results/validation/figure2_metadata.json`.

`scripts/generate_figure3.py` consumes only `results/points/d3_d3brs_point.json`. It recomputes both weak-axis angles and writes them, the two eigenvalues, the intended interpretation, and the PDF SHA-256 to `results/validation/figure3_metadata.json`.

## Integrity validation

`SHA256SUMS` covers every packaged file except itself. The sibling archive checksum covers the complete `.tar.gz`. Extracted-archive verification must reproduce the candidate's checksum pass.
