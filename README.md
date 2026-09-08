# DESY3 CPL posterior-gain reproducibility package v1.1.1

This provenance synchronization patch reproduces the manuscript-facing numerical analysis for `paper1_Y3_submission_finaljournal_v3.tex` (SHA-256 `ef0a8f631010a9c1d632dc5d1e87934de3896ecface2c8294259d4444fda4271`). Its scientific payload is byte-identical to public release v1.1.0; only version, manuscript-provenance, audit, and integrity metadata are updated. It does not replace or modify v1.1.0 or v1.0.0.

The package contains analysis code, machine-readable outputs, figures, provenance records, and integrity hashes. It intentionally contains no raw DES chain files.

## Scientific scope

- `results/points/`: point estimates for BS -> BRS, BR -> BRS, and D3 -> D3+BRS.
- `results/bootstrap/`: current directional bootstrap realizations and summaries.
- `results/robustness/`: hard-prior and 68% HPD robustness outputs used in manuscript Table 3.
- `results/validation/`: the standalone 12-check synthetic validation, metadata for all three figure generators, published-covariance check, and full v3 release audit.
- `results/production/`: preserved production outputs, separated from regenerated audit outputs.
- `figures/`: the three manuscript-facing covariance-geometry figures. The original historical plotting scripts were not retained. `scripts/generate_figure1.py` and `scripts/generate_figure2.py` are explicitly reconstructed replacement generators that reproduce the validated scientific geometry from packaged point-result JSON files; they are not represented as the original plotters. Figure 3 explicitly identifies dominant `v1` and subdominant `v2`.

The historical external-only BR, BS, BRS, and PBRS chains are identified by filename and SHA-256 in `input/external_chain_SHA256SUMS.txt`; they are not redistributed. The public D3 and D3+BRS chains are identified in `input/public_d3_chain_SHA256SUMS.txt`, but are also omitted from this compact package. Supply verified local copies to rerun chain-level calculations.

## Quick validation

Python 3.10 or newer is recommended.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/validate_synthetic_gaussian.py
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/desy3_mplconfig \
  .venv/bin/python scripts/generate_figure1.py
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/desy3_mplconfig \
  .venv/bin/python scripts/generate_figure2.py
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/desy3_mplconfig \
  .venv/bin/python scripts/generate_figure3.py
.venv/bin/python scripts/verify_v3_release.py
shasum -a 256 -c SHA256SUMS
```

The final verifier checks all manuscript-facing point values, bootstrap summaries, Table 3 robustness values, public-chain reproduction checks, the 12 synthetic checks, all three figure metadata/hash records, Figure 3 semantics, input provenance records, and the absence of redistributed historical chains.

## Reproducing from chains

Use `scripts/paper1_y3_reproduce.py` for point, bootstrap, and robustness calculations. The transition-specific scripts preserve the adopted production definitions:

```bash
.venv/bin/python scripts/paper1_y3_reproduce.py point REF_CHAIN NEW_CHAIN \
  --label BS_to_BRS --output results/points/bs_brs_point.json

.venv/bin/python scripts/bootstrap_directional_diagnostics.py \
  BR_CHAIN BRS_CHAIN --label BR_to_BRS --nboot 1000 --seed 20260906 \
  --output-prefix results/bootstrap/br_brs_directional_bootstrap

.venv/bin/python scripts/hard_prior_directional_robustness.py \
  BS_CHAIN BRS_CHAIN --label BS_to_BRS \
  --output-prefix results/robustness/bs_brs_hard_prior_directional
```

See `docs/NUMERICAL_PROVENANCE.md`, `docs/CHAIN_PROVENANCE.md`, `docs/REPRODUCIBILITY_MAP.md`, and `docs/ENVIRONMENT.md` for exact inputs, seeds, definitions, environment, and output mappings.

## Interpretation guardrails

Generalized covariance gain is a reference-dependent multiplicative precision comparison, not an additive Fisher contribution. FoM gain does not determine directional complementarity, and covariance contraction is distinct from posterior-center motion. D3 -> D3+BRS adds the full BRS block and is not an RSD-only test. Published DES Y6 FoM gains do not determine this directional decomposition.

For D3 -> D3+BRS, the dominant direction is `v1` with Lambda1 about 37.22 and a 50.82-degree angle from the weak D3 axis. It is not weak-axis aligned. The subdominant `v2`, with Lambda2 about 7.06, lies about 2.23 degrees from the weak D3 axis.

## Release status

Version 1.1.1 is a provenance synchronization patch for the final copy-edited journal-submission manuscript. No scientific results, numerical outputs, figures, tables, analysis scripts, bootstrap products, validation results, or conclusions changed relative to v1.1.0. Versioned releases are archived under the Zenodo concept DOI `10.5281/zenodo.22342577`; the published v1.1.0 and v1.0.0 releases remain unchanged.
