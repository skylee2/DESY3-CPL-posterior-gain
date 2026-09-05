# DES Y3 CPL Posterior-Gain Reproducibility Package

## Purpose

This repository accompanies the paper “Growth-Sensitive Information Preferentially Tightens the Distance-Probe Weak Direction in CPL Dark Energy: A Posterior-Level Test with Historical DES Year 3 Chains.” It preserves the code, compact numerical outputs, figures, and provenance records used to audit the covariance-pair generalized-gain analysis. It does not contain the manuscript or the DES chain data.

## Scientific headline

For the historical external-only BS to BRS comparison, the dominant incremental precision gain is approximately $\Lambda_1=3.15$ and lies approximately $4.63$ degrees from the pre-existing BS weak covariance axis.

## Contents

- `scripts/`: the clean documented reproducer and independent printed-covariance validator.
- `scripts/production/`: unchanged original bootstrap production scripts retained as provenance objects.
- `results/`: the archived v07.2 bootstrap realizations and summary, plus the deterministic printed-covariance output.
- `figures/`: the three final figure PDFs copied byte-for-byte from the frozen submission package.
- `input/`: external-input instructions and verified filenames/SHA256 identities; no chain data.
- `docs/`: manuscript-to-code mapping, numerical and chain provenance, and validation details.
- `MANIFEST.md`, `SHA256SUMS`, `VERSION`, `CHANGELOG.md`, `CITATION.cff`, `LICENSE`, and `NOTICE`: release metadata and integrity records.

## External inputs

DES chains are not redistributed. The four historical external-only inputs expected by the production records are:

- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.br_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.bs_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.brs_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.pbrs_w0wa_realy3dat.txt`

Verify authorized local copies against [`input/external_chain_SHA256SUMS.txt`](input/external_chain_SHA256SUMS.txt). Public DES Y3 D3 and D3+BRS chains are likewise not bundled; supply their local paths explicitly to the clean reproducer.

## Quick validation

From the package root:

```bash
python scripts/check_published_covariances.py
```

This deterministic check requires no chain data and uses only the six-decimal covariance matrices printed in the paper.

## Full-chain reproduction

The canonical interface is `scripts/paper1_y3_reproduce.py`; no additional wrapper is needed. Inspect all commands with:

```bash
python scripts/paper1_y3_reproduce.py --help
```

For the central point estimate:

```bash
python scripts/paper1_y3_reproduce.py point \
  /path/to/des_y3_chains/chain_...bs_w0wa_realy3dat.txt \
  /path/to/des_y3_chains/chain_...brs_w0wa_realy3dat.txt \
  --label "BS -> BRS" \
  --output /path/to/output/bs_to_brs_point.json
```

The same `point` command accepts the BR/BRS pair and the public D3/D3+BRS pair. The script never downloads data and never assumes chains are included.

## Bootstrap reproduction

The clean implementation can run the declared multiplier bootstrap:

```bash
python scripts/paper1_y3_reproduce.py bootstrap \
  /path/to/des_y3_chains/chain_...bs_w0wa_realy3dat.txt \
  /path/to/des_y3_chains/chain_...brs_w0wa_realy3dat.txt \
  --nboot 1000 --seed 20260904 \
  --csv /path/to/output/bs_brs_bootstrap.csv \
  --output /path/to/output/bs_brs_bootstrap.json
```

The unchanged script `scripts/production/bs_brs_bootstrap_v072.py` and its archived CSV/summary document the actual manuscript-production bootstrap. Exact percentile values depend on the RNG seed, stream, and retained-row convention, so the clean implementation is not represented as byte-identical to the production realization stream.

## Reproducibility map

See [`docs/REPRODUCIBILITY_MAP.md`](docs/REPRODUCIBILITY_MAP.md).

## Provenance

See [`docs/CHAIN_PROVENANCE.md`](docs/CHAIN_PROVENANCE.md) and [`docs/NUMERICAL_PROVENANCE.md`](docs/NUMERICAL_PROVENANCE.md). The raw chains are deliberately excluded.

## Software requirements

- Python 3.10 or newer
- NumPy
- SciPy

Install the declared minimum dependencies with `python -m pip install -r requirements.txt` in an isolated environment.

## Citation

Use `CITATION.cff` to cite this software package and cite the accompanying paper. The DOI is to be added after archival; a repository URL will be added once one exists.

## License

The MIT license applies to original author-created code and documentation within the scope explained in `NOTICE`. It does not license or redistribute DES chain data or relicense third-party survey products.

## Version

Version `1.0.0`, repository candidate dated 2026-09-05. This package has not been published or uploaded.
