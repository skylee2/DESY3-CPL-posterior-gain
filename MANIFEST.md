# Package manifest — v1.1.0

This manifest describes release candidate v1.1.0 only.

The package contains 75 regular files. `SHA256SUMS` is the exhaustive file-level inventory for the other 74 files and excludes only itself.

- `README.md`: scope, quick validation, interpretation guardrails, and release status.
- `RELEASE_AUDIT_V3.md`: v1.0.0-to-v3 discrepancy table and actions.
- `docs/`: numerical provenance, chain provenance, validation details, and manuscript-to-output map.
- `input/`: external and public input checksum records; no raw chains.
- `scripts/`: clean reproduction, adopted production, validation, figure-generation, and release-audit code. `generate_figure1.py` and `generate_figure2.py` are reconstructed replacement generators; the original historical plotting scripts were not retained. `generate_figure3.py` is the corrected Figure 3 replacement generator.
- `results/points/`: regenerated point estimates.
- `results/bootstrap/`: regenerated transition-level bootstrap outputs.
- `results/robustness/`: regenerated hard-prior and HPD outputs.
- `results/validation/`: standalone and release-wide validation artifacts, including `figure1_metadata.json`, `figure2_metadata.json`, and `figure3_metadata.json`.
- `results/production/`: separated adopted production artifacts.
- `figures/`: manuscript-facing PDF figures.
- `SHA256SUMS`: file-level integrity manifest, excluding itself.

`LICENSE`, `NOTICE`, `CITATION.cff`, `CHANGELOG.md`, `VERSION`, `requirements.txt`, and `requirements-lock.txt` provide release metadata and environment requirements.
