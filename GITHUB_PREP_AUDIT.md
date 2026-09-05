# GitHub Preparation Audit

## Repository candidate

`release/DESY3-CPL-posterior-gain/`, prepared locally from the audited v1.0.0 release candidate. No remote repository or release service was contacted.

## Version

`1.0.0`

## Public files included

The candidate includes the intended `scripts/`, `results/`, `figures/`, `input/`, and `docs/` trees, together with `README.md`, `CITATION.cff`, `LICENSE`, `NOTICE`, `CHANGELOG.md`, `MANIFEST.md`, `SHA256SUMS`, `VERSION`, `requirements.txt`, `.gitignore`, the prior release audit, and this preparation audit.

## Files deliberately excluded

Raw DES chains, private correspondence, email screenshots, local audit scratch files, temporary files, caches, bytecode, operating-system metadata, LaTeX auxiliary files, local environments, and the frozen manuscript are excluded.

## Raw DES chain exclusion

No `chain_*.txt` data file is present. The filename-and-SHA256 provenance records in `input/` remain included. The `.gitignore` explicitly rejects `chain_*.txt` without excluding the checksum metadata files or required archived CSV output.

## Privacy scan

A recursive text/binary scan found no author-home absolute path, private email address, private correspondence, credential, secret, private key, or email screenshot. File-type and filename checks found no unexpected binary or temporary artifact. No immutable production script contains an absolute local path, and the clean reproducer accepts explicit input and output paths.

## CITATION.cff validation

`CITATION.cff` parses as YAML and declares CFF 1.2.0, the required package title, Seokcheon Lee as author, and version 1.0.0. It contains no invented DOI, journal, volume, page, or repository URL. The manuscript title is retained as an otherwise unadorned preferred citation.

## README status

The README uses neutral repository-candidate wording, makes no publication claim, contains no private local path or nonexistent URL, and marks the DOI as to be added after archival.

## Deterministic validation

`python scripts/check_published_covariances.py` passed using the existing isolated release-validation Python environment. No bootstrap calculation was rerun. The clean reproducer help command, Python syntax checks, and archived JSON syntax checks also passed.

## Git status

All intended public files are to be staged for the single initial commit. The post-commit working tree must be clean, with no configured remote.

## Initial commit

One local root commit is to be created with subject `Initial v1.0.0 reproducibility package`. No tag is to be created.

## Remaining placeholders

The archival DOI and public repository URL remain unset. The DOI is explicitly marked for addition after archival; the repository URL can be added after remote creation.

## Final verdict

A. Ready to create remote GitHub repository
