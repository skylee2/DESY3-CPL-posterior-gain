# Release manifest

| Path | Purpose | Origin | Status |
|---|---|---|---|
| `.gitignore` | Excludes local artifacts, environments, and raw DES chain files | RELEASE METADATA | permanent |
| `README.md` | Reader-facing purpose, contents, inputs, and commands | DOCUMENTATION | permanent |
| `CITATION.cff` | Machine-readable software and preferred-paper citation | RELEASE METADATA | permanent |
| `LICENSE` | MIT license for original author-created software/documentation | RELEASE METADATA | permanent |
| `NOTICE` | License scope and external-data/product boundaries | RELEASE METADATA | permanent |
| `CHANGELOG.md` | v1.0.0 local staging history | RELEASE METADATA | permanent |
| `MANIFEST.md` | Complete permanent-file inventory | RELEASE METADATA | permanent |
| `RELEASE_AUDIT.md` | Construction, coverage, validation, and privacy audit | DOCUMENTATION | permanent |
| `GITHUB_PREP_AUDIT.md` | Local GitHub-repository preparation and validation audit | DOCUMENTATION | permanent |
| `SHA256SUMS` | SHA256 integrity list for all other permanent files | RELEASE METADATA | permanent; self-excluded |
| `VERSION` | Machine-readable release version | RELEASE METADATA | permanent |
| `requirements.txt` | Minimum dependencies used by the Python scripts | RELEASE METADATA | permanent |
| `scripts/paper1_y3_reproduce.py` | Canonical clean chain-level reproducer | AUTHOR CLEAN REPRODUCER | unchanged copy; syntax/help passed |
| `scripts/check_published_covariances.py` | Chain-free printed-covariance validator | AUTHOR CLEAN REPRODUCER | unchanged copy; syntax/runtime passed |
| `scripts/production/bs_brs_bootstrap_v071.py` | Historical v07.1 production bootstrap | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/bs_brs_bootstrap_v072.py` | Historical v07.2 production bootstrap used for archived results | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/__init__.py` | Package marker for preserved D3 production modules | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/config.py` | Chain configuration used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/generalized_modes.py` | Generalized-mode calculations used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/information_modes.py` | Covariance/FoM/ordinary-mode calculations used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/io_chains.py` | Weighted chain reader used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/posterior_shift.py` | Posterior-shift calculation used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/robustness.py` | Hard-boundary, HPD, local, and KDE estimators used by D3 driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/run_y3_stage1.py` | Original D3/D3+BRS validation-production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy; not rerun |
| `scripts/production/d3_stage1/analysis/validate_cpl.py` | Pair-analysis routine called by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `scripts/production/d3_stage1/analysis/weighted_stats.py` | Weighted statistics used by D3 production driver | ORIGINAL PRODUCTION FILE | immutable byte-identical copy |
| `results/bs_brs_bootstrap_v072.csv` | 1,000 v07.2 bootstrap realizations | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/bs_brs_bootstrap_summary_v072.txt` | v07.2 point/bootstrap summary | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/check_published_covariances_output.txt` | Deterministic rounded-covariance validation output | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/d3_stage1/core_validation.json` | Original D3/D3+BRS full-posterior point output | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/d3_stage1/stage1_validation.json` | Original D3 robustness/KDE/400-bootstrap output | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/d3_stage1/stage1_summary.md` | Reader-facing original D3 production summary | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/d3_stage1/public_d3_closure.json` | Public D3 chain-ingestion/ESS/covariance closure | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy |
| `results/d3_stage1/local_ess_bootstrap.json` | Archived 20,000-realization very-local diagnostic output | ARCHIVED NUMERICAL OUTPUT | immutable byte-identical copy; not rerun |
| `figures/fig_BS_to_BRS_geometry_v07.pdf` | Final manuscript Figure 1 | FINAL PAPER FIGURE | matches frozen submission by SHA256 |
| `figures/fig1_BR_to_BRS_geometry_v06.pdf` | Final manuscript Figure 2 | FINAL PAPER FIGURE | matches frozen submission by SHA256 |
| `figures/fig2_D3_to_D3BRS_geometry_v06.pdf` | Final manuscript Figure 3 | FINAL PAPER FIGURE | matches frozen submission by SHA256 |
| `input/README.md` | External-input acquisition/verification guidance | DOCUMENTATION | no data included |
| `input/external_chain_SHA256SUMS.txt` | Verified filenames and hashes for four excluded external chains | PROVENANCE RECORD | unchanged copy |
| `input/public_d3_chain_SHA256SUMS.txt` | Verified filenames and hashes for excluded public D3/D3+BRS chains | PROVENANCE RECORD | locally computed from unambiguous analysis inputs |
| `docs/REPRODUCIBILITY_MAP.md` | Manuscript Table/Figure/Appendix to code map | DOCUMENTATION | permanent |
| `docs/NUMERICAL_PROVENANCE.md` | Result-level scripts, outputs, values, and classifications | DOCUMENTATION | permanent |
| `docs/CHAIN_PROVENANCE.md` | Conservative historical-chain provenance and exclusion policy | DOCUMENTATION | permanent |
| `docs/VALIDATION.md` | Printed-covariance, closure, and production/clean validation | DOCUMENTATION | permanent |

Total permanent files: 44, including `SHA256SUMS` itself.
