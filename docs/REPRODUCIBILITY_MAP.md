# Reproducibility map

| Manuscript-facing item | Script | Machine-readable output |
|---|---|---|
| BS -> BRS point values | `scripts/paper1_y3_reproduce.py point` | `results/points/bs_brs_point.json` |
| BR -> BRS point values | same | `results/points/br_brs_point.json` |
| D3 -> D3+BRS point values | same | `results/points/d3_d3brs_point.json` |
| Reverse BRS -> D3+BRS appendix check | same | `results/points/brs_d3brs_point.json` |
| BS bootstrap | `scripts/production/bs_brs_bootstrap_v072.py` | `results/bootstrap/bs_brs_bootstrap_v072.csv` and summary |
| BR bootstrap | `scripts/bootstrap_directional_diagnostics.py` | `results/bootstrap/br_brs_directional_bootstrap.{csv,json}` |
| D3 bootstrap | same | `results/bootstrap/d3_d3brs_directional_bootstrap.{csv,json}` |
| BS and BR hard-prior checks | `scripts/hard_prior_directional_robustness.py` | `results/robustness/*_hard_prior_directional.{csv,json}` |
| D3 hard-prior and HPD checks | `scripts/paper1_y3_reproduce.py robustness` | `results/robustness/d3_d3brs_robustness.json` |
| D3 public-chain closure | `scripts/production/d3_stage1/analysis/run_y3_stage1.py` | `results/production/d3_stage1/*.json` |
| Synthetic Gaussian validation | `scripts/validate_synthetic_gaussian.py` | `results/validation/synthetic_gaussian_validation.json` |
| Figure 1: BS -> BRS geometry | `scripts/generate_figure1.py` (reconstructed replacement; original historical plotter not retained) | `figures/fig_BS_to_BRS_geometry_v07.pdf` plus `results/validation/figure1_metadata.json` |
| Figure 2: BR -> BRS geometry | `scripts/generate_figure2.py` (reconstructed replacement; original historical plotter not retained) | `figures/fig1_BR_to_BRS_geometry_v06.pdf` plus `results/validation/figure2_metadata.json` |
| Figure 3 | `scripts/generate_figure3.py` | PDF plus `results/validation/figure3_metadata.json` |
| Complete v3 audit, including Tables 1--3 and supporting appendix values | `scripts/verify_v3_release.py` | `results/validation/v3_release_audit.json` |

The files under `results/production/` preserve adopted production artifacts. The parallel `points`, `bootstrap`, `robustness`, and `validation` directories hold cleanly regenerated, audit-facing products. Internal version strings such as `v072` identify the adopted method implementation; they are not manuscript-version claims.

Input filename fields in machine-readable results are portable basename identifiers, not paths relative to the release directory. Raw chains remain outside the package; callers supply their verified local locations as command-line arguments.

The Figure 1 and Figure 2 replacements consume only `results/points/bs_brs_point.json` and `results/points/br_brs_point.json`, respectively. They recompute the displayed native-coordinate angle from the packaged vectors, draw the 68% and 95% Gaussian-equivalent covariance ellipses from the packaged covariance matrices, and do not read or redistribute raw chain files. They reproduce the validated scientific geometry and manuscript-facing numerical labels, not the byte-level rendering history of the unavailable original plotters.
