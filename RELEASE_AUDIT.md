# DES Y3 CPL Posterior-Gain v1.0.0 Release Audit

## 1. Source directories

- Frozen journal package (read-only source for final figures): `paper1_Y3_submission_package/`
- Reproducibility artifacts (read-only source for code/results): `paper1_Y3_reproducibility/`
- Historical source package (read-only duplicate/provenance inspection only): `paper1_Y3_draft_v07_source_package/`
- Verified external checksum record: project-root `external_chain_SHA256SUMS.txt`
- Project production code searched read-only: DESY6 `analysis/`, `notebooks/`, `runs/`, and local Git-visible history
- Archived D3 outputs searched read-only: `results/y3_validation_d3_to_d3_brs/` and `results/pre_stage2_validation/`
- Public D3 inputs hashed in place, not copied: `data/des_y3_public/`

The frozen manuscript/PDF, original scripts, original outputs, and any external chain data were not modified.

### Pre-copy source inventory

| Source file | Purpose | Paper result | Copy? | Reason |
|---|---|---|---|---|
| `paper1_Y3_reproducibility/paper1_y3_reproduce.py` | Clean documented chain-level calculations | point spectra, means/covariances, ESS/FoM, bootstrap, hard-boundary/HPD checks | yes | canonical public reproducer |
| `paper1_Y3_reproducibility/check_published_covariances.py` | Independent rounded-covariance calculation | BS/BRS gains, weak axis, angle, FoMs, closure | yes | chain-free validation utility |
| `paper1_Y3_reproducibility/check_published_covariances_output.txt` | Archived deterministic output | same | yes | validation checkpoint |
| `paper1_Y3_reproducibility/production/bs_brs_bootstrap_v071.py` | Historical earlier production bootstrap | BS/BRS bootstrap development provenance | yes | immutable production provenance |
| `paper1_Y3_reproducibility/production/bs_brs_bootstrap_v072.py` | Final archived production bootstrap | Table 2 and BS/BRS point values | yes | immutable production provenance |
| `paper1_Y3_reproducibility/production/bs_brs_bootstrap_v072.csv` | 1,000 realization output | Table 2 distributions/fractions/correlation | yes | archived numerical output |
| `paper1_Y3_reproducibility/production/bs_brs_bootstrap_summary_v072.txt` | Aggregate production summary | Table 2 and headline BS/BRS result | yes | archived numerical output |
| `paper1_Y3_reproducibility/README_DESY3_repro.md` | Internal precursor documentation | general usage/provenance notes | no | replaced by public staging README/docs; contains historical manuscript-version language |
| submission Figure 1 PDF | final BS/BRS geometry | Figure 1 | yes | authoritative frozen figure |
| submission Figure 2 PDF | final BR/BRS geometry | Figure 2 | yes | authoritative frozen figure |
| submission Figure 3 PDF | final D3/D3+BRS geometry | Figure 3 | yes | authoritative frozen figure |
| `external_chain_SHA256SUMS.txt` | exact external filename/hash identities | Appendix C/input provenance | yes | verified checksum record, contains no chain data |
| `analysis/run_y3_stage1.py` and nine directly required `analysis/` modules | original D3/D3+BRS production/validation pipeline | D3 point result, Table 3, D3 bootstrap and diagnostics | yes | exact production driver/dependencies established by matching archived outputs |
| `results/y3_validation_d3_to_d3_brs/core_validation.json` | D3/D3+BRS full-posterior result | D3 section | yes | original numerical production artifact |
| `results/y3_validation_d3_to_d3_brs/stage1_validation.json` | D3 robustness/KDE/bootstrap payload | Table 3 and D3 validation | yes | original numerical production artifact; not rerun |
| `results/y3_validation_d3_to_d3_brs/stage1_summary.md` | human-readable D3 production summary | D3 result/Table 3 | yes | archived summary matches manuscript values |
| `results/pre_stage2_validation/public_d3_closure.json` | public-chain ingestion and covariance closure | D3 ESS/covariance | yes | supporting archived output |
| `results/pre_stage2_validation/local_ess_bootstrap.json` | very-local conditional diagnostic | validation discussion | yes | supporting archived output; ensemble not rerun |
| `data/des_y3_public/d3_w0wa_nla_realy3dat.txt` and `d3_brs_w0wa_nla_realy3dat.txt` | exact public D3 analysis inputs | D3 result/Table 3 | checksum record only | raw chains remain external; paths and numerical identity match archived production JSON |
| four duplicate bootstrap files in `paper1_Y3_draft_v07_source_package/` | historical duplicate scripts/results | Table 2 | no | all four are byte-identical to selected reproducibility sources |
| historical package audits/README files | working-history notes | none | no | internal notes are outside public release scope |
| historical/frozen manuscript PDFs and TeX | manuscript sources/builds | paper | no | release describes the paper; frozen submission remains separate |
| any raw `chain_*.txt` | external DES posterior samples | chain-level results | no | redistribution prohibited |

No exact figure-generation script, complete historical Table 1 summary, BR/BRS original production output/bootstrap artifact, or Appendix B original production artifact was found. The clean reproducer regenerated the frozen Table 1 and BR/BRS point values, so these are historical-production provenance gaps rather than public reproducibility failures. D3 point/Table 3/400-bootstrap production code and outputs were located and archived.

## 2. Reference UDM release inspected

Reference inspected read-only: `UDM2502_Revised_Manuscript/release/UDMCLASS_patch_v1.0.0` under the Projects tree. Its release philosophy informed this staging layout: concise top-level README, clear separation of scripts/inputs/results/figures/docs, immutable provenance artifacts, explicit third-party data exclusions, machine-readable citation metadata, MIT license plus scope notice, version/changelog, complete manifest, and SHA256 coverage. No UDM scientific content, DOI, repository URL, Git metadata, or release artifact was copied.

## 3. Files included

The release contains 42 permanent files. The complete path/purpose/origin/status inventory is in `MANIFEST.md`. Newly added artifacts comprise the exact D3 Stage-1 production driver and directly imported analysis modules, five archived D3 numerical outputs, and one public-D3 two-chain checksum record.

## 4. Files deliberately excluded

The following raw DES historical external-only chains are explicitly excluded:

- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.br_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.bs_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.brs_w0wa_realy3dat.txt`
- `chain_2pt_NG_final_2ptunblind_02_26_21_wnz_maglim_covupdate.fits.scales-ml_3x2pt_8_6_0.5_v0.40.ini.pbrs_w0wa_realy3dat.txt`

Public D3/D3+BRS chain data are also excluded. Manuscript files, working audits, historical drafts, LaTeX auxiliaries, raw data, notebooks, caches, logs, private correspondence, screenshots, and duplicate historical copies are not included.

## 5. Production provenance

The selected v07.1 and v07.2 production scripts and the v07.2 CSV/summary are exact-byte copies. Source and release SHA256 values match:

| Artifact | SHA256 |
|---|---|
| `bs_brs_bootstrap_v071.py` | `ce0ceba96f7e0f369b81cd2241267aa97dc078440376c2ed76687e721af4a11d` |
| `bs_brs_bootstrap_v072.py` | `234159278c893db2f307e19dfff95158863ae820ce2971821f66e3e73e621df6` |
| `bs_brs_bootstrap_v072.csv` | `81bdc7e1d6ddde51eca9b279eec8179470ba802cee1083caf1a6ba01b1ffe8fb` |
| `bs_brs_bootstrap_summary_v072.txt` | `be834bedd3364121bbd3afdcb9d091c22d3e74e7545f5f26ff903150422fbaa4` |

The corresponding four files in the historical source package were also checked and are byte-identical duplicates; they were not copied twice.

The D3 production driver and its directly imported modules were copied unchanged under `scripts/production/d3_stage1/analysis/`. Their source/release SHA256 values are identical. Archived `core_validation.json`, `stage1_validation.json`, and `stage1_summary.md` match the manuscript D3 covariance, ESS, generalized gains, angles, hard-boundary/HPD values, and 400-bootstrap summary. The public D3 closure and very-local diagnostic outputs were also preserved. None of the stored bootstrap ensembles was rerun.

The author reports completed sign-off on the external-chain provenance: Jessica Muir confirmed analysis use/supply, release-label filename interpretation, matching common BS/BRS priors, and the documented $A_s,n_s$ treatment. The author also manually recomputed and confirmed all four external-chain hashes. This provenance item is RESOLVED.

## 6. Manuscript-to-code reproducibility coverage

| Item | Coverage |
|---|---|
| Table 1 external-chain summaries | CLEAN REPRODUCER COVERS RESULT; historical production summary not located (non-blocking) |
| Table 2 central BS/BRS bootstrap | COMPLETE |
| Table 3 D3 robustness | ORIGINAL PRODUCTION ARTIFACT AND SCRIPT FOUND; clean reproducer also covers result |
| Table 4 historical metadata | PROVENANCE RESOLVED; author sign-off complete |
| Figure 1 BS/BRS geometry | GENERATOR NOT LOCATED; numerical content reproducible and final PDF verified (non-blocking) |
| Figure 2 BR/BRS geometry | GENERATOR NOT LOCATED; numerical content reproducible and final PDF verified (non-blocking) |
| Figure 3 D3/D3+BRS geometry | GENERATOR NOT LOCATED; production numerical artifact retained and final PDF verified (non-blocking) |
| BS/BRS headline point gain and angle | COMPLETE |
| BS weak-axis printed-covariance result | VALIDATION ONLY (plus clean chain interface) |
| BR/BRS point result | CLEAN REPRODUCER COVERS RESULT; original production artifact not located (non-blocking) |
| D3/D3+BRS point result | ORIGINAL PRODUCTION ARTIFACT AND SCRIPT FOUND; clean reproduction passed |
| BR/BRS exact production bootstrap | PROVENANCE GAP (non-blocking to frozen point-result reproducibility) |
| D3/D3+BRS exact production bootstrap | ORIGINAL PRODUCTION ARTIFACT AND SCRIPT FOUND; ensemble not rerun |
| Appendix A determinant closure | VALIDATION ONLY |
| Appendix B reciprocal diagnostic | PARTIAL |
| Appendix C filename/hash identity | COMPLETE for external and public-D3 checksum records; author sign-off resolved |
| Appendix D generalized-eigenvalue calculation | VALIDATION ONLY |

Detailed commands and limitations are recorded in `docs/REPRODUCIBILITY_MAP.md` and `docs/NUMERICAL_PROVENANCE.md`.

## 7. Syntax/runtime checks

- `python -m py_compile` passed for the four existing Python scripts and all ten newly copied D3 production-module files.
- `--help` passed for the clean reproducer and each of its `point`, `bootstrap`, and `robustness` subcommands.
- `--help` passed for both immutable production scripts.
- No script implements `--check-only`; none was added because the existing clean subcommand interface is sufficient.
- Runtime checks used a temporary isolated Python 3.14 environment with NumPy 2.5.2 and SciPy 1.18.1 installed from `requirements.txt`. This temporary environment is not part of the release.
- No bootstrap or robustness ensemble was rerun. During this targeted closure, only clean deterministic `point` commands were run against the unambiguously identified local BR/BS/BRS and public D3/D3+BRS chains. All completed successfully.

## 8. Numerical validation

`python scripts/check_published_covariances.py` completed successfully. Its stdout matched `results/check_published_covariances_output.txt` byte-for-byte and returned:

- $\Lambda_1=3.150001102925091$
- $\Lambda_2=1.0067755437366135$
- $\Delta\phi_{\rm weak}=4.6329895023056835$ deg
- $e_{\rm weak}^{\rm BS}=(-0.08613667,0.99628333)$, up to sign
- $\mathrm{FoM}_{\rm BS}=33.40694396953772$
- $\mathrm{FoM}_{\rm BRS}=59.49198095418582$

Determinant closure passed: the generalized-eigenvalue product and covariance-determinant ratio agree to floating-point precision. An independent read-only calculation from the archived 1,000-row CSV reproduced the 16/50/84 summaries, 56.9%, 91.6%, both zero exceedance counts, and correlation $-0.204985...$.

Clean point reproduction additionally returned BS/BRS $(3.14989043,1.00667765)$, BR/BRS $(6.73433387,1.49721969)$, and D3/D3+BRS $(37.22034260,7.06184552)$, with the frozen Table 1/D3 ESS, means, covariances, FoMs, and determinant closure. This confirms public reproducibility of the frozen point results independently of incomplete historical production provenance.

## 9. Figure identity checks

All release figures are byte-identical to the frozen submission-package figures:

| Figure | SHA256 |
|---|---|
| Figure 1, `fig_BS_to_BRS_geometry_v07.pdf` | `3e26fa7d6b7e2ecda53b4f64c134341cdec84117733b34881d72ddd70ac7ad10` |
| Figure 2, `fig1_BR_to_BRS_geometry_v06.pdf` | `470a1a253f16ced2f5eb0bc524b2c93ef4a5ca46f9a1ae525e423e3a264fcd4f` |
| Figure 3, `fig2_D3_to_D3BRS_geometry_v06.pdf` | `932f2e623df6e421045039677b876ef55aed0156ab3480a226035e3321e4fdc9` |

The figures were not regenerated or altered.

Targeted generator search used exact output filenames, output-path strings, manuscript covariance values, figure/axis terminology, plotting calls, script comments, project notebooks, and local Git-visible filenames. The only project plotting module found, `analysis/make_figures.py`, explicitly permits development-mock output only, writes different filenames, and does not contain the paper covariance inputs. It was rejected as a generator candidate. Figures 1--3 are each classified **GENERATOR NOT LOCATED**.

## 10. Privacy/public-release audit

The staging tree contains no raw `chain_*.txt` file, private email address, private correspondence text, screenshot, absolute `/Users/...` path, account path, API key/token, password, private key, Codex log, shell history, cache, bytecode, or personal working note. Production scripts and archived JSON contain only portable relative data paths; immutable provenance copies were not modified.

## 11. License/citation metadata

The MIT license and scope notice follow the philosophy of the prior author release. The notice limits the license claim to original author-created code/documentation and explicitly excludes chain-data redistribution and relicensing of DES/DESI/Pantheon or other survey products. `CITATION.cff` contains author, version, local staging date, and an unpublished preferred paper title, with no invented DOI, journal metadata, or repository URL.

## 12. SHA256 integrity

Top-level `SHA256SUMS` covers every permanent release file except itself using stable relative paths. A clean `shasum -a 256 -c SHA256SUMS` verification passed after package assembly. No temporary audit logs or caches are listed.

## 13. Remaining gaps

- Exact figure-generation scripts were not located; final PDFs are preserved and verified and their numerical content is reproducible, but exact historical plot regeneration is unavailable.
- A complete historical Table 1 production summary and original BR/BRS point/bootstrap artifact were not located; the clean reproducer regenerates the frozen point/table values.
- The exact original production artifact for Appendix B was not located; the clean generic point interface covers the calculation.
- The standalone generator for `local_ess_bootstrap.json` was not identified, although the immutable output and D3 production summary are retained.

All remaining gaps concern historical production provenance. They are transparently disclosed and are **NON-BLOCKING** because the frozen scientific point results are independently reproducible, central BS/BRS and D3 validation outputs are archived, chain identities are fixed by SHA256, and no unavailable file is misrepresented as original production code.

## 14. Final verdict

A. Ready for GitHub/Zenodo preparation
