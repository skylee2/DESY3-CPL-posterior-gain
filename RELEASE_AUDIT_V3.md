# Release-package audit against manuscript v3

Audit date: 2026-09-08. Authoritative source: `paper1_Y3_submission_finaljournal_v3.tex`, SHA-256 `56919dc92f8466f3c2b23ffb8c860e49b8d833f97c73ed913e1c85893d79e340`.

## Discrepancies and actions

| Item | Current v3 value or requirement | v1.0.0 package state | Source in old package | Action in v1.1.0 candidate |
|---|---|---|---|---|
| Primary BS -> BRS point result | (3.1499, 1.0067), det 3.1709, r_gain 1.0062, angle 4.632 deg, D_ref about 0.25 | Numerical values agreed, but documentation covered only part of the current v3 audit | `docs/NUMERICAL_PROVENANCE.md`; BS production summary | Regenerated point JSON and verified every field and determinant closure |
| BS bootstrap | Current v3 intervals; 56.9%, 91.6%, 100%, 100% | Adopted v07.2 realization existed and agreed | `results/bs_brs_bootstrap_*` | Regenerated from verified inputs; retained adopted production artifact separately |
| BR -> BRS point result | (6.7343, 1.4972), det 10.0828, r_gain 1.1721, angle 3.345 deg, D_ref 1.382 | Point values reproduced, but exact bootstrap provenance was explicitly marked missing | old `docs/NUMERICAL_PROVENANCE.md` | Added point JSON and full 1000-realization bootstrap CSV/JSON/text |
| BR bootstrap | Lambda medians 6.718 and 1.500; r_gain 1.174; angle 3.35 deg | No production artifact included | old provenance “SOURCE NOT LOCATED” entry | Regenerated and archived exact machine-readable products |
| Figure 1 and Figure 2 generation | Validated BS -> BRS and BR -> BRS covariance geometry, means, weak axes, generalized directions, gains, and native-coordinate angles | Final PDFs were retained, but the original historical plotting scripts were not | final PDF artifacts and point-analysis outputs | Added explicitly reconstructed replacement generators driven only by packaged point-result JSON; no claim is made that they are the original plotters |
| D3 -> D3+BRS point result | (37.2203, 7.06185), det 262.8443, r_gain 1.32560, dominant angle 50.82 deg, D_ref 0.699 | Point values agreed | D3 Stage-1 JSON | Regenerated clean point JSON and reran original Stage-1 pipeline |
| D3 directional interpretation | Dominant v1 is 50.82 deg from weak axis; subdominant v2 is 2.23 deg away | Figure title said “weak-axis alignment”; legend called Lambda2 about 7.06 the near-weak-axis direction and Lambda1 about 37.23 the “other direction” | `figures/fig2_D3_to_D3BRS_geometry_v06.pdf` | Regenerated same-named PDF with explicit dominant-v1 and subdominant-v2 labels and metadata |
| D3 bootstrap | Median dominant angle 50.84 deg; 4.0% on weak side of 45 deg | Aggregate spectrum existed, but no packaged directional-realization output established the v3 interpretation | D3 Stage-1 summary/JSON | Added 400-realization directional CSV/JSON/text and verified 16/400 |
| Table 3 robustness | All BS, BR, and D3 hard-prior rows plus D3 68% HPD row | D3 values existed; complete BS/BR machine-readable Table 3 artifacts were absent | D3 Stage-1 outputs; documentation gap | Added BS and BR CSV/JSON and regenerated D3 JSON; verifier checks all rows |
| Synthetic validation | Analytic (10, 5), determinant 50, r_gain 169/97; exact recovered values; 12/12 checks | Current standalone validation and machine-readable output absent | no corresponding v1.0.0 artifact | Added documented seed, construction, tolerances, named tests, and exact JSON; earlier undocumented validation attempts are generically marked as superseded |
| Public-chain checks | z_p 0.24435; Q_0.05(w0) -1.39847; closure about 1.1e-15 | Values existed in D3 production JSON but were not part of a release-wide v3 gate | D3 Stage-1 outputs | Reran and incorporated into the v3 release verifier |
| Other manuscript numbers | Table 1 summaries and covariances; D3 vectors, axis ratios, angle matrix, local diagnostic; reverse BRS -> D3+BRS check | Values were dispersed across point and D3 production outputs | old clean reproducer and D3 Stage-1 outputs | Regenerated the reverse comparison and added explicit release-wide checks for all supporting values |
| Covariance check documentation | Current v3 manuscript | Script named obsolete `paper1_Y3_v09_final.tex` | `scripts/check_published_covariances.py` | Updated documentation string only; calculation unchanged |
| Dependencies | NumPy, SciPy, pandas, Matplotlib | `requirements.txt` omitted pandas and Matplotlib | `requirements.txt` | Added runtime dependencies needed by the preserved D3 pipeline and figure generator |
| Release metadata | Current title and version | Title contained “Historical”; version 1.0.0 | README and citation metadata | Updated candidate metadata to v1.1.0 and current manuscript title |
| Integrity products | Hashes for rebuilt candidate and archive | Existing hashes/tarball described the older package | `SHA256SUMS`; v1.0.0 tarball | Rebuilt candidate from verified sources; generated new hashes and a distinct archive |

## Files deliberately not carried forward

- `scripts/production/bs_brs_bootstrap_v071.py`: superseded by the adopted v07.2 implementation.
- old v1.0.0 README, provenance, validation, reproducibility, manifest, release-note, checksum, and citation files: they describe the older package state.
- old Figure 3: its directional labels are inconsistent with v3.
- all raw historical DES external-only chains: excluded under the release policy.

## Version recommendation

Use **v1.1.0**. The update is minor-level because it adds new supported analysis surfaces (standalone validation, BR/D3 directional bootstraps, full robustness products, a Figure 3 generator, and a release-wide verification gate) while preserving the analysis interface and scientific conclusions. Although one figure correction alone might be patch-level, the aggregate reproducibility expansion is larger than a metadata or bug-fix-only v1.0.1 release.

This candidate is not published to GitHub or Zenodo.
