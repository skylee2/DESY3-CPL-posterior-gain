# Numerical provenance

Classifications used below are `PRODUCTION OUTPUT`, `CLEAN REPRODUCTION`, `INDEPENDENT VALIDATION`, and `TEXTUAL CHECK ONLY`. Where the precise generating artifact was not found in the designated source directories, this is stated explicitly.

| Result | Inputs | Script | Archived output | Paper value | Classification |
|---|---|---|---|---|---|
| BS $\to$ BRS point spectrum | BS, BRS | `scripts/production/bs_brs_bootstrap_v072.py`; independently `scripts/paper1_y3_reproduce.py point` | `results/bs_brs_bootstrap_summary_v072.txt` | $(\Lambda_1,\Lambda_2)=(3.1499,1.0067)$ | PRODUCTION OUTPUT; CLEAN REPRODUCTION |
| BS $\to$ BRS incremental gains | BS, BRS | same | same | $(\delta_1,\delta_2)=(2.1499,0.0067)$ | CLEAN REPRODUCTION; TEXTUAL CHECK ONLY for rounded paper values |
| BS $\to$ BRS effective-rank statistic | BS, BRS | same | same | $r_{\rm gain}\simeq1.0062$ | PRODUCTION OUTPUT; CLEAN REPRODUCTION |
| BS dominant-gain/weak-axis angle | BS, BRS | same | same | $\Delta\phi_{\rm weak}\simeq4.632$ deg | PRODUCTION OUTPUT; CLEAN REPRODUCTION |
| BS/BRS FoMs | BS, BRS | `scripts/paper1_y3_reproduce.py point` | printed-covariance check output gives the rounded-matrix check | $33.41\to59.49$ | CLEAN REPRODUCTION; INDEPENDENT VALIDATION |
| BS weak covariance axis | printed $C_{\rm BS}$, or BS chain | `scripts/check_published_covariances.py`; clean reproducer | `results/check_published_covariances_output.txt` | $e_{\rm weak}^{\rm BS}\simeq(-0.086,0.996)$, up to sign | INDEPENDENT VALIDATION; CLEAN REPRODUCTION |
| BS $\to$ BRS bootstrap $\Lambda_1$ | BS, BRS | `scripts/production/bs_brs_bootstrap_v072.py` | CSV and summary | $3.142^{+0.112}_{-0.101}$ | PRODUCTION OUTPUT |
| BS $\to$ BRS bootstrap $\Lambda_2$ | BS, BRS | same | same | $1.006^{+0.034}_{-0.033}$ | PRODUCTION OUTPUT |
| $P(\Lambda_2\ge1)$ | BS, BRS | same | same | 56.9% | PRODUCTION OUTPUT |
| Bootstrap weak-axis angle | BS, BRS | same | same | $4.65^{+0.27}_{-0.26}$ deg | PRODUCTION OUTPUT |
| Angular stability | BS, BRS | same | CSV and summary | 91.6% within 5 deg; 0/1000 reaching 10 deg; 0/1000 crossing 45 deg | PRODUCTION OUTPUT |
| Bootstrap correlation | BS, BRS | same | same | $\mathrm{Corr}(\Lambda_1,\Delta\phi_{\rm weak})=-0.205$ | PRODUCTION OUTPUT |
| BR $\to$ BRS point spectrum | BR, BRS | `scripts/paper1_y3_reproduce.py point` | no original production output located | $(6.7343,1.4972)$ | CLEAN REPRODUCTION verified locally; historical production provenance gap is non-blocking |
| D3 $\to$ D3+BRS point spectrum | public D3, D3+BRS | original `scripts/production/d3_stage1/analysis/run_y3_stage1.py`; clean point command | `results/d3_stage1/core_validation.json`, `stage1_validation.json`, `stage1_summary.md` | $(37.2203,7.06185)$ | PRODUCTION OUTPUT; ORIGINAL SCRIPT FOUND; CLEAN REPRODUCTION |
| External-chain sample counts, ESS, means, covariances, FoMs, condition numbers | BR, BS, BRS | `scripts/paper1_y3_reproduce.py point` | BS/BRS summary and printed-covariance check; no complete historical Table 1 output | Table 1 values | CLEAN REPRODUCTION verified locally; partial historical production provenance |
| D3/D3+BRS ESS and point diagnostics | public D3, D3+BRS | original D3 Stage-1 driver; clean point command | D3 Stage-1 JSON/summary and `public_d3_closure.json` | manuscript values including $N_{\rm eff}=5881.15,3114.64$ | PRODUCTION OUTPUT; CLEAN REPRODUCTION |
| D3 hard-boundary and 68% HPD checks | D3, D3+BRS | original D3 Stage-1 driver; clean robustness command | `stage1_validation.json`, `stage1_summary.md` | Table 3 values | PRODUCTION OUTPUT; ORIGINAL SCRIPT FOUND; CLEAN REPRODUCTION |
| D3 local-ESS conditional diagnostic | D3, D3+BRS | historical generator not separately identified | `results/d3_stage1/local_ess_bootstrap.json`; summarized in `stage1_summary.md` | manuscript local-ESS/interval values | PRODUCTION OUTPUT retained; exact standalone generator not located |
| BR/BRS bootstrap intervals | BR, BRS | generic clean bootstrap implementation | none located | manuscript validation interval | PROVENANCE GAP — SOURCE NOT LOCATED for exact production script/output; non-blocking to point-result reproduction |
| D3/D3+BRS bootstrap intervals | public D3, D3+BRS | original D3 Stage-1 driver | `stage1_validation.json`, `stage1_summary.md` | manuscript 400-realization interval | PRODUCTION OUTPUT; ORIGINAL SCRIPT FOUND; bootstrap not rerun during gap closure |

The production v07.2 summary records 9,026 retained BS samples and 10,093 retained BRS samples for its bootstrap ingestion. Table 1 reports stored chain rows and effective sample sizes separately; these quantities should not be conflated.

The archived v07.2 CSV contains 1,000 realizations plus its header. Its summary records seed `20260904`. The clean reproducer implements the declared statistical prescription, but its output is not asserted to be byte-for-byte identical to the archived stream because exact results can depend on RNG/row-retention details.

During targeted gap closure, deterministic clean `point` commands were run (no bootstrap). They reproduced BS/BRS $(3.14989043,1.00667765)$, BR/BRS $(6.73433387,1.49721969)$, and D3/D3+BRS $(37.22034260,7.06184552)$ together with the frozen ESS, means, covariances, FoMs, and determinant closures. This establishes public numerical reproducibility independently of whether every historical ad hoc production summary survives.
