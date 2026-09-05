# Validation

## A. Printed-covariance independent check

Run:

```bash
python scripts/check_published_covariances.py
```

Using only the six-decimal $C_{\rm BS}$ and $C_{\rm BRS}$ matrices printed in the manuscript, the independent utility returns approximately:

- $\Lambda_1=3.1500011$
- $\Lambda_2=1.0067755$
- $\Delta\phi_{\rm weak}=4.63299$ deg
- $e_{\rm weak}^{\rm BS}=(-0.08614,0.99628)$, up to sign
- $\mathrm{FoM}_{\rm BS}=33.4069$
- $\mathrm{FoM}_{\rm BRS}=59.4920$

The tiny difference from full-chain values is expected because the manuscript covariances are printed to six decimal places.

## B. Determinant closure

For the printed-covariance check,

```text
Lambda_1 * Lambda_2       = 3.1713440731683407
det(C_BS) / det(C_BRS)    = 3.171344073168344
(FoM_BRS / FoM_BS)^2      = 3.171344073168344 (to displayed precision)
```

Thus the determinant identity

```text
Lambda_1 Lambda_2 = det(C_ref)/det(C_new) = (FoM_new/FoM_ref)^2
```

closes to floating-point precision.

## C. Production versus clean reproducer

The scripts in `scripts/production/` document what was actually used during manuscript development and are retained unchanged, including their historical names and defaults. `scripts/paper1_y3_reproduce.py` is a separate clean, documented implementation of weighted chain ingestion, covariance construction, generalized gains, alignment angles, multiplier bootstrap, and robustness cuts.

Both forms are retained intentionally: production files establish provenance, while the clean reproducer provides the public-facing interface. The archived v07.2 CSV and summary are the authoritative realization-level and aggregate production outputs included here. The clean bootstrap is not asserted to reproduce their bytes because the exact RNG stream can depend on row-retention details even when the declared seed and statistical prescription agree.
