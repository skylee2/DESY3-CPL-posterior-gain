# Chain provenance and redistribution policy

## Historical external-only chains

The BR, BS, BRS, and PBRS products are historical DES external-only chains. Their complete filenames and SHA-256 hashes are in `input/external_chain_SHA256SUMS.txt`.

The headers identify PolyChord with 500 live points and 60 repeats and a shared prior file, `extensions_fiducial/extensions_ini_files/priors_ml_datasettings.ini`. BS has seven varied parameters with fixed `n_s = 0.97` and `A_s = 2.19e-9`. BRS has nine varied parameters and samples `n_s` and `A_s` over the header-recorded prior ranges. PBRS is tracked for provenance although it is not required for the three directional comparisons in this package.

These raw chains are not publicly redistributed in this candidate. The checksum records permit an authorized local copy to be verified without conflating historical external-only inputs with public D3 products.

## Public D3 products

The public inputs are:

- `d3_w0wa_nla_realy3dat.txt` for DES Y3 3x2pt;
- `d3_brs_w0wa_nla_realy3dat.txt` for DES Y3 3x2pt+BRS.

Their hashes are in `input/public_d3_chain_SHA256SUMS.txt`. D3 -> D3+BRS adds the full BRS external block and is not an RSD-only test.

Machine-readable results store these two basenames as portable input identifiers. They are not paths into a parent development repository. To reproduce an analysis, obtain the public files, verify their hashes, and pass their actual local paths explicitly on the command line; generated output will retain only the basename identifiers.

## Packaged data policy

Only checksums, derived numerical products, scripts, figures, and documentation are included. No filename matching a raw historical `chain_*.txt` is present. The release-wide verifier enforces this policy.
