# External chain inputs

No DES chain is distributed in this package. Users with authorized or public local copies must provide chain paths explicitly to `scripts/paper1_y3_reproduce.py` or to the archived production scripts.

The exact filenames and SHA256 identities of the four historical external-only chains are recorded in `external_chain_SHA256SUMS.txt`. Verify a directory without copying data into this package, for example:

```bash
cd /path/to/des_y3_chains
shasum -a 256 -c /path/to/DESY3_CPL_posterior_gain_v1.0.0/input/external_chain_SHA256SUMS.txt
```

The clean reproducer also accepts explicit paths to the public DES Y3 D3 and D3+BRS products used for the projected-LSS comparison. Those products are not bundled here.
