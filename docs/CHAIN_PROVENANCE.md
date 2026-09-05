# Chain provenance

The historical external-only chains used for the central comparisons were run for the DES Year 3 extensions analysis. They were supplied directly by Jessica Muir from the DES analysis team. DES regards chains run for that Y3 paper as public products, although these particular external-only chains were not all uploaded to the public web release. The author has reviewed and signed off on this provenance statement.

Filename interpretation follows the DES Y3 extensions release labeling scheme: `br`, `bs`, `brs`, and `pbrs` identify the corresponding historical external-probe combinations. Exact local file identity is fixed by the filename/SHA256 pairs in `../input/external_chain_SHA256SUMS.txt`. The author manually recomputed all four hashes and confirmed exact agreement with that record.

For the BS/BRS comparison, parameters common to both runs use the same priors according to the DES-team clarification. The BS run fixes $A_s$ and $n_s$, while the BRS run varies and marginalizes over them. The author has signed off on this interpretation. This document records the clarification without reproducing private correspondence.

The exact public D3 inputs were identified unambiguously from the paths recorded in the archived D3 production JSON and by agreement of their sample counts, ESS values, means, covariances, and generalized gains with the frozen manuscript. Their filename/SHA256 pairs are recorded in `../input/public_d3_chain_SHA256SUMS.txt`.

Raw chain files are deliberately **not redistributed** in this release. No private email address, screenshot, or extended quotation from private correspondence is included. Private-provenance author sign-off is **RESOLVED**.
