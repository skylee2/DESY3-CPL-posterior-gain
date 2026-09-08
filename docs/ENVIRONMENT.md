# Verified environment

The full candidate audit was run with:

- Python 3.14.6;
- Darwin 25.5.0 arm64;
- the exact packages in `requirements-lock.txt`.

`requirements.txt` gives supported minimum versions. `requirements-lock.txt` records the fresh environment used for the archived regeneration and validation. Matplotlib was run with the noninteractive `Agg` backend for deterministic PDF generation.
