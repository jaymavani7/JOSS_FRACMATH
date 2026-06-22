# FRACMATH Smoke Checks

This folder contains lightweight checks intended for reviewers before running the full MATLAB, Abaqus, and plotting workflows.

The smoke check verifies that key paper assets, benchmark input files, result files, and documentation files are present and readable. It does not replace the full reproducibility workflow in `REPRODUCIBILITY.md`.

From the repository root, run:

```matlab
addpath('tests')
run_smoke_checks
```

Expected result:

```text
FRACMATH smoke checks passed.
```

The full benchmark simulations remain documented in `REPRODUCIBILITY.md`.
