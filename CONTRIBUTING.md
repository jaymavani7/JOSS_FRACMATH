# Contributing to FRACMATH

Thank you for your interest in improving FRACMATH. This repository is prepared for open review and for reproducible research use.

## Reporting issues

Please use GitHub Issues to report problems, unclear documentation, failed reproduction steps, or suggested improvements. Include:

- the operating system,
- MATLAB, Abaqus, Python, and compiler versions when relevant,
- the folder and script you ran,
- the exact command or MATLAB call,
- the full error message or unexpected output.

## Asking for support

For reviewer or user support, open a GitHub Issue with a short title and the workflow name, for example `3PB MATLAB solver`, `Abaqus UMAT`, `Nooru-Mohamed`, `torsion`, or `JOSS paper`.

## Proposing changes

Small documentation fixes can be proposed directly in a pull request. For solver, model, or benchmark changes, please open an issue first so the expected behavior, validation data, and reproducibility impact can be discussed.

## Development expectations

- Keep benchmark input files and generated outputs traceable to the script that created them.
- Document any new material parameters, meshes, or solver options.
- Prefer focused commits that separate source-code changes from generated-result changes.
- Run `tests/run_smoke_checks.m` before submitting a pull request.

## Authorship and citation

Contributions will be acknowledged according to their role and scope. If you use FRACMATH in research, please cite the software metadata in `CITATION.cff` and the associated JOSS paper after publication.
