# 3D torsion benchmark

This folder contains the 3D notched-beam torsion benchmark used in the FRACMATH paper.

## Folder map

| Path | Purpose |
| --- | --- |
| `working/` | Main MATLAB torsion solver, mesh files, and output folders |
| `visulize mesh/` | Mesh visualization cases and Abaqus model files |

## Main workflow

Use the instructions in `working/README.md`.

Short version:

1. Open MATLAB.
2. Change the current folder to `Torsion/working`.
3. Run `run_torsion`.
4. Inspect `out_torsion_LIVE_ONLY_OLIVER/`.

## Main outputs

| Path | Contents |
| --- | --- |
| `working/out_torsion_LIVE_ONLY_OLIVER/` | Main Oliver-bandwidth torsion outputs |
| `visulize mesh/` | Mesh visualization figures and input files |
