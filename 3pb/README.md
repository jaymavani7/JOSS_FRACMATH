# 2D three-point bending benchmark

This folder contains the main 2D notched three-point bending validation case for the FRACMAT paper.

The benchmark compares:

1. A vectorized MATLAB continuum damage mechanics solver.
2. An Abaqus/Standard model using an Oliver-matched UMAT.
3. Python scripts that create the final MATLAB-vs-Abaqus comparison figures.

## Folder map

| Path | Purpose |
| --- | --- |
| `matlab/` | MATLAB 2D solver, mesh input files, and MATLAB results |
| `abaqus/` | Abaqus model generation, UMAT, ODB extraction, and Abaqus results |
| `comparison/` | Scripts and figures comparing MATLAB and Abaqus results |

## Recommended run order

1. Run the MATLAB solver from `3pb/matlab/`.
2. Run the Abaqus/UMAT workflow from `3pb/abaqus/`.
3. Run the comparison scripts from `3pb/comparison/`.

Detailed instructions are in each subfolder README and in the root `REPRODUCIBILITY.md`.

## Main outputs

| Output | Location |
| --- | --- |
| MATLAB load-CMOD curve | `matlab/Gregoire_3PB/results/matlab_load_cmod.csv` |
| Abaqus load-CMOD curve | `abaqus/Gregoire_3PB/results/abaqus_load_cmod.csv` |
| MATLAB timing log | `matlab/Gregoire_3PB/results/matlab_timing.txt` |
| Abaqus timing log | `abaqus/Gregoire_3PB/results/abaqus_timing.txt` |
| Comparison summary | `comparison/comparison_summary.md` |
| Final comparison figures | `comparison/*.png` and `comparison/*.pdf` |

The Abaqus ODB file is large and is stored through Git LFS.
