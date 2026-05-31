# Abaqus/UMAT 2D 3PB workflow

This folder contains the Abaqus/Standard validation workflow for the same notched 3PB benchmark used by the MATLAB solver.

## Main command

Run from this folder:

```text
3pb/abaqus
```

Command:

```bash
abaqus cae noGUI=run_3pb_abaqus_OLIVER_T3_FAST.py
```

## What the workflow does

1. Builds the 2D CPS3 Abaqus model.
2. Writes the Oliver T3 gradient table used by the UMAT.
3. Compiles and runs `cdm_umat_2d_OLIVER_T3_FAST.for`.
4. Extracts load-CMOD and damage information from the ODB.
5. Writes result CSV files and damage figures.

## Important files

| File or folder | Purpose |
| --- | --- |
| `run_3pb_abaqus_OLIVER_T3_FAST.py` | Main Abaqus build-run-extract-plot script |
| `cdm_umat_2d_OLIVER_T3_FAST.for` | Abaqus UMAT source code |
| `extract_damage.py` | Damage extraction helper |
| `extract_peak_omega.py` | Peak-damage extraction helper |
| `Gregoire_3PB/` | Abaqus job folder, ODB, logs, and extracted results |
| `Gregoire_3PB/results/` | Final Abaqus CSV files and figures |

## Optional controls

Set these before running Abaqus if needed:

```bash
set ABQ_CPUS=4
set ABQ_FIELD_FREQ=100
set ABQ_AUTO_PLOT=1
```

Meaning:

- `ABQ_CPUS`: number of Abaqus CPUs/domains.
- `ABQ_FIELD_FREQ`: field-output frequency.
- `ABQ_AUTO_PLOT`: set to `1` to plot immediately after the run.

## Main outputs

| File | Meaning |
| --- | --- |
| `Gregoire_3PB/results/abaqus_load_cmod.csv` | Abaqus CMOD and load response |
| `Gregoire_3PB/results/abaqus_timing.txt` | Abaqus timing summary |
| `Gregoire_3PB/results/abaqus_load_cmod_fig.png` | Abaqus load-CMOD plot |
| `Gregoire_3PB/results/abaqus_fig_damage_peak.png` | Damage at peak load |
| `Gregoire_3PB/results/abaqus_fig_damage_postpeak.png` | Damage after peak |
| `Gregoire_3PB/results/abaqus_fig_damage_last_step.png` | Final damage state |
| `Gregoire_3PB/Gregoire_3PB.odb` | Abaqus output database, stored through Git LFS |

## Requirements

- Abaqus/CAE and Abaqus/Standard.
- A Fortran compiler configured for Abaqus UMAT compilation.
- Python with `numpy` and `matplotlib` if plotting outside Abaqus Python.
