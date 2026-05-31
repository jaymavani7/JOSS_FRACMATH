# 3PB comparison figures

This folder contains the scripts and stored outputs used to compare the MATLAB and Abaqus/UMAT 3PB simulations.

## Commands

Run from this folder:

```text
3pb/comparison
```

Commands:

```bash
python plot_comparison.py
python plot_image_comparison.py
```

## Inputs

The scripts read the stored result files from:

- `../matlab/Gregoire_3PB/results/`
- `../abaqus/Gregoire_3PB/results/`

## Outputs

| File | Purpose |
| --- | --- |
| `comparison_summary.md` | Markdown and LaTeX summary table |
| `comparison_load_cmod.png` | MATLAB-vs-Abaqus load-CMOD comparison |
| `comparison_load_cmod.pdf` | PDF version of the load-CMOD comparison |
| `runtime_comparison.png` | Runtime comparison figure |
| `runtime_comparison.pdf` | PDF version of the runtime comparison |
| `performance_breakdown.png` | Combined performance summary figure |
| `performance_breakdown.pdf` | PDF version of the performance summary |
| `time_comparison_bar.png` | Bar chart comparing time |

## Stored comparison values

The stored summary reports:

- MATLAB peak load: 3.90 kN.
- Abaqus + UMAT peak load: 3.84 kN.
- CMOD at peak: approximately 0.026 mm for both solvers.
- MATLAB wall-clock time: 97.18 s.
- Abaqus wall-clock time: 1336.00 s.
