# 3PB comparison figures

This folder contains the regenerated MATLAB and Abaqus/UMAT 3PB comparison figures used in the JOSS paper.

## Commands

Run from this folder:

```bash
python plot_comparison.py
```

## Inputs

- `../matlab/Gregoire_3PB/results/`
- `../abaqus/Gregoire_3PB/results/`

## Outputs

| File | Purpose |
| --- | --- |
| `comparison_summary.md` | Stored MATLAB--Abaqus summary values |
| `comparison_load_cmod.png` | Load-CMOD comparison |
| `runtime_comparison.png` | Runtime comparison |
| `performance_breakdown.png` | Combined performance breakdown |
| `comparison_load_cmod_panels.png` | Side-by-side response and timing panels |

## Stored comparison values

- MATLAB peak load: 3.64 kN.
- Abaqus + UMAT peak load: 3.61 kN.
- CMOD at peak: MATLAB 0.022811 mm; Abaqus 0.022485 mm.
- MATLAB solver wall-clock time: 547.58 s.
- Abaqus submit-to-done wall-clock time: 1996.25 s.
- Abaqus internal `.msg/.dat` wall-clock time: 1968.00 s.
- Abaqus Load-CMOD rows: 10001.
