# Torsion working folder

This folder contains the main MATLAB workflow for the 3D torsion benchmark.

## Main file to run

```matlab
run_torsion
```

Run it from this folder:

```text
Torsion/working
```

## Step-by-step

1. Open MATLAB.
2. Change the current folder to `Torsion/working`.
3. Run:

   ```matlab
   run_torsion
   ```

4. Inspect:

   ```text
   out_torsion_LIVE_ONLY_OLIVER
   ```

## Required inputs

The solver expects these mesh files in the current folder:

- `Job-1_nodes.txt`
- `Job-1_elements.txt`
- `Job-1_left_nodes.txt`
- `Job-1_right_nodes.txt`

## Main outputs

| Output pattern | Meaning |
| --- | --- |
| `out_torsion_LIVE_ONLY_OLIVER/*T_theta_CMOD_static.csv` | torsion response table |
| `out_torsion_LIVE_ONLY_OLIVER/*torque_theta*.png` | torque-rotation curves |
| `out_torsion_LIVE_ONLY_OLIVER/*load_theta*.png` | load-rotation curves |
| `out_torsion_LIVE_ONLY_OLIVER/*PEAK_DAMAGE_VECTOR.csv` | peak-damage data |
| `out_torsion_LIVE_ONLY_OLIVER/*PEAK_LOAD_DAMAGE_iso.png` | peak-damage view |
| `out_torsion_LIVE_ONLY_OLIVER/*animation.mp4` | saved animation |

## Notes

- The script uses a fully vectorized TET4 damage solver.
- The crack-band length uses the direction-dependent Oliver bandwidth computed from the maximum principal strain direction.
- Units are N, mm, and MPa.
