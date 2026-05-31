# MATLAB 2D 3PB solver

This folder contains the MATLAB implementation of the 2D notched three-point bending continuum damage benchmark.

## Main file to run

```matlab
solver_main_3pb
```

Run it from this folder:

```text
3pb/matlab
```

## Step-by-step

1. Open MATLAB.
2. Change the current folder to `3pb/matlab`.
3. Run:

   ```matlab
   solver_main_3pb
   ```

4. Wait for the solver to finish the displacement-control load steps.
5. Inspect the result folder:

   ```text
   Gregoire_3PB/results
   ```

## Required inputs

The solver expects these files in `Gregoire_3PB/`:

- `nodes.txt`
- `elements.txt`
- `top_nodes.txt`
- `left_nodes.txt`
- `right_nodes.txt`
- `cmod1.txt`
- `cmod2.txt`

## Main outputs

| File | Meaning |
| --- | --- |
| `Gregoire_3PB/results/matlab_load_cmod.csv` | CMOD and load response |
| `Gregoire_3PB/results/matlab_timing.txt` | Timing and peak-load information |
| `Gregoire_3PB/results/fig_mesh.png` | Mesh and boundary-condition figure |
| `Gregoire_3PB/results/matlab_load_cmod_fig.png` | Load-CMOD curve |
| `Gregoire_3PB/results/fig_damage_peak.png` | Damage at peak load |
| `Gregoire_3PB/results/fig_damage_postpeak.png` | Damage after peak |
| `Gregoire_3PB/results/fig_damage_last_step.png` | Final damage state |
| `Gregoire_3PB/results/simulation_video.mp4` | Simulation animation |

## Notes

- The solver uses the modified von Mises equivalent strain, exponential softening, and direction-dependent Oliver crack-band regularization.
- The default material and solver parameters are defined near the top of `solver_main_3pb.m`.
- Keep the working folder at `3pb/matlab` because paths are relative.
