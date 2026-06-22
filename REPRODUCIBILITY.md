# Reproducibility guide

This guide lists the reviewer-facing steps for reproducing the simulations and comparison figures stored in this repository.

Run each workflow from the folder listed in that section. The scripts use relative paths.

## 0. Clone and prepare the repository

```bash
git clone https://github.com/Jaykumar9033/FRACMATH.git
cd FRACMATH
git lfs install
git lfs pull
python -m pip install -r requirements.txt
```

Required external software:

- MATLAB for all MATLAB workflows.
- Abaqus/CAE and Abaqus/Standard for the Abaqus validation workflow.
- A Fortran compiler configured with Abaqus for the UMAT.
- Python with `numpy` and `matplotlib` for plotting.

Optional fast smoke check from the repository root:

```matlab
addpath('tests')
run_smoke_checks
```

This check verifies that key paper assets, benchmark input files, generated result files, and documentation files are present and readable. It is not a substitute for the full workflows below.

## 1. Reproduce the 2D 3PB MATLAB simulation

Working folder:

```text
3pb/matlab
```

MATLAB command:

```matlab
solver_main_3pb
```

The solver reads:

- `Gregoire_3PB/nodes.txt`
- `Gregoire_3PB/elements.txt`
- `Gregoire_3PB/top_nodes.txt`
- `Gregoire_3PB/left_nodes.txt`
- `Gregoire_3PB/right_nodes.txt`
- `Gregoire_3PB/cmod1.txt`
- `Gregoire_3PB/cmod2.txt`

The main outputs are written to:

```text
3pb/matlab/Gregoire_3PB/results
```

Expected reviewer files:

- `matlab_load_cmod.csv`
- `matlab_timing.txt`
- `matlab_load_cmod_fig.png`
- `fig_damage_peak.png`
- `fig_damage_postpeak.png`
- `fig_damage_last_step.png`
- `simulation_video.mp4`

## 2. Reproduce the 2D 3PB Abaqus/UMAT simulation

Working folder:

```text
3pb/abaqus
```

Main command:

```bash
abaqus cae noGUI=run_3pb_abaqus_OLIVER_T3_FAST.py
```

Optional environment variables:

```bash
set ABQ_CPUS=4
set ABQ_FIELD_FREQ=100
set ABQ_AUTO_PLOT=1
```

Important files:

- `run_3pb_abaqus_OLIVER_T3_FAST.py`: builds the model, writes the Oliver T3 bandwidth table, runs Abaqus, extracts response data, and plots results.
- `cdm_umat_2d_OLIVER_T3_FAST.for`: Abaqus UMAT with the same modified von Mises damage model and Oliver direction-dependent bandwidth used in MATLAB.
- `Gregoire_3PB/Gregoire_3PB.odb`: Abaqus output database, tracked with Git LFS.

The main outputs are written to:

```text
3pb/abaqus/Gregoire_3PB/results
```

Expected reviewer files:

- `abaqus_load_cmod.csv`
- `abaqus_timing.txt`
- `abaqus_load_cmod_fig.png`
- `abaqus_fig_damage_peak.png`
- `abaqus_fig_damage_postpeak.png`
- `abaqus_fig_damage_last_step.png`

## 3. Recreate the MATLAB-vs-Abaqus comparison figures

Working folder:

```text
3pb/comparison
```

Commands:

```bash
python plot_comparison.py
python plot_image_comparison.py
```

Inputs are read from the MATLAB and Abaqus result folders. Outputs are written into `3pb/comparison/`.

Expected reviewer files:

- `comparison_summary.md`
- `comparison_load_cmod.png`
- `runtime_comparison.png`
- `performance_breakdown.png`
- `time_comparison_bar.png`

The stored summary reports peak load, CMOD at peak, wall-clock time, and total process time for MATLAB and Abaqus.

## 4. Reproduce the 3D Nooru-Mohamed benchmark

Working folder:

```text
Noor mohammad/Mesh
```

Optional mesh and boundary-condition visualization:

```matlab
visulizaiton
```

Damage simulation command:

```matlab
opts = struct();
opts.nIncr = 900;
opts.load_path = '4c';
opts.snapshot_stride = 1;
opts.live_damage_thresh = 0.005;
opts.save_damage_thresh = 0.95;
opts.show_live = true;
opts.live_stride = 1;
opts.show_mesh = true;
opts.save_show_mesh = true;
opts.damage_colormap = 'turbo';
opts.damage_clim_mode = 'visible';
opts.bandwidth_method = 'oliver';
damage_static('Job-1', opts);
```

The solver reads files such as:

- `Job-1_nodes.txt`
- `Job-1_elements.txt`
- `Job-1_top_nodes.txt`
- `Job-1_bottom_nodes.txt`
- `Job-1_left_nodes.txt`
- `Job-1_right_nodes.txt`
- `Job-1_BCs.txt`

The main outputs are written to:

```text
Noor mohammad/Mesh/out_NR_vectorized_LIVE_damage_mesh
```

Additional stored snapshots are in:

```text
```

## 5. Reproduce the 3D torsion benchmark

Working folder:

```text
Torsion/working
```

MATLAB command:

```matlab
run_torsion
```

The solver reads:

- `Job-1_nodes.txt`
- `Job-1_elements.txt`
- `Job-1_left_nodes.txt`
- `Job-1_right_nodes.txt`

The main outputs are written to:

```text
Torsion/working/out_torsion_LIVE_ONLY_OLIVER
```

Expected reviewer files include load/torque curves, CMOD tables, damage snapshots, and an animation MP4.

## 6. Rebuild the theory manual

Working folder:

```text
doc
```

The PDF is already stored as:

```text
doc/theory_manual.pdf
```

To rebuild from source with a local LaTeX installation:

```bash
pdflatex theory_manual.tex
pdflatex theory_manual.tex
```
