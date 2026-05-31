# FRACMATH JOSS Repository

This repository contains finite element fracture and damage simulations used for the JOSS project. It includes MATLAB solvers, Abaqus comparison workflows, generated plots, and GitHub Pages documentation for the main example cases.

## Repository layout

- `3PB/` - three-point bending benchmark with MATLAB, Abaqus, and comparison results.
- `Noor mohammad/` - Nooru-Mohamed 3D mesh, boundary condition visualization, and live damage simulation files.
- `Torsion/` - 3D torsion damage simulation and output visualizations.
- `index.html` - GitHub Pages landing page for the repository.
- `assets/css/style.css` - styling for the GitHub Pages site.

## Software requirements

Use the tools that match the workflow you want to run:

- MATLAB R2020b or newer is recommended for the `.m` solvers and visualization scripts.
- Python 3.9 or newer with `numpy`, `pandas`, and `matplotlib` for plotting/comparison scripts.
- Abaqus/CAE with Abaqus/Standard for the Abaqus 3PB workflow.
- A Fortran compiler configured with Abaqus for UMAT runs.

## Main files to run

### 1. Three-point bending MATLAB solver

Run this file from the `3PB/Matlab` folder:

```matlab
solver_main_3pb
```

Expected working folder:

```text
3PB/Matlab
```

The solver reads mesh files from `3PB/Matlab/Gregoire_3PB/` and writes plots, CSV data, timing information, and video output to:

```text
3PB/Matlab/Gregoire_3PB/results/
```

### 2. Three-point bending Abaqus pipeline

Run this command from the `3PB/abaqus` folder:

```bash
abaqus cae noGUI=RUN_3PB_ABAQUS_FULL_MATCH_MATLAB.py
```

Expected working folder:

```text
3PB/abaqus
```

This workflow builds the Abaqus model, runs the UMAT job, extracts Load-CMOD data, and writes Abaqus result plots to:

```text
3PB/abaqus/abaqus_results/
```

Required UMAT file:

```text
3PB/abaqus/scm_umat_2d_OLIVER_MATCH_MATLAB.for
```

### 3. Three-point bending comparison plots

Run this file from the comparison plot folder:

```bash
python make_plots.py
```

Expected working folder:

```text
3PB/Comparison/comparison_plots
```

The script reads MATLAB and Abaqus result files from:

```text
3PB/Comparison/Matlab_result/
3PB/Comparison/abaqus_result/
```

It writes publication-style comparison plots and a summary file into:

```text
3PB/Comparison/comparison_plots/
```

### 4. Nooru-Mohamed 3D damage simulation

Run this MATLAB file from the `Noor mohammad/Mesh` folder. The file is `damage_static_NR_vectorized_LIVE_damage_different_colours.m`:

```matlab
opts = struct();
opts.nIncr = 900;
opts.load_path = '4c';
opts.snapshot_stride = 1;
opts.show_live = true;
opts.show_mesh = true;
opts.bandwidth_method = 'oliver';
damage_static_NR_vectorized_LIVE_damage_different_colours('Job-1', opts);
```

Expected working folder:

```text
Noor mohammad/Mesh
```

Required input files include:

```text
Job-1_nodes.txt
Job-1_elements.txt
Job-1_top_nodes.txt
Job-1_bottom_nodes.txt
Job-1_left_nodes.txt
Job-1_right_nodes.txt
```

### 5. Nooru-Mohamed mesh and boundary condition visualization

Run this MATLAB file from the `Noor mohammad/Mesh` folder. The file is `visulizaiton.m`:

```matlab
visulizaiton
```

It creates the 3D mesh and boundary condition figures, including:

```text
nooru_mesh_3D.png
nooru_BC_2D_no_dimensions.png
```

### 6. Torsion damage simulation

Run this MATLAB file from the `Torsion/working` folder. The file is `run_torsion_stress_strain_static_fast_modvm_multiview_video.m`:

```matlab
run_torsion_stress_strain_static_fast_modvm_multiview_video
```

Expected working folder:

```text
Torsion/working
```

Required input files include:

```text
Job-1_nodes.txt
Job-1_elements.txt
Job-1_left_nodes.txt
Job-1_right_nodes.txt
Job-1_left_elements.txt
Job-1_right_elements.txt
```

The output folder is:

```text
Torsion/working/out_torsion_LIVE_ONLY_OLIVER/
```

## Suggested run order

1. Run the MATLAB 3PB solver: `3PB/Matlab/solver_main_3pb.m`.
2. Run the Abaqus 3PB pipeline: `3PB/abaqus/RUN_3PB_ABAQUS_FULL_MATCH_MATLAB.py`.
3. Copy or confirm the latest MATLAB/Abaqus outputs in `3PB/Comparison/Matlab_result/` and `3PB/Comparison/abaqus_result/`.
4. Generate comparison plots with `3PB/Comparison/comparison_plots/make_plots.py`.
5. Run the Nooru-Mohamed and torsion examples from their own working folders as needed.

## GitHub Pages

This repository includes a static GitHub Pages site in the repository root.

To publish it:

1. Push the repository to GitHub.
2. Open the repository on GitHub: `https://github.com/jaymavani7/JOSS_FRACMATH`.
3. Go to **Settings** -> **Pages**.
4. Under **Build and deployment**, choose **Deploy from a branch**.
5. Select branch `main` and folder `/root`.
6. Save the settings.

After GitHub builds the site, it should be available at:

```text
https://jaymavani7.github.io/JOSS_FRACMATH/
```

## Notes

- Run MATLAB scripts from the folders listed above because several scripts read input files using relative paths.
- Large generated files such as `.odb`, `.mp4`, Abaqus scratch files, and increment image sequences can make the repository large. Consider using Git LFS or excluding temporary solver files before final publication.
- The repository currently contains generated results so the GitHub Pages site can show figures without rerunning every simulation.


