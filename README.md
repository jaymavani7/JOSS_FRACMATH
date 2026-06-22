# FRACMATH paper repository

This repository contains the code, input data, generated results, and documentation for the JOSS paper:

**FRACMATH: A vectorized MATLAB framework for continuum damage mechanics with crack-band regularization**

Authors: Jaykumar Mavani and Madura Pathirage, Department of Civil, Construction, and Environmental Engineering, University of New Mexico.

The GitHub repository is `Jaykumar9033/FRACMATH`, and the software described in the paper is `FRACMATH`.

## What is included

| Path | Contents | Main purpose |
| --- | --- | --- |
| `3pb/` | MATLAB solver, Abaqus/UMAT workflow, results, and comparison plots | Main 2D notched three-point bending validation case |
| `Noor mohammad/` | 3D Nooru-Mohamed mesh, MATLAB solver, and saved damage outputs | 3D mixed-mode validation case |
| `Torsion/` | 3D torsion mesh, MATLAB solver, visualization scripts, and outputs | 3D notched beam torsion validation case |
| `doc/` | Theory manual in LaTeX and PDF form | Formulation, derivations, and solver notes |
| `paper/` | JOSS manuscript, bibliography, and figures | Submission paper package |
| `tests/` | Lightweight MATLAB smoke checks | Fast repository completeness check |
| `MANIFEST.md` | File map for code, inputs, and generated outputs | Helps reviewers locate each item |
| `REPRODUCIBILITY.md` | One-by-one analysis instructions | Main reviewer run guide |
| `CONTRIBUTING.md` | Issue, support, and contribution guidance | Open-source workflow guidance |
| `CHANGELOG.md` | Version and submission-preparation notes | Release history |
| `CITATION.cff` | Citation metadata | Software citation |
| `requirements.txt` | Python plotting dependencies | Used by comparison and plotting scripts |

Large Abaqus output database files (`*.odb`) are tracked with Git LFS.

## Main code entry points

| Analysis | Folder to open first | File or command to run |
| --- | --- | --- |
| 2D 3PB MATLAB solver | `3pb/matlab/` | `solver_main_3pb` in MATLAB |
| 2D 3PB Abaqus + UMAT | `3pb/abaqus/` | `abaqus cae noGUI=run_3pb_abaqus_OLIVER_T3_FAST.py` |
| 2D 3PB comparison plots | `3pb/comparison/` | `python plot_comparison.py` and `python plot_image_comparison.py` |
| 3D Nooru-Mohamed mesh view | `Noor mohammad/Mesh/` | `visulizaiton` in MATLAB |
| 3D Nooru-Mohamed damage solve | `Noor mohammad/Mesh/` | `damage_static('Job-1', opts)` in MATLAB |
| 3D torsion damage solve | `Torsion/working/` | `run_torsion` in MATLAB |
| Theory manual | `doc/` | `theory_manual.pdf`, or build `theory_manual.tex` with LaTeX |

Run scripts from the folders listed above. Several scripts use relative paths to find mesh files, boundary-condition files, and result folders.

## Quick setup for reviewers

1. Install Git LFS before cloning or pulling this repository.

   ```bash
   git lfs install
   git lfs pull
   ```

2. Install the Python plotting dependencies.

   ```bash
   python -m pip install -r requirements.txt
   ```

3. Use MATLAB for the MATLAB workflows.

   Recommended: MATLAB R2020b or newer. The scripts use vectorized matrix operations and sparse linear solves.

4. Use Abaqus/CAE and Abaqus/Standard for the Abaqus workflow.

   The Abaqus 3PB run requires a Fortran compiler configured for Abaqus UMAT compilation.

5. Read `REPRODUCIBILITY.md` for the complete step-by-step run order.

6. Optionally run the lightweight MATLAB smoke check from the repository root.

   ```matlab
   addpath('tests')
   run_smoke_checks
   ```

## Analysis overview

### 1. 2D notched three-point bending benchmark

The 3PB benchmark is the main validation case. It compares a vectorized MATLAB continuum damage mechanics solver with an Abaqus/Standard UMAT implementation using the same CPS3 mesh, material parameters, displacement-control loading, and Oliver direction-dependent crack-band regularization.

Important folders:

- `3pb/matlab/Gregoire_3PB/`: MATLAB input mesh files and MATLAB-generated results.
- `3pb/abaqus/Gregoire_3PB/`: Abaqus model files, ODB output, extracted CSV data, and Abaqus-generated figures.
- `3pb/comparison/`: final MATLAB-vs-Abaqus comparison figures and summary table.

### 2. 3D Nooru-Mohamed benchmark

The Nooru-Mohamed case demonstrates the 3D tetrahedral MATLAB damage solver under mixed-mode loading. It includes mesh files, boundary-condition node sets, the vectorized Newton-Raphson solver, live-damage visualization, and saved damage geometry outputs.

Important folders:

- `Noor mohammad/Mesh/`: main mesh, solver, visualization script, and primary output folder.

### 3. 3D torsion benchmark

The torsion case demonstrates the 3D solver for a notched beam under torsional loading. It includes mesh files, solver script, output curves, snapshots, and an animation.

Important folder:

- `Torsion/working/`: main solver folder and result folders.

## Expected key outputs

| Workflow | Representative outputs |
| --- | --- |
| 3PB MATLAB | `3pb/matlab/Gregoire_3PB/results/matlab_load_cmod.csv`, damage figures, timing log |
| 3PB Abaqus | `3pb/abaqus/Gregoire_3PB/results/abaqus_load_cmod.csv`, damage figures, timing log, `Gregoire_3PB.odb` |
| 3PB comparison | `3pb/comparison/comparison_summary.md`, `comparison_load_cmod.png`, `runtime_comparison.png`, `performance_breakdown.png` |
| Nooru-Mohamed | `Noor mohammad/Mesh/out_NR_vectorized_LIVE_damage_mesh/`, final state, result CSV, damage PNGs |
| Torsion | `Torsion/working/out_torsion_LIVE_ONLY_OLIVER/`, load/torque curves, snapshots, animation |

## Software requirements

Use the tools required by the workflow you want to reproduce:

- MATLAB R2020b or newer recommended.
- Abaqus/CAE and Abaqus/Standard for the Abaqus workflows.
- A Fortran compiler configured with Abaqus for UMAT compilation.
- Python 3.9+ with `numpy` and `matplotlib` for comparison plots outside Abaqus.
- A LaTeX distribution if rebuilding `doc/theory_manual.pdf` from source.
- Git LFS for cloning the repository with Abaqus `.odb` files.

Install Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

## Reproducibility notes

- Generated results are included so reviewers can inspect figures and tables without rerunning every expensive simulation.
- The Abaqus output database is large and stored through Git LFS.
- MATLAB and Abaqus scripts should be run from their own benchmark folders because input and output paths are relative.
- If MATLAB prompts about function names not matching file names, run the command using the file name shown in the tables above.

## Citation

If you use this repository, cite it using `CITATION.cff`. Update `CITATION.cff` with the JOSS DOI after publication.

## License

This repository is released under the MIT License. See `LICENSE`.
