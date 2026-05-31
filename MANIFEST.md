# Repository manifest

This file lists the main source-code, input-data, and generated-result locations for the FRACMAT JOSS repository.

## Root files

| File | Purpose |
| --- | --- |
| `README.md` | Main repository overview and quick start |
| `REPRODUCIBILITY.md` | Step-by-step run guide for each analysis |
| `CITATION.cff` | Citation metadata for the software archive |
| `LICENSE` | MIT License |
| `requirements.txt` | Python plotting dependencies |
| `.gitattributes` | Git LFS tracking for large Abaqus `.odb` files |
| `.gitignore` | Local scratch-file ignore rules |

## 2D three-point bending benchmark

| Path | Type | Description |
| --- | --- | --- |
| `3pb/README.md` | documentation | Overview of the complete 3PB benchmark |
| `3pb/matlab/solver_main_3pb.m` | source | Main MATLAB 2D continuum damage solver |
| `3pb/matlab/make_3pb_inp.py` | source | Helper for Abaqus input/export workflows |
| `3pb/matlab/export.m` | source | MATLAB export helper |
| `3pb/matlab/Gregoire_3PB/*.txt` | input | Mesh, node sets, CMOD nodes, and boundary-condition data |
| `3pb/matlab/Gregoire_3PB/results/` | output | MATLAB curves, timing log, figures, and video |
| `3pb/abaqus/run_3pb_abaqus_OLIVER_T3_FAST.py` | source | Abaqus build-run-extract-plot workflow |
| `3pb/abaqus/cdm_umat_2d_OLIVER_T3_FAST.for` | source | Abaqus UMAT implementation |
| `3pb/abaqus/extract_damage.py` | source | Abaqus ODB damage extraction helper |
| `3pb/abaqus/extract_peak_omega.py` | source | Peak-damage extraction helper |
| `3pb/abaqus/Gregoire_3PB/` | input/output | Abaqus model files, ODB, logs, extracted results |
| `3pb/comparison/plot_comparison.py` | source | Load-CMOD and runtime comparison plotting |
| `3pb/comparison/plot_image_comparison.py` | source | Image-panel comparison plotting |
| `3pb/comparison/comparison_summary.md` | output | MATLAB-vs-Abaqus summary table |

## 3D Nooru-Mohamed benchmark

| Path | Type | Description |
| --- | --- | --- |
| `Noor mohammad/README.md` | documentation | Overview of the Nooru-Mohamed benchmark |
| `Noor mohammad/Mesh/README.md` | documentation | One-by-one run instructions |
| `Noor mohammad/Mesh/damage_static_NR_vectorized_LIVE_damage_different_colours.m` | source | Main 3D vectorized MATLAB damage solver |
| `Noor mohammad/Mesh/visulizaiton.m` | source | Mesh and boundary-condition visualization |
| `Noor mohammad/Mesh/Job-1_*.txt` | input | Mesh and node-set files |
| `Noor mohammad/Mesh/out_NR_vectorized_LIVE_damage_mesh/` | output | Main solver output folder |
| `Noor mohammad/out_ultra/` | output | Additional stored damage snapshots and curve data |

## 3D torsion benchmark

| Path | Type | Description |
| --- | --- | --- |
| `Torsion/README.md` | documentation | Overview of the torsion benchmark |
| `Torsion/working/README.md` | documentation | One-by-one run instructions |
| `Torsion/working/run_torsion_stress_strain_static_fast_modvm_multiview_video.m` | source | Main 3D torsion MATLAB solver |
| `Torsion/working/visulize.m` | source | Visualization helper |
| `Torsion/working/Job-1_*.txt` | input | Mesh and node-set files |
| `Torsion/working/out_torsion_LIVE_ONLY_OLIVER/` | output | Main torsion output folder |
| `Torsion/visulize mesh/` | source/input | Mesh visualization cases and Abaqus model files |

## Documentation

| Path | Type | Description |
| --- | --- | --- |
| `doc/README.md` | documentation | Theory manual instructions |
| `doc/theory_manual.tex` | source | LaTeX source for the theory manual |
| `doc/theory_manual.pdf` | output | Compiled theory manual |
