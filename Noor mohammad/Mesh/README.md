# Nooru-Mohamed mesh and solver

This folder contains the main MATLAB workflow for the 3D Nooru-Mohamed benchmark.

## Main files

| File | Purpose |
| --- | --- |
| `damage_static.m` | Main vectorized 3D CDM damage solver |
| `visulizaiton.m` | Mesh and boundary-condition visualization |
| `Job-1_nodes.txt` | Node coordinates |
| `Job-1_elements.txt` | Tetrahedral element connectivity |
| `Job-1_BCs.txt` | Boundary-condition data |
| `Job-1_top_nodes.txt`, `Job-1_bottom_nodes.txt`, `Job-1_left_nodes.txt`, `Job-1_right_nodes.txt` | Boundary node sets |

## Step 1: visualize mesh and boundary conditions

Open MATLAB and change the current folder to:

```text
Noor mohammad/Mesh
```

Run:

```matlab
visulizaiton
```

## Step 2: run the damage analysis

Use this MATLAB command block:

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

## Outputs

The default output folder is:

```text
out_NR_vectorized_LIVE_damage_mesh
```

Representative outputs include:

- `Job-1_fast_NR_results.csv`
- `Job-1_fast_NR_final_state.mat`
- `Job-1_FINAL_damage_geometry_with_mesh.png`
- increment damage images in `increment_damaged_geometry_png/`

## Notes

- The solver uses TET4 elements, modified von Mises equivalent strain, exponential softening, and an Oliver direction-dependent crack-band bandwidth.
- The folder name contains a space, so use quotes if navigating from a terminal.
