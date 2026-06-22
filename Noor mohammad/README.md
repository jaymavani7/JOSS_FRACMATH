# 3D Nooru-Mohamed benchmark

This folder contains the 3D Nooru-Mohamed validation case used in the FRACMATH paper.

## Folder map

| Path | Purpose |
| --- | --- |
| `Mesh/` | Main mesh files, MATLAB solver, visualization script, and primary output folder |

## Main workflow

Use the instructions in `Mesh/README.md`.

Short version:

1. Open MATLAB.
2. Change the current folder to `Noor mohammad/Mesh`.
3. Optionally run `visulizaiton` to inspect the mesh and boundary conditions.
4. Run `damage_static('Job-1', opts)` with the options shown in `Mesh/README.md`.

## Main outputs

| Path | Contents |
| --- | --- |
| `Mesh/out_NR_vectorized_LIVE_damage_mesh/` | Primary result CSV, final state, and damage images |
