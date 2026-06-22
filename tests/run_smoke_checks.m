function run_smoke_checks()

root = fileparts(fileparts(mfilename('fullpath')));

requiredFiles = {
    'README.md'
    'LICENSE'
    'CITATION.cff'
    'CONTRIBUTING.md'
    'CHANGELOG.md'
    'REPRODUCIBILITY.md'
    'paper/paper.md'
    'paper/paper.bib'
    'paper/images/fig_mesh.png'
    'paper/images/fig_b1_mesh_abq.png'
    'paper/images/fig_b1_damage.png'
    'paper/images/fig_b1_results.png'
    'paper/images/fig_b1_compact.png'
    'paper/images/fig_b2_mesh.png'
    'paper/images/nooru_damage_evolution_3x3.png'
    'paper/images/nooru_damage_evolution_selected.png'
    'paper/images/fig_b3_mesh.png'
    'paper/images/fig_b3_damage_evolution.png'
    'paper/images/torsion_max_damage.png'
    'paper/images/torsion_peak_damage_iso.png'
    'paper/images/load_cmod_comparison.png'
    'paper/images/nooru_BC_2D.png'
    'doc/theory_manual.pdf'
    '3pb/matlab/solver_main_3pb.m'
    '3pb/matlab/Gregoire_3PB/nodes.txt'
    '3pb/matlab/Gregoire_3PB/elements.txt'
    '3pb/matlab/Gregoire_3PB/results/matlab_load_cmod.csv'
    '3pb/abaqus/cdm_umat_2d_OLIVER_T3_FAST.for'
    '3pb/comparison/plot_comparison.py'
    'Noor mohammad/Mesh/damage_static.m'
    'Noor mohammad/Mesh/Job-1_nodes.txt'
    'Torsion/working/run_torsion.m'
    'Torsion/working/Job-1_nodes.txt'
    };

for i = 1:numel(requiredFiles)
    path = fullfile(root, requiredFiles{i});
    assert(isfile(path), 'Missing required file: %s', requiredFiles{i});
end

nodes = readmatrix(fullfile(root, '3pb/matlab/Gregoire_3PB/nodes.txt'));
elements = readmatrix(fullfile(root, '3pb/matlab/Gregoire_3PB/elements.txt'));
loadCmod = readmatrix(fullfile(root, '3pb/matlab/Gregoire_3PB/results/matlab_load_cmod.csv'), 'NumHeaderLines', 1);

assert(size(nodes, 1) > 0 && size(nodes, 2) >= 3, '3PB node table has an unexpected shape.');
assert(size(elements, 1) > 0 && size(elements, 2) >= 4, '3PB element table has an unexpected shape.');
assert(size(loadCmod, 1) > 0 && size(loadCmod, 2) >= 2, '3PB load-CMOD table has an unexpected shape.');
assert(max(loadCmod(:, 2)) > 0, '3PB load-CMOD table does not contain positive load values.');

fprintf('FRACMATH smoke checks passed.\n');
end
