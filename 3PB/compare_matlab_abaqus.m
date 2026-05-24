function compare_matlab_abaqus()
% COMPARE_MATLAB_ABAQUS  Overlay MATLAB and Abaqus load-CMOD curves
% and write a side-by-side timing + memory comparison.
%
%   Expects (in ./Gregoire_3PB/results/):
%       matlab_load_cmod.csv     written by solver_main_3pb
%       matlab_timing.txt        written by solver_main_3pb
%       abaqus_load_cmod.csv     written by run_3pb_abaqus.py
%       abaqus_timing.txt        written by run_3pb_abaqus.py
%
%   Outputs:
%       fig_compare_load_cmod.png/pdf
%       comparison_table.txt
%
%   If only the MATLAB files are present, the script still writes a
%   single-curve plot and a one-column timing report.
%
%   Authors: [Name to be added]

clc;
fprintf('==== MATLAB vs Abaqus comparison ====\n');

res_dir = fullfile('Gregoire_3PB', 'results');
if ~exist(res_dir, 'dir')
    error('Folder %s not found.', res_dir);
end

% --- read load-CMOD CSVs ----------------------------------------------
m_csv = fullfile(res_dir, 'matlab_load_cmod.csv');
a_csv = fullfile(res_dir, 'abaqus_load_cmod.csv');

m_data = []; a_data = [];
if exist(m_csv, 'file'), m_data = readmatrix(m_csv); end
if exist(a_csv, 'file'), a_data = readmatrix(a_csv); end

if isempty(m_data) && isempty(a_data)
    error(['Neither matlab_load_cmod.csv nor abaqus_load_cmod.csv ' ...
           'is present in %s.'], res_dir);
end

% --- read timing logs -------------------------------------------------
m_tim = read_timing(fullfile(res_dir, 'matlab_timing.txt'));
a_tim = read_timing(fullfile(res_dir, 'abaqus_timing.txt'));

% =====================================================================
% Plot: overlay
% =====================================================================
fh = figure('Color','w','Position',[100 100 580 420], 'Visible','off');
ax = axes('Parent', fh); hold(ax, 'on'); grid(ax,'on');

leg_h = [];
if ~isempty(m_data)
    h = plot(ax, m_data(:,1), m_data(:,2)/1000, '-', ...
             'Color', [0.10 0.35 0.75], 'LineWidth', 2.0, ...
             'DisplayName', 'MATLAB (this work)');
    leg_h(end+1) = h;
end
if ~isempty(a_data)
    % subsample so the markers don't crowd
    n = size(a_data, 1);
    idx = unique(round(linspace(1, n, min(40, n))));
    h = plot(ax, a_data(idx,1), a_data(idx,2)/1000, 'o', ...
             'Color', [0.90 0.40 0.10], 'MarkerSize', 6, ...
             'MarkerFaceColor','w', 'LineWidth', 1.2, ...
             'DisplayName', 'Abaqus + UMAT');
    leg_h(end+1) = h;
end

xlabel(ax, 'CMOD [mm]', 'FontSize', 11);
ylabel(ax, 'Load [kN]', 'FontSize', 11);
title(ax,  'MATLAB vs Abaqus + UMAT: load vs CMOD');
legend(leg_h, 'Location','northeast');

xmax = 0; ymax = 0;
if ~isempty(m_data); xmax = max(xmax, max(m_data(:,1))); ymax = max(ymax, max(m_data(:,2))); end
if ~isempty(a_data); xmax = max(xmax, max(a_data(:,1))); ymax = max(ymax, max(a_data(:,2))); end
set(ax, 'XLim', [0, xmax*1.05], 'YLim', [0, ymax/1000 * 1.18]);

basepath = fullfile(res_dir, 'fig_compare_load_cmod');
print(fh, [basepath '.png'], '-dpng', '-r200');
print(fh, [basepath '.pdf'], '-dpdf', '-painters');
close(fh);
fprintf('  wrote %s.png + .pdf\n', basepath);

% =====================================================================
% Timing + memory table
% =====================================================================
out_path = fullfile(res_dir, 'comparison_table.txt');
fid = fopen(out_path, 'w');
fprintf(fid, 'MATLAB vs Abaqus + UMAT -- 3PB (Gregoire D=100, a/D=0.2)\n');
fprintf(fid, '%s\n', repmat('=', 1, 60));
fprintf(fid, '%-22s %14s %14s %10s\n', 'Quantity', 'MATLAB', 'Abaqus', 'Ratio');
fprintf(fid, '%s\n', repmat('-', 1, 60));

rows = {
    {'Peak load (N)',       'peak_load',  '%.2f'};
    {'CMOD@peak (mm)',      'cmod_peak',  '%.4f'};
    {'Wall-clock (s)',      'wall',       '%.2f'};
    {'Peak RAM (MB)',       'ram_mb',     '%.1f'};
    {'DOFs',                'dofs',       '%d'};
    {'Load steps',          'steps',      '%d'};
};

for k = 1:numel(rows)
    label = rows{k}{1}; key = rows{k}{2}; fmt = rows{k}{3};
    mv = get_field(m_tim, key);
    av = get_field(a_tim, key);
    mstr = if_num_format(mv, fmt);
    astr = if_num_format(av, fmt);
    rstr = '---';
    if isfinite(mv) && isfinite(av) && mv ~= 0
        rstr = sprintf('%.2fx', av/mv);
    end
    fprintf(fid, '%-22s %14s %14s %10s\n', label, mstr, astr, rstr);
end

fprintf(fid, '%s\n', repmat('-', 1, 60));
fprintf(fid, '\nMATLAB breakdown:\n');
fprintf(fid, '  assembly: %.2f s\n', get_field(m_tim, 't_asm'));
fprintf(fid, '  damage:   %.2f s\n', get_field(m_tim, 't_dam'));
fprintf(fid, '  solve:    %.2f s\n', get_field(m_tim, 't_solve'));
fclose(fid);
fprintf('  wrote %s\n', out_path);

% --- echo to console too -----
type(out_path);

end


% =====================================================================
function s = read_timing(path)
% Parse the timing .txt written by solver_main_3pb or
% run_3pb_abaqus.py. Returns a struct with peak_load, cmod_peak,
% wall, ram_mb, dofs, steps, t_asm, t_dam, t_solve (NaN if missing).

    s = struct('peak_load', NaN, 'cmod_peak', NaN, 'wall', NaN, ...
                'ram_mb', NaN, 'dofs', NaN, 'steps', NaN, ...
                't_asm', NaN, 't_dam', NaN, 't_solve', NaN);
    if ~exist(path, 'file'), return; end
    fid = fopen(path, 'r');
    while ~feof(fid)
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        ln_low = lower(ln);

        if contains(ln_low, 'peak load')
            v = pull_number(ln); if ~isnan(v); s.peak_load = v; end
        elseif contains(ln_low, 'cmod@peak') || ...
               contains(ln_low, 'cmod at peak')
            v = pull_number(ln); if ~isnan(v); s.cmod_peak = v; end
        elseif contains(ln_low, 'wall-clock') || ...
               contains(ln_low, 'wallclock')
            v = pull_number(ln); if ~isnan(v); s.wall = v; end
        elseif contains(ln_low, 'assembly')
            v = pull_number(ln); if ~isnan(v); s.t_asm = v; end
        elseif contains(ln_low, 'damage')
            v = pull_number(ln); if ~isnan(v); s.t_dam = v; end
        elseif contains(ln_low, 'solve')
            v = pull_number(ln); if ~isnan(v); s.t_solve = v; end
        elseif contains(ln_low, 'peak ram')
            v = pull_number(ln); if ~isnan(v); s.ram_mb = v; end
        elseif contains(ln_low, 'mesh') && contains(ln_low, 'dof')
            % "Mesh: 2691 CPS3, 2850 DOFs"
            tok = regexp(ln, '(\d+)\s*DOFs', 'tokens', 'once');
            if ~isempty(tok), s.dofs = str2double(tok{1}); end
        elseif contains(ln_low, 'load steps')
            v = pull_number(ln); if ~isnan(v); s.steps = v; end
        end
    end
    fclose(fid);
end


function v = pull_number(ln)
% Return the first numeric token on the line (ignoring % values).
    parts = regexp(ln, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    if isempty(parts), v = NaN; return; end
    v = str2double(parts{1});
end


function v = get_field(s, key)
    if isfield(s, key), v = s.(key); else, v = NaN; end
end


function s = if_num_format(v, fmt)
    if isnan(v)
        s = '[--]';
    elseif contains(fmt, '%d')
        s = sprintf(fmt, round(v));
    else
        s = sprintf(fmt, v);
    end
end