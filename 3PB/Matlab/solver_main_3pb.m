function solver_main_3pb()
% SOLVER_MAIN_3PB  Single-file MATLAB solver for the notched 3PB case.
%
%   Reads the .txt mesh files written by export_3pb.m (or make_mesh_3pb.py),
%   runs the vectorized CDM solver (modified von Mises eq strain, exponential
%   softening, crack-band regularization, modified Newton-Raphson with
%   secant tangent refactored once per load step), saves all plots
%   and the timing log.
%
%   Requires:
%       Gregoire_3PB/nodes.txt, elements.txt, top_nodes.txt,
%       left_nodes.txt, right_nodes.txt, cmod1.txt, cmod2.txt
%
%   Outputs (in Gregoire_3PB/results/):
%       matlab_load_cmod.csv              CMOD[mm], Load[N]
%       matlab_timing.txt                 peak load + breakdown + RAM
%       fig_mesh.png/pdf                  mesh + BCs + load (visualization)
%       matlab_load_cmod_fig.png/pdf      Load vs CMOD curve
%       fig_damage_peak.png/pdf           damage contour at peak load
%       fig_damage_postpeak.png/pdf       damage contour late post-peak
%       simulation_video.mp4              Video of the simulation process

clc;
fprintf('==== MATLAB CDM 3PB solver ====\n');

% --- where the .txt files are ------------------------------------------
case_dir = 'Gregoire_3PB';
if ~exist(case_dir, 'dir')
    error('Folder %s not found. Run export_3pb first to generate mesh files.', ...
          case_dir);
end

res_dir = fullfile(case_dir, 'results');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end

% --- material + solver parameters --------------------------------------
p.E         = 37000;       % MPa
p.nu        = 0.20;
p.t         = 50;          % mm  (thickness)
p.ft        = 3.50;        % MPa
p.fc        = 35.0;        % MPa
p.GF        = 0.090;       % N/mm
p.OMEGA_MAX = 1 - 1e-12;
p.eps0      = p.ft / p.E;

p.max_disp  = -0.2;        % mm  total midspan deflection
p.num_steps = 10000;
p.tol       = 1e-6;
p.max_iter  = 30;
p.cmod_limit = 0.5;

live_every   = 5;          % refresh live window every N steps

% =====================================================================
% 1. Load the mesh
% =====================================================================
[nodes, elems, dof] = load_mesh(case_dir);
nN = size(nodes, 1);
nE = size(elems, 1);
fprintf('  mesh: %d nodes, %d CPS3 elements, %d DOFs\n', nN, nE, 2*nN);

% --- save static mesh figure -------------------------------------------
fig_mesh(nodes, elems, dof, res_dir);

% =====================================================================
% 2. Pre-compute element data
% =====================================================================
[B_all, area_v, hband_v, dof_mat] = precompute_T3(nodes, elems);
D_el = (p.E / (1 - p.nu^2)) * [1 p.nu 0; p.nu 1 0; 0 0 (1-p.nu)/2];
DB   = pagemtimes(D_el, B_all);
Ke0  = pagemtimes(permute(B_all,[2 1 3]), DB) .* ...
       reshape(area_v * p.t, 1, 1, []);
[II, JJ] = sparse_indices(dof_mat);
ef_e = max(p.eps0/2 + p.GF ./ (hband_v * p.ft), p.eps0 + 1e-12);

% =====================================================================
% 3. Open the LIVE figure & Setup Video Writer
% =====================================================================
live = open_live_fig(nodes, elems);

% Setup Video Writer
video_path = fullfile(res_dir, 'simulation_video.mp4');
v_writer = VideoWriter(video_path, 'MPEG-4');
v_writer.FrameRate = 10; % Frames per second (adjust to speed up/slow down)
v_writer.Quality = 95;
open(v_writer);
fprintf('  Saving video to: %s\n', video_path);

% =====================================================================
% 4. Solve (vectorised modified Newton-Raphson, secant tangent)
% =====================================================================
omega = zeros(nE, 1);
kappa = zeros(nE, 1);
u     = zeros(2*nN, 1);

CMOD = []; F = []; Disp = [];
t_asm = 0; t_dam = 0; t_solve = 0;

peak_load_so_far = 0;
snap_peak = struct('u',[],'omega',[],'load',0);
snap_pp   = struct('u',[],'omega',[],'load',0);

ram_before = ram_bytes();
wall0 = tic;

fprintf('  step    u(mm)   load(kN)    CMOD(mm)   dmax    iter\n');

for step = 1:p.num_steps
    u_tgt = step * (p.max_disp / p.num_steps);
    u(dof.prescribed) = u_tgt;

    % --- assemble secant stiffness at start of step --------------------
    tic;
    K = assemble_K(Ke0, omega, II, JJ, nN);
    t_asm = t_asm + toc;

    r = zeros(2*nN, 1);
    for it = 1:p.max_iter
        tic;
        r = K * u;
        r_free = r(dof.free);
        if it == 1, nrm0 = max(norm(r), 1e-12); end
        if norm(r_free) / nrm0 < p.tol, break; end
        du = K(dof.free, dof.free) \ (-r_free);
        u(dof.free)       = u(dof.free) + du;
        u(dof.prescribed) = u_tgt;
        t_solve = t_solve + toc;
    end

    % --- damage update -------------------------------------------------
    tic;
    [omega, kappa] = damage_update(u, B_all, dof_mat, kappa, omega, ef_e, p);
    t_dam = t_dam + toc;

    F_now = -sum(r(dof.prescribed));
    C_now = max(0, mean(u(2*dof.cmod2-1)) - mean(u(2*dof.cmod1-1)));
    d_max = max(omega);

    CMOD(end+1) = C_now;   %#ok<AGROW>
    F(end+1)    = F_now;   %#ok<AGROW>
    Disp(end+1) = -u_tgt;  %#ok<AGROW>

    % --- snapshots for final figures -----------------------------------
    if F_now > peak_load_so_far
        peak_load_so_far  = F_now;
        snap_peak.u       = u;
        snap_peak.omega   = omega;
        snap_peak.load    = F_now;
    end
    if C_now > 0.3 && isempty(snap_pp.u)
        snap_pp.u     = u;
        snap_pp.omega = omega;
        snap_pp.load  = F_now;
    end

    % --- console progress ----------------------------------------------
    if mod(step, max(1, round(p.num_steps/20))) == 0 || step <= 3
        fprintf('  %4d  %7.4f  %9.3f  %10.5f  %6.3f  %4d\n', ...
                step, -u_tgt, F_now/1000, C_now, d_max, it);
    end

    % --- update live figure AND WRITE VIDEO FRAME -----------------------
    if mod(step, live_every) == 0 || step == 1 || step == p.num_steps
        update_live_fig(live, nodes, u, omega, CMOD, F, step, ...
                          F_now, C_now, d_max);
        
        % Capture the current figure and write to video
        frame = getframe(live.fh);
        writeVideo(v_writer, frame);
    end

    if C_now >= p.cmod_limit
        fprintf('  CMOD limit reached at step %d\n', step);
        break;
    end
end

% Close the video writer
close(v_writer);
fprintf('  Video saved successfully.\n');

t_total = toc(wall0);
ram_after  = ram_bytes();
peak_ram_MB = max(0, (ram_after - ram_before)) / 2^20;
if peak_ram_MB == 0
    try, m = memory; peak_ram_MB = m.MemUsedMATLAB/2^20; catch; peak_ram_MB = NaN; end
end

if isempty(snap_pp.u)
    snap_pp.u = u; snap_pp.omega = omega; snap_pp.load = F(end);
end

% =====================================================================
% 5. Save results
% =====================================================================
[pk_load, ip] = max(F);
pk_cmod = CMOD(ip);

csv_path = fullfile(res_dir, 'matlab_load_cmod.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, '# cmod[mm], load[N]\n');
for i = 1:numel(F)
    fprintf(fid, '%.6e, %.6e\n', CMOD(i), F(i));
end
fclose(fid);
fprintf('  wrote %s\n', csv_path);

t_path = fullfile(res_dir, 'matlab_timing.txt');
fid = fopen(t_path, 'w');
fprintf(fid, 'MATLAB 3PB solver\n');
fprintf(fid, 'Peak load:    %.2f N\n',    pk_load);
fprintf(fid, 'CMOD@peak:    %.6f mm\n',   pk_cmod);
fprintf(fid, 'Wall-clock:   %.2f s\n',    t_total);
fprintf(fid, '  assembly:   %.2f s  (%.1f%%)\n', t_asm,   100*t_asm/t_total);
fprintf(fid, '  damage:     %.2f s  (%.1f%%)\n', t_dam,   100*t_dam/t_total);
fprintf(fid, '  solve:      %.2f s  (%.1f%%)\n', t_solve, 100*t_solve/t_total);
fprintf(fid, 'Peak RAM:     %.1f MB\n',   peak_ram_MB);
fprintf(fid, 'Mesh:         %d CPS3, %d DOFs\n', nE, 2*nN);
fprintf(fid, 'Load steps:   %d\n',        numel(F));
fclose(fid);
fprintf('  wrote %s\n', t_path);

% =====================================================================
% 6. Final static figures
% =====================================================================
fig_load_cmod(CMOD, F, pk_load, pk_cmod, res_dir);

fig_damage(nodes, elems, snap_peak.u, snap_peak.omega, ...
           sprintf('Peak load: %.2f kN', snap_peak.load/1000), ...
           fullfile(res_dir, 'fig_damage_peak'));

fig_damage(nodes, elems, snap_pp.u, snap_pp.omega, ...
           sprintf('Post-peak: %.2f kN', snap_pp.load/1000), ...
           fullfile(res_dir, 'fig_damage_postpeak'));

fprintf('\n==== summary ====\n');
fprintf('  peak load : %.2f N at CMOD = %.4f mm\n', pk_load, pk_cmod);
fprintf('  wall-clock: %.1f s  (asm %.1f, dam %.1f, solve %.1f)\n', ...
        t_total, t_asm, t_dam, t_solve);
fprintf('  peak RAM  : %.1f MB\n', peak_ram_MB);
fprintf('  output -> %s/\n', res_dir);

end % solver_main_3pb


% =====================================================================
%                       LIVE FIGURE HELPERS
% =====================================================================
function live = open_live_fig(nodes, elems)
    fh = figure('Name','CDM 3PB – Live', ...
                'Color','w', ...
                'Position',[60 60 1100 420], ...
                'NumberTitle','off');

    ax1 = subplot(1,2,1,'Parent',fh);
    ph  = patch('Parent',ax1, ...
                'Faces',elems, 'Vertices',nodes, ...
                'FaceVertexCData', zeros(size(nodes,1),1), ...
                'FaceColor','interp', 'EdgeColor','none');
    clim(ax1, [0 1]);
    colormap(ax1, crack_cmap());
    cb1 = colorbar(ax1);
    cb1.Label.String = '\omega  (damage)';
    axis(ax1,'equal','tight');
    xlabel(ax1,'x [mm]'); ylabel(ax1,'y [mm]');
    title(ax1,'Damage contour  (deformed \times scale)', 'FontSize',10);
    set(ax1,'XLim',[min(nodes(:,1))-10, max(nodes(:,1))+10], ...
            'YLim',[-5, max(nodes(:,2))+5]);

    ax2 = subplot(1,2,2,'Parent',fh);
    hold(ax2,'on');
    grid(ax2,'on');
    lh  = plot(ax2, NaN, NaN, '-', 'Color',[0.10 0.35 0.75], 'LineWidth',2);
    mh  = plot(ax2, NaN, NaN, 'o', 'Color','r', ...
               'MarkerFaceColor','r', 'MarkerSize',8);
    xlabel(ax2,'CMOD [mm]'); ylabel(ax2,'Load [kN]');
    title(ax2,'Load vs CMOD  (live)', 'FontSize',10);
    
    % Set initial live limit here
    xlim(ax2,[0 0.35]); 
    ylim(ax2,[0 6]);

    live.fh  = fh;
    live.ax1 = ax1;
    live.ph  = ph;
    live.ax2 = ax2;
    live.lh  = lh;
    live.mh  = mh;
    live.nodes = nodes;

    drawnow limitrate;
end

function update_live_fig(live, nodes, u, omega, CMOD, F, step, ...
                          F_now, C_now, d_max)
    if ~ishandle(live.fh), return; end

    nN = size(nodes,1);

    Uxy   = reshape(u, 2, []).';
    umax  = max(abs(Uxy(:)));
    scale = min(1, max(5, 3 / max(umax, 1e-12)));
    deformed = nodes + scale * Uxy;

    omega_node = elem2node_avg(omega, live.ph.Faces, nN);

    set(live.ph, 'Vertices', deformed, 'FaceVertexCData', omega_node);
    set(live.lh, 'XData', CMOD,  'YData', F/1000);
    set(live.mh, 'XData', C_now, 'YData', F_now/1000);

    % Only auto-scale the Y-axis (Load), leave X-axis locked at 0.35
    if max(F/1000) * 1.15 > live.ax2.YLim(2)
        ylim(live.ax2, [0, max(F/1000)*1.25]);
    end

    live.fh.Name = sprintf( ...
        'CDM 3PB  |  step %d  |  load %.2f kN  |  CMOD %.4f mm  |  dmax %.3f', ...
        step, F_now/1000, C_now, d_max);

    drawnow limitrate;
end

function omega_node = elem2node_avg(omega_e, elems, nN)
    acc = zeros(nN, 1);
    cnt = zeros(nN, 1);
    for k = 1:3
        acc = acc + accumarray(elems(:,k), omega_e, [nN 1]);
        cnt = cnt + accumarray(elems(:,k), ones(size(elems,1),1), [nN 1]);
    end
    omega_node = acc ./ max(cnt, 1);
end


% =====================================================================
%                    MESH LOADING
% =====================================================================
function [nodes, elems, dof] = load_mesh(d)
    raw_n = load(fullfile(d, 'nodes.txt'));
    raw_e = load(fullfile(d, 'elements.txt'));

    map = zeros(max(raw_n(:,1)), 1);
    map(raw_n(:,1)) = 1:size(raw_n,1);

    nodes = raw_n(:, 2:3);
    elems = [map(raw_e(:,2)), map(raw_e(:,3)), map(raw_e(:,4))];

    top   = map(load_id(d,'top_nodes.txt'));
    left  = map(load_id(d,'left_nodes.txt'));
    right = map(load_id(d,'right_nodes.txt'));
    c1    = map(load_id(d,'cmod1.txt'));
    c2    = map(load_id(d,'cmod2.txt'));

    nN = size(nodes,1);
    fix = [];
    for n = left(:)';  fix = [fix; 2*n-1; 2*n]; end  %#ok<AGROW>
    for n = right(:)'; fix = [fix; 2*n];        end  %#ok<AGROW>
    pres    = 2 * top(:);
    all_dof = (1:2*nN)';

    dof.fixed      = unique(fix);
    dof.prescribed = pres;
    dof.free       = setdiff(all_dof, unique([dof.fixed; pres]));
    dof.cmod1      = c1;
    dof.cmod2      = c2;
end


function v = load_id(d, fname)
    f = fullfile(d, fname);
    v = [];
    if exist(f,'file')
        try, v = load(f); catch; end
        v = v(:);
    end
end


% =====================================================================
%                    ELEMENT PRE-COMPUTATION
% =====================================================================
function [B_all, area_v, hband_v, dof_mat] = precompute_T3(nodes, elems)
    x1 = nodes(elems(:,1),1);  y1 = nodes(elems(:,1),2);
    x2 = nodes(elems(:,2),1);  y2 = nodes(elems(:,2),2);
    x3 = nodes(elems(:,3),1);  y3 = nodes(elems(:,3),2);

    area_v = 0.5*abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));
    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;
    inv2A = 1 ./ (2*area_v);

    nE = numel(area_v);
    B_all = zeros(3, 6, nE);
    B_all(1,1,:) = b1.*inv2A;  B_all(1,3,:) = b2.*inv2A;  B_all(1,5,:) = b3.*inv2A;
    B_all(2,2,:) = c1.*inv2A;  B_all(2,4,:) = c2.*inv2A;  B_all(2,6,:) = c3.*inv2A;
    B_all(3,1,:) = c1.*inv2A;  B_all(3,2,:) = b1.*inv2A;  B_all(3,3,:) = c2.*inv2A;
    B_all(3,4,:) = b2.*inv2A;  B_all(3,5,:) = c3.*inv2A;  B_all(3,6,:) = b3.*inv2A;

    hband_v = sqrt(area_v);

    dof_mat = [2*elems(:,1)-1, 2*elems(:,1), ...
               2*elems(:,2)-1, 2*elems(:,2), ...
               2*elems(:,3)-1, 2*elems(:,3)];
end


function [II, JJ] = sparse_indices(dof_mat)
    [lr, lc] = ndgrid(1:6, 1:6);
    II = dof_mat(:, lr(:));
    JJ = dof_mat(:, lc(:));
end


function K = assemble_K(Ke0, omega, II, JJ, nN)
    nE = size(Ke0,3);
    Ke = Ke0 .* reshape(1 - omega, 1, 1, nE);
    V  = reshape(Ke, 36, nE).';
    K  = sparse(II(:), JJ(:), V(:), 2*nN, 2*nN);
end


% =====================================================================
%                        DAMAGE UPDATE
% =====================================================================
function [omega_new, kappa_new] = damage_update(u, B_all, dof_mat, ...
                                                kappa_old, omega_old, ...
                                                ef_e, p)
    nE   = size(B_all,3);
    u_e  = reshape(u(dof_mat).', 6, 1, nE);
    strain = squeeze(pagemtimes(B_all, u_e)).';

    ex = strain(:,1); ey = strain(:,2); gxy = strain(:,3);
    me  = (ex+ey)/2;
    rad = sqrt(((ex-ey)/2).^2 + (gxy/2).^2);
    e1  = me + rad;  e2 = me - rad;
    e3  = -(p.nu/(1-p.nu)) .* (e1+e2);
    I1  = e1+e2+e3;
    J2  = 0.5*((e1-e2).^2 + (e2-e3).^2 + (e3-e1).^2);

    k  = p.fc/p.ft;
    a1 = (k-1)/(2*k*(1-2*p.nu));
    a2 = 1/(2*k);
    a3 = ((k-1)/(1-2*p.nu))^2;
    a4 = 12*k/(1+p.nu)^2;
    eq_s = a1*I1 + a2*sqrt(max(a3*I1.^2 + a4*J2, 0));
    eq_s = max(eq_s, 0);

    kappa_new = max(kappa_old, eq_s);
    omega_new = zeros(nE,1);
    m = kappa_new > p.eps0;
    if any(m)
        km  = kappa_new(m);
        den = max(ef_e(m) - p.eps0, 1e-15);
        omega_new(m) = 1 - (p.eps0 ./ km) .* exp(-(km - p.eps0) ./ den);
    end
    omega_new = max(omega_new, omega_old);
    omega_new = min(max(omega_new,0), p.OMEGA_MAX);
    bad = ~isfinite(omega_new);
    omega_new(bad) = omega_old(bad);
end


% =====================================================================
%          PUBLICATION-QUALITY STATIC FIGURES
% =====================================================================

% ---------------------------------------------------------------------
%  DAMAGE / CRACK FIGURE
% ---------------------------------------------------------------------
function fig_damage(nodes, elems, u, omega, label, basepath)
% Publication-quality crack visualization without dimension lines.

    % ---- deformed coordinates ----------------------------------------
    Uxy   = reshape(u, 2, []).';
    umax  = max(abs(Uxy(:)));
    scale = min(1, max(30, 2.0 / max(umax, 1e-12)));
    deformed = nodes + scale * Uxy;

    % ---- specimen geometry -------------------------------------------
    x_min = min(nodes(:,1));  x_max = max(nodes(:,1));
    y_min = min(nodes(:,2));  y_max = max(nodes(:,2));

    % ---- figure setup ------------------------------------------------
    fh = figure('Color','w', ...
                'Position',[100 100 900 400], ...
                'Visible','off', ...
                'PaperUnits','centimeters', ...
                'PaperSize', [18 8], ...
                'PaperPosition',[0 0 18 8]);

    % Adjusted position to use more space now that dimensions are gone
    ax = axes('Parent',fh, 'Units','normalized', 'Position',[0.10 0.15 0.75 0.75]);
    hold(ax,'on');

    % ---- layer 1: full mesh in light grey ----------------------------
    patch('Parent',ax, ...
          'Faces',    elems, ...
          'Vertices', deformed, ...
          'FaceColor',[0.92 0.92 0.93], ...
          'EdgeColor',[0.75 0.77 0.80], ...
          'LineWidth',0.15);

    % ---- layer 2: crack zone colored by damage ----------------------
    THRESH    = 0.01;
    crack_idx = find(omega > THRESH);
    if ~isempty(crack_idx)
        patch('Parent',ax, ...
              'Faces',    elems(crack_idx,:), ...
              'Vertices', deformed, ...
              'CData',    omega(crack_idx), ...
              'FaceColor','flat', ...
              'EdgeColor','none');
    end

    % ---- colormap & colorbar ----------------------------------------
    colormap(ax, crack_cmap());
    clim(ax, [0 1]);
    cb = colorbar(ax, 'Location','eastoutside');
    cb.Label.String  = '\omega  (damage)';
    cb.Label.FontSize = 10;

    % ---- axes styling ------------------------------------------------
    axis(ax,'equal','tight');
    ax.Box      = 'on';
    ax.LineWidth = 0.8;
    ax.FontSize  = 10;
    
    % Tighten limits now that dimension arrows are removed
    xlim(ax, [x_min - 5, x_max + 5]);
    ylim(ax, [y_min - 5, y_max + 5]);

    xlabel(ax, '$x$ [mm]', 'Interpreter','latex');
    ylabel(ax, '$y$ [mm]', 'Interpreter','latex');

    % ---- label box (Top-Left) ----------------------------------------
    text(ax, 0.03, 0.94, label, ...
         'Units','normalized', ...
         'Interpreter','none', ...
         'FontSize', 11, 'FontWeight', 'bold', ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'BackgroundColor', 'w', 'EdgeColor', [0.2 0.2 0.2], ...
         'Margin', 6, 'LineWidth', 0.8);

    % ---- export ------------------------------------------------------
    save_fig_hq(fh, basepath);
    close(fh);
end

% ---------------------------------------------------------------------
%  LOAD vs CMOD FIGURE
% ---------------------------------------------------------------------
function fig_load_cmod(CMOD, F, pk_load, pk_cmod, res_dir)

    fh = figure('Color','w', ...
                'Position',[100 100 500 390], ...
                'Visible','off', ...
                'PaperUnits','centimeters', ...
                'PaperSize',[8.8 7.0], ...
                'PaperPosition',[0 0 8.8 7.0]);

    ax = axes('Parent',fh, ...
              'Units','normalized','Position',[0.14 0.14 0.82 0.78]);
    hold(ax,'on');
    grid(ax,'on');
    ax.GridColor      = [0.80 0.80 0.80];
    ax.GridLineStyle  = ':';
    ax.FontSize       = 9;
    ax.LineWidth      = 0.7;
    ax.Box            = 'on';

    % shaded fill under curve
    fill(ax, [CMOD, fliplr(CMOD)], [F/1000, zeros(1,numel(F))], ...
         [0.08 0.30 0.72], 'FaceAlpha',0.10, 'EdgeColor','none');

    % main curve
    plot(ax, CMOD, F/1000, '-', ...
         'Color',[0.08 0.30 0.72], 'LineWidth',2.0);

    % peak star marker
    plot(ax, pk_cmod, pk_load/1000, 'pentagram', ...
         'MarkerSize',13, ...
         'MarkerFaceColor',[0.98 0.82 0.0], ...
         'MarkerEdgeColor',[0.40 0.28 0.0], ...
         'LineWidth',1.0);

    % peak label — positioned to avoid overlap within the 0.35 limit
    lbl_x  = pk_cmod + 0.35 * 0.03; 
    lbl_y  = pk_load/1000 * 1.04;
    ha_str = 'left';
    if lbl_x > 0.35 * 0.72
        lbl_x  = pk_cmod - 0.35 * 0.03;
        ha_str = 'right';
    end
    text(ax, lbl_x, lbl_y, ...
         sprintf('$P_{\\rm peak} = %.2f$ kN\nCMOD $= %.4f$ mm', ...
                 pk_load/1000, pk_cmod), ...
         'Interpreter','latex','FontSize',8.5,'FontWeight','bold', ...
         'Color',[0.30 0.20 0.0], ...
         'HorizontalAlignment',ha_str,'VerticalAlignment','bottom', ...
         'BackgroundColor','w','EdgeColor',[0.60 0.50 0.20], ...
         'Margin',3,'LineWidth',0.5);

    xlabel(ax, 'CMOD [mm]',  'Interpreter','latex','FontSize',10);
    ylabel(ax, 'Load [kN]',  'Interpreter','latex','FontSize',10);
    title(ax,  'Load--CMOD response', 'FontSize',10,'FontWeight','bold');

    % Set exact x-axis limit here
    xlim(ax, [0, 0.35]);
    ylim(ax, [0, pk_load/1000 * 1.28]);

    save_fig_hq(fh, fullfile(res_dir, 'matlab_load_cmod_fig'));
    close(fh);
end

% ---------------------------------------------------------------------
%  MESH FIGURE
% ---------------------------------------------------------------------
function fig_mesh(nodes, elems, dof, res_dir)

    x_min = min(nodes(:,1));  x_max = max(nodes(:,1));
    y_min = min(nodes(:,2));  y_max = max(nodes(:,2));
    span  = x_max - x_min;
    ht    = y_max - y_min;

    % 1. Increase the figure width to create physical space for the label
    fh = figure('Color','w', ...
                'Position',[100 100 800 320], ... 
                'Visible','off', ...
                'PaperUnits','centimeters', ...
                'PaperSize',[18.0 7.2], ...
                'PaperPosition',[0 0 18.0 7.2]);

    % 2. Constrain the axes area so it doesn't expand into the label space
    ax = axes('Parent',fh,'Units','normalized','Position',[0.08 0.18 0.80 0.74]);
    hold(ax,'on');

    % --- Mesh Visualization ---
    patch('Parent',ax,'Faces',elems,'Vertices',nodes, ...
          'FaceColor',[0.93 0.93 0.96], ...
          'EdgeColor',[0.50 0.58 0.72],'LineWidth',0.20);

    % --- Supports ---
    n_fix = ceil(dof.fixed/2);
    cnt   = accumarray(n_fix(:), 1, [size(nodes,1),1]);
    for ni = find(cnt >= 2)'
        plot_support(ax, nodes(ni,1), y_min, 6, 'pin');
    end
    for ni = find(cnt == 1)'
        plot_support(ax, nodes(ni,1), y_min, 6, 'roller');
    end

    % --- Load Arrow ---
    load_nodes = unique(ceil(dof.prescribed/2));
    xC = mean(nodes(load_nodes,1));
    yT = max(nodes(load_nodes,2));
    quiver(ax, xC, yT+22, 0, -18, 0, ...
           'Color',[0.0 0.52 0.38],'LineWidth',2.2, ...
           'MaxHeadSize',0.55,'AutoScale','off');
    text(ax, xC+9, yT+15, '$P$', ...
         'Interpreter','latex','FontSize',12, ...
         'Color',[0.0 0.52 0.38],'FontWeight','bold');

    % --- Dimensions ---
    % Span (Bottom)
    ann_y = y_min - 13;
    draw_dim_arrow(ax, x_min, ann_y, x_max, ann_y, [0.22 0.22 0.22], 'horizontal');
    text(ax, (x_min+x_max)/2, ann_y - 3, ...
         sprintf('$L_{\\rm span} = %.0f$ mm', span), ...
         'Interpreter','latex','FontSize',8, ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'Color',[0.20 0.20 0.20]);

    % Height (Right-side) - FIXED POSITIONING
    rr_arrow = x_max + span*0.04; 
    rr_text  = rr_arrow + span*0.03;
    
    draw_dim_arrow(ax, rr_arrow, y_min, rr_arrow, y_max, [0.22 0.22 0.22], 'vertical');
    
    text(ax, rr_text, (y_min+y_max)/2, ...
         sprintf('$H = %.0f$ mm', ht), ...
         'Interpreter','latex','FontSize',9,'Color',[0.20 0.20 0.20], ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'Rotation',90);

    % --- Limits ---
    axis(ax,'equal');
    ax.Box      = 'on';
    ax.LineWidth = 0.7;
    ax.FontSize  = 8;
    ax.TickDir   = 'out';
    
    % Ensure limits contain all nodes AND all annotation positions
    xlim(ax,[x_min-18, x_max+45]); 
    ylim(ax,[y_min-25, y_max+30]);

    xlabel(ax,'$x$ [mm]','Interpreter','latex','FontSize',9);
    ylabel(ax,'$y$ [mm]','Interpreter','latex','FontSize',9);
    title(ax,'FE mesh — boundary conditions and load point', ...
          'FontSize',10,'FontWeight','bold');

    save_fig_hq(fh, fullfile(res_dir,'fig_mesh'));
    close(fh);
end

% =====================================================================
%  COLORMAP  — blue -> cyan -> green -> yellow -> red  (crack style)
% =====================================================================
function c = crack_cmap()
% Spectral crack colormap matching reference publication images.
% Undamaged background stays cool blue-grey; crack tip is dark red.
    n     = 256;
    stops = [0.84 0.88 0.95;   % 0.00  light blue-grey (undamaged)
             0.18 0.42 0.86;   % 0.10  blue
             0.05 0.72 0.88;   % 0.30  cyan
             0.18 0.80 0.32;   % 0.50  green
             0.96 0.90 0.08;   % 0.70  yellow
             0.98 0.44 0.04;   % 0.85  orange
             0.82 0.04 0.04];  % 1.00  dark red
    pos   = [0, 0.10, 0.30, 0.50, 0.70, 0.85, 1.00];
    t     = linspace(0,1,n)';
    c     = zeros(n,3);
    for ch = 1:3
        c(:,ch) = interp1(pos, stops(:,ch), t, 'pchip');
    end
    c = min(max(c,0),1);
end

% Alias kept for any legacy calls
function c = damage_cmap(), c = crack_cmap(); end


% =====================================================================
%  SHARED HELPERS
% =====================================================================
function draw_dim_arrow(ax, x1, y1, x2, y2, col, dir)
% Double-headed dimension arrow with proper directional triangular markers
    plot(ax, [x1 x2], [y1 y2], '-', 'Color',col, 'LineWidth',1.0);
    
    if strcmp(dir, 'horizontal')
        plot(ax, x1, y1, '<', 'Color',col, 'MarkerFaceColor',col, 'MarkerSize',5);
        plot(ax, x2, y2, '>', 'Color',col, 'MarkerFaceColor',col, 'MarkerSize',5);
    elseif strcmp(dir, 'vertical')
        plot(ax, x1, y1, 'v', 'Color',col, 'MarkerFaceColor',col, 'MarkerSize',5);
        plot(ax, x2, y2, '^', 'Color',col, 'MarkerFaceColor',col, 'MarkerSize',5);
    end
end


function plot_support(ax, x, y, s, kind)
    patch('Parent',ax, ...
          'XData',[x-s, x+s, x], 'YData',[y-s, y-s, y], ...
          'FaceColor',[0.32 0.32 0.38],'EdgeColor','k','LineWidth',0.7);
    if strcmp(kind,'roller')
        % small circles under roller
        th = linspace(0,2*pi,24);
        r  = s * 0.45;
        for dx = [-s*0.6, s*0.6]
            fill(ax, x+dx + r*cos(th), y-s-r + r*sin(th), ...
                 [0.50 0.50 0.55],'EdgeColor','k','LineWidth',0.5);
        end
        plot(ax, [x-s*1.4, x+s*1.4],[y-s-2*r-0.5, y-s-2*r-0.5], ...
             '-k','LineWidth',1.0);
    end
end


function save_fig_hq(fh, basepath)
% Export 600 dpi PNG and vector-quality PDF (painters renderer).
    print(fh, [basepath '.png'], '-dpng',     '-r600');
    print(fh, [basepath '.pdf'], '-dpdf',     '-painters', '-r600');
end


% =====================================================================
%                        SYSTEM HELPERS
% =====================================================================
function r = ram_bytes()
    r = 0;
    if isunix
        try
            [~,txt] = system(sprintf('grep VmPeak /proc/%d/status', feature('getpid')));
            tok = regexp(txt,'(\d+)','tokens','once');
            if ~isempty(tok), r = str2double(tok{1})*1024; end
        catch
        end
    elseif ispc
        try, m = memory; r = m.MemUsedMATLAB; catch; end
    end
end