function solver_main_3pb_OLIVER_bandwidth_vectorized()

clc;
fprintf('==== MATLAB CDM 3PB solver: vectorized modified Newton-Raphson ====\n');

maxNumCompThreads(1);

case_dir = 'Gregoire_3PB';
if ~exist(case_dir, 'dir')
    error('Folder %s not found. Run export_3pb first to generate mesh files.', ...
          case_dir);
end

res_dir = fullfile(case_dir, 'results');
if ~exist(res_dir, 'dir'); mkdir(res_dir); end

p.E         = 37000;
p.nu        = 0.20;
p.t         = 50;
p.ft        = 3.50;
p.fc        = 35.0;
p.GF        = 0.090;
p.OMEGA_MAX = 1 - 1e-12;
p.eps0      = p.ft / p.E;

p.max_disp  = -0.2;
p.num_steps = 10000;
p.tol       = 1e-6;
p.max_iter  = 30;
p.cmod_limit = 0.5;

live_every   = 5;

[nodes, elems, dof] = load_mesh(case_dir);
nN = size(nodes, 1);
nE = size(elems, 1);
fprintf('  mesh: %d nodes, %d CPS3 elements, %d DOFs\n', nN, nE, 2*nN);

fig_mesh(nodes, elems, dof, res_dir);

[B_all, area_v, gradN_all, dof_mat] = precompute_T3(nodes, elems);
D_el = (p.E / (1 - p.nu^2)) * [1 p.nu 0; p.nu 1 0; 0 0 (1-p.nu)/2];
DB   = pagemtimes(D_el, B_all);
Ke0  = pagemtimes(permute(B_all,[2 1 3]), DB) .* ...
       reshape(area_v * p.t, 1, 1, []);
[II, JJ] = sparse_indices(dof_mat);

live = open_live_fig(nodes, elems);

video_path = fullfile(res_dir, 'simulation_video.mp4');
v_writer = VideoWriter(video_path, 'MPEG-4');
v_writer.FrameRate = 10;
v_writer.Quality = 95;
open(v_writer);
fprintf('  Saving video to: %s\n', video_path);

omega = zeros(nE, 1);
kappa = zeros(nE, 1);
u     = zeros(2*nN, 1);

CMOD  = zeros(p.num_steps, 1);
F     = zeros(p.num_steps, 1);
Disp  = zeros(p.num_steps, 1);
Hmean = zeros(p.num_steps, 1);
Hmin  = zeros(p.num_steps, 1);
Hmax  = zeros(p.num_steps, 1);
steps_done = 0;

K_next = [];

t_asm = 0; t_dam = 0; t_solve = 0; t_viz = 0;

peak_load_so_far = 0;
snap_peak = struct('u',[],'omega',[],'load',0);
snap_pp   = struct('u',[],'omega',[],'load',0);

ram_before = ram_bytes();
wall0 = tic;

fprintf('  step    u(mm)   load(kN)    CMOD(mm)   dmax    iter\n');

for step = 1:p.num_steps
    u_tgt = step * (p.max_disp / p.num_steps);
    u(dof.prescribed) = u_tgt;

    tic;
    if isempty(K_next)
        K = assemble_K(Ke0, omega, II, JJ, nN);
    else
        K = K_next;
        K_next = [];
    end

    Kff  = K(dof.free, dof.free);
    Kfac = factor_free_stiffness(Kff);
    t_asm = t_asm + toc;

    r = zeros(2*nN, 1);
    for it = 1:p.max_iter
        tic;
        r = K * u;
        r_free = r(dof.free);
        if it == 1
            nrm0 = max(norm(r_free), 1e-12);
        end
        if norm(r_free) / nrm0 < p.tol
            t_solve = t_solve + toc;
            break;
        end
        du = -(Kfac \ r_free);
        u(dof.free)       = u(dof.free) + du;
        u(dof.prescribed) = u_tgt;
        t_solve = t_solve + toc;
    end

    tic;
    [omega, kappa, h_oliver] = damage_update(u, B_all, gradN_all, dof_mat, kappa, omega, p);
    t_dam = t_dam + toc;

    tic;
    K_next  = assemble_K(Ke0, omega, II, JJ, nN);
    r_react = K_next * u;
    t_asm = t_asm + toc;

    F_now = -sum(r_react(dof.prescribed));
    C_now = max(0, mean(u(2*dof.cmod2-1)) - mean(u(2*dof.cmod1-1)));
    d_max = max(omega);

    steps_done   = step;
    CMOD(step)   = C_now;
    F(step)      = F_now;
    Disp(step)   = -u_tgt;
    Hmean(step)  = mean(h_oliver);
    Hmin(step)   = min(h_oliver);
    Hmax(step)   = max(h_oliver);

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

    if mod(step, max(1, round(p.num_steps/20))) == 0 || step <= 3
        fprintf('  %4d  %7.4f  %9.3f  %10.5f  %6.3f  %4d\n', ...
                step, -u_tgt, F_now/1000, C_now, d_max, it);
    end

    if mod(step, live_every) == 0 || step == 1 || step == p.num_steps
        tv = tic;
        update_live_fig(live, nodes, u, omega, CMOD(1:step).', F(1:step).', step, ...
                          F_now, C_now, d_max);
        frame = getframe(live.fh);
        writeVideo(v_writer, frame);
        t_viz = t_viz + toc(tv);
    end

    if C_now >= p.cmod_limit
        fprintf('  CMOD limit reached at step %d\n', step);
        break;
    end
end

t_total  = toc(wall0);
t_solver = t_total - t_viz;

CMOD  = CMOD(1:steps_done).';
F     = F(1:steps_done).';
Disp  = Disp(1:steps_done).';
Hmean = Hmean(1:steps_done).';
Hmin  = Hmin(1:steps_done).';
Hmax  = Hmax(1:steps_done).';

close(v_writer);
fprintf('  Video saved successfully.\n');
ram_after  = ram_bytes();
peak_ram_MB = max(0, (ram_after - ram_before)) / 2^20;
if peak_ram_MB == 0
    try, m = memory; peak_ram_MB = m.MemUsedMATLAB/2^20; catch; peak_ram_MB = NaN; end
end

if isempty(snap_pp.u)
    snap_pp.u = u; snap_pp.omega = omega; snap_pp.load = F(end);
end

[pk_load, ip] = max(F);
pk_cmod = CMOD(ip);

csv_path = fullfile(res_dir, 'matlab_load_cmod.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, '# cmod[mm], load[N]\n');
fprintf(fid, '%.6e, %.6e\n', [CMOD(:), F(:)].');
fclose(fid);
fprintf('  wrote %s\n', csv_path);

hcsv_path = fullfile(res_dir, 'matlab_oliver_bandwidth_history.csv');
fid = fopen(hcsv_path, 'w');
fprintf(fid, '# step, mean_h_oliver[mm], min_h_oliver[mm], max_h_oliver[mm]\n');
fprintf(fid, '%d, %.6e, %.6e, %.6e\n', [(1:numel(Hmean)).', Hmean(:), Hmin(:), Hmax(:)].');
fclose(fid);
fprintf('  wrote %s\n', hcsv_path);

t_path = fullfile(res_dir, 'matlab_timing.txt');
fid = fopen(t_path, 'w');
fprintf(fid, 'MATLAB 3PB solver\n');
fprintf(fid, 'Peak load:    %.2f N\n',    pk_load);
fprintf(fid, 'CMOD@peak:    %.6f mm\n',   pk_cmod);
fprintf(fid, 'Solver wall-clock: %.2f s   <-- compare to Abaqus\n', t_solver);
fprintf(fid, 'End-to-end:        %.2f s   (incl. live fig + video)\n', t_total);
fprintf(fid, 'Visualization:     %.2f s\n', t_viz);
fprintf(fid, '  assembly:   %.2f s  (%.1f%%)\n', t_asm,   100*t_asm/t_solver);
fprintf(fid, '  damage:     %.2f s  (%.1f%%)\n', t_dam,   100*t_dam/t_solver);
fprintf(fid, '  solve:      %.2f s  (%.1f%%)\n', t_solve, 100*t_solve/t_solver);
fprintf(fid, 'Peak RAM:     %.1f MB\n',   peak_ram_MB);
fprintf(fid, 'Mesh:         %d CPS3, %d DOFs\n', nE, 2*nN);
fprintf(fid, 'Load steps:   %d\n',        numel(F));
if ~isempty(Hmean)
    fprintf(fid, 'Oliver h mean final: %.6e mm\n', Hmean(end));
    fprintf(fid, 'Oliver h min final:  %.6e mm\n', Hmin(end));
    fprintf(fid, 'Oliver h max final:  %.6e mm\n', Hmax(end));
end
fclose(fid);
fprintf('  wrote %s\n', t_path);

fig_load_cmod(CMOD, F, pk_load, pk_cmod, res_dir);

fig_damage(nodes, elems, snap_peak.u, snap_peak.omega, ...
           sprintf('Peak load: %.2f kN', snap_peak.load/1000), ...
           fullfile(res_dir, 'fig_damage_peak'));

fig_damage(nodes, elems, snap_pp.u, snap_pp.omega, ...
           sprintf('Post-peak: %.2f kN', snap_pp.load/1000), ...
           fullfile(res_dir, 'fig_damage_postpeak'));

fig_damage(nodes, elems, u, omega, ...
           sprintf('Last step: %.2f kN', F(end)/1000), ...
           fullfile(res_dir, 'fig_damage_last_step'));

fprintf('\n==== summary ====\n');
fprintf('  peak load : %.2f N at CMOD = %.4f mm\n', pk_load, pk_cmod);
fprintf('  wall-clock: %.1f s  (asm %.1f, dam %.1f, solve %.1f)\n', ...
        t_total, t_asm, t_dam, t_solve);
fprintf('  peak RAM  : %.1f MB\n', peak_ram_MB);
fprintf('  output -> %s/\n', res_dir);

end

function live = open_live_fig(nodes, elems)

    FULL_DAMAGE_THRESH = 0.99;

    fh = figure('Name','CDM 3PB – Live', ...
                'Color','w', ...
                'Position',[60 60 1100 420], ...
                'NumberTitle','off');

    ax1 = subplot(1,2,1,'Parent',fh);
    hold(ax1,'on');

    ph_bg = patch('Parent',ax1, ...
                  'Faces',elems, ...
                  'Vertices',nodes, ...
                  'FaceColor',[0.92 0.92 0.93], ...
                  'EdgeColor',[0.78 0.80 0.82], ...
                  'LineWidth',0.10);

    ph_fd = patch('Parent',ax1, ...
                  'Faces',zeros(0,3), ...
                  'Vertices',nodes, ...
                  'FaceColor',[0.0 0.0 0.0], ...
                  'EdgeColor',[0.0 0.0 0.0], ...
                  'LineWidth',0.20);

    axis(ax1,'equal','tight');
    xlabel(ax1,'x [mm]');
    ylabel(ax1,'y [mm]');
    title(ax1,sprintf('Fully damaged elements only  (\omega \ge %.2f)', FULL_DAMAGE_THRESH), ...
          'FontSize',10);
    set(ax1,'XLim',[min(nodes(:,1))-10, max(nodes(:,1))+10], ...
            'YLim',[-5, max(nodes(:,2))+5]);

    ax2 = subplot(1,2,2,'Parent',fh);
    hold(ax2,'on');
    grid(ax2,'on');
    lh  = plot(ax2, NaN, NaN, '-', 'Color',[0.10 0.35 0.75], 'LineWidth',2);
    mh  = plot(ax2, NaN, NaN, 'o', 'Color','r', ...
               'MarkerFaceColor','r', 'MarkerSize',8);
    xlabel(ax2,'CMOD [mm]');
    ylabel(ax2,'Load [kN]');
    title(ax2,'Load vs CMOD  (live)', 'FontSize',10);

    xlim(ax2,[0 0.35]);
    ylim(ax2,[0 6]);

    live.fh    = fh;
    live.ax1   = ax1;
    live.ax2   = ax2;
    live.ph_bg = ph_bg;
    live.ph_fd = ph_fd;
    live.lh    = lh;
    live.mh    = mh;
    live.nodes = nodes;
    live.elems = elems;
    live.full_damage_thresh = FULL_DAMAGE_THRESH;

    drawnow limitrate;
end

function update_live_fig(live, nodes, u, omega, CMOD, F, step, ...
                          F_now, C_now, d_max)

    if ~ishandle(live.fh), return; end

    THRESH_FULL = live.full_damage_thresh;
    full_idx = find(omega >= THRESH_FULL);

    if isempty(full_idx)
        set(live.ph_fd, 'Faces', zeros(0,3), 'Vertices', nodes);
    else
        set(live.ph_fd, 'Faces', live.elems(full_idx,:), 'Vertices', nodes);
    end

    set(live.lh, 'XData', CMOD,  'YData', F/1000);
    set(live.mh, 'XData', C_now, 'YData', F_now/1000);

    if max(F/1000) * 1.15 > live.ax2.YLim(2)
        ylim(live.ax2, [0, max(F/1000)*1.25]);
    end

    live.fh.Name = sprintf( ...
        'CDM 3PB | step %d | load %.2f kN | CMOD %.4f mm | dmax %.3f | omega>=%.2f elems %d', ...
        step, F_now/1000, C_now, d_max, THRESH_FULL, numel(full_idx));

    drawnow limitrate;
end

function omega_node = elem2node_avg(omega_e, elems, nN)

    ids = elems(:);
    val = repelem(omega_e(:), 3);

    acc = accumarray(ids, val, [nN 1], @sum, 0);
    cnt = accumarray(ids, 1,   [nN 1], @sum, 0);

    omega_node = acc ./ max(cnt, 1);
end

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

    fix_left  = reshape([2*left(:)-1, 2*left(:)].', [], 1);
    fix_right = 2*right(:);
    fix       = [fix_left; fix_right];

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

function [B_all, area_v, gradN_all, dof_mat] = precompute_T3(nodes, elems)
    x1 = nodes(elems(:,1),1);  y1 = nodes(elems(:,1),2);
    x2 = nodes(elems(:,2),1);  y2 = nodes(elems(:,2),2);
    x3 = nodes(elems(:,3),1);  y3 = nodes(elems(:,3),2);

    area_signed = 0.5*((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));
    area_v = abs(area_signed);
    if any(area_v <= 1e-14)
        bad = find(area_v <= 1e-14, 1, 'first');
        error('Zero or near-zero CPS3 area at element %d.', bad);
    end

    b1 = y2-y3; b2 = y3-y1; b3 = y1-y2;
    c1 = x3-x2; c2 = x1-x3; c3 = x2-x1;
    inv2A_signed = 1 ./ (2*area_signed);

    g1x = b1 .* inv2A_signed;  g1y = c1 .* inv2A_signed;
    g2x = b2 .* inv2A_signed;  g2y = c2 .* inv2A_signed;
    g3x = b3 .* inv2A_signed;  g3y = c3 .* inv2A_signed;

    nE = numel(area_v);
    B_all = zeros(3, 6, nE);
    B_all(1,1,:) = g1x;  B_all(1,3,:) = g2x;  B_all(1,5,:) = g3x;
    B_all(2,2,:) = g1y;  B_all(2,4,:) = g2y;  B_all(2,6,:) = g3y;
    B_all(3,1,:) = g1y;  B_all(3,2,:) = g1x;  B_all(3,3,:) = g2y;
    B_all(3,4,:) = g2x;  B_all(3,5,:) = g3y;  B_all(3,6,:) = g3x;

    gradN_all = [g1x, g1y, g2x, g2y, g3x, g3y];

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

function Kfac = factor_free_stiffness(Kff)

    try
        Kfac = decomposition(Kff, 'chol');
    catch
        Kfac = decomposition(Kff, 'lu');
    end
end

function [omega_new, kappa_new, h_oliver] = damage_update(u, B_all, gradN_all, dof_mat, ...
                                                kappa_old, omega_old, p)
    nE   = size(B_all,3);
    u_e  = reshape(u(dof_mat).', 6, 1, nE);
    strain = squeeze(pagemtimes(B_all, u_e)).';

    ex = strain(:,1); ey = strain(:,2); gxy = strain(:,3);
    me  = (ex+ey)/2;
    rad = sqrt(((ex-ey)/2).^2 + (gxy/2).^2);
    e1  = me + rad;  e2 = me - rad;

    theta_p = 0.5 * atan2(gxy, ex - ey);
    nx = cos(theta_p);
    ny = sin(theta_p);

    iso = abs(ex-ey) + abs(gxy) < 1e-18;
    nx(iso) = 1.0;
    ny(iso) = 0.0;

    h_oliver = oliver_bandwidth_T3(gradN_all, nx, ny);

    ef_e = max(p.eps0/2 + p.GF ./ (h_oliver * p.ft), p.eps0 + 1e-12);

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

function h = oliver_bandwidth_T3(gradN_all, nx, ny)

    g1x = gradN_all(:,1); g1y = gradN_all(:,2);
    g2x = gradN_all(:,3); g2y = gradN_all(:,4);
    g3x = gradN_all(:,5); g3y = gradN_all(:,6);

    den = abs(g1x.*nx + g1y.*ny) + ...
          abs(g2x.*nx + g2y.*ny) + ...
          abs(g3x.*nx + g3y.*ny);

    h = 2.0 ./ max(den, 1e-14);
    h = max(h, 1e-12);
end

function fig_damage(nodes, elems, u, omega, label, basepath)

    deformed = nodes;

    x_min = min(nodes(:,1));  x_max = max(nodes(:,1));
    y_min = min(nodes(:,2));  y_max = max(nodes(:,2));

    THRESH_FULL = 0.99;
    full_idx    = find(omega >= THRESH_FULL);

    fh = figure('Color','w', ...
                'Position',[100 100 1380 420], ...
                'Visible','off', ...
                'Renderer','painters', ...
                'InvertHardcopy','off', ...
                'PaperUnits','centimeters', ...
                'PaperSize',[25.5 7.8], ...
                'PaperPosition',[0 0 25.5 7.8]);

    ax = axes('Parent',fh, ...
              'Units','normalized', ...
              'Position',[0.065 0.18 0.760 0.73]);
    hold(ax,'on');

    patch('Parent',ax, ...
          'Faces',elems, ...
          'Vertices',deformed, ...
          'FaceColor',[0.955 0.955 0.965], ...
          'EdgeColor',[0.86 0.875 0.89], ...
          'LineWidth',0.045);

    if ~isempty(full_idx)
        patch('Parent',ax, ...
              'Faces',elems(full_idx,:), ...
              'Vertices',deformed, ...
              'FaceColor',[0.00 0.00 0.00], ...
              'EdgeColor',[0.00 0.00 0.00], ...
              'LineWidth',0.12);
    end

    axis(ax,'equal');
    ax.Box       = 'off';
    ax.LineWidth = 0.75;
    ax.FontSize  = 10;
    ax.TickDir   = 'out';
    ax.Layer     = 'top';

    xlim(ax, [x_min - 5, x_max + 5]);
    ylim(ax, [y_min - 5, y_max + 5]);

    xlabel(ax, '$x$ [mm]', 'Interpreter','latex', 'FontSize',11);
    ylabel(ax, '$y$ [mm]', 'Interpreter','latex', 'FontSize',11);
    title(ax, 'Fully damaged elements only', ...
          'Interpreter','latex', ...
          'FontSize',13, ...
          'FontWeight','bold');

    axp = axes('Parent',fh, ...
               'Units','normalized', ...
               'Position',[0.850 0.19 0.135 0.70]);
    hold(axp,'on');
    axis(axp,[0 1 0 1]);
    axis(axp,'off');

    text(axp,0.00,0.98,'Legend', ...
         'FontSize',11.5, ...
         'FontWeight','bold', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'Interpreter','none', ...
         'Clipping','off');

    plot(axp,[0.02 0.22],[0.83 0.83],'-', ...
         'Color',[0.86 0.875 0.89], ...
         'LineWidth',8);
    text(axp,0.30,0.83,'FE mesh', ...
         'FontSize',10.5, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle', ...
         'Interpreter','none', ...
         'Clipping','off');

    plot(axp,[0.02 0.22],[0.70 0.70],'-', ...
         'Color',[0.00 0.00 0.00], ...
         'LineWidth',8);
    text(axp,0.30,0.70,sprintf('$\\omega \\ge %.2f$',THRESH_FULL), ...
         'FontSize',10.5, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','middle', ...
         'Interpreter','latex', ...
         'Clipping','off');

    text(axp,0.00,0.51,'Result', ...
         'FontSize',11.0, ...
         'FontWeight','bold', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'Interpreter','none', ...
         'Clipping','off');

    text(axp,0.00,0.40,label, ...
         'FontSize',10.2, ...
         'FontWeight','bold', ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'Interpreter','none', ...
         'Clipping','off');

    text(axp,0.00,0.29,sprintf('Threshold: \\omega >= %.2f',THRESH_FULL), ...
         'FontSize',9.8, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'Interpreter','tex', ...
         'Clipping','off');

    text(axp,0.00,0.19,sprintf('Fully damaged: %d elements',numel(full_idx)), ...
         'FontSize',9.8, ...
         'HorizontalAlignment','left', ...
         'VerticalAlignment','top', ...
         'Interpreter','none', ...
         'Clipping','off');

    save_fig_hq(fh, basepath);
    close(fh);
end

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

    fill(ax, [CMOD, fliplr(CMOD)], [F/1000, zeros(1,numel(F))], ...
         [0.08 0.30 0.72], 'FaceAlpha',0.10, 'EdgeColor','none');

    plot(ax, CMOD, F/1000, '-', ...
         'Color',[0.08 0.30 0.72], 'LineWidth',2.0);

    plot(ax, pk_cmod, pk_load/1000, 'pentagram', ...
         'MarkerSize',13, ...
         'MarkerFaceColor',[0.98 0.82 0.0], ...
         'MarkerEdgeColor',[0.40 0.28 0.0], ...
         'LineWidth',1.0);

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

    xlim(ax, [0, 0.35]);
    ylim(ax, [0, pk_load/1000 * 1.28]);

    save_fig_hq(fh, fullfile(res_dir, 'matlab_load_cmod_fig'));
    close(fh);
end

function fig_mesh(nodes, elems, dof, res_dir)

    x_min = min(nodes(:,1));  x_max = max(nodes(:,1));
    y_min = min(nodes(:,2));  y_max = max(nodes(:,2));
    span  = x_max - x_min;
    ht    = y_max - y_min;

    fh = figure('Color','w', ...
                'Position',[100 100 800 320], ...
                'Visible','off', ...
                'PaperUnits','centimeters', ...
                'PaperSize',[18.0 7.2], ...
                'PaperPosition',[0 0 18.0 7.2]);

    ax = axes('Parent',fh,'Units','normalized','Position',[0.08 0.18 0.80 0.74]);
    hold(ax,'on');

    patch('Parent',ax,'Faces',elems,'Vertices',nodes, ...
          'FaceColor',[0.93 0.93 0.96], ...
          'EdgeColor',[0.50 0.58 0.72],'LineWidth',0.20);

    n_fix = ceil(dof.fixed/2);
    cnt   = accumarray(n_fix(:), 1, [size(nodes,1),1]);
    for ni = find(cnt >= 2)'
        plot_support(ax, nodes(ni,1), y_min, 6, 'pin');
    end
    for ni = find(cnt == 1)'
        plot_support(ax, nodes(ni,1), y_min, 6, 'roller');
    end

    load_nodes = unique(ceil(dof.prescribed/2));
    xC = mean(nodes(load_nodes,1));
    yT = max(nodes(load_nodes,2));
    quiver(ax, xC, yT+22, 0, -18, 0, ...
           'Color',[0.0 0.52 0.38],'LineWidth',2.2, ...
           'MaxHeadSize',0.55,'AutoScale','off');
    text(ax, xC+9, yT+15, '$P$', ...
         'Interpreter','latex','FontSize',12, ...
         'Color',[0.0 0.52 0.38],'FontWeight','bold');

    ann_y = y_min - 13;
    draw_dim_arrow(ax, x_min, ann_y, x_max, ann_y, [0.22 0.22 0.22], 'horizontal');
    text(ax, (x_min+x_max)/2, ann_y - 3, ...
         sprintf('$L_{\\rm span} = %.0f$ mm', span), ...
         'Interpreter','latex','FontSize',8, ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'Color',[0.20 0.20 0.20]);

    rr_arrow = x_max + span*0.04;
    rr_text  = rr_arrow + span*0.03;

    draw_dim_arrow(ax, rr_arrow, y_min, rr_arrow, y_max, [0.22 0.22 0.22], 'vertical');

    text(ax, rr_text, (y_min+y_max)/2, ...
         sprintf('$H = %.0f$ mm', ht), ...
         'Interpreter','latex','FontSize',9,'Color',[0.20 0.20 0.20], ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'Rotation',90);

    axis(ax,'equal');
    ax.Box      = 'on';
    ax.LineWidth = 0.7;
    ax.FontSize  = 8;
    ax.TickDir   = 'out';

    xlim(ax,[x_min-18, x_max+45]);
    ylim(ax,[y_min-25, y_max+30]);

    xlabel(ax,'$x$ [mm]','Interpreter','latex','FontSize',9);
    ylabel(ax,'$y$ [mm]','Interpreter','latex','FontSize',9);
    title(ax,'FE mesh — boundary conditions and load point', ...
          'FontSize',10,'FontWeight','bold');

    save_fig_hq(fh, fullfile(res_dir,'fig_mesh'));
    close(fh);
end

function c = crack_cmap()

    n     = 256;
    stops = [0.84 0.88 0.95;
             0.18 0.42 0.86;
             0.05 0.72 0.88;
             0.18 0.80 0.32;
             0.96 0.90 0.08;
             0.98 0.44 0.04;
             0.82 0.04 0.04];
    pos   = [0, 0.10, 0.30, 0.50, 0.70, 0.85, 1.00];
    t     = linspace(0,1,n)';
    c     = zeros(n,3);
    for ch = 1:3
        c(:,ch) = interp1(pos, stops(:,ch), t, 'pchip');
    end
    c = min(max(c,0),1);
end

function c = damage_cmap(), c = crack_cmap(); end

function draw_dim_arrow(ax, x1, y1, x2, y2, col, dir)

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

    set(fh,'Color','w','InvertHardcopy','off');

    try
        exportgraphics(fh, [basepath '.png'], ...
                       'Resolution',600, ...
                       'BackgroundColor','white');
        exportgraphics(fh, [basepath '.pdf'], ...
                       'ContentType','vector', ...
                       'BackgroundColor','white');
    catch

        set(fh,'PaperPositionMode','auto');
        print(fh, [basepath '.png'], '-dpng', '-r600');
        print(fh, [basepath '.pdf'], '-dpdf', '-painters', '-r600');
    end
end

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
