function damage_static_NR_ultra(prefix, opts)
    % Static solve via FSOLVE on free DOFs (secant damage tangent) - CPU only.
    %
    % Mesh inputs required in working folder:
    %   <prefix>_nodes.txt, _elements.txt, _top_nodes.txt, _bottom_nodes.txt,
    %   _left_nodes.txt, _right_nodes.txt
    %
    % Live visualization (during run):
    %   Tiles 1-3 : 3D internal-crack views from three angles (deformed mesh)
    %   Tile 4    : Combined Ps vs delta_s and P vs delta
    %
    % Publication-quality outputs (after run) -- matches solver_main_3pb style:
    %   <prefix>_fig_mesh.png/pdf            FE mesh + boundary conditions (3-view)
    %   <prefix>_load_disp.png/pdf           Load-displacement curve with peak markers
    %   <prefix>_damage_peak.png/pdf         Damage at peak load (3-view)
    %   <prefix>_damage_postpeak.png/pdf     Damage post-peak (3-view)
    %
    % Other outputs:
    %   <prefix>_curve_4c.csv                Numerical results
    %   <prefix>_simulation.mp4              Run video
    %   <prefix>_final_rotation.mp4          Final rotation inspection video
    %   <prefix>_crack_view_<deg>_deg.png    Snapshots every 45 degrees

    if nargin < 1 || isempty(prefix), prefix = 'Job-1'; end
    if nargin < 2 || isempty(opts),   opts   = struct();   end
    get_opt = @(f,def) local_get(opts,f,def);

    % ===== Material (Nooru-Mohamed paper values) =====
    E      = get_opt('E',    29000.0);   % Young's modulus [MPa]
    nu     = get_opt('nu',   0.20);      % Poisson ratio
    GF     = get_opt('GF',   0.11);      % Fracture energy [N/mm]  (Gf = 110 J/m2)
    ft     = get_opt('ft',   3.00);      % Tensile strength [MPa]
    k_tc   = get_opt('k',    10.0);      % Tension/compression ratio
    kappa0 = get_opt('kappa0', ft/E);    % Damage threshold (default = ft/E)
    j0     = kappa0;

    % ===== Loading controls =====
    nIncr     = get_opt('nIncr',       500);
    Uy_end    = get_opt('Uy_end',      0.50);  % Total normal (vertical) displacement [mm]
    gamma     = get_opt('gamma',       0.60);  % Shear/normal coupling: delta_s = gamma * delta
    load_path = get_opt('load_path',  '4c');   % '4c' or 'shear_force'
    Fx_tot    = get_opt('Fx_tot',      1.0e4); % Applied shear force for 'shear_force' path [N]
    tol       = get_opt('tol',         1e-6);
    maxEval   = get_opt('maxFunEvals', 4000);
    maxIter   = get_opt('maxIter',     400);

    % ===== 3D Visualization Controls =====
    def_scale    = get_opt('def_scale', 15.0);   % Scale factor for 3D deformation visual
    crack_thresh = get_opt('crack_thresh', 0.80);% Show elements with damage > 80%
    hull_alpha   = get_opt('hull_alpha', 0.15);  % Transparency of undamaged outer mesh

    % ===== Output Directory =====
    out_dir = 'out_ultra';
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    % ===== Mesh & node sets =====
    nd  = readmatrix([prefix '_nodes.txt']); ids = nd(:,1);
    p   = nd(:,2:4);  np = size(p,1);
    ed  = readmatrix([prefix '_elements.txt']);
    [~,T] = ismember(ed(:,2:5), ids);  ne = size(T,1);

    Ntop   = uint32(readmatrix([prefix '_top_nodes.txt']));
    Nbot   = uint32(readmatrix([prefix '_bottom_nodes.txt']));
    Nleft  = uint32(readmatrix([prefix '_left_nodes.txt']));
    Nright = uint32(readmatrix([prefix '_right_nodes.txt']));

    % ===== Precompute (B, V, he) and unit Ke =====
    [B3, V3, he3, ~] = precompute_TET4_vector(p, T);
    [Dunit, ~, ~]    = iso3D_D(1.0, nu);
    Ke_unit3 = pagemtimes(permute(B3,[2 1 3]), pagemtimes(repmat(Dunit,[1 1 ne]), B3));
    Ke_unit3 = bsxfun(@times, Ke_unit3, reshape(V3,1,1,ne));   % 12x12xne

    % ===== Element DOFs =====
    ndof = 3*np;
    EDOF = zeros(12, ne, 'uint32');
    for e = 1:ne
        n4 = T(e,:);
        EDOF(:,e) = uint32([3*n4(1)-2:3*n4(1), 3*n4(2)-2:3*n4(2), ...
                            3*n4(3)-2:3*n4(3), 3*n4(4)-2:3*n4(4)]).';
    end
    idx_glob = double(EDOF(:));

    % ===== Prebuild stiffness triplets (sparsity pattern) =====
    nnz_per_e = 12*12;
    I_trip    = zeros(ne*nnz_per_e, 1, 'uint32');
    J_trip    = zeros(ne*nnz_per_e, 1, 'uint32');
    Val_unit  = zeros(ne*nnz_per_e, 1);
    eid_rep   = zeros(ne*nnz_per_e, 1, 'uint32');

    base = 0;
    for e = 1:ne
        dofs = double(EDOF(:,e));
        [jj,ii] = meshgrid(dofs, dofs);
        blk  = Ke_unit3(:,:,e);
        kN   = numel(blk);
        I_trip(base+(1:kN)) = uint32(ii(:));
        J_trip(base+(1:kN)) = uint32(jj(:));
        Val_unit(base+(1:kN)) = blk(:);
        eid_rep(base+(1:kN))  = uint32(e);
        base = base + kN;
    end
    Kpat = sparse(double(I_trip), double(J_trip), double(Val_unit~=0), ndof, ndof);

    % ===== Damage state variables =====
    U  = zeros(ndof,1);
    be = (E * j0 ./ GF) .* he3(:);        % exponential softening parameter per element
    kappa = j0 * ones(ne,1);              % irreversible history, updated each increment

    % ===== Boundary condition helpers =====
    dof         = @(nid,comp) 3*(double(nid)-1)+comp;
    fixY_bottom = dof(Nbot, 2);
    presc_top_uy= dof(Ntop, 2);

    xb = p(double(Nbot),1); [~,ib] = min(abs(xb - mean(xb)));
    fixX_pin = dof(Nbot(ib), 1);
    [~, iPinZ] = min(p(:,3));  fixZ_pin = 3*(iPinZ-1)+3;

    base_fix = unique([fixY_bottom; fixX_pin; fixZ_pin]);
    Lx = dof(Nleft,  1);
    Rx = dof(Nright, 1);

    % External forces template (for 'shear_force' path only)
    Fext_base = sparse(ndof,1);
    if strcmpi(load_path,'shear_force')
        if ~isempty(Nleft),  Fext_base(Lx) = +Fx_tot / max(numel(Nleft),1);  end
        if ~isempty(Nright), Fext_base(Rx) = -Fx_tot / max(numel(Nright),1); end
    end

    % ===== FSOLVE options =====
    base_opts = optimoptions('fsolve', ...
        'SpecifyObjectiveGradient', true, ...
        'FunctionTolerance',        tol, ...
        'StepTolerance',            max(tol^1.5, 1e-12), ...
        'MaxFunctionEvaluations',   maxEval, ...
        'MaxIterations',            maxIter, ...
        'Display',                  'off');

    % ===== Results arrays =====
    delta_s_res = zeros(nIncr,1);   % shear displacement delta_s [mm]
    delta_res   = zeros(nIncr,1);   % normal displacement delta    [mm]
    Ps_res      = zeros(nIncr,1);   % shear force Ps               [N]
    P_res       = zeros(nIncr,1);   % normal force P               [N]

    % ===== Snapshots for publication-quality damage figures =====
    peak_metric_so_far = 0;
    snap_peak = struct('U',[], 'De',[], 'P',0, 'Ps',0, 'delta',0, 'deltas',0, 'step',0);
    snap_pp   = struct('U',[], 'De',[], 'P',0, 'Ps',0, 'delta',0, 'deltas',0, 'step',0);

    % ===== Precompute Faces for 3D Visualization =====
    fprintf('Extracting external faces for live 3D visualization...\n');
    allF = [T(:,[1 2 3]); T(:,[1 2 4]); T(:,[1 3 4]); T(:,[2 3 4])];
    owner_per_face = repmat((1:ne)', 4, 1);

    [sF, ~] = sort(allF, 2);
    [uF, ~, iC] = unique(sF, 'rows');
    counts = accumarray(iC, 1);
    extFaces = uF(counts == 1, :);

    % ===== Figure: 2x2 Layout for Multiple Views and Combined Plot =====
    fig = figure('Color','w','Name','Nooru-Mohamed path 4c: Multi-View Live', ...
                 'Position',[50 50 1400 900]);
    tl  = tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');

    % Setup Arrays to hold handles for the 3 different 3D views
    view_angles = [45 30; -45 30; 0 90]; % [Iso 1, Iso 2, Top View]
    ax3d   = gobjects(3,1);
    hHull  = gobjects(3,1);
    hCrack = gobjects(3,1);

    % Initialize the 3 visualizer tiles
    for v = 1:3
        ax3d(v) = nexttile(tl, v);
        hold(ax3d(v),'on'); axis(ax3d(v),'equal', 'vis3d');
        view(ax3d(v), view_angles(v,:));
        colormap(ax3d(v), crack_cmap());   % publication crack colormap

        hHull(v) = patch(ax3d(v), 'Faces', extFaces, 'Vertices', p, ...
                      'FaceColor', [0.85 0.86 0.90], 'EdgeColor', 'none', ...
                      'FaceAlpha', hull_alpha, 'AmbientStrength', 0.6);

        hCrack(v) = patch(ax3d(v), 'Faces', [], 'Vertices', p, ...
                       'FaceVertexCData', [], 'FaceColor', 'flat', ...
                       'EdgeColor', 'none', 'AmbientStrength', 0.8);

        camlight(ax3d(v), 'headlight'); lighting(ax3d(v), 'gouraud');
        title(ax3d(v), sprintf('View %d  (\\omega > %.2f)', v, crack_thresh), ...
              'FontSize', 10, 'FontWeight', 'bold');
        try, clim(ax3d(v),[0 1]); catch, caxis(ax3d(v),[0 1]); end %#ok<NOSEMI>
        grid(ax3d(v), 'on'); box(ax3d(v), 'on');
    end
    colorbar(ax3d(3), 'Location', 'eastoutside');

    % Tile 4 - Combined Load vs Displacement Plot
    ax4 = nexttile(tl, 4);
    hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
    ax4.GridColor = [0.80 0.80 0.80]; ax4.GridLineStyle = ':';
    xlabel(ax4,'Displacement [mm]','FontSize',10);
    ylabel(ax4,'Load [N]','FontSize',10);
    title(ax4,'Load - Displacement (live)','FontSize',10,'FontWeight','bold');

    hPs = plot(ax4, NaN, NaN, '-', 'Color',[0.82 0.10 0.10], 'LineWidth', 2, ...
               'DisplayName', 'Shear (P_s vs \delta_s)');
    hP  = plot(ax4, NaN, NaN, '-', 'Color',[0.10 0.30 0.75], 'LineWidth', 2, ...
               'DisplayName', 'Normal (P vs \delta)');
    hPsMark = plot(ax4, NaN, NaN, 'o', 'Color',[0.82 0.10 0.10], ...
                   'MarkerFaceColor',[0.82 0.10 0.10],'MarkerSize',7, ...
                   'HandleVisibility','off');
    hPMark  = plot(ax4, NaN, NaN, 'o', 'Color',[0.10 0.30 0.75], ...
                   'MarkerFaceColor',[0.10 0.30 0.75],'MarkerSize',7, ...
                   'HandleVisibility','off');

    yline(ax4, 0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    ylim(ax4, [-5000 30000]);
    legend(ax4, 'Location', 'northwest', 'FontSize', 9, 'Box','off');
    ytickformat(ax4,'%d');

    stride = max(1, round(nIncr/120));

    % ===== Initialize Simulation Video Writer =====
    vid_sim_name = fullfile(out_dir, [prefix '_simulation.mp4']);
    vid_sim = VideoWriter(vid_sim_name, 'MPEG-4');
    vid_sim.FrameRate = 15;
    open(vid_sim);
    fprintf('Recording simulation video to: %s\n', vid_sim_name);

    % ===== Incremental loading loop =====
    fprintf('Starting Nooru-Mohamed Static Implicit Solver...\n');
    for s = 1:nIncr
        lambda = s / nIncr;
        Uy_t   = lambda * Uy_end;

        % Path 4c: proportional shear + normal
        if strcmpi(load_path,'4c')
            deltas = gamma * Uy_t;
            uxL_t  = +0.5 * deltas;
            uxR_t  = -0.5 * deltas;
        else
            uxL_t = NaN; uxR_t = NaN;
        end

        % Dirichlet DOF set for this increment
        dir_nodes = [presc_top_uy; base_fix];
        if strcmpi(load_path,'4c')
            dir_nodes = [dir_nodes; Lx; Rx];
        end
        dir_nodes = unique(dir_nodes);
        isDir = false(ndof,1);  isDir(dir_nodes) = true;
        free  = find(~isDir);

        % External force vector
        if strcmpi(load_path,'shear_force')
            Fext = lambda * Fext_base;
        else
            Fext = sparse(ndof,1);
        end

        % Initial guess - previous converged state with updated BCs
        Utrial = U;
        Utrial(presc_top_uy) = Uy_t;
        if strcmpi(load_path,'4c')
            if ~isempty(Lx), Utrial(Lx) = uxL_t; end
            if ~isempty(Rx), Utrial(Rx) = uxR_t; end
        end
        Utrial(base_fix) = 0;
        x0 = Utrial(free);

        % Jacobian sparsity pattern restricted to free DOFs
        Jpat_ff = spones(Kpat(free, free));
        if isempty(Jpat_ff)
            error('Empty Jacobian pattern: no free DOFs - check boundary conditions.');
        end
        fopts = optimoptions(base_opts, 'JacobPattern', Jpat_ff);

        fun = @(x) residual_and_jac( ...
            x, Utrial, free, dir_nodes, presc_top_uy, Lx, Rx, base_fix, ...
            Uy_t, uxL_t, uxR_t, B3, EDOF, Ke_unit3, idx_glob, ...
            E, nu, k_tc, j0, be, kappa, Fext, ...
            I_trip, J_trip, eid_rep, Val_unit );

        [x_sol, ~, exitflag, output] = fsolve(fun, x0, fopts);
        if exitflag <= 0
            warning('FSOLVE did not fully converge at increment %d (flag=%d, iters=%d).', ...
                    s, exitflag, output.iterations);
        end

        % Accept converged solution
        U(free)          = x_sol;
        U(presc_top_uy)  = Uy_t;
        if strcmpi(load_path,'4c')
            if ~isempty(Lx), U(Lx) = uxL_t; end
            if ~isempty(Rx), U(Rx) = uxR_t; end
        end
        U(base_fix) = 0;

        % ---- Update damage history after convergence ----
        ue12_bc = reshape(U(EDOF), 12, 1, ne);
        epsv_bc = pagemtimes(B3, ue12_bc);
        eeq_bc  = reshape(eqv_strain_modified_vm_vec(epsv_bc, nu, k_tc), [], 1);

        kappa = max(kappa, eeq_bc);     % irreversible update

        De_bc        = zeros(ne,1);
        act_bc       = (kappa >= j0);
        De_bc(act_bc)= 1 - (j0 ./ kappa(act_bc)) .* exp(-be(act_bc) .* (kappa(act_bc)-j0));
        De_bc        = min(max(De_bc, 0), 0.999999);
        s_e_bc       = E .* (1 - De_bc);

        % Internal force assembly
        fe12_bc = pagemtimes(Ke_unit3, ue12_bc);
        fe12_bc = bsxfun(@times, fe12_bc, reshape(s_e_bc,1,1,ne));
        Fint_bc = assemble_accum(idx_glob, fe12_bc, ndof);

        % ---- Separate Ps (shear) and P (normal) ----
        P_normal  = sum(Fint_bc(presc_top_uy));

        if strcmpi(load_path,'4c')
            P_shear = abs(sum(Fint_bc(Lx)));
        else
            P_shear = sum(Fint_bc(Lx));
        end

        % Store for plotting
        if strcmpi(load_path,'4c')
            delta_s_res(s) = uxL_t - uxR_t;
        else
            delta_s_res(s) = mean(U(Lx)) - mean(U(Rx));
        end
        delta_res(s) = mean(U(presc_top_uy));
        Ps_res(s)    = P_shear;
        P_res(s)     = P_normal;

        % ---- Snapshot capture (peak metric = max(|Ps|, |P|)) ----
        peak_metric_now = max(abs(P_shear), abs(P_normal));
        if peak_metric_now > peak_metric_so_far
            peak_metric_so_far = peak_metric_now;
            snap_peak.U      = U;
            snap_peak.De     = De_bc;
            snap_peak.P      = P_normal;
            snap_peak.Ps     = P_shear;
            snap_peak.delta  = delta_res(s);
            snap_peak.deltas = delta_s_res(s);
            snap_peak.step   = s;
        end
        % Post-peak snapshot: once load drops to ~50% of peak metric and lambda > 0.5
        if lambda > 0.5 && peak_metric_now < 0.5*peak_metric_so_far && isempty(snap_pp.U)
            snap_pp.U      = U;
            snap_pp.De     = De_bc;
            snap_pp.P      = P_normal;
            snap_pp.Ps     = P_shear;
            snap_pp.delta  = delta_res(s);
            snap_pp.deltas = delta_s_res(s);
            snap_pp.step   = s;
        end

        % ---- Live plot update ----
        if mod(s,stride)==0 || s==1 || s==nIncr
            % Update curves
            set(hPs,'XData', delta_s_res(1:s), 'YData', Ps_res(1:s));
            set(hP, 'XData', delta_res(1:s),   'YData', P_res(1:s));
            % Update current-point markers
            set(hPsMark,'XData', delta_s_res(s), 'YData', Ps_res(s));
            set(hPMark, 'XData', delta_res(s),   'YData', P_res(s));

            % Update 3D Deformation Geometry
            Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
            P_def = p + def_scale * [Ux, Uy, Uz];

            % Filter internal Cracks
            is_cracked = De_bc > crack_thresh;
            crack_faces = []; crack_colors = [];
            if any(is_cracked)
                crack_mask = repmat(is_cracked, 4, 1);
                crack_faces = allF(crack_mask, :);
                crack_colors = De_bc(owner_per_face(crack_mask));
            end

            % Broadcast geometry to all 3 viewpoints
            for v = 1:3
                set(hHull(v), 'Vertices', P_def);
                if any(is_cracked)
                    set(hCrack(v), 'Vertices', P_def, ...
                                   'Faces', crack_faces, ...
                                   'FaceVertexCData', crack_colors);
                else
                    set(hCrack(v), 'Vertices', P_def, 'Faces', []);
                end
            end

            drawnow;
            frame = getframe(fig);
            writeVideo(vid_sim, frame);
        end
    end

    close(vid_sim);

    % Fallback if post-peak snapshot was never reached
    if isempty(snap_pp.U)
        snap_pp.U      = U;
        snap_pp.De     = De_bc;
        snap_pp.P      = P_res(end);
        snap_pp.Ps     = Ps_res(end);
        snap_pp.delta  = delta_res(end);
        snap_pp.deltas = delta_s_res(end);
        snap_pp.step   = nIncr;
    end

    % ===== Save Data Results =====
    out_csv = fullfile(out_dir, [prefix '_curve_4c.csv']);
    writetable( table(delta_s_res, delta_res, Ps_res, P_res, ...
                      'VariableNames', {'delta_s_mm','delta_mm','Ps_N','P_N'}), ...
                out_csv );
    fprintf('Done. Results saved to %s\n', out_csv);

    % =====================================================================
    % ===== FINAL SLOW ROTATION, VIDEO, AND PHOTO EXPORT =====
    % =====================================================================
    fprintf('Starting final slow rotation to save inspection video and photos...\n');
    vid_rot_name = fullfile(out_dir, [prefix '_final_rotation.mp4']);
    vid_rot = VideoWriter(vid_rot_name, 'MPEG-4');
    vid_rot.FrameRate = 30;
    open(vid_rot);

    for ang = 1:2:360
        camorbit(ax3d(1), 2, 0, 'data', [0 0 1]);
        camorbit(ax3d(2), 2, 0, 'data', [0 0 1]);
        camorbit(ax3d(3), 2, 0, 'data', [0 0 1]);

        drawnow;
        frame = getframe(fig);
        writeVideo(vid_rot, frame);

        if mod(ang, 45) == 1 || ang == 360
            img_name = fullfile(out_dir, sprintf('%s_crack_view_%03d_deg.png', prefix, ang));
            exportgraphics(fig, img_name, 'Resolution', 300);
        end
    end

    close(vid_rot);
    fprintf('Media export complete! Videos and Photos are in the %s/ folder.\n', out_dir);

    % =====================================================================
    % ===== PUBLICATION-QUALITY STATIC FIGURES (solver_main_3pb style) ====
    % =====================================================================
    fprintf('\n==== Generating publication-quality figures ====\n');

    bc.Ntop   = double(Ntop);
    bc.Nbot   = double(Nbot);
    bc.Nleft  = double(Nleft);
    bc.Nright = double(Nright);

    % 1) FE mesh + boundary conditions (3-view)
    fig_mesh_3d(p, extFaces, bc, view_angles, ...
                fullfile(out_dir, [prefix '_fig_mesh']));

    % 2) Load - displacement curve with peak markers + shaded fills
    fig_load_disp(delta_s_res, Ps_res, delta_res, P_res, ...
                  fullfile(out_dir, [prefix '_load_disp']));

    % 3) Damage at peak (3-view)
    fig_damage_3d(p, extFaces, allF, owner_per_face, ...
                  snap_peak.U, snap_peak.De, def_scale, view_angles, ...
                  sprintf('Peak  P = %.2f kN,  P_s = %.2f kN  (\\delta = %.4f mm)', ...
                          snap_peak.P/1000, snap_peak.Ps/1000, snap_peak.delta), ...
                  fullfile(out_dir, [prefix '_damage_peak']));

    % 4) Damage post-peak (3-view)
    fig_damage_3d(p, extFaces, allF, owner_per_face, ...
                  snap_pp.U, snap_pp.De, def_scale, view_angles, ...
                  sprintf('Post-peak  P = %.2f kN,  P_s = %.2f kN  (\\delta = %.4f mm)', ...
                          snap_pp.P/1000, snap_pp.Ps/1000, snap_pp.delta), ...
                  fullfile(out_dir, [prefix '_damage_postpeak']));

    fprintf('Publication figures saved to %s/\n', out_dir);
    fprintf('  %s_fig_mesh.png/pdf\n', prefix);
    fprintf('  %s_load_disp.png/pdf\n', prefix);
    fprintf('  %s_damage_peak.png/pdf\n', prefix);
    fprintf('  %s_damage_postpeak.png/pdf\n', prefix);
end


% =========================================================================
%  Residual & Secant Jacobian
% =========================================================================
function [r_free, J_free] = residual_and_jac(x_free, Uin, free, dir_nodes, ...
                                              presc_top_uy, Lx, Rx, base_fix, ...
                                              Uy_t, uxL_t, uxR_t, ...
                                              B3, EDOF, Ke_unit3, idx_glob, ...
                                              E, nu, k_tc, j0, be, kappa_history, ...
                                              Fext, I_trip, J_trip, eid_rep, Val_unit_cached)
    ndof = numel(Uin);
    ne   = size(EDOF,2);

    % Enforce Dirichlet BCs on full displacement vector
    U = Uin;
    U(free)          = x_free;
    U(presc_top_uy)  = Uy_t;
    if ~isnan(uxL_t) && ~isempty(Lx), U(Lx) = uxL_t; end
    if ~isnan(uxR_t) && ~isempty(Rx), U(Rx) = uxR_t; end
    U(base_fix) = 0;

    % Element strains
    ue12 = reshape(U(EDOF), 12, 1, ne);
    epsv = pagemtimes(B3, ue12);

    % Equivalent strain
    eeq = reshape(eqv_strain_modified_vm_vec(epsv, nu, k_tc), [], 1);

    % Irreversible local history (frozen at start of increment - secant approach)
    kappa_loc = max(kappa_history, eeq);

    % Damage
    De    = zeros(ne,1);
    act   = (kappa_loc >= j0);
    De(act) = 1 - (j0 ./ kappa_loc(act)) .* exp(-be(act) .* (kappa_loc(act)-j0));
    De    = min(max(De, 0), 0.999999);
    s_e   = E .* (1 - De);

    % Internal forces
    fe12 = pagemtimes(Ke_unit3, ue12);
    fe12 = bsxfun(@times, fe12, reshape(s_e,1,1,ne));
    Fint = assemble_accum(idx_glob, fe12, ndof);

    % Residual on free DOFs
    R = Fext - Fint;
    R(dir_nodes) = 0;
    r_free = R(free);

    % Secant tangent stiffness
    sval   = Val_unit_cached .* s_e(eid_rep);
    K      = sparse(double(I_trip), double(J_trip), double(sval), ndof, ndof);
    J_free = -K(free, free);
end


% =========================================================================
%  Helper functions
% =========================================================================
function val = local_get(s, field, def)
    if ~isfield(s,field) || isempty(s.(field))
        val = def;
    else
        val = s.(field);
    end
end

% ------------------------------------------------------------------
function [D, lambda, G] = iso3D_D(E, nu)
    G      = E / (2*(1+nu));
    lambda = E*nu / ((1+nu)*(1-2*nu));
    D = [ lambda+2*G, lambda,      lambda,      0, 0, 0;
          lambda,     lambda+2*G,  lambda,      0, 0, 0;
          lambda,     lambda,      lambda+2*G,  0, 0, 0;
          0,          0,           0,           G, 0, 0;
          0,          0,           0,           0, G, 0;
          0,          0,           0,           0, 0, G ];
end

% ------------------------------------------------------------------
function [B3, V3, he3, hmin] = precompute_TET4_vector(nodes, tets)
    ne   = size(tets,1);
    B3   = zeros(6,12,ne);
    V3   = zeros(1,1,ne);
    he3  = zeros(ne,1);
    hmin = inf;

    edges = [1 2;1 3;1 4;2 3;2 4;3 4];

    for e = 1:ne
        idx = tets(e,:);  Xe = nodes(idx,:);
        x1 = Xe(1,:).';  x2 = Xe(2,:).';
        x3 = Xe(3,:).';  x4 = Xe(4,:).';

        Jac   = [x2-x1, x3-x1, x4-x1];
        detJ  = det(Jac);
        V     = abs(detJ)/6;
        if V <= 1e-16, error('Non-positive element volume at e=%d',e); end

        V3(1,1,e) = V;
        JTinv     = (Jac.') \ eye(3);
        g2 = JTinv(:,1);  g3 = JTinv(:,2);
        g4 = JTinv(:,3);  g1 = -(g2+g3+g4);
        grads = [g1 g2 g3 g4];

        B = zeros(6,12);
        for a = 1:4
            dNx = grads(1,a);  dNy = grads(2,a);  dNz = grads(3,a);
            c   = (3*(a-1)+1):(3*a);
            B(1,c) = [dNx  0    0  ];
            B(2,c) = [0    dNy  0  ];
            B(3,c) = [0    0    dNz];
            B(4,c) = [dNy  dNx  0  ];
            B(5,c) = [0    dNz  dNy];
            B(6,c) = [dNz  0    dNx];
        end
        B3(:,:,e) = B;

        he_vol   = (12*V)^(1/3);
        he_min_e = inf;
        for k = 1:6
            i = edges(k,1);  j = edges(k,2);
            lij = norm(Xe(i,:) - Xe(j,:));
            if lij < he_min_e, he_min_e = lij; end
            if lij < hmin,     hmin     = lij; end
        end
        he3(e,1) = min(he_vol, he_min_e);
    end
end

% ------------------------------------------------------------------
function eeq = eqv_strain_modified_vm_vec(epsv, nu, k)
    exx = squeeze(epsv(1,1,:));  eyy = squeeze(epsv(2,1,:));
    ezz = squeeze(epsv(3,1,:));
    gxy = squeeze(epsv(4,1,:));  gyz = squeeze(epsv(5,1,:));
    gzx = squeeze(epsv(6,1,:));
    exy = 0.5*gxy;  eyz = 0.5*gyz;  ezx = 0.5*gzx;

    I1  = exx + eyy + ezz;
    sxx = exx - I1/3;  syy = eyy - I1/3;  szz = ezz - I1/3;
    J2  = 0.5*(sxx.^2 + syy.^2 + szz.^2 + 2*(exy.^2 + eyz.^2 + ezx.^2));

    a1 = (k-1) ./ (1-2*nu);
    a2 = 12*k  ./ (1+nu).^2;

    term1 = (k-1) ./ (2*k*(1-2*nu)) .* I1;
    term2 = (1/(2*k)) .* sqrt( a1.^2 .* I1.^2 + a2.*J2 );
    eeq   = term1 + term2;
    eeq   = reshape(eeq, 1, 1, []);
end

% ------------------------------------------------------------------
function Fglob = assemble_accum(idx_glob, vals_12x1xne, ndof)
    subs  = double(idx_glob);
    vals  = vals_12x1xne(:);
    Fglob = accumarray(subs, vals, [ndof,1], @sum, 0);
end


% =========================================================================
%  PUBLICATION-QUALITY FIGURE HELPERS (matches solver_main_3pb)
% =========================================================================

% ---------------------------------------------------------------------
%  CRACK COLORMAP -- blue-grey -> blue -> cyan -> green -> yellow -> red
% ---------------------------------------------------------------------
function c = crack_cmap()
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


% ---------------------------------------------------------------------
%  3D MESH FIGURE  -- 3-view layout with BC highlights and load arrows
% ---------------------------------------------------------------------
function fig_mesh_3d(p, extFaces, bc, view_angles, basepath)
    fh = figure('Color','w','Position',[100 100 1400 480],'Visible','off', ...
                'PaperUnits','centimeters','PaperSize',[36 9.5], ...
                'PaperPosition',[0 0 36 9.5]);
    tl = tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');
    title(tl,'FE mesh -- boundary conditions and load configuration', ...
          'FontWeight','bold','FontSize',12);

    xr = [min(p(:,1)) max(p(:,1))];
    yr = [min(p(:,2)) max(p(:,2))];
    zr = [min(p(:,3)) max(p(:,3))];
    Lx_geom = max(xr(2)-xr(1), 1);

    for v = 1:3
        ax = nexttile(tl);
        hold(ax,'on'); axis(ax,'equal','vis3d');
        view(ax, view_angles(v,:));

        % Outer hull
        patch(ax, 'Faces',extFaces, 'Vertices',p, ...
              'FaceColor',[0.93 0.94 0.97], ...
              'EdgeColor',[0.55 0.62 0.78], ...
              'LineWidth',0.2, 'FaceAlpha',0.35);

        % Bottom support nodes (fixed Y) - grey triangles
        if ~isempty(bc.Nbot)
            scatter3(ax, p(bc.Nbot,1), p(bc.Nbot,2), p(bc.Nbot,3), ...
                     22, [0.30 0.30 0.35], 'v', 'filled', ...
                     'MarkerEdgeColor','k', 'LineWidth',0.4, ...
                     'DisplayName','Fixed (bottom)');
        end

        % Top loaded nodes - green
        if ~isempty(bc.Ntop)
            scatter3(ax, p(bc.Ntop,1), p(bc.Ntop,2), p(bc.Ntop,3), ...
                     22, [0.0 0.55 0.30], '^', 'filled', ...
                     'MarkerEdgeColor','k', 'LineWidth',0.4, ...
                     'DisplayName','Loaded (top, P)');
        end

        % Left lip nodes - red
        if ~isempty(bc.Nleft)
            scatter3(ax, p(bc.Nleft,1), p(bc.Nleft,2), p(bc.Nleft,3), ...
                     18, [0.82 0.10 0.10], 'o', 'filled', ...
                     'MarkerEdgeColor','k', 'LineWidth',0.3, ...
                     'DisplayName','Left lip (+P_s)');
        end

        % Right lip nodes - blue
        if ~isempty(bc.Nright)
            scatter3(ax, p(bc.Nright,1), p(bc.Nright,2), p(bc.Nright,3), ...
                     18, [0.10 0.30 0.75], 'o', 'filled', ...
                     'MarkerEdgeColor','k', 'LineWidth',0.3, ...
                     'DisplayName','Right lip (-P_s)');
        end

        % Load arrow at top centroid
        if ~isempty(bc.Ntop)
            cT = mean(p(bc.Ntop,:),1);
            arrL = 0.20 * (yr(2)-yr(1) + 1e-9);
            quiver3(ax, cT(1), cT(2)+arrL, cT(3), 0, -arrL*0.85, 0, 0, ...
                    'Color',[0.0 0.55 0.30],'LineWidth',2.2, ...
                    'MaxHeadSize',0.6,'AutoScale','off', ...
                    'HandleVisibility','off');
            text(ax, cT(1), cT(2)+arrL*1.1, cT(3), 'P', ...
                 'Interpreter','tex','FontSize',12,'FontWeight','bold', ...
                 'Color',[0.0 0.55 0.30], ...
                 'HorizontalAlignment','center');
        end

        % Shear arrows on lips
        if ~isempty(bc.Nleft) && ~isempty(bc.Nright)
            cL = mean(p(bc.Nleft,:),1);
            cR = mean(p(bc.Nright,:),1);
            arrS = 0.15 * Lx_geom;
            quiver3(ax, cL(1)-arrS, cL(2), cL(3), arrS*0.85, 0, 0, 0, ...
                    'Color',[0.82 0.10 0.10],'LineWidth',1.8, ...
                    'MaxHeadSize',0.6,'AutoScale','off', ...
                    'HandleVisibility','off');
            quiver3(ax, cR(1)+arrS, cR(2), cR(3), -arrS*0.85, 0, 0, 0, ...
                    'Color',[0.10 0.30 0.75],'LineWidth',1.8, ...
                    'MaxHeadSize',0.6,'AutoScale','off', ...
                    'HandleVisibility','off');
            text(ax, cL(1)-arrS, cL(2), cL(3), 'P_s', ...
                 'Interpreter','tex','FontSize',10,'FontWeight','bold', ...
                 'Color',[0.82 0.10 0.10], ...
                 'HorizontalAlignment','right');
            text(ax, cR(1)+arrS, cR(2), cR(3), 'P_s', ...
                 'Interpreter','tex','FontSize',10,'FontWeight','bold', ...
                 'Color',[0.10 0.30 0.75], ...
                 'HorizontalAlignment','left');
        end

        grid(ax,'on'); box(ax,'on');
        ax.FontSize = 8; ax.LineWidth = 0.6;
        xlabel(ax,'$x$ [mm]','Interpreter','latex');
        ylabel(ax,'$y$ [mm]','Interpreter','latex');
        zlabel(ax,'$z$ [mm]','Interpreter','latex');
        title(ax, sprintf('View %d  (az=%d^{\\circ}, el=%d^{\\circ})', ...
                          v, view_angles(v,1), view_angles(v,2)), ...
              'FontSize',9,'FontWeight','bold');

        if v == 1
            legend(ax,'Location','northeastoutside','FontSize',7,'Box','off');
        end
    end

    save_fig_hq(fh, basepath);
    close(fh);
end


% ---------------------------------------------------------------------
%  LOAD - DISPLACEMENT FIGURE  -- shaded fills, peak markers, labels
% ---------------------------------------------------------------------
function fig_load_disp(delta_s, Ps, delta, P, basepath)
    fh = figure('Color','w','Position',[100 100 700 520],'Visible','off', ...
                'PaperUnits','centimeters','PaperSize',[14 10.5], ...
                'PaperPosition',[0 0 14 10.5]);
    ax = axes('Parent',fh,'Units','normalized','Position',[0.12 0.13 0.83 0.79]);
    hold(ax,'on'); grid(ax,'on');
    ax.GridColor=[0.80 0.80 0.80]; ax.GridLineStyle=':';
    ax.FontSize=10; ax.LineWidth=0.7; ax.Box='on';

    % Shaded fills under each curve
    fill(ax, [delta_s(:); flipud(delta_s(:))], ...
             [Ps(:);      zeros(size(Ps(:)))], ...
         [0.82 0.10 0.10],'FaceAlpha',0.10,'EdgeColor','none', ...
         'HandleVisibility','off');
    fill(ax, [delta(:);   flipud(delta(:))], ...
             [P(:);       zeros(size(P(:)))], ...
         [0.10 0.30 0.75],'FaceAlpha',0.10,'EdgeColor','none', ...
         'HandleVisibility','off');

    % Curves
    hPs = plot(ax, delta_s, Ps, '-','Color',[0.82 0.10 0.10],'LineWidth',2.0, ...
               'DisplayName','Shear  (P_s vs \delta_s)');
    hP  = plot(ax, delta,   P,  '-','Color',[0.10 0.30 0.75],'LineWidth',2.0, ...
               'DisplayName','Normal (P vs \delta)');

    % Peak markers
    [pk_Ps, iPs] = max(abs(Ps)); pk_Ps_signed = Ps(iPs);
    [pk_P,  iP ] = max(abs(P));  pk_P_signed  = P(iP);
    plot(ax, delta_s(iPs), pk_Ps_signed, 'pentagram', ...
         'MarkerSize',13, ...
         'MarkerFaceColor',[0.98 0.82 0.0], ...
         'MarkerEdgeColor',[0.40 0.10 0.10],'LineWidth',1.0, ...
         'HandleVisibility','off');
    plot(ax, delta(iP), pk_P_signed, 'pentagram', ...
         'MarkerSize',13, ...
         'MarkerFaceColor',[0.98 0.82 0.0], ...
         'MarkerEdgeColor',[0.10 0.10 0.40],'LineWidth',1.0, ...
         'HandleVisibility','off');

    % Peak label box for Ps
    xr_data = max([delta_s(:); delta(:)]);
    text(ax, delta_s(iPs) + xr_data*0.02, pk_Ps_signed*1.05, ...
         sprintf('$P_{s,\\rm peak} = %.2f$ kN\n$\\delta_s = %.4f$ mm', ...
                 pk_Ps/1000, delta_s(iPs)), ...
         'Interpreter','latex','FontSize',8.5,'FontWeight','bold', ...
         'Color',[0.40 0.10 0.10], ...
         'HorizontalAlignment','left','VerticalAlignment','bottom', ...
         'BackgroundColor','w','EdgeColor',[0.60 0.20 0.20], ...
         'Margin',3,'LineWidth',0.5);

    % Peak label box for P
    text(ax, delta(iP) + xr_data*0.02, pk_P_signed*1.05, ...
         sprintf('$P_{\\rm peak} = %.2f$ kN\n$\\delta = %.4f$ mm', ...
                 pk_P/1000, delta(iP)), ...
         'Interpreter','latex','FontSize',8.5,'FontWeight','bold', ...
         'Color',[0.10 0.10 0.40], ...
         'HorizontalAlignment','left','VerticalAlignment','bottom', ...
         'BackgroundColor','w','EdgeColor',[0.20 0.20 0.60], ...
         'Margin',3,'LineWidth',0.5);

    yline(ax, 0, 'k--','LineWidth',0.8,'HandleVisibility','off');

    xlabel(ax,'Displacement [mm]','Interpreter','latex','FontSize',11);
    ylabel(ax,'Load [N]','Interpreter','latex','FontSize',11);
    title(ax,'Load -- displacement response (Nooru-Mohamed path 4c)', ...
          'FontSize',11,'FontWeight','bold');
    legend(ax,[hPs,hP],'Location','northwest','FontSize',9,'Box','off');

    y_top = max(abs([Ps(:); P(:)]))*1.30;
    y_bot = min([0; Ps(:); P(:)])*1.20;
    ylim(ax,[y_bot y_top]);
    xlim(ax,[0 xr_data*1.05]);

    save_fig_hq(fh, basepath);
    close(fh);
end


% ---------------------------------------------------------------------
%  3D DAMAGE FIGURE  -- 3-view layout with hull + crack overlay
% ---------------------------------------------------------------------
function fig_damage_3d(p, extFaces, allF, owner_per_face, U, De, def_scale, ...
                       view_angles, label, basepath)
    if isempty(U) || isempty(De)
        return;
    end

    % Deformed coordinates
    Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
    P_def = p + def_scale * [Ux, Uy, Uz];

    % Filter cracked faces (light threshold so we see process zone too)
    THRESH = 0.05;
    is_cracked = De > THRESH;
    crack_faces = []; crack_colors = [];
    if any(is_cracked)
        crack_mask  = repmat(is_cracked, 4, 1);
        crack_faces = allF(crack_mask, :);
        crack_colors = De(owner_per_face(crack_mask));
    end

    fh = figure('Color','w','Position',[100 100 1500 480],'Visible','off', ...
                'PaperUnits','centimeters','PaperSize',[38 9.5], ...
                'PaperPosition',[0 0 38 9.5]);
    tl = tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');
    title(tl, label, ...
          'Interpreter','tex','FontWeight','bold','FontSize',11);

    cm = crack_cmap();

    for v = 1:3
        ax = nexttile(tl);
        hold(ax,'on'); axis(ax,'equal','vis3d');
        view(ax, view_angles(v,:));

        % Layer 1: full mesh in light grey (deformed)
        patch(ax,'Faces',extFaces, 'Vertices',P_def, ...
              'FaceColor',[0.92 0.92 0.93], ...
              'EdgeColor','none', ...
              'FaceAlpha',0.18, ...
              'AmbientStrength',0.6);

        % Layer 2: cracked elements colored by damage
        if ~isempty(crack_faces)
            patch(ax,'Faces',crack_faces, 'Vertices',P_def, ...
                  'FaceVertexCData',crack_colors, ...
                  'FaceColor','flat','EdgeColor','none', ...
                  'AmbientStrength',0.8);
        end

        camlight(ax,'headlight'); lighting(ax,'gouraud');
        colormap(ax, cm);
        try, clim(ax,[0 1]); catch, caxis(ax,[0 1]); end %#ok<NOSEMI>

        grid(ax,'on'); box(ax,'on');
        ax.FontSize = 8; ax.LineWidth = 0.6;
        xlabel(ax,'$x$ [mm]','Interpreter','latex');
        ylabel(ax,'$y$ [mm]','Interpreter','latex');
        zlabel(ax,'$z$ [mm]','Interpreter','latex');
        title(ax, sprintf('View %d  (az=%d^{\\circ}, el=%d^{\\circ})', ...
                          v, view_angles(v,1), view_angles(v,2)), ...
              'FontSize',9,'FontWeight','bold');

        if v == 3
            cb = colorbar(ax,'Location','eastoutside');
            cb.Label.String = '\omega  (damage)';
            cb.Label.FontSize = 9;
        end
    end

    save_fig_hq(fh, basepath);
    close(fh);
end


% ---------------------------------------------------------------------
%  SAVE FIGURE -- 600 dpi PNG + vector-quality PDF
% ---------------------------------------------------------------------
function save_fig_hq(fh, basepath)
    print(fh, [basepath '.png'], '-dpng', '-r600');
    try
        print(fh, [basepath '.pdf'], '-dpdf', '-painters','-r600');
    catch
        try
            exportgraphics(fh, [basepath '.pdf'], 'ContentType','vector');
        catch
        end
    end
end