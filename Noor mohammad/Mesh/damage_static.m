function damage_static(prefix, opts)

if nargin < 1 || isempty(prefix), prefix = 'Job-1'; end
if nargin < 2 || isempty(opts),   opts   = struct();  end
get_opt = @(f,def) local_get(opts,f,def);

E      = get_opt('E',    29000.0);
nu     = get_opt('nu',   0.20);
GF     = get_opt('GF',   0.11);
ft     = get_opt('ft',   3.00);
k_tc   = get_opt('k',    10.0);
kappa0 = get_opt('kappa0', ft/E);
j0     = kappa0;

nIncr       = get_opt('nIncr',       900);
Uy_end      = get_opt('Uy_end',      0.50);
gamma       = get_opt('gamma',       0.60);
load_path   = get_opt('load_path',  '4c');
Fx_tot      = get_opt('Fx_tot',      1.0e4);
tol         = get_opt('tol',         1e-6);
step_tol    = get_opt('step_tol',    1e-10);
maxIter     = get_opt('maxIter',     40);
use_line_search = get_opt('use_line_search', false);
maxLineSearch   = get_opt('maxLineSearch',   6);
print_stride    = get_opt('print_stride',    10);

def_scale       = get_opt('def_scale',       15.0);
full_damage_thr = get_opt('full_damage_thresh', 0.95);

live_damage_thr = get_opt('live_damage_thresh', get_opt('plot_damage_thresh', 0.005));
save_damage_thr = get_opt('save_damage_thresh', full_damage_thr);

save_all_increment_geometry = get_opt('save_all_increment_geometry', true);
snapshot_stride = max(1, round(get_opt('snapshot_stride', 1)));
save_when_no_full_damage = get_opt('save_when_no_full_damage', false);

show_live       = get_opt('show_live', true);
live_stride     = max(1, round(get_opt('live_stride', 1)));
image_res       = get_opt('image_res', 180);
view_angle      = get_opt('view_angle', [0 90]);
damage_clim     = get_opt('damage_clim', [0.0 1.0]);
out_dir         = get_opt('out_dir', 'out_NR_vectorized_LIVE_damage_mesh');

show_mesh       = get_opt('show_mesh', true);
save_show_mesh  = get_opt('save_show_mesh', true);
mesh_alpha      = get_opt('mesh_alpha', 0.10);
mesh_edge_color = get_opt('mesh_edge_color', [0.55 0.58 0.62]);
mesh_face_color = get_opt('mesh_face_color', [0.88 0.90 0.94]);
mesh_line_width = get_opt('mesh_line_width', 0.12);

damage_edge       = get_opt('damage_edge', true);
damage_edge_color = get_opt('damage_edge_color', [0.20 0.20 0.20]);
damage_edge_alpha = get_opt('damage_edge_alpha', 0.25);
damage_line_width = get_opt('damage_line_width', 0.08);

damage_colormap = get_opt('damage_colormap', 'turbo');
damage_clim_mode = lower(string(get_opt('damage_clim_mode', 'visible')));
use_discrete_damage_colors = get_opt('use_discrete_damage_colors', false);
n_damage_bands = max(2, round(get_opt('n_damage_bands', 12)));

bandwidth_method = lower(string(get_opt('bandwidth_method', 'oliver')));

if ~exist(out_dir,'dir'), mkdir(out_dir); end
inc_dir = fullfile(out_dir, 'increment_damaged_geometry_png');
if save_all_increment_geometry && ~exist(inc_dir,'dir'), mkdir(inc_dir); end

is4c = strcmpi(load_path,'4c');

nd = readmatrix([prefix '_nodes.txt']);
ids = nd(:,1);
p   = nd(:,2:4);
np  = size(p,1);

ed = readmatrix([prefix '_elements.txt']);
[okElem,T] = ismember(ed(:,2:5), ids);
if any(~okElem(:))
    error('Some element connectivity node IDs were not found in %s_nodes.txt.', prefix);
end
T  = double(T);
ne = size(T,1);

Ntop   = read_nodeset_local([prefix '_top_nodes.txt'],    ids, 'top');
Nbot   = read_nodeset_local([prefix '_bottom_nodes.txt'], ids, 'bottom');
Nleft  = read_nodeset_local([prefix '_left_nodes.txt'],   ids, 'left');
Nright = read_nodeset_local([prefix '_right_nodes.txt'],  ids, 'right');

fprintf('\nFAST VECTORIZED MODIFIED NEWTON DAMAGE SOLVER\n');
fprintf('Nodes = %d | Elements = %d | Increments = %d\n', np, ne, nIncr);
fprintf('Live damage threshold: omega >= %.4f\n', live_damage_thr);
fprintf('Saved PNG damage threshold: omega >= %.4f\n', save_damage_thr);

fprintf('Precomputing vectorized TET4 B matrices...\n');
[B3, V3, gradN3, valid] = precompute_TET4_vectorized(p,T);
if any(~valid)
    warning('%d invalid or near-zero tetrahedra detected.', nnz(~valid));
end

[Dunit, ~, ~] = iso3D_D(1.0, nu);
Ke_unit3 = pagemtimes(permute(B3,[2 1 3]), pagemtimes(repmat(Dunit,[1 1 ne]), B3));
Ke_unit3 = bsxfun(@times, Ke_unit3, reshape(V3,1,1,ne));

ndof = 3*np;
EDOF = build_edof(T);
idx_glob = EDOF(:);

fprintf('Prebuilding global sparse triplet pattern...\n');
[I_trip, J_trip, Val_unit, eid_rep] = build_triplet_pattern(EDOF, Ke_unit3);

U = zeros(ndof,1);
kappa = j0 * ones(ne,1);
h_band_bc = zeros(ne,1);

dof = @(nid,comp) 3*(double(nid)-1)+comp;

fixY_bottom = dof(Nbot, 2);
presc_top_uy = dof(Ntop, 2);

xb = p(double(Nbot),1);
[~,ib] = min(abs(xb - mean(xb)));
fixX_pin = dof(Nbot(ib), 1);

[~,iPinZ] = min(p(:,3));
fixZ_pin = 3*(iPinZ-1)+3;

Lx = dof(Nleft,  1);
Rx = dof(Nright, 1);

if is4c
    base_fix = unique([fixY_bottom; fixZ_pin]);
else
    base_fix = unique([fixY_bottom; fixX_pin; fixZ_pin]);
end

dir_nodes = [presc_top_uy; base_fix];
if is4c
    dir_nodes = [dir_nodes; Lx; Rx];
end
dir_nodes = unique(dir_nodes(:));

isDir = false(ndof,1);
isDir(dir_nodes) = true;
free = find(~isDir);

fprintf('Free DOFs = %d | Prescribed/fixed DOFs = %d\n', numel(free), numel(dir_nodes));

Fext_base = sparse(ndof,1);
if strcmpi(load_path,'shear_force')
    if ~isempty(Nleft),  Fext_base(Lx) = +Fx_tot / max(numel(Nleft),1);  end
    if ~isempty(Nright), Fext_base(Rx) = -Fx_tot / max(numel(Nright),1); end
end

delta_s_res = zeros(nIncr,1);
delta_res   = zeros(nIncr,1);
Ps_res      = zeros(nIncr,1);
P_res       = zeros(nIncr,1);
maxD_res    = zeros(nIncr,1);
meanD_res   = zeros(nIncr,1);
nD005_res   = zeros(nIncr,1);
nD050_res   = zeros(nIncr,1);
nD950_res   = zeros(nIncr,1);
nPlot_res   = zeros(nIncr,1);
iter_res    = zeros(nIncr,1);
relres_res  = zeros(nIncr,1);
mean_h_res  = zeros(nIncr,1);
min_h_res   = zeros(nIncr,1);
max_h_res   = zeros(nIncr,1);
saved_png   = cell(nIncr,1);

fprintf('Extracting faces for damaged-geometry images...\n');
allF = [T(:,[1 2 3]); T(:,[1 2 4]); T(:,[1 3 4]); T(:,[2 3 4])];
owner_per_face = repmat((1:ne)',4,1);
extFaces = external_hull(T);

if show_live
    vis_state = 'on';
else
    vis_state = 'off';
end

fig = figure('Color','w', ...
    'Name','Live damaged geometry', ...
    'Visible',vis_state, ...
    'Position',[80 80 950 820]);

ax = axes(fig);
hold(ax,'on');
axis(ax,'equal');
axis(ax,'vis3d');
box(ax,'on');
grid(ax,'on');
colormap(ax, select_damage_cmap(damage_colormap, n_damage_bands, use_discrete_damage_colors));
try
    clim(ax, damage_clim);
catch
    caxis(ax, damage_clim);
end

hHull = patch(ax, ...
    'Faces', extFaces, ...
    'Vertices', p, ...
    'FaceColor', mesh_face_color, ...
    'EdgeColor', mesh_edge_color, ...
    'LineWidth', mesh_line_width, ...
    'FaceAlpha', mesh_alpha, ...
    'Visible', ternary_vis(show_live && show_mesh));
try_set(hHull,'EdgeAlpha',0.45);

if damage_edge
    dam_edge_color = damage_edge_color;
else
    dam_edge_color = 'none';
end

hDam = patch(ax, ...
    'Faces', zeros(0,3), ...
    'Vertices', p, ...
    'FaceVertexCData', zeros(0,1), ...
    'FaceColor','flat', ...
    'EdgeColor', dam_edge_color, ...
    'LineWidth', damage_line_width, ...
    'FaceAlpha',1.0, ...
    'AmbientStrength',0.85);
try_set(hDam,'EdgeAlpha',damage_edge_alpha);

cb = colorbar(ax,'Location','eastoutside');
ylabel(cb,'Damage \omega');
xlabel(ax,'x [mm]');
ylabel(ax,'y [mm]');
zlabel(ax,'z [mm]');
view(ax, view_angle);
camproj(ax,'orthographic');
camup(ax,[0 1 0]);
lighting(ax,'gouraud');

xyzMin = min(p,[],1);
xyzMax = max(p,[],1);
span = max(xyzMax - xyzMin);
extraDisp = def_scale * max(abs([Uy_end, gamma*Uy_end, Fx_tot*0]));
margin = 0.08*span + extraDisp;
xlim(ax,[xyzMin(1)-margin, xyzMax(1)+margin]);
ylim(ax,[xyzMin(2)-margin, xyzMax(2)+margin]);
zlim(ax,[xyzMin(3)-margin, xyzMax(3)+margin]);
daspect(ax,[1 1 1]);
axis(ax,'manual');

fprintf('Starting modified Newton solver...\n');

tic_total = tic;

for s = 1:nIncr
    lambda_load = s/nIncr;
    Uy_t = lambda_load * Uy_end;

    if is4c
        deltas = gamma * Uy_t;
        uxL_t = +0.5*deltas;
        uxR_t = -0.5*deltas;
    else
        uxL_t = NaN;
        uxR_t = NaN;
    end

    if strcmpi(load_path,'shear_force')
        Fext = lambda_load * Fext_base;
    else
        Fext = sparse(ndof,1);
    end

    Utrial = U;
    Utrial(presc_top_uy) = Uy_t;
    Utrial(base_fix) = 0.0;
    if is4c
        if ~isempty(Lx), Utrial(Lx) = uxL_t; end
        if ~isempty(Rx), Utrial(Rx) = uxR_t; end
    end

    x = Utrial(free);
    relres = Inf;
    converged = false;

    for it = 1:maxIter
        [R, Jff, Fint_it] = residual_and_jac_NR( ...
            x, Utrial, free, B3, EDOF, Ke_unit3, idx_glob, ...
            E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, Fext, ...
            I_trip, J_trip, eid_rep, Val_unit, ndof);

        scaleR = max([1.0, double(full(norm(Fext(free),inf))), double(full(norm(Fint_it(free),inf)))]);
        relres = double(full(norm(R,inf))) / scaleR;

        if relres <= tol
            converged = true;
            break;
        end

        dx = -Jff \ R;

        if any(~isfinite(dx))
            warning('Non-finite Newton correction at increment %d, iteration %d.', s, it);
            break;
        end

        if use_line_search
            norm0 = double(full(norm(R,inf)));
            alpha = 1.0;
            accepted = false;
            x_best = x + dx;
            best_norm = Inf;

            for ls = 1:maxLineSearch
                x_try = x + alpha*dx;
                Rtry = residual_only_NR( ...
                    x_try, Utrial, free, B3, EDOF, Ke_unit3, idx_glob, ...
                    E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, Fext, ndof);

                ntry = double(full(norm(Rtry,inf)));
                if ntry < best_norm
                    best_norm = ntry;
                    x_best = x_try;
                end
                if ntry <= (1 - 1e-4*alpha)*norm0 || ntry < tol*scaleR
                    x = x_try;
                    accepted = true;
                    break;
                end
                alpha = 0.5*alpha;
            end

            if ~accepted
                x = x_best;
            end
        else
            x = x + dx;
        end

        if double(full(norm(dx,inf))) <= step_tol * max(1.0, double(full(norm(x,inf))))
            converged = true;
            break;
        end
    end

    if ~converged
        warning('Newton did not fully converge at increment %d. it=%d, relres=%.3e', s, it, relres);
    end

    U(free) = x;
    U(presc_top_uy) = Uy_t;
    U(base_fix) = 0.0;
    if is4c
        if ~isempty(Lx), U(Lx) = uxL_t; end
        if ~isempty(Rx), U(Rx) = uxR_t; end
    end

    [De_bc, ~, Fint_bc, eeq_bc, h_band_bc] = damage_state_and_internal_force( ...
        U, B3, EDOF, Ke_unit3, idx_glob, E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, ndof);

    kappa = max(kappa, eeq_bc);

    [De_bc, ~, Fint_bc, ~, h_band_bc] = damage_state_and_internal_force( ...
        U, B3, EDOF, Ke_unit3, idx_glob, E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, ndof);

    Freac = Fint_bc - Fext;
    P_normal = double(full(sum(Freac(presc_top_uy))));

    if is4c
        if isempty(Lx)
            P_shear = 0.0;
        else
            P_shear = double(full(abs(sum(Freac(Lx)))));
        end
    else
        if isempty(Lx)
            P_shear = 0.0;
        else
            P_shear = double(full(sum(Fint_bc(Lx))));
        end
    end

    if is4c
        delta_s_res(s) = uxL_t - uxR_t;
    else
        delta_s_res(s) = safe_mean(U(Lx)) - safe_mean(U(Rx));
    end

    delta_res(s) = safe_mean(U(presc_top_uy));
    Ps_res(s) = P_shear;
    P_res(s) = P_normal;
    maxD_res(s) = double(max(De_bc));
    meanD_res(s) = double(mean(De_bc));
    nD005_res(s) = double(nnz(De_bc >= 0.005));
    nD050_res(s) = double(nnz(De_bc >= 0.050));
    nD950_res(s) = double(nnz(De_bc >= 0.950));
    nPlot_res(s) = double(nnz(De_bc >= live_damage_thr));
    iter_res(s) = double(it);
    relres_res(s) = double(relres);
    mean_h_res(s) = double(mean(h_band_bc));
    min_h_res(s)  = double(min(h_band_bc));
    max_h_res(s)  = double(max(h_band_bc));

    do_save = save_all_increment_geometry && (mod(s,snapshot_stride)==0 || s==1 || s==nIncr);
    if do_save
        [did_save, png_name] = save_damage_geometry_snapshot( ...
            fig, ax, hHull, hDam, p, U, def_scale, allF, owner_per_face, De_bc, ...
            save_damage_thr, s, nIncr, lambda_load, inc_dir, prefix, image_res, ...
            save_when_no_full_damage, save_show_mesh, damage_clim, damage_clim_mode);

        if did_save
            saved_png{s} = png_name;
        end
    end

    do_live = show_live && (mod(s,live_stride)==0 || s==1 || s==nIncr);
    if do_live
        update_damage_geometry_live( ...
            fig, ax, hHull, hDam, p, U, def_scale, allF, owner_per_face, De_bc, ...
            live_damage_thr, s, nIncr, lambda_load, maxD_res(s), nPlot_res(s), show_mesh, damage_clim, damage_clim_mode);
        drawnow limitrate;
    end

    if mod(s,print_stride)==0 || s==1 || s==nIncr
        fprintf('inc %4d/%4d | it=%2d | relres=%.2e | maxD=%.4f | Nlive(omega>=%.3f)=%d | Nfull(omega>=0.950)=%d | hOliver(mean)=%.3e | Ps=%.3e | P=%.3e\n', ...
            s, nIncr, it, relres, maxD_res(s), live_damage_thr, nPlot_res(s), nD950_res(s), mean_h_res(s), P_shear, P_normal);
    end
end

elapsed = toc(tic_total);

final_png = fullfile(out_dir, [prefix '_FINAL_damage_geometry_with_mesh.png']);
[~, ~] = save_damage_geometry_snapshot( ...
    fig, ax, hHull, hDam, p, U, def_scale, allF, owner_per_face, De_bc, ...
    save_damage_thr, nIncr, nIncr, 1.0, out_dir, [prefix '_FINAL'], image_res, true, save_show_mesh, damage_clim, damage_clim_mode);

helper_final = fullfile(out_dir, sprintf('%s_FINAL_inc_%04d_damageMesh_omega_ge_%s.png', ...
    prefix, nIncr, fmt_underscore(save_damage_thr)));
if exist(helper_final,'file')
    if exist(final_png,'file'), delete(final_png); end
    movefile(helper_final, final_png);
else
    exportgraphics(fig, final_png, 'Resolution', image_res);
end

final_fig = fullfile(out_dir, [prefix '_FINAL_damage_geometry_with_mesh.fig']);
savefig(fig, final_fig);

out_csv = fullfile(out_dir, [prefix '_fast_NR_results.csv']);
results = table(delta_s_res, delta_res, Ps_res, P_res, maxD_res, meanD_res, ...
                nD005_res, nD050_res, nD950_res, nPlot_res, ...
                iter_res, relres_res, mean_h_res, min_h_res, max_h_res, saved_png, ...
    'VariableNames', {'delta_s_mm','delta_mm','Ps_N','P_N', ...
                      'max_damage','mean_damage', ...
                      'nElem_D_ge_0p005','nElem_D_ge_0p050','nElem_D_ge_0p950', ...
                      'nElem_D_ge_liveThreshold','newton_iterations','relative_residual', ...
                      'mean_h_oliver_mm','min_h_oliver_mm','max_h_oliver_mm', ...
                      'saved_png'});

writetable(results, out_csv);

out_mat = fullfile(out_dir, [prefix '_fast_NR_final_state.mat']);
save(out_mat, 'U','kappa','De_bc','delta_s_res','delta_res','Ps_res','P_res', ...
              'maxD_res','meanD_res','nD005_res','nD050_res','nD950_res', ...
              'nPlot_res','iter_res','relres_res','mean_h_res','min_h_res','max_h_res','h_band_bc', ...
              'live_damage_thr','save_damage_thr','full_damage_thr', ...
              'elapsed','damage_colormap','damage_clim_mode','use_discrete_damage_colors','n_damage_bands', ...
              'bandwidth_method','-v7.3');

fprintf('\nDONE.\n');
fprintf('Elapsed time: %.2f s\n', elapsed);
fprintf('Saved results CSV: %s\n', out_csv);
fprintf('Saved final damaged geometry/mesh PNG: %s\n', final_png);
fprintf('Saved final MATLAB figure: %s\n', final_fig);
fprintf('Saved final state MAT: %s\n', out_mat);
fprintf('Increment PNG folder: %s\n', inc_dir);

if ~show_live
    close(fig);
end

end

function [R,Jff,Fint] = residual_and_jac_NR( ...
    x, Utrial, free, B3, EDOF, Ke_unit3, idx_glob, ...
    E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, Fext, ...
    I_trip, J_trip, eid_rep, Val_unit, ndof)

Uloc = Utrial;
Uloc(free) = x;

[~, s_e, Fint, ~] = damage_state_and_internal_force( ...
    Uloc, B3, EDOF, Ke_unit3, idx_glob, E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, ndof);

Rfull = Fint - Fext;
R = Rfull(free);

vals = Val_unit .* s_e(eid_rep);
K = sparse(I_trip, J_trip, vals, ndof, ndof);
K = 0.5*(K + K.');

Jff = K(free,free);

if ~isempty(Jff)
    dmax = double(full(max(abs(diag(Jff)))));
    if isempty(dmax) || dmax <= 0 || ~isfinite(dmax)
        dmax = 1.0;
    end
    Jff = Jff + speye(size(Jff))*1e-12*dmax;
end

end

function R = residual_only_NR( ...
    x, Utrial, free, B3, EDOF, Ke_unit3, idx_glob, ...
    E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, Fext, ndof)

Uloc = Utrial;
Uloc(free) = x;

[~, ~, Fint, ~] = damage_state_and_internal_force( ...
    Uloc, B3, EDOF, Ke_unit3, idx_glob, E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa, ndof);

Rfull = Fint - Fext;
R = Rfull(free);

end

function [De, s_e, Fint, eeq, h_band, crack_n] = damage_state_and_internal_force( ...
    U, B3, EDOF, Ke_unit3, idx_glob, E, nu, k_tc, j0, GF, gradN3, bandwidth_method, kappa_old, ndof)

ne = size(EDOF,2);

ue12 = reshape(U(EDOF), 12, 1, ne);
epsv = pagemtimes(B3, ue12);

[eeq, crack_n] = eqv_strain_modified_vm_vec(epsv, nu, k_tc);
eeq = reshape(eeq, [], 1);

h_band = oliver_bandwidth_TET4(gradN3, crack_n, bandwidth_method);
be = (E * j0 ./ GF) .* h_band(:);

kappa_trial = max(kappa_old, eeq);

De = zeros(ne,1);
active = kappa_trial >= j0;
De(active) = 1 - (j0 ./ kappa_trial(active)) .* exp(-be(active) .* (kappa_trial(active)-j0));
De = min(max(De,0),0.999999);

s_e = E .* max(1 - De, 1e-8);

fe12 = pagemtimes(Ke_unit3, ue12);
fe12 = bsxfun(@times, fe12, reshape(s_e,1,1,ne));

Fint = assemble_accum(idx_glob, fe12, ndof);

end

function [B_all,V_el,gradN_all,valid] = precompute_TET4_vectorized(p,T)

ne = size(T,1);

X1 = p(T(:,1),:);
X2 = p(T(:,2),:);
X3 = p(T(:,3),:);
X4 = p(T(:,4),:);

a = X2 - X1;
b = X3 - X1;
c = X4 - X1;

detJ = dot(a, cross(b,c,2), 2);
V_el = abs(detJ)/6;

valid = isfinite(detJ) & abs(detJ) > eps & V_el > eps;

safeDet = detJ;
bad = abs(safeDet) <= eps | ~isfinite(safeDet);
safeDet(bad) = eps;

g2 = cross(b,c,2) ./ safeDet;
g3 = cross(c,a,2) ./ safeDet;
g4 = cross(a,b,2) ./ safeDet;
g1 = -(g2 + g3 + g4);

gx = [g1(:,1), g2(:,1), g3(:,1), g4(:,1)];
gy = [g1(:,2), g2(:,2), g3(:,2), g4(:,2)];
gz = [g1(:,3), g2(:,3), g3(:,3), g4(:,3)];

B_all = zeros(6,12,ne);

for aNode = 1:4
    c0 = 3*(aNode-1);
    gxv = reshape(gx(:,aNode),1,1,ne);
    gyv = reshape(gy(:,aNode),1,1,ne);
    gzv = reshape(gz(:,aNode),1,1,ne);

    B_all(1,c0+1,:) = gxv;
    B_all(2,c0+2,:) = gyv;
    B_all(3,c0+3,:) = gzv;

    B_all(4,c0+1,:) = gyv;
    B_all(4,c0+2,:) = gxv;

    B_all(5,c0+2,:) = gzv;
    B_all(5,c0+3,:) = gyv;

    B_all(6,c0+1,:) = gzv;
    B_all(6,c0+3,:) = gxv;
end

V_el(~valid) = max(V_el(~valid), eps);

gradN_all = zeros(ne,4,3);
gradN_all(:,:,1) = gx;
gradN_all(:,:,2) = gy;
gradN_all(:,:,3) = gz;

end

function h_band = oliver_bandwidth_TET4(gradN_all, crack_n, bandwidth_method)

if nargin < 3 || isempty(bandwidth_method)
    bandwidth_method = "oliver";
end

bandwidth_method = lower(string(bandwidth_method));
if bandwidth_method ~= "oliver"
    warning('Unknown bandwidth_method "%s". Using Oliver directional bandwidth.', char(bandwidth_method));
end

crack_n = double(crack_n);
normal_norm = sqrt(sum(crack_n.^2,2));
bad_normal = ~isfinite(normal_norm) | normal_norm <= 1e-14;
normal_norm(bad_normal) = 1.0;
crack_n = crack_n ./ normal_norm;
crack_n(bad_normal,:) = repmat([1 0 0], nnz(bad_normal), 1);

gx = gradN_all(:,:,1);
gy = gradN_all(:,:,2);
gz = gradN_all(:,:,3);

proj = gx .* crack_n(:,1) + gy .* crack_n(:,2) + gz .* crack_n(:,3);
denom = sum(abs(proj),2);

h_band = 2.0 ./ max(denom, eps);

end

function EDOF = build_edof(T)

ne = size(T,1);
EDOF = zeros(12,ne);

for a = 1:4
    rows = 3*(a-1) + (1:3);
    n = T(:,a).';
    EDOF(rows(1),:) = 3*n - 2;
    EDOF(rows(2),:) = 3*n - 1;
    EDOF(rows(3),:) = 3*n;
end

end

function [I_trip, J_trip, Val_unit, eid_rep] = build_triplet_pattern(EDOF, Ke_unit3)

ne = size(EDOF,2);

I_blk = repmat(EDOF, 12, 1);
J_blk = repelem(EDOF, 12, 1);

I_trip = I_blk(:);
J_trip = J_blk(:);

Val_pages = reshape(Ke_unit3, 144, ne);
Val_unit = Val_pages(:);

eid_rep = repelem((1:ne)', 144);

end

function [D,lambda,mu] = iso3D_D(E,nu)

mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

D = [lambda+2*mu, lambda,      lambda,      0,  0,  0;
     lambda,      lambda+2*mu, lambda,      0,  0,  0;
     lambda,      lambda,      lambda+2*mu, 0,  0,  0;
     0,           0,           0,           mu, 0,  0;
     0,           0,           0,           0,  mu, 0;
     0,           0,           0,           0,  0,  mu];

end

function [eeq, crack_n] = eqv_strain_modified_vm_vec(epsv,nu,k)

exx = reshape(epsv(1,1,:),[],1);
eyy = reshape(epsv(2,1,:),[],1);
ezz = reshape(epsv(3,1,:),[],1);
gxy = reshape(epsv(4,1,:),[],1);
gyz = reshape(epsv(5,1,:),[],1);
gxz = reshape(epsv(6,1,:),[],1);

exy = 0.5*gxy;
eyz = 0.5*gyz;
exz = 0.5*gxz;

I1 = exx + eyy + ezz;
em = I1/3;

dxx = exx - em;
dyy = eyy - em;
dzz = ezz - em;

J2 = 0.5*(dxx.^2 + dyy.^2 + dzz.^2 + 2*(exy.^2 + eyz.^2 + exz.^2));

denom = max(abs(1-2*nu),1e-12);
term1 = ((k-1)/(2*k*denom)) .* I1;

rad = (((k-1)/denom).*I1).^2 + (12*k/((1+nu)^2)).*J2;
rad = max(rad,0);

term2 = (1/(2*k)) .* sqrt(rad);

eeq = max(0, term1 + term2);

[~, crack_n] = max_principal_strain_direction_vec(exx, eyy, ezz, exy, eyz, exz);

end

function [lambda1, n1] = max_principal_strain_direction_vec(exx, eyy, ezz, exy, eyz, exz)

ne = numel(exx);

p1 = exy.^2 + eyz.^2 + exz.^2;
q  = (exx + eyy + ezz) / 3;

p2 = (exx-q).^2 + (eyy-q).^2 + (ezz-q).^2 + 2*p1;
p  = sqrt(max(p2,0) / 6);

p_safe = max(p, eps);
B11 = (exx - q) ./ p_safe;
B22 = (eyy - q) ./ p_safe;
B33 = (ezz - q) ./ p_safe;
B12 = exy ./ p_safe;
B23 = eyz ./ p_safe;
B13 = exz ./ p_safe;

r = 0.5 * (B11.*(B22.*B33 - B23.^2) ...
        - B12.*(B12.*B33 - B23.*B13) ...
        + B13.*(B12.*B23 - B22.*B13));
r = min(max(r,-1),1);
phi = acos(r) / 3;

lambda1 = q + 2*p.*cos(phi);

diag_case = p1 < 1e-28;
spherical = p < 1e-28;

a11 = exx - lambda1;
a22 = eyy - lambda1;
a33 = ezz - lambda1;

row1 = [a11, exy, exz];
row2 = [exy, a22, eyz];
row3 = [exz, eyz, a33];

c12 = cross(row1,row2,2);
c13 = cross(row1,row3,2);
c23 = cross(row2,row3,2);

n12 = sum(c12.^2,2);
n13 = sum(c13.^2,2);
n23 = sum(c23.^2,2);

n1 = c12;
use13 = n13 > n12 & n13 >= n23;
use23 = n23 > n12 & n23 > n13;
n1(use13,:) = c13(use13,:);
n1(use23,:) = c23(use23,:);

nrm = sqrt(sum(n1.^2,2));
bad = ~isfinite(nrm) | nrm < 1e-20 | diag_case | spherical;

if any(bad)
    vals = [exx, eyy, ezz];
    [~, imax] = max(vals, [], 2);
    n_fallback = zeros(ne,3);
    n_fallback(imax == 1,1) = 1;
    n_fallback(imax == 2,2) = 1;
    n_fallback(imax == 3,3) = 1;
    if any(spherical)
        n_fallback(spherical,:) = repmat([1 0 0], nnz(spherical), 1);
    end
    n1(bad,:) = n_fallback(bad,:);
    nrm(bad) = sqrt(sum(n1(bad,:).^2,2));
end

n1 = n1 ./ max(nrm,eps);

end

function F = assemble_accum(idx_glob, fe12, ndof)

F = accumarray(idx_glob(:), fe12(:), [ndof,1], @sum, 0, true);

end

function [did_save, png_name] = save_damage_geometry_snapshot( ...
    fig, ax, hHull, hDam, p, U, def_scale, allF, owner_per_face, De, ...
    plot_thr, inc, nIncr, lambda_load, out_folder, prefix, image_res, save_when_empty, show_mesh_in_saved, damage_clim, damage_clim_mode)

if nargin < 19 || isempty(show_mesh_in_saved)
    show_mesh_in_saved = true;
end
if nargin < 20 || isempty(damage_clim)
    damage_clim = [0 1];
end
if nargin < 21 || isempty(damage_clim_mode)
    damage_clim_mode = "visible";
end

face_damage = De(owner_per_face);
keep_faces = face_damage >= plot_thr;
nfaces = nnz(keep_faces);

did_save = false;
png_name = '';

if nfaces == 0 && ~save_when_empty
    return;
end

oldHullVis = '';
oldClim = get_current_clim(ax);

Ux = U(1:3:end);
Uy = U(2:3:end);
Uz = U(3:3:end);
P_def = p + def_scale*[Ux,Uy,Uz];

if nargin >= 3 && ~isempty(hHull) && isgraphics(hHull)
    oldHullVis = get(hHull,'Visible');
    set(hHull,'Vertices',P_def);
    if show_mesh_in_saved
        set(hHull,'Visible','on');
    else
        set(hHull,'Visible','off');
    end
end

if nfaces > 0
    crack_faces = allF(keep_faces,:);
    crack_colors = face_damage(keep_faces);

    set(hDam, ...
        'Vertices', P_def, ...
        'Faces', crack_faces, ...
        'FaceVertexCData', crack_colors, ...
        'Visible','on');
else
    set(hDam, ...
        'Vertices', P_def, ...
        'Faces', zeros(0,3), ...
        'FaceVertexCData', zeros(0,1), ...
        'Visible','off');
end

if nfaces > 0
    apply_damage_color_limits(ax, face_damage(keep_faces), plot_thr, damage_clim, damage_clim_mode);
else
    apply_damage_color_limits(ax, [], plot_thr, damage_clim, damage_clim_mode);
end

title(ax, sprintf('Damage with mesh | inc=%d/%d | lambda=%.4f | damaged faces=%d | omega >= %.3f', ...
    inc, nIncr, lambda_load, nfaces, plot_thr), ...
    'FontSize',12,'FontWeight','bold','Interpreter','none');

view(ax,2);
camproj(ax,'orthographic');
camup(ax,[0 1 0]);
drawnow;

png_name = fullfile(out_folder, sprintf('%s_inc_%04d_damageMesh_omega_ge_%s.png', ...
    prefix, inc, fmt_underscore(plot_thr)));

exportgraphics(fig, png_name, 'Resolution', image_res);

if ~isempty(oldHullVis) && isgraphics(hHull)
    set(hHull,'Visible',oldHullVis);
end

if ~isempty(oldClim)
    set_current_clim(ax, oldClim);
end

did_save = true;

end

function update_damage_geometry_live( ...
    fig, ax, hHull, hDam, p, U, def_scale, allF, owner_per_face, De, ...
    live_thr, inc, nIncr, lambda_load, maxD, nLive, show_mesh_live, damage_clim, damage_clim_mode)

if nargin < 17 || isempty(show_mesh_live)
    show_mesh_live = true;
end
if nargin < 18 || isempty(damage_clim)
    damage_clim = [0 1];
end
if nargin < 19 || isempty(damage_clim_mode)
    damage_clim_mode = "visible";
end

if ~isgraphics(fig) || ~isgraphics(ax) || ~isgraphics(hDam)
    return;
end

face_damage = De(owner_per_face);
keep_faces = face_damage >= live_thr;
nfaces = nnz(keep_faces);

Ux = U(1:3:end);
Uy = U(2:3:end);
Uz = U(3:3:end);
P_def = p + def_scale*[Ux,Uy,Uz];

if ~isempty(hHull) && isgraphics(hHull)
    set(hHull, 'Vertices', P_def, 'Visible', ternary_vis(show_mesh_live));
end

if nfaces > 0
    set(hDam, ...
        'Vertices', P_def, ...
        'Faces', allF(keep_faces,:), ...
        'FaceVertexCData', face_damage(keep_faces), ...
        'Visible','on');
else
    set(hDam, ...
        'Vertices', P_def, ...
        'Faces', zeros(0,3), ...
        'FaceVertexCData', zeros(0,1), ...
        'Visible','off');
end

if nfaces > 0
    apply_damage_color_limits(ax, face_damage(keep_faces), live_thr, damage_clim, damage_clim_mode);
else
    apply_damage_color_limits(ax, [], live_thr, damage_clim, damage_clim_mode);
end

title(ax, sprintf('LIVE damage + mesh | inc=%d/%d | lambda=%.4f | maxD=%.4f | N(omega >= %.3f)=%d | faces=%d', ...
    inc, nIncr, lambda_load, maxD, live_thr, nLive, nfaces), ...
    'FontSize',12,'FontWeight','bold','Interpreter','none');

view(ax,2);
camproj(ax,'orthographic');
camup(ax,[0 1 0]);

end

function extFaces = external_hull(T)

allF = [T(:,[1 2 3]); T(:,[1 2 4]); T(:,[1 3 4]); T(:,[2 3 4])];
sF = sort(allF,2);
[uf,~,ic] = unique(sF,'rows');
cnt = accumarray(ic,1);
extFaces = uf(cnt == 1,:);

end

function try_set(h, prop, val)

if isempty(h) || ~isgraphics(h)
    return;
end

try
    set(h, prop, val);
catch

end

end

function v = ternary_vis(tf)

if tf
    v = 'on';
else
    v = 'off';
end

end

function nodes = read_nodeset_local(file_name, ids, set_name)

raw = readmatrix(file_name);
raw = raw(:);
raw = raw(isfinite(raw));

[tf,loc] = ismember(raw, ids);
if any(~tf)
    warning('%d nodes in %s set were not found in nodes file and were ignored.', nnz(~tf), set_name);
end

nodes = uint32(loc(tf));
nodes = nodes(nodes > 0);

end

function m = safe_mean(x)

if isempty(x)
    m = NaN;
else
    m = double(full(mean(x)));
end

end

function c = full_damage_cmap()

n = 256;
stops = [0.15 0.20 0.85;
         0.00 0.70 0.90;
         0.12 0.80 0.30;
         0.98 0.90 0.05;
         0.98 0.45 0.03;
         0.82 0.04 0.04];

pos = [0, 0.25, 0.50, 0.75, 0.90, 1.00];

t = linspace(0,1,n).';
c = zeros(n,3);

for ch = 1:3
    c(:,ch) = interp1(pos, stops(:,ch), t, 'pchip');
end

c = min(max(c,0),1);

end

function c = select_damage_cmap(name, n_bands, make_discrete)

if nargin < 1 || isempty(name)
    name = 'turbo';
end
if nargin < 2 || isempty(n_bands)
    n_bands = 12;
end
if nargin < 3 || isempty(make_discrete)
    make_discrete = false;
end

name = lower(string(name));
n = 256;

switch name
    case "turbo"
        c = builtin_colormap_safe('turbo', n);
    case "parula"
        c = builtin_colormap_safe('parula', n);
    case "jet"
        c = builtin_colormap_safe('jet', n);
    case "hot"
        c = builtin_colormap_safe('hot', n);
    case "cool"
        c = builtin_colormap_safe('cool', n);
    case "spring"
        c = builtin_colormap_safe('spring', n);
    case "autumn"
        c = builtin_colormap_safe('autumn', n);
    case "winter"
        c = builtin_colormap_safe('winter', n);
    case "fire"
        c = interp_colormap([ ...
            0.02 0.02 0.08;
            0.15 0.05 0.40;
            0.80 0.05 0.05;
            1.00 0.55 0.00;
            1.00 0.95 0.05], n);
    case "yellow_red"
        c = interp_colormap([ ...
            1.00 1.00 0.75;
            1.00 0.85 0.05;
            1.00 0.45 0.00;
            0.85 0.05 0.02;
            0.35 0.00 0.00], n);
    case "blue_green_yellow_red"
        c = interp_colormap([ ...
            0.05 0.10 0.75;
            0.00 0.70 0.85;
            0.10 0.75 0.25;
            1.00 0.92 0.05;
            1.00 0.42 0.00;
            0.80 0.00 0.00], n);
    case "damage_bands"
        c = interp_colormap([ ...
            0.05 0.15 0.95;
            0.00 0.70 1.00;
            0.00 0.75 0.20;
            0.95 0.95 0.05;
            1.00 0.50 0.00;
            0.90 0.00 0.00;
            0.35 0.00 0.00], max(n_bands,2));
        return;
    otherwise
        warning('Unknown damage_colormap "%s". Using turbo.', char(name));
        c = builtin_colormap_safe('turbo', n);
end

if make_discrete
    c = discretize_colormap(c, n_bands);
end

end

function c = builtin_colormap_safe(name, n)

try
    f = str2func(char(name));
    c = f(n);
catch
    if strcmpi(char(name),'turbo')
        c = jet(n);
    else
        c = full_damage_cmap();
    end
end

end

function c = interp_colormap(stops, n)

m = size(stops,1);
pos = linspace(0,1,m);
t = linspace(0,1,n).';
c = zeros(n,3);

for ch = 1:3
    c(:,ch) = interp1(pos, stops(:,ch), t, 'pchip');
end

c = min(max(c,0),1);

end

function c = discretize_colormap(c0, n_bands)

n_bands = max(2, round(n_bands));
idx = round(linspace(1, size(c0,1), n_bands));
c = c0(idx,:);

end

function apply_damage_color_limits(ax, visible_damage, threshold, default_clim, mode)

mode = lower(string(mode));

switch mode
    case "full"
        lim = default_clim;

    case "visible"
        lim = [threshold, default_clim(2)];

    case "auto"
        visible_damage = visible_damage(isfinite(visible_damage));
        if isempty(visible_damage)
            lim = default_clim;
        else
            dmin = min(visible_damage);
            dmax = max(visible_damage);
            if abs(dmax - dmin) < 1e-8
                dmin = max(0, dmin - 0.025);
                dmax = min(1, dmax + 0.025);
            end
            lim = [dmin, dmax];
        end

    otherwise
        lim = [threshold, default_clim(2)];
end

if numel(lim) ~= 2 || any(~isfinite(lim)) || lim(2) <= lim(1)
    lim = [0 1];
end

lim(1) = max(0, lim(1));
lim(2) = min(1, lim(2));
if lim(2) <= lim(1)
    lim = [0 1];
end

set_current_clim(ax, lim);

end

function lim = get_current_clim(ax)

try
    lim = clim(ax);
catch
    lim = caxis(ax);
end

end

function set_current_clim(ax, lim)

try
    clim(ax, lim);
catch
    caxis(ax, lim);
end

end

function val = local_get(opts, field, default)

if isstruct(opts) && isfield(opts,field) && ~isempty(opts.(field))
    val = opts.(field);
else
    val = default;
end

end

function s = fmt_underscore(x)

s = sprintf('%.4e', x);
s = strrep(s, '+', '');
s = strrep(s, '-', 'm');
s = strrep(s, '.', 'p');

end
