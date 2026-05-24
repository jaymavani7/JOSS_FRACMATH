% =========================================================================
% ENHANCED LIVE VISUALIZATION IMPLICIT CDM SOLVER (MODIFIED NEWTON-RAPHSON)
% Crack Band Model + Oliver's Method + Multi-Angle Crack View + Video Export
% =========================================================================
clc; clear; close all;

%% 1. MATERIAL AND SOLVER SETTINGS
E       = 35000;        % Young's Modulus (MPa)
nu      = 0.20;         % Poisson's ratio
ft      = 3.0;          % Tensile strength (MPa)
GF      = 0.08;         % Fracture energy (N/mm)
k_tc    = 10.0;         % Compressive/Tensile strength ratio
kappa0  = ft / E;       % Damage threshold

% Control settings
Uy_end      = 0.30;     % Total prescribed displacement (mm)
n_steps     = 1000;     % Number of load increments
max_iter    = 50;       % Max iterations per load step for mNR
tol         = 1e-4;     % Convergence tolerance (Force norm)
load_dir    = 1;        % 1=X, 2=Y, 3=Z direction for pulling

% Visualization settings
def_scale       = 15.0; % Scale factor for 3D deformation visual
crack_thresh    = 0.9999; % Only show internal elements with damage > 80%
hull_alpha      = 0.12; % Transparency of the undamaged outer mesh

% ---- OUTPUT SETTINGS ----
output_folder   = 'CDM_Output';     % Folder for saved images & video
save_video      = true;             % Save MP4 video of simulation
video_fps       = 15;               % Frames per second in output video
% Snapshot at these fraction-of-max-load thresholds (e.g. 0.25 = 25% of peak)
snap_load_fracs = [0.25, 0.50, 0.75, 0.90, 1.00];

if ~exist(output_folder, 'dir'), mkdir(output_folder); end

%% 2. LOAD MESH AND NODE FILES
fprintf('Loading mesh and node sets...\n');
nodes       = readmatrix('Job-3_nodes.txt'); 
elements    = readmatrix('Job-3_elements.txt');
N_left      = readmatrix('Job-3_left_nodes.txt');
N_right     = readmatrix('Job-3_right_nodes.txt');
N_cmod1     = readmatrix('Job-3_CMOD1_nodes.txt');
N_cmod2     = readmatrix('Job-3_CMOD2_nodes.txt');

p  = nodes(:, 2:4);
t  = elements(:, 2:5);
np = size(p, 1);
ne = size(t, 1);
ndof = 3 * np;

%% 3. PRECOMPUTATION (FULLY VECTORIZED)
fprintf('Precomputing B-matrices, volumes, and h_e (Oliver method)...\n');

edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
he = zeros(ne, 1);
for i = 1:6
    edge_len = vecnorm(p(t(:, edges(i,1)), :) - p(t(:, edges(i,2)), :), 2, 2);
    he = max(he, edge_len);
end

eps_f = GF ./ (he .* ft) + 0.5 * kappa0;
eps_f = max(eps_f, kappa0 * 1.01); 

B3 = zeros(6, 12, ne);
V  = zeros(ne, 1);
for e = 1:ne
    idx = t(e, :); Xe = p(idx, :);
    x1 = Xe(1,:)'; x2 = Xe(2,:)'; x3 = Xe(3,:)'; x4 = Xe(4,:)';
    J  = [x2-x1, x3-x1, x4-x1];
    V(e) = abs(det(J))/6;
    JTinv = (J') \ eye(3);
    g2 = JTinv(:,1); g3 = JTinv(:,2); g4 = JTinv(:,3); g1 = -(g2+g3+g4);
    grads = [g1 g2 g3 g4];
    B = zeros(6,12);
    for a = 1:4
        c = (3*(a-1)+1):(3*a);
        B(1,c)=[grads(1,a) 0 0]; B(2,c)=[0 grads(2,a) 0]; B(3,c)=[0 0 grads(3,a)];
        B(4,c)=[grads(2,a) grads(1,a) 0]; B(5,c)=[0 grads(3,a) grads(2,a)]; B(6,c)=[grads(3,a) 0 grads(1,a)];
    end
    B3(:,:,e) = B;
end

G = E / (2*(1+nu)); lambda = E*nu / ((1+nu)*(1-2*nu));
D0 = [lambda+2*G, lambda, lambda, 0, 0, 0;
      lambda, lambda+2*G, lambda, 0, 0, 0;
      lambda, lambda, lambda+2*G, 0, 0, 0;
      0, 0, 0, G, 0, 0;
      0, 0, 0, 0, G, 0;
      0, 0, 0, 0, 0, G];

Ke_elastic = pagemtimes(permute(B3,[2 1 3]), pagemtimes(repmat(D0,[1 1 ne]), B3));
Ke_elastic = bsxfun(@times, Ke_elastic, reshape(V, 1, 1, ne));

EDOF = zeros(12, ne);
for e = 1:ne  
    n4 = t(e,:);
    EDOF(:,e) = [3*n4(1)-2:3*n4(1), 3*n4(2)-2:3*n4(2), 3*n4(3)-2:3*n4(3), 3*n4(4)-2:3*n4(4)]';
end
I_idx = repmat(EDOF, 12, 1);
J_idx = repelem(EDOF, 12, 1);

%% 4. BOUNDARY CONDITIONS
dof = @(nid, comp) 3*(nid-1) + comp;

fix_left    = [dof(N_left,1); dof(N_left,2); dof(N_left,3)];
presc_right = dof(N_right, load_dir);
fix_right   = setdiff([dof(N_right,1); dof(N_right,2); dof(N_right,3)], presc_right); 

dir_dofs  = unique([fix_left; presc_right; fix_right]);
free_dofs = setdiff(1:ndof, dir_dofs)';

%% 5. SINGLE CRACK COLOUR
% Fully-damaged elements are shown in this one solid colour only.
% Change the RGB values here to pick a different colour.
crack_color = [1.00, 0.18, 0.05];   % vivid red  (try [1 0.6 0] for orange)

%% 6. SURFACE DATA FOR OUTER HULL
allF = [t(:,[1 2 3]); t(:,[1 2 4]); t(:,[1 3 4]); t(:,[2 3 4])];
owner_per_face = repmat((1:ne)', 4, 1);

[sF, ~] = sort(allF, 2);
[uF, ~, iC] = unique(sF, 'rows');
counts = accumarray(iC, 1);
extFaces = uF(counts == 1, :);

%% 7. FIGURE LAYOUT  (4 views + Load-CMOD panel)
fprintf('Setting up multi-angle live figure...\n');

fig = figure('Color',[0.08 0.08 0.12], 'Name','Live Crack Simulation - Multi-Angle', ...
             'Units','normalized', 'Position',[0.02 0.04 0.96 0.92]);

tl = tiledlayout(fig, 2, 3, 'Padding','tight', 'TileSpacing','compact');

% Camera angles: [azimuth, elevation]
view_angles = [-30 25;   % Isometric
                0   0;   % Front
               90   0;   % Side
                0  90];  % Top
view_names  = {'Isometric', 'Front (Y-Z)', 'Side (X-Z)', 'Top (X-Y)'};

bg   = [0.08 0.08 0.12];
grid_col = [0.25 0.25 0.30];

axV = gobjects(4,1);
hHullArr  = gobjects(4,1);
hCrackArr = gobjects(4,1);

for v = 1:4
    axV(v) = nexttile(tl, v);
    set(axV(v), 'Color', bg, 'XColor', grid_col, 'YColor', grid_col, 'ZColor', grid_col, ...
                'GridColor', grid_col, 'GridAlpha', 0.4);
    hold(axV(v), 'on'); axis(axV(v), 'equal', 'vis3d');

    % Outer ghost hull — undamaged body
    hHullArr(v) = patch(axV(v), 'Faces', extFaces, 'Vertices', p, ...
                  'FaceColor', [0.55 0.60 0.65], 'EdgeColor', 'none', ...
                  'FaceAlpha', hull_alpha, 'AmbientStrength', 0.5);

    % Crack patch — single solid colour, no colormap needed
    hCrackArr(v) = patch(axV(v), 'Faces', [], 'Vertices', p, ...
                   'FaceColor', crack_color, ...
                   'EdgeColor', 'none', 'AmbientStrength', 0.85, 'SpecularStrength', 0.5);

    view(axV(v), view_angles(v,1), view_angles(v,2));
    camlight(axV(v), 'headlight'); lighting(axV(v), 'gouraud');

    title(axV(v), view_names{v}, 'Color','w', 'FontSize', 10, 'FontWeight','bold');
    grid(axV(v), 'on'); box(axV(v), 'on');
end

% Load-CMOD panel (spans bottom 2 tiles)
axG = nexttile(tl, 5, [1 2]);
set(axG, 'Color', [0.10 0.10 0.15], 'XColor','w', 'YColor','w', ...
         'GridColor',[0.3 0.3 0.4], 'GridAlpha',0.5, 'FontSize',10);
hold(axG,'on'); grid(axG,'on'); box(axG,'on');
xlabel(axG, 'Crack Mouth Opening Displacement  (CMOD)  [mm]', 'Color','w', 'FontSize',11);
ylabel(axG, 'Total Reaction Force  [N]',                       'Color','w', 'FontSize',11);
title(axG,  'Load vs. CMOD', 'Color','w', 'FontSize',12, 'FontWeight','bold');
hLine = plot(axG, 0, 0, '-', 'LineWidth', 2.5, 'Color',[0.25 0.75 1.0]);
% Snapshot markers
hSnapPts = plot(axG, nan, nan, 'o', 'MarkerSize',9, 'MarkerFaceColor',[1 0.4 0.1], ...
                'MarkerEdgeColor','w', 'LineStyle','none');

% Info text overlay on isometric view
hInfo = text(axV(1), 0, 0, 0, '', 'Color','w', 'FontSize', 9, ...
             'BackgroundColor',[0 0 0 0.5], 'Margin', 4, ...
             'VerticalAlignment','top', 'HorizontalAlignment','left', ...
             'Units','normalized', 'Position',[0.02 0.97]);

drawnow;

%% 8. VIDEO WRITER SETUP
if save_video
    vid_path = fullfile(output_folder, 'CDM_simulation.mp4');
    vidObj = VideoWriter(vid_path, 'MPEG-4');
    vidObj.FrameRate = video_fps;
    vidObj.Quality   = 95;
    open(vidObj);
    fprintf('Video writer opened: %s\n', vid_path);
end

%% 9. SNAPSHOT STATE TRACKING
snap_taken   = false(size(snap_load_fracs));
snap_cmods   = [];
snap_loads   = [];
peak_load    = 0;   % updated each step; used for fraction-based snapshots

%% 10. MODIFIED NEWTON-RAPHSON SOLVER
fprintf('Starting Modified Newton-Raphson Implicit Solver...\n');

U       = zeros(ndof, 1);
kappa   = kappa0 * ones(ne, 1);
damage  = zeros(ne, 1);

history_load = zeros(n_steps, 1);
history_cmod = zeros(n_steps, 1);
dU_step = Uy_end / n_steps;

for step = 1:n_steps
    
    U(presc_right) = U(presc_right) + dU_step;
    
    Ke_sec = bsxfun(@times, Ke_elastic, reshape(1 - damage, 1, 1, ne));
    K_global = sparse(I_idx(:), J_idx(:), Ke_sec(:), ndof, ndof);
    K_ff = K_global(free_dofs, free_dofs);
    L_solver = decomposition(K_ff, 'chol', 'lower');
    
    iter = 0; res_norm = 1.0;
    
    while res_norm > tol && iter < max_iter
        iter = iter + 1;
        ue12 = reshape(U(EDOF), 12, 1, ne);
        epsv = pagemtimes(B3, ue12);
        
        exx=squeeze(epsv(1,1,:)); eyy=squeeze(epsv(2,1,:)); ezz=squeeze(epsv(3,1,:));
        gxy=squeeze(epsv(4,1,:)); gyz=squeeze(epsv(5,1,:)); gzx=squeeze(epsv(6,1,:));
        I1 = exx + eyy + ezz; 
        J2 = 0.5*((exx-I1/3).^2 + (eyy-I1/3).^2 + (ezz-I1/3).^2 + ...
                   0.5*(gxy.^2 + gyz.^2 + gzx.^2));
        
        a1 = (k_tc-1)/(1-2*nu); a2 = 12*k_tc/(1+nu)^2;
        eeq = (k_tc-1)/(2*k_tc*(1-2*nu)).*I1 + ...
              (1/(2*k_tc)).*sqrt( (a1^2).*I1.^2 + a2.*J2 );
        
        kappa = max(kappa, eeq);
        act   = (kappa > kappa0);
        
        denom = max(eps_f(act) - kappa0, 1e-18);
        damage(act) = 1.0 - (kappa0 ./ kappa(act)) .* ...
                      exp(-(kappa(act) - kappa0) ./ denom);
        damage = min(max(damage, 0), 0.999999999999999);
        
        se   = bsxfun(@times, pagemtimes(repmat(D0,[1 1 ne]), epsv), reshape(1-damage,1,1,ne));
        fe12 = bsxfun(@times, pagemtimes(permute(B3,[2 1 3]), se), reshape(V,1,1,ne));
        F_int = accumarray(EDOF(:), fe12(:), [ndof, 1]);
        
        R      = -F_int; 
        R_free = R(free_dofs);
        res_norm = norm(R_free);
        if res_norm < tol, break; end
        
        dU_free = L_solver \ R_free;
        U(free_dofs) = U(free_dofs) + dU_free;
    end
    
    Total_Load = sum(F_int(presc_right));
    u_cmod1 = U(dof(N_cmod1(1), 1:3));
    u_cmod2 = U(dof(N_cmod2(1), 1:3));
    CMOD_val = norm(u_cmod1 - u_cmod2);
    
    history_load(step) = Total_Load;
    history_cmod(step) = CMOD_val;
    peak_load = max(peak_load, Total_Load);
    
    fprintf('Step %3d/%d | Iters: %2d | Load: %8.2f N | CMOD: %.4f mm | MaxDmg: %.4f\n', ...
        step, n_steps, iter, Total_Load, CMOD_val, max(damage));

    % ---- BUILD DEFORMED GEOMETRY ----
    Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
    P_def = p + def_scale * [Ux, Uy, Uz];
    
    is_cracked = damage > crack_thresh;
    crack_mask = repmat(is_cracked, 4, 1);
    
    if any(is_cracked)
        crack_faces = allF(crack_mask, :);
    end
    
    % ---- UPDATE ALL 4 VIEWS ----
    for v = 1:4
        set(hHullArr(v),  'Vertices', P_def);
        if any(is_cracked)
            % Single solid colour — no CData needed
            set(hCrackArr(v), 'Vertices', P_def, 'Faces', crack_faces);
        end
    end
    
    % ---- UPDATE LOAD-CMOD CURVE ----
    set(hLine, 'XData', history_cmod(1:step), 'YData', history_load(1:step));
    
    % ---- INFO TEXT ----
    set(hInfo, 'String', sprintf( ...
        'Step %d/%d\nIter: %d\nLoad: %.1f N\nCMOD: %.4f mm\nMax d: %.4f', ...
        step, n_steps, iter, Total_Load, CMOD_val, max(damage)));
    
    drawnow;

    % ---- FRACTION-BASED SNAPSHOTS ----
    if peak_load > 0
        for si = 1:numel(snap_load_fracs)
            if ~snap_taken(si) && (Total_Load >= snap_load_fracs(si) * peak_load)
                snap_taken(si) = true;
                snap_cmods(end+1) = CMOD_val;  %#ok<SAGROW>
                snap_loads(end+1) = Total_Load; %#ok<SAGROW>
                set(hSnapPts, 'XData', snap_cmods, 'YData', snap_loads);
                
                fname = fullfile(output_folder, ...
                    sprintf('snap_step%04d_load%.0f_CMOD%.4f.png', step, Total_Load, CMOD_val));
                exportgraphics(fig, fname, 'Resolution', 200);
                fprintf('  >> Snapshot saved: %s\n', fname);
            end
        end
    end
    
    % ---- WRITE VIDEO FRAME ----
    if save_video
        frame = getframe(fig);
        writeVideo(vidObj, frame);
    end

    % ---- STOP CONDITION ----
    if CMOD_val >= 0.1
        fprintf('\nCMOD reached %.4f mm at step %d. Stopping.\n', CMOD_val, step);
        history_load = history_load(1:step);
        history_cmod = history_cmod(1:step);
        break;
    end
end

%% 11. CLOSE VIDEO AND SAVE FINAL MULTI-ANGLE FIGURE
if save_video
    close(vidObj);
    fprintf('Video saved: %s\n', fullfile(output_folder, 'CDM_simulation.mp4'));
end

% Final high-res multi-view export
final_fig_path = fullfile(output_folder, 'final_multiview_crack.png');
exportgraphics(fig, final_fig_path, 'Resolution', 300);
fprintf('Final multi-view image saved: %s\n', final_fig_path);

%% 12. POST-PROCESSING: INDIVIDUAL ANGLE SNAPSHOTS (clean, no UI chrome)
fprintf('\nGenerating clean per-angle final snapshots...\n');

angle_labels = {'isometric','front','side','top'};
for v = 1:4
    figS = figure('Color',[0.08 0.08 0.12], 'Visible','off', ...
                  'Units','pixels', 'Position',[0 0 900 700]);
    axS  = axes(figS, 'Color',[0.08 0.08 0.12], ...
                'XColor',[0.3 0.3 0.4], 'YColor',[0.3 0.3 0.4], 'ZColor',[0.3 0.3 0.4]);
    hold(axS,'on'); axis(axS,'equal','vis3d');

    % Ghost outer hull
    patch(axS, 'Faces', extFaces, 'Vertices', P_def, ...
          'FaceColor',[0.55 0.60 0.65], 'EdgeColor','none', ...
          'FaceAlpha', hull_alpha, 'AmbientStrength',0.5);

    % Fully-damaged zone — single solid colour, no colormap/colorbar
    if any(is_cracked)
        patch(axS, 'Faces', crack_faces, 'Vertices', P_def, ...
              'FaceColor', crack_color, ...
              'EdgeColor','none', 'AmbientStrength',0.85, 'SpecularStrength',0.5);
    end

    view(axS, view_angles(v,1), view_angles(v,2));
    camlight(axS,'headlight'); lighting(axS,'gouraud');
    grid(axS,'on'); box(axS,'on');
    title(axS, sprintf('Fully Damaged Zone  —  %s view', view_names{v}), ...
          'Color','w', 'FontSize',13, 'FontWeight','bold');

    savename = fullfile(output_folder, sprintf('final_%s.png', angle_labels{v}));
    exportgraphics(figS, savename, 'Resolution', 300);
    close(figS);
    fprintf('  Saved: %s\n', savename);
end

fprintf('\n=== Simulation Complete. All outputs in: %s ===\n', output_folder);