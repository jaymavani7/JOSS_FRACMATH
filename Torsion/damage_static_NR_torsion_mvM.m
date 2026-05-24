% =========================================================================
% LIVE VISUALIZATION IMPLICIT CDM SOLVER (MODIFIED NEWTON-RAPHSON)
% Crack Band Model + Oliver's Method + Internal Crack Visualization
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
crack_thresh    = 0.80; % Only show internal elements with damage > 80%
hull_alpha      = 0.15; % Transparency of the undamaged outer mesh (0 to 1)

%% 2. LOAD MESH AND NODE FILES
fprintf('Loading mesh and node sets...\n');
nodes       = readmatrix('Job-2_nodes.txt'); 
elements    = readmatrix('Job-2_elements.txt');
N_left      = readmatrix('Job-2_left_nodes.txt');
N_right     = readmatrix('Job-2_right_nodes.txt');
N_cmod1     = readmatrix('Job-2_CMOD1_nodes.txt');
N_cmod2     = readmatrix('Job-2_CMOD2_nodes.txt');

p = nodes(:, 2:4);       % Coordinates
t = elements(:, 2:5);    % Connectivity
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
    J = [x2-x1, x3-x1, x4-x1];
    V(e) = abs(det(J))/6;
    
    JTinv = (J') \ eye(3);
    g2 = JTinv(:,1); g3 = JTinv(:,2); g4 = JTinv(:,3); g1 = -(g2+g3+g4);
    grads = [g1 g2 g3 g4];
    
    B = zeros(6,12);
    for a=1:4
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

fix_left = [dof(N_left,1); dof(N_left,2); dof(N_left,3)];
presc_right = dof(N_right, load_dir);
fix_right   = setdiff([dof(N_right,1); dof(N_right,2); dof(N_right,3)], presc_right); 

dir_dofs  = unique([fix_left; presc_right; fix_right]);
free_dofs = setdiff(1:ndof, dir_dofs)';

%% 5. LIVE VISUALIZATION SETUP (INTERNAL CRACK VIEW)
fprintf('Setting up live graphics...\n');

allF = [t(:,[1 2 3]); t(:,[1 2 4]); t(:,[1 3 4]); t(:,[2 3 4])];
owner_per_face = repmat((1:ne)', 4, 1);

[sF, ~] = sort(allF, 2);
[uF, ~, iC] = unique(sF, 'rows');
counts = accumarray(iC, 1);
extFaces = uF(counts == 1, :);

fig = figure('Color','w','Name','Live Crack Simulation','Units','normalized','Position',[0.1 0.1 0.8 0.6]);
tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

axD = nexttile(tl, 1); 
hold(axD, 'on'); axis(axD, 'equal', 'vis3d');
colormap(axD, flipud(hot)); 

hHull = patch(axD, 'Faces', extFaces, 'Vertices', p, ...
              'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', ...
              'FaceAlpha', hull_alpha, 'AmbientStrength', 0.6);

hCrack = patch(axD, 'Faces', [], 'Vertices', p, ...
               'FaceVertexCData', [], 'FaceColor', 'flat', ...
               'EdgeColor', 'none', 'AmbientStrength', 0.8);

view(axD, [-30 30]); 
camlight(axD, 'headlight'); lighting(axD, 'gouraud');
title(axD, sprintf('Internal Damage Map (d > %.2f)', crack_thresh));
caxis(axD, [0 1]); 
colorbar(axD, 'Location', 'westoutside');
grid(axD, 'on'); box(axD, 'on');

axG = nexttile(tl, 2); 
hold(axG, 'on'); grid(axG, 'on'); box(axG, 'on');
xlabel(axG, 'Crack Mouth Opening Displacement (CMOD) [mm]', 'FontSize', 11);
ylabel(axG, 'Total Reaction Force [N]', 'FontSize', 11);
title(axG, 'Live Load vs. CMOD');
hLine = plot(axG, 0, 0, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);

drawnow; 

%% 6. MODIFIED NEWTON-RAPHSON SOLVER
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
        J2 = 0.5*((exx-I1/3).^2 + (eyy-I1/3).^2 + (ezz-I1/3).^2 + 0.5*(gxy.^2 + gyz.^2 + gzx.^2));
        
        a1 = (k_tc-1)/(1-2*nu); a2 = 12*k_tc/(1+nu)^2;
        eeq = (k_tc-1)/(2*k_tc*(1-2*nu)).*I1 + (1/(2*k_tc)).*sqrt( (a1^2).*I1.^2 + a2.*J2 );
        
        kappa = max(kappa, eeq);
        act = (kappa > kappa0);
        
        denom = max(eps_f(act) - kappa0, 1e-18);
        damage(act) = 1.0 - (kappa0 ./ kappa(act)) .* exp(-(kappa(act) - kappa0) ./ denom);
        damage = min(max(damage, 0), 0.999999999999999);
        
        se = bsxfun(@times, pagemtimes(repmat(D0, [1 1 ne]), epsv), reshape(1-damage, 1, 1, ne));
        fe12 = bsxfun(@times, pagemtimes(permute(B3,[2 1 3]), se), reshape(V, 1, 1, ne));
        F_int = accumarray(EDOF(:), fe12(:), [ndof, 1]);
        
        R = -F_int; 
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
    
    fprintf('Step %3d/%d | Iters: %2d | Load: %8.2f N | CMOD: %.4f mm\n', ...
        step, n_steps, iter, Total_Load, CMOD_val);

    % --- STOP CONDITION ---
    if CMOD_val >= 0.1
        fprintf('\nCMOD reached %.4f mm at step %d. Stopping.\n', CMOD_val, step);
        history_load = history_load(1:step);
        history_cmod = history_cmod(1:step);
        
        % Final plot update before exit
        if isvalid(fig)
            set(hLine, 'XData', history_cmod, 'YData', history_load);
            Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
            P_def = p + def_scale * [Ux, Uy, Uz];
            set(hHull, 'Vertices', P_def);
            is_cracked = damage > crack_thresh;
            if any(is_cracked)
                crack_mask = repmat(is_cracked, 4, 1); 
                crack_faces = allF(crack_mask, :);
                crack_colors = damage(owner_per_face(crack_mask));
                set(hCrack, 'Vertices', P_def, 'Faces', crack_faces, 'FaceVertexCData', crack_colors);
            end
            drawnow;
        end
        break;
    end
    
    % --- LIVE VISUALIZATION UPDATE ---
    if isvalid(fig)
        set(hLine, 'XData', history_cmod(1:step), 'YData', history_load(1:step));
        
        Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
        P_def = p + def_scale * [Ux, Uy, Uz];
        
        set(hHull, 'Vertices', P_def);
        
        is_cracked = damage > crack_thresh;
        if any(is_cracked)
            crack_mask = repmat(is_cracked, 4, 1); 
            crack_faces = allF(crack_mask, :);
            crack_colors = damage(owner_per_face(crack_mask));
            set(hCrack, 'Vertices', P_def, 'Faces', crack_faces, 'FaceVertexCData', crack_colors);
        end
        
        drawnow; 
    end
end
fprintf('\nSimulation Complete.\n');