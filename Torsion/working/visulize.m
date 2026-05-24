%% =========================================================
%  visualize_mesh_BC.m
%  Visualise a C3D4 (4-node tetrahedral) FE mesh and its
%  boundary conditions (left / right node and element sets).
%
%  Files expected in the same folder (or update paths below):
%    Job-1_nodes.txt          -> nodeID  x  y  z
%    Job-1_elements.txt       -> elemID  n1 n2 n3 n4
%    Job-1_left_nodes.txt     -> list of left-BC node IDs
%    Job-1_right_nodes.txt    -> list of right-BC node IDs
%    Job-1_left_elements.txt  -> list of left-BC element IDs
%    Job-1_right_elements.txt -> list of right-BC element IDs
%% =========================================================

clearvars; close all; clc;

%% ── 1. LOAD DATA ─────────────────────────────────────────
fprintf('Loading nodes ...\n');
raw_nodes    = readmatrix('Job-1_nodes.txt');          % [N x 4]
nodeID       = raw_nodes(:,1);
coords       = raw_nodes(:,2:4);                       % [N x 3] (x,y,z)

fprintf('Loading elements ...\n');
raw_elems    = readmatrix('Job-1_elements.txt');       % [E x 5]
elemID       = raw_elems(:,1);
conn_id      = raw_elems(:,2:5);                       % [E x 4] node IDs

fprintf('Loading boundary condition sets ...\n');
left_nodeIDs  = readmatrix('Job-1_left_nodes.txt');
right_nodeIDs = readmatrix('Job-1_right_nodes.txt');
left_elemIDs  = readmatrix('Job-1_left_elements.txt');
right_elemIDs = readmatrix('Job-1_right_elements.txt');

fprintf('  Nodes     : %d\n', size(coords,1));
fprintf('  Elements  : %d\n', size(conn_id,1));
fprintf('  Left  BCs : %d nodes | %d elements\n', ...
        numel(left_nodeIDs),  numel(left_elemIDs));
fprintf('  Right BCs : %d nodes | %d elements\n', ...
        numel(right_nodeIDs), numel(right_elemIDs));

%% ── 2. BUILD NODE-ID → ROW-INDEX MAP (vectorised) ───────
% Node IDs may be non-contiguous; use a sparse lookup vector.
maxID   = max(nodeID);
nodeMap = zeros(maxID, 1, 'uint32');
nodeMap(nodeID) = uint32(1:numel(nodeID));

%% ── 3. MAP CONNECTIVITY IDs → ROW INDICES (vectorised) ──
% conn_id contains raw node IDs; convert to row indices in 'coords'.
conn = nodeMap(conn_id);                               % [E x 4], row indices

%% ── 4. EXTRACT SURFACE FACES (vectorised) ────────────────
% Each C3D4 tet has 4 triangular faces:
%   f1: [n1 n2 n3]   f2: [n1 n2 n4]
%   f3: [n1 n3 n4]   f4: [n2 n3 n4]
%
% A face on the free surface is shared by exactly ONE element.
% Identify those faces by sorting each face's node triplet and
% using 'unique' with occurrence counting.

nE    = size(conn, 1);
face_local = [1 2 3;   % face 1
              1 2 4;   % face 2
              1 3 4;   % face 3
              2 3 4];  % face 4   [4 x 3]

% Build all 4*nE faces: [4*nE x 3] global row indices
all_faces = reshape(conn(:, face_local'), [], 3);
%   conn(:, face_local') -> [nE x 12], reshape to [4*nE x 3]

% Sort each row so face [a b c] == face [b a c] etc.
faces_sorted = sort(all_faces, 2);                     % [4*nE x 3]

% Find rows that appear exactly once -> boundary / surface faces
[~, ia, ic]  = unique(faces_sorted, 'rows', 'stable');
counts        = accumarray(ic, 1);                     % occurrence count
surf_mask     = counts(ic) == 1;                       % logical [4*nE x 1]
surf_faces    = all_faces(surf_mask, :);               % [nSurf x 3]

fprintf('  Surface triangles: %d\n', size(surf_faces,1));

%% ── 5. MAP BC IDs → ROW INDICES ─────────────────────────
left_node_idx  = nodeMap(left_nodeIDs);
right_node_idx = nodeMap(right_nodeIDs);

% Map BC element IDs to element row positions
maxEID = max(elemID);
elemMap = zeros(maxEID, 1, 'uint32');
elemMap(elemID) = uint32(1:numel(elemID));

left_elem_rows  = elemMap(left_elemIDs);
right_elem_rows = elemMap(right_elemIDs);

%% ── 6. EXTRACT BC ELEMENT FACES (for highlighting) ──────
% Build faces for left-BC elements only
left_conn  = conn(left_elem_rows,  :);                 % [nLe x 4]
right_conn = conn(right_elem_rows, :);                 % [nRe x 4]

left_faces  = reshape(left_conn(:, face_local'),  [], 3);
right_faces = reshape(right_conn(:, face_local'), [], 3);

% Keep only faces that are on the surface (appear once globally)
left_sorted  = sort(left_faces,  2);
right_sorted = sort(right_faces, 2);

% Use ismember to intersect with the global surface face set
[left_is_surf]  = ismember(left_sorted,  faces_sorted(surf_mask,:), 'rows');
[right_is_surf] = ismember(right_sorted, faces_sorted(surf_mask,:), 'rows');

left_surf_faces  = left_faces(left_is_surf,  :);
right_surf_faces = right_faces(right_is_surf, :);

%% ── 7. COLOUR MAP FOR SURFACE NORMALS (optional depth cue) ─
% Compute approximate z-centroid for each surface face for shading
x = coords(:,1);  y = coords(:,2);  z = coords(:,3);

surf_z = mean(z(surf_faces), 2);                       % [nSurf x 1]

%% ── 8. VISUALISE ─────────────────────────────────────────
fig = figure('Name','FE Mesh + BCs', 'Color','w', ...
             'Position',[100 100 1200 700]);

%% 8a. Full mesh surface (grey, semi-transparent)
ax1 = subplot(1,2,1);
p_mesh = patch('Faces',    surf_faces, ...
               'Vertices', coords,     ...
               'FaceVertexCData', surf_z, ...
               'FaceColor','flat',     ...
               'EdgeColor','none',     ...
               'FaceAlpha', 0.85);
colormap(ax1, gray(256));
axis equal tight;  view(3);
camlight headlight;  lighting gouraud;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Full Mesh (Surface)', 'FontSize',13);
grid on;

%% 8b. Mesh + BC highlighting
ax2 = subplot(1,2,2);

% Base mesh (grey, very transparent)
patch('Faces',    surf_faces, ...
      'Vertices', coords,     ...
      'FaceColor',[0.78 0.82 0.86], ...
      'EdgeColor','none',     ...
      'FaceAlpha', 0.25);
hold on;

% Left-BC surface faces (blue)
if ~isempty(left_surf_faces)
    patch('Faces',    left_surf_faces, ...
          'Vertices', coords,          ...
          'FaceColor',[0.13 0.47 0.71],...  % blue
          'EdgeColor','none',          ...
          'FaceAlpha', 0.90,           ...
          'DisplayName','Left BC faces');
end

% Right-BC surface faces (red)
if ~isempty(right_surf_faces)
    patch('Faces',    right_surf_faces, ...
          'Vertices', coords,           ...
          'FaceColor',[0.84 0.15 0.16], ...  % red
          'EdgeColor','none',           ...
          'FaceAlpha', 0.90,            ...
          'DisplayName','Right BC faces');
end

% Left-BC nodes (blue filled circles)
scatter3(x(left_node_idx), y(left_node_idx), z(left_node_idx), ...
         28, [0.13 0.47 0.71], 'filled',  ...
         'DisplayName','Left BC nodes', ...
         'MarkerEdgeColor','k', 'LineWidth',0.3);

% Right-BC nodes (red filled circles)
scatter3(x(right_node_idx), y(right_node_idx), z(right_node_idx), ...
         28, [0.84 0.15 0.16], 'filled',  ...
         'DisplayName','Right BC nodes', ...
         'MarkerEdgeColor','k', 'LineWidth',0.3);

axis equal tight;  view(3);
camlight headlight;  lighting gouraud;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Mesh + Boundary Conditions', 'FontSize',13);
legend('Location','best','FontSize',9);
grid on;

%% ── 9. SEPARATE CLOSE-UP FIGURE FOR EACH BC SET ─────────
figure('Name','BC Close-up', 'Color','w', 'Position',[150 150 1200 500]);

%% Left BC
ax3 = subplot(1,2,1);
patch('Faces',    surf_faces, ...
      'Vertices', coords,     ...
      'FaceColor',[0.88 0.90 0.92], ...
      'EdgeColor','none',     ...
      'FaceAlpha', 0.25);
hold on;
if ~isempty(left_surf_faces)
    patch('Faces',    left_surf_faces, ...
          'Vertices', coords,          ...
          'FaceColor',[0.13 0.47 0.71],...
          'EdgeColor','none',          ...
          'FaceAlpha', 0.90);
end
scatter3(x(left_node_idx), y(left_node_idx), z(left_node_idx), ...
         40, [0.13 0.47 0.71], 'filled', ...
         'MarkerEdgeColor','k','LineWidth',0.4);

% Zoom to left BC region
pad = 20;
xl = [min(x(left_node_idx))-pad, max(x(left_node_idx))+pad];
yl = [min(y(left_node_idx))-pad, max(y(left_node_idx))+pad];
zl = [min(z(left_node_idx))-pad, max(z(left_node_idx))+pad];
xlim(xl); ylim(yl); zlim(zl);
axis equal;  view(3);
camlight headlight; lighting gouraud;
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Left BC  (%d nodes, %d elem faces)', ...
      numel(left_node_idx), size(left_surf_faces,1)), 'FontSize',12);
grid on;

%% Right BC
ax4 = subplot(1,2,2);
patch('Faces',    surf_faces, ...
      'Vertices', coords,     ...
      'FaceColor',[0.88 0.90 0.92], ...
      'EdgeColor','none',     ...
      'FaceAlpha', 0.25);
hold on;
if ~isempty(right_surf_faces)
    patch('Faces',    right_surf_faces, ...
          'Vertices', coords,           ...
          'FaceColor',[0.84 0.15 0.16], ...
          'EdgeColor','none',           ...
          'FaceAlpha', 0.90);
end
scatter3(x(right_node_idx), y(right_node_idx), z(right_node_idx), ...
         40, [0.84 0.15 0.16], 'filled', ...
         'MarkerEdgeColor','k','LineWidth',0.4);

xr = [min(x(right_node_idx))-pad, max(x(right_node_idx))+pad];
yr = [min(y(right_node_idx))-pad, max(y(right_node_idx))+pad];
zr = [min(z(right_node_idx))-pad, max(z(right_node_idx))+pad];
xlim(xr); ylim(yr); zlim(zr);
axis equal;  view(3);
camlight headlight; lighting gouraud;
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Right BC  (%d nodes, %d elem faces)', ...
      numel(right_node_idx), size(right_surf_faces,1)), 'FontSize',12);
grid on;

%% ── 10. PRINT SUMMARY ────────────────────────────────────
fprintf('\n── Mesh Summary ─────────────────────────────────\n');
fprintf('  Total nodes     : %d\n',  numel(nodeID));
fprintf('  Total elements  : %d\n',  numel(elemID));
fprintf('  Surface faces   : %d\n',  size(surf_faces,1));
fprintf('  Bounding box\n');
fprintf('    X : [%.3f, %.3f]\n', min(x), max(x));
fprintf('    Y : [%.3f, %.3f]\n', min(y), max(y));
fprintf('    Z : [%.3f, %.3f]\n', min(z), max(z));
fprintf('\n── Left BC ──────────────────────────────────────\n');
fprintf('  Nodes    : %d\n', numel(left_node_idx));
fprintf('  Elements : %d\n', numel(left_elem_rows));
fprintf('  Surf elem faces highlighted: %d\n', size(left_surf_faces,1));
fprintf('\n── Right BC ─────────────────────────────────────\n');
fprintf('  Nodes    : %d\n', numel(right_node_idx));
fprintf('  Elements : %d\n', numel(right_elem_rows));
fprintf('  Surf elem faces highlighted: %d\n', size(right_surf_faces,1));
fprintf('─────────────────────────────────────────────────\n');