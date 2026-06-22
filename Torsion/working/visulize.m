clearvars; close all; clc;

fprintf('Loading nodes ...\n');
raw_nodes    = readmatrix('Job-1_nodes.txt');
nodeID       = raw_nodes(:,1);
coords       = raw_nodes(:,2:4);

fprintf('Loading elements ...\n');
raw_elems    = readmatrix('Job-1_elements.txt');
elemID       = raw_elems(:,1);
conn_id      = raw_elems(:,2:5);

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

maxID   = max(nodeID);
nodeMap = zeros(maxID, 1, 'uint32');
nodeMap(nodeID) = uint32(1:numel(nodeID));

conn = nodeMap(conn_id);

nE    = size(conn, 1);
face_local = [1 2 3;
              1 2 4;
              1 3 4;
              2 3 4];

all_faces = reshape(conn(:, face_local'), [], 3);

faces_sorted = sort(all_faces, 2);

[~, ia, ic]  = unique(faces_sorted, 'rows', 'stable');
counts        = accumarray(ic, 1);
surf_mask     = counts(ic) == 1;
surf_faces    = all_faces(surf_mask, :);

fprintf('  Surface triangles: %d\n', size(surf_faces,1));

left_node_idx  = nodeMap(left_nodeIDs);
right_node_idx = nodeMap(right_nodeIDs);

maxEID = max(elemID);
elemMap = zeros(maxEID, 1, 'uint32');
elemMap(elemID) = uint32(1:numel(elemID));

left_elem_rows  = elemMap(left_elemIDs);
right_elem_rows = elemMap(right_elemIDs);

left_conn  = conn(left_elem_rows,  :);
right_conn = conn(right_elem_rows, :);

left_faces  = reshape(left_conn(:, face_local'),  [], 3);
right_faces = reshape(right_conn(:, face_local'), [], 3);

left_sorted  = sort(left_faces,  2);
right_sorted = sort(right_faces, 2);

[left_is_surf]  = ismember(left_sorted,  faces_sorted(surf_mask,:), 'rows');
[right_is_surf] = ismember(right_sorted, faces_sorted(surf_mask,:), 'rows');

left_surf_faces  = left_faces(left_is_surf,  :);
right_surf_faces = right_faces(right_is_surf, :);

x = coords(:,1);  y = coords(:,2);  z = coords(:,3);

surf_z = mean(z(surf_faces), 2);

fig = figure('Name','FE Mesh + BCs', 'Color','w', ...
             'Position',[100 100 1200 700]);

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

ax2 = subplot(1,2,2);

patch('Faces',    surf_faces, ...
      'Vertices', coords,     ...
      'FaceColor',[0.78 0.82 0.86], ...
      'EdgeColor','none',     ...
      'FaceAlpha', 0.25);
hold on;

if ~isempty(left_surf_faces)
    patch('Faces',    left_surf_faces, ...
          'Vertices', coords,          ...
          'FaceColor',[0.13 0.47 0.71],...
          'EdgeColor','none',          ...
          'FaceAlpha', 0.90,           ...
          'DisplayName','Left BC faces');
end

if ~isempty(right_surf_faces)
    patch('Faces',    right_surf_faces, ...
          'Vertices', coords,           ...
          'FaceColor',[0.84 0.15 0.16], ...
          'EdgeColor','none',           ...
          'FaceAlpha', 0.90,            ...
          'DisplayName','Right BC faces');
end

scatter3(x(left_node_idx), y(left_node_idx), z(left_node_idx), ...
         28, [0.13 0.47 0.71], 'filled',  ...
         'DisplayName','Left BC nodes', ...
         'MarkerEdgeColor','k', 'LineWidth',0.3);

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

figure('Name','BC Close-up', 'Color','w', 'Position',[150 150 1200 500]);

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