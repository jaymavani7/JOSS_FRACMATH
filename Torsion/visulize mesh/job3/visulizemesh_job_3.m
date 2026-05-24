%% =========================================================
%  visualize_mesh_Job3_WITH_DIMENSIONS_AND_CORRECT_LOAD.m
%  Visualises FE mesh + BC nodes + dimensions + load at node 67
% =========================================================

clear; clc; close all;

%% ── 0. FILE PATHS ───────────────────────────────────────────────────────
fprintf('Select the NODES file (Job-3_nodes.txt) ...\n');
[fn, fp] = uigetfile('*.txt','Select NODES file (Job-3_nodes.txt)');
if isequal(fn,0), error('No node file selected. Aborting.'); end
nodesFile = fullfile(fp, fn);

fprintf('Select the ELEMENTS file (Job-3_elements.txt) ...\n');
[fn, fp] = uigetfile('*.txt','Select ELEMENTS file (Job-3_elements.txt)');
if isequal(fn,0), error('No element file selected. Aborting.'); end
elemFile = fullfile(fp, fn);

%% ── 1. READ NODES ───────────────────────────────────────────────────────
fprintf('Reading nodes ... ');
nodeData = readmatrix(nodesFile, 'FileType','text');

nodeID = nodeData(:,1);
X      = nodeData(:,2);
Y      = nodeData(:,3);
Z      = nodeData(:,4);

fprintf('%d nodes loaded.\n', numel(nodeID));

maxID  = max(nodeID);
id2row = zeros(maxID,1,'int32');
id2row(nodeID) = 1:numel(nodeID);

%% ── 2. READ ELEMENTS ────────────────────────────────────────────────────
fprintf('Reading elements ... ');
elemData = readmatrix(elemFile, 'FileType','text');

% C3D10: first 4 are corner nodes
E     = elemData(:,2:5);
nElem = size(E,1);

fprintf('%d elements loaded.\n', nElem);

%% ── 3. EXTRACT OUTER SURFACE FACES ──────────────────────────────────────
fprintf('Finding outer surface faces ... ');

faceConn = [1 2 3;
            1 2 4;
            1 3 4;
            2 3 4];

nFacesTotal = nElem * 4;
allFaces    = zeros(nFacesTotal, 3, 'int32');

for f = 1:4
    rows = (f-1)*nElem + 1 : f*nElem;
    allFaces(rows,:) = int32(E(:, faceConn(f,:)));
end

allFaces = sort(allFaces, 2);

[sortedF, ia] = sortrows(allFaces);
isDup = all(diff(sortedF,1,1) == 0, 2);

dupMask = false(nFacesTotal,1);
dupMask([find(isDup); find(isDup)+1]) = true;

outerFaces = allFaces(ia(~dupMask), :);

fprintf('%d outer triangular faces.\n', size(outerFaces,1));

%% ── 4. DEFINE BC NODES AND LOAD NODE ────────────────────────────────────
% According to your correction:
%   node 67 is the LOAD POINT
% so do not treat node 67 as just a generic plotted load at center.
%
% Keep only support / BC nodes here:
BC.names  = {'left\_left','left\_right','right\_right'};
BC.nodeID = [14, 45, 64];
BC.colors = [0.85 0.12 0.08;   % red
             0.06 0.53 0.85;   % blue
             0.10 0.68 0.10];  % green

% Load node
LOAD.nodeID = 67;
LOAD.name   = 'Load point';
LOAD.color  = [0.00 0.45 0.25];   % dark green
LOAD.dir    = [0 -1 0];           % downward in -Y direction

%% ── 5. GEOMETRY DIMENSIONS ──────────────────────────────────────────────
xmin = min(X); xmax = max(X);
ymin = min(Y); ymax = max(Y);
zmin = min(Z); zmax = max(Z);

Lx = xmax - xmin;
Ly = ymax - ymin;
Lz = zmax - zmin;

padX = 0.12 * max(Lx, eps);
padY = 0.12 * max(Ly, eps);
padZ = 0.12 * max(Lz, eps);

%% ── 6. CREATE FIGURE ────────────────────────────────────────────────────
fig = figure('Name','Job-3 | FE Mesh + Boundary Conditions + Dimensions + Load', ...
             'Color','w', ...
             'NumberTitle','off', ...
             'Units','normalized', ...
             'Position',[0.05 0.05 0.88 0.88]);

ax = axes('Parent', fig);
hold(ax, 'on');

%% ── 7. PLOT MESH SURFACE ────────────────────────────────────────────────
patch(ax, ...
    'Faces',    double(outerFaces), ...
    'Vertices', [X, Y, Z], ...
    'FaceColor',[0.72 0.72 0.76], ...
    'EdgeColor',[0.30 0.30 0.35], ...
    'FaceAlpha', 0.18, ...
    'EdgeAlpha', 0.12, ...
    'LineWidth', 0.35, ...
    'DisplayName','Mesh surface');

%% ── 8. PLOT BC NODES ────────────────────────────────────────────────────
for k = 1:numel(BC.nodeID)

    nid = BC.nodeID(k);

    if nid > maxID || id2row(nid) == 0
        warning('BC node %d not found in node list.', nid);
        continue;
    end

    r   = id2row(nid);
    lbl = strrep(BC.names{k}, '\_', '_');

    scatter3(ax, X(r), Y(r), Z(r), 220, ...
             BC.colors(k,:), ...
             'filled', ...
             'MarkerEdgeColor','k', ...
             'LineWidth',1.4, ...
             'DisplayName', sprintf('%s  node %d', lbl, nid));

    text(ax, X(r), Y(r), Z(r), ...
         sprintf('  %s\n  node %d', lbl, nid), ...
         'FontSize', 10, ...
         'FontWeight','bold', ...
         'Color', BC.colors(k,:));
end

%% ── 9. PLOT LOAD NODE 67 ────────────────────────────────────────────────
if LOAD.nodeID > maxID || id2row(LOAD.nodeID) == 0
    error('Load node %d not found in node list.', LOAD.nodeID);
end

rL = id2row(LOAD.nodeID);
xL = X(rL);
yL = Y(rL);
zL = Z(rL);

% Show the load node as a special marker
scatter3(ax, xL, yL, zL, 280, ...
         LOAD.color, ...
         'filled', ...
         'd', ...
         'MarkerEdgeColor','k', ...
         'LineWidth',1.4, ...
         'DisplayName', sprintf('Load node %d', LOAD.nodeID));

text(ax, xL, yL, zL, ...
     sprintf('  load node %d', LOAD.nodeID), ...
     'FontSize', 10, ...
     'FontWeight','bold', ...
     'Color', LOAD.color);

%% ── 10. DRAW LOAD ARROW EXACTLY AT NODE 67 ──────────────────────────────
loadDir = LOAD.dir(:).';
loadDir = loadDir / norm(loadDir);

% Choose arrow length based on specimen size
arrowLength = 0.22 * max([Lx Ly Lz]);

% Start point is above the load node, end point is at the node
startPt = [xL, yL, zL] - loadDir * arrowLength;

quiver3(ax, startPt(1), startPt(2), startPt(3), ...
        loadDir(1)*arrowLength, ...
        loadDir(2)*arrowLength, ...
        loadDir(3)*arrowLength, ...
        0, ...
        'Color', LOAD.color, ...
        'LineWidth', 2.6, ...
        'MaxHeadSize', 0.8, ...
        'DisplayName', 'Applied load P');

text(ax, startPt(1), startPt(2), startPt(3), ...
     '  P', ...
     'Color', LOAD.color, ...
     'FontSize', 18, ...
     'FontWeight', 'bold');

%% ── 11. ADD DIMENSION LINES ─────────────────────────────────────────────
dimColor = [0.05 0.05 0.05];

% X dimension
p1x = [xmin, ymin - padY, zmin - padZ];
p2x = [xmax, ymin - padY, zmin - padZ];
drawDimLine3D(ax, p1x, p2x, sprintf('L_x = %.3g mm', Lx), dimColor);

% Y dimension
p1y = [xmax + padX, ymin, zmin - padZ];
p2y = [xmax + padX, ymax, zmin - padZ];
drawDimLine3D(ax, p1y, p2y, sprintf('H_y = %.3g mm', Ly), dimColor);

% Z dimension
p1z = [xmax + padX, ymax + padY, zmin];
p2z = [xmax + padX, ymax + padY, zmax];
drawDimLine3D(ax, p1z, p2z, sprintf('B_z = %.3g mm', Lz), dimColor);

% Extension lines
plot3(ax, [xmin xmin], [ymin ymin-padY], [zmin zmin-padZ], '--', ...
      'Color', dimColor, 'LineWidth', 0.8, 'HandleVisibility','off');
plot3(ax, [xmax xmax], [ymin ymin-padY], [zmin zmin-padZ], '--', ...
      'Color', dimColor, 'LineWidth', 0.8, 'HandleVisibility','off');

plot3(ax, [xmax xmax+padX], [ymin ymin], [zmin zmin-padZ], '--', ...
      'Color', dimColor, 'LineWidth', 0.8, 'HandleVisibility','off');
plot3(ax, [xmax xmax+padX], [ymax ymax], [zmin zmin-padZ], '--', ...
      'Color', dimColor, 'LineWidth', 0.8, 'HandleVisibility','off');

%% ── 12. AXIS / VIEW SETTINGS ────────────────────────────────────────────
axis(ax,'equal');
grid(ax,'on');
box(ax,'on');

xlabel(ax, 'X [mm]', 'FontSize', 12);
ylabel(ax, 'Y [mm]', 'FontSize', 12);
zlabel(ax, 'Z [mm]', 'FontSize', 12);

title(ax, ...
    { 'Job-3 - FE Mesh + Boundary Conditions + Dimensions', ...
      sprintf('%d nodes | %d C3D10 elements | Lx = %.3g mm | Ly = %.3g mm | Lz = %.3g mm', ...
      numel(nodeID), nElem, Lx, Ly, Lz) }, ...
    'FontSize', 13, ...
    'FontWeight','bold');

legend(ax, 'Location','best', 'FontSize',10);

% Better limits so dimensions and arrow are visible
xlim(ax, [xmin - 0.20*Lx, xmax + 0.30*Lx]);
ylim(ax, [ymin - 0.28*Ly, ymax + 0.35*Ly]);
zlim(ax, [zmin - 0.30*max(Lz,1), zmax + 0.30*max(Lz,1)]);

view(ax, 35, 20);

%% ── 13. VIEW BUTTONS ────────────────────────────────────────────────────
uicontrol('Style','pushbutton','String','Isometric', ...
    'Units','normalized','Position',[0.01 0.20 0.07 0.04], ...
    'Callback',@(~,~) view(ax,35,20));

uicontrol('Style','pushbutton','String','Top XY', ...
    'Units','normalized','Position',[0.01 0.15 0.07 0.04], ...
    'Callback',@(~,~) view(ax,0,90));

uicontrol('Style','pushbutton','String','Front XZ', ...
    'Units','normalized','Position',[0.01 0.10 0.07 0.04], ...
    'Callback',@(~,~) view(ax,0,0));

uicontrol('Style','pushbutton','String','Side YZ', ...
    'Units','normalized','Position',[0.01 0.05 0.07 0.04], ...
    'Callback',@(~,~) view(ax,-90,0));

%% ── 14. SAVE FIGURE ─────────────────────────────────────────────────────
outPNG = fullfile(pwd, 'Job3_mesh_BC_dimensions_correctLoadNode67.png');
outFIG = fullfile(pwd, 'Job3_mesh_BC_dimensions_correctLoadNode67.fig');

exportgraphics(fig, outPNG, 'Resolution', 300);
savefig(fig, outFIG);

fprintf('\nSaved:\n');
fprintf('  PNG: %s\n', outPNG);
fprintf('  FIG: %s\n', outFIG);
fprintf('\nDone.\n');

%% ========================================================================
% LOCAL FUNCTION
%% ========================================================================
function drawDimLine3D(ax, p1, p2, labelText, colorVal)

    p1 = p1(:).';
    p2 = p2(:).';

    v = p2 - p1;
    L = norm(v);

    if L < eps
        return;
    end

    e = v / L;

    plot3(ax, [p1(1) p2(1)], ...
              [p1(2) p2(2)], ...
              [p1(3) p2(3)], ...
              '-', ...
              'Color', colorVal, ...
              'LineWidth', 1.8, ...
              'HandleVisibility','off');

    ah = 0.08 * L;

    quiver3(ax, p1(1), p1(2), p1(3), ...
            ah*e(1), ah*e(2), ah*e(3), ...
            0, ...
            'Color', colorVal, ...
            'LineWidth', 1.8, ...
            'MaxHeadSize', 1.2, ...
            'HandleVisibility','off');

    quiver3(ax, p2(1), p2(2), p2(3), ...
            -ah*e(1), -ah*e(2), -ah*e(3), ...
            0, ...
            'Color', colorVal, ...
            'LineWidth', 1.8, ...
            'MaxHeadSize', 1.2, ...
            'HandleVisibility','off');

    pm = 0.5 * (p1 + p2);

    text(ax, pm(1), pm(2), pm(3), ...
         ['  ' labelText '  '], ...
         'FontSize', 11, ...
         'FontWeight','bold', ...
         'Color', colorVal, ...
         'BackgroundColor','w', ...
         'Margin', 2, ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','middle');
end