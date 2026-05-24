function plot_nooru_mesh_BC
% PLOT_NOORU_MESH_BC  Publication-quality mesh + BC figures for the
% Nooru-Mohamed double-edge-notched specimen (3D TET4 mesh).
%
% Reads the Job-1 text files, extracts the exterior surface of the
% tetrahedral mesh, colours the four grip/boundary node sets, and renders:
%
%   Fig.1 -> oblique 3D view
%   Fig.2 -> clean 2D front-view BC schematic WITHOUT dimension labels
%
% Outputs:
%   nooru_mesh_3D.pdf
%   nooru_mesh_3D.png
%   nooru_BC_2D_no_dimensions.pdf
%   nooru_BC_2D_no_dimensions.png
%
% Connectivity assumed: [elemID n1 n2 n3 n4]
% Nodes assumed:        [nodeID x y z]
%
% Just run:
%   >> plot_nooru_mesh_BC

% ======================================================================
% CONFIG
% ======================================================================
F.nodes  = 'Job-1_nodes.txt';
F.elems  = 'Job-1_elements.txt';
F.top    = 'Job-1_top_nodes.txt';
F.bottom = 'Job-1_bottom_nodes.txt';
F.left   = 'Job-1_left_nodes.txt';
F.right  = 'Job-1_right_nodes.txt';

C.bulk   = [0.86 0.86 0.89];     % bulk surface fill
C.edge   = [0.42 0.42 0.48];     % mesh edge colour
C.top    = [0.84 0.15 0.16];     % red
C.bottom = [0.12 0.47 0.71];     % blue
C.left   = [0.17 0.63 0.17];     % green
C.right  = [0.58 0.40 0.74];     % purple

view3D     = [-58 16];           % [azimuth elevation]
exportFigs = true;               % save figures
% ======================================================================

%% ---------------------------------------------------------------------
% READ MESH
% ----------------------------------------------------------------------
Nd = readNumeric(F.nodes);        % [id x y z]
El = readNumeric(F.elems);        % [id n1 n2 n3 n4]

nodeID = Nd(:,1);
V      = Nd(:,2:4);
conn   = El(:,2:5);

maxID = max(nodeID);
id2row = zeros(maxID,1);
id2row(nodeID) = 1:numel(nodeID);

T = id2row(conn);
nNode = size(V,1);

%% ---------------------------------------------------------------------
% EXTERIOR SURFACE OF TETRAHEDRAL MESH
% ----------------------------------------------------------------------
Fall = [T(:,[1 2 3]); ...
        T(:,[1 2 4]); ...
        T(:,[1 3 4]); ...
        T(:,[2 3 4])];

[sF,~] = sort(Fall,2);
[uF,~,ic] = unique(sF,'rows');
cnt  = accumarray(ic,1);

surf = uF(cnt==1,:);              % exterior triangular faces

%% ---------------------------------------------------------------------
% GRIP NODE SETS
% ----------------------------------------------------------------------
sets  = {'top','bottom','left','right'};
cols  = {C.top,C.bottom,C.left,C.right};
inSet = false(nNode,numel(sets));

for k = 1:numel(sets)
    ids = readNumeric(F.(sets{k}));
    inSet(id2row(ids(:,1)),k) = true;
end

faceInSet = false(size(surf,1),numel(sets));

for k = 1:numel(sets)
    m = inSet(:,k);
    faceInSet(:,k) = m(surf(:,1)) & m(surf(:,2)) & m(surf(:,3));
end

bulkFace = ~any(faceInSet,2);

%% ---------------------------------------------------------------------
% GEOMETRY LIMITS
% ----------------------------------------------------------------------
xr = [min(V(:,1)) max(V(:,1))];
yr = [min(V(:,2)) max(V(:,2))];
zr = [min(V(:,3)) max(V(:,3))];

Lspec = xr(2) - xr(1);
Hspec = yr(2) - yr(1);

%% =====================================================================
% FIGURE 1: 3D VIEW
% ======================================================================
Vp = V(:,[1 3 2]);                % plot as x, z, y

f1 = figure('Color','w','Units','centimeters','Position',[2 2 17 18]);
ax1 = axes(f1);
hold(ax1,'on');

patch(ax1, ...
    'Faces',surf(bulkFace,:), ...
    'Vertices',Vp, ...
    'FaceColor',C.bulk, ...
    'EdgeColor',C.edge, ...
    'LineWidth',0.15, ...
    'FaceAlpha',1);

for k = 1:numel(sets)
    patch(ax1, ...
        'Faces',surf(faceInSet(:,k),:), ...
        'Vertices',Vp, ...
        'FaceColor',cols{k}, ...
        'EdgeColor',[0.15 0.15 0.15], ...
        'LineWidth',0.15, ...
        'FaceAlpha',0.97);
end

axis(ax1,'equal');
axis(ax1,'tight');
view(ax1,view3D);

camlight(ax1,'headlight');
lighting(ax1,'gouraud');
material(ax1,'dull');

grid(ax1,'off');
box(ax1,'off');

xlabel(ax1,'x (mm)');
ylabel(ax1,'z (mm)');
zlabel(ax1,'y (mm)');

set(ax1,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

% 3D load arrows
zf = zr(2);
Larrow3D = 0.16*Hspec;

% Top tension arrows
xa = linspace(xr(1),xr(2),6);

quiver3(ax1, ...
    xa, ...
    zf*ones(size(xa)), ...
    (yr(2)+0.04*Larrow3D)*ones(size(xa)), ...
    zeros(size(xa)), ...
    zeros(size(xa)), ...
    Larrow3D*ones(size(xa)), ...
    0, ...
    'Color',C.top, ...
    'LineWidth',1.4, ...
    'MaxHeadSize',0.6);

text(ax1, ...
    xr(2), zf, yr(2)+Larrow3D*1.5, ...
    '  P  (tension)', ...
    'Color',C.top, ...
    'FontSize',11, ...
    'FontWeight','bold');

% Left shear arrows
ya = linspace(0,yr(2),6);

quiver3(ax1, ...
    (xr(1)-1.05*Larrow3D)*ones(size(ya)), ...
    zf*ones(size(ya)), ...
    ya, ...
    Larrow3D*ones(size(ya)), ...
    zeros(size(ya)), ...
    zeros(size(ya)), ...
    0, ...
    'Color',C.left, ...
    'LineWidth',1.4, ...
    'MaxHeadSize',0.6);

text(ax1, ...
    xr(1)-1.6*Larrow3D, zf, yr(2), ...
    'P_s ', ...
    'HorizontalAlignment','right', ...
    'Color',C.left, ...
    'FontSize',11, ...
    'FontWeight','bold');

% Legend
hL = gobjects(1,4);

for k = 1:numel(sets)
    hL(k) = patch(ax1,NaN,NaN,cols{k}, ...
        'EdgeColor','k', ...
        'LineWidth',0.3);
end

legend(ax1,hL, ...
    {'Top grip','Bottom grip','Left grip','Right grip'}, ...
    'Location','northeastoutside', ...
    'FontSize',10, ...
    'Box','off');

title(ax1, ...
    'Nooru-Mohamed specimen: mesh and boundary grips', ...
    'FontWeight','bold', ...
    'FontSize',12);

%% =====================================================================
% FIGURE 2: 2D FRONT VIEW WITHOUT DIMENSIONS
% ======================================================================
z = V(:,3);
tol = 1e-6;
isZ0 = abs(z) < tol;

front = surf(isZ0(surf(:,1)) & ...
             isZ0(surf(:,2)) & ...
             isZ0(surf(:,3)), :);

f2 = figure('Color','w','Units','centimeters','Position',[2 2 16 14]);
ax2 = axes(f2);
hold(ax2,'on');

patch(ax2, ...
    'Faces',front, ...
    'Vertices',V(:,1:2), ...
    'FaceColor',[0.95 0.95 0.97], ...
    'EdgeColor',[0.55 0.55 0.6], ...
    'LineWidth',0.18);

% Highlight grip edges
drawGripEdge(ax2, V, front, inSet(:,1), C.top,    2.6); % top
drawGripEdge(ax2, V, front, inSet(:,2), C.bottom, 2.6); % bottom
drawGripEdge(ax2, V, front, inSet(:,3), C.left,   2.6); % left
drawGripEdge(ax2, V, front, inSet(:,4), C.right,  2.6); % right

% Fixed supports
hatchEdge(ax2, [xr(1) yr(1); xr(2) yr(1)], 'down',  C.bottom);
hatchEdge(ax2, [xr(2) yr(1); xr(2) yr(1)+0.5*Hspec], 'right', C.right);

% Load arrows
La = 0.12*Hspec;

% Top tension arrows
xa = linspace(xr(1),xr(2),7);

quiver(ax2, ...
    xa, ...
    (yr(2)+0.04*Hspec)*ones(size(xa)), ...
    zeros(size(xa)), ...
    0.18*Hspec*ones(size(xa)), ...
    0, ...
    'Color',C.top, ...
    'LineWidth',1.3, ...
    'MaxHeadSize',0.5);

text(ax2, ...
    mean(xr), yr(2)+0.32*Hspec, ...
    'P  (tension)', ...
    'Color',C.top, ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',11);

% Left shear arrows
ya = linspace(yr(1)+0.2*Hspec,yr(2),6);

quiver(ax2, ...
    (xr(1)-0.25*Lspec)*ones(size(ya)), ...
    ya, ...
    0.22*Lspec*ones(size(ya)), ...
    zeros(size(ya)), ...
    0, ...
    'Color',C.left, ...
    'LineWidth',1.3, ...
    'MaxHeadSize',0.5);

text(ax2, ...
    xr(1)-0.36*Lspec, yr(1)+0.65*Hspec, ...
    'P_s', ...
    'Color',C.left, ...
    'Rotation',90, ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',11);

% Notch labels only
text(ax2, ...
    xr(1)+0.22*Lspec, yr(1)+0.11*Hspec, ...
    '\leftarrow notch', ...
    'Color',[0.2 0.2 0.2], ...
    'FontSize',9, ...
    'BackgroundColor','w', ...
    'Margin',1);

text(ax2, ...
    xr(2)-0.08*Lspec, yr(1)-0.14*Hspec, ...
    'notch \rightarrow', ...
    'Color',[0.2 0.2 0.2], ...
    'FontSize',9, ...
    'HorizontalAlignment','right', ...
    'BackgroundColor','w', ...
    'Margin',1);

axis(ax2,'equal');

% Clean axis limits without dimension space
xlim(ax2,[xr(1)-0.35*Lspec, xr(2)+0.25*Lspec]);
ylim(ax2,[yr(1)-0.20*Hspec, yr(2)+0.40*Hspec]);

xlabel(ax2,'x (mm)');
ylabel(ax2,'y (mm)');

set(ax2,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

box(ax2,'on');
grid(ax2,'off');

title(ax2, ...
    'Boundary conditions', ...
    'FontWeight','bold', ...
    'FontSize',12);

%% ---------------------------------------------------------------------
% EXPORT
% ----------------------------------------------------------------------
if exportFigs
    exportgraphics(f1,'nooru_mesh_3D.pdf','ContentType','vector');
    exportgraphics(f1,'nooru_mesh_3D.png','Resolution',300);

    exportgraphics(f2,'nooru_BC_2D_no_dimensions.pdf','ContentType','vector');
    exportgraphics(f2,'nooru_BC_2D_no_dimensions.png','Resolution',300);

    fprintf('Saved:\n');
    fprintf('  nooru_mesh_3D.pdf\n');
    fprintf('  nooru_mesh_3D.png\n');
    fprintf('  nooru_BC_2D_no_dimensions.pdf\n');
    fprintf('  nooru_BC_2D_no_dimensions.png\n');
end

fprintf('Nodes: %d | Elements: %d | Surface tris: %d\n', ...
    nNode, size(T,1), size(surf,1));

fprintf('Geometry size from mesh:\n');
fprintf('  x length = %.3f mm\n',Lspec);
fprintf('  y height = %.3f mm\n',Hspec);
fprintf('  z thickness = %.3f mm\n',zr(2)-zr(1));

end

% ======================================================================
% HELPER FUNCTIONS
% ======================================================================

function A = readNumeric(fname)
% Robust whitespace-delimited numeric read.

if exist(fname,'file') ~= 2
    error('File not found: %s. Run from the folder containing the txt files.',fname);
end

try
    A = readmatrix(fname);
catch
    A = load(fname);
end

A = A(all(~isnan(A),2),:);
end

% ----------------------------------------------------------------------
function drawGripEdge(ax, V, faces, nodeMask, col, lw)
% Draw boundary edges of the front mesh whose endpoints are in a node set.

E = [faces(:,[1 2]); ...
     faces(:,[2 3]); ...
     faces(:,[3 1])];

Es = sort(E,2);
[uE,~,ic] = unique(Es,'rows');
bnd = uE(accumarray(ic,1)==1,:);

sel = nodeMask(bnd(:,1)) & nodeMask(bnd(:,2));
bnd = bnd(sel,:);

X = [V(bnd(:,1),1), V(bnd(:,2),1)]';
Y = [V(bnd(:,1),2), V(bnd(:,2),2)]';

line(ax, X, Y, ...
    'Color',col, ...
    'LineWidth',lw);
end

% ----------------------------------------------------------------------
function hatchEdge(ax, seg, dir, col)
% Simple fixed-support hatching along a segment.

p1 = seg(1,:);
p2 = seg(2,:);

nTicks = 14;
t = linspace(0,1,nTicks)';

P = p1 + t.*(p2-p1);

axisScale = max(range(xlim(ax)),range(ylim(ax)));
L = 0.018*axisScale;

switch dir
    case 'down'
        d = [-1 -1];
    case 'right'
        d = [ 1 -1];
    otherwise
        d = [ 1  1];
end

d = d./norm(d)*L*2.2;

for i = 1:nTicks
    line(ax, ...
        [P(i,1), P(i,1)+d(1)], ...
        [P(i,2), P(i,2)+d(2)], ...
        'Color',col, ...
        'LineWidth',0.8);
end

line(ax, ...
    [p1(1), p2(1)], ...
    [p1(2), p2(2)], ...
    'Color',col, ...
    'LineWidth',1.6);
end