clear; clc; close all;

dataDir = uigetdir(pwd, 'Select folder containing Job-2 output files');
if isequal(dataDir, 0), error('No folder selected. Aborting.'); end

prefix = 'Job-2';

nodesFile  = fullfile(dataDir, [prefix '_nodes.txt']);
elemFile   = fullfile(dataDir, [prefix '_elements.txt']);
leftFile   = fullfile(dataDir, [prefix '_left_nodes.txt']);
rightFile  = fullfile(dataDir, [prefix '_right_nodes.txt']);
cmod1File  = fullfile(dataDir, [prefix '_CMOD1_nodes.txt']);
cmod2File  = fullfile(dataDir, [prefix '_CMOD2_nodes.txt']);

fprintf('Reading nodes ... ');
nodeData = readmatrix(nodesFile, 'FileType','text');
nodeID   = nodeData(:,1);
X        = nodeData(:,2);
Y        = nodeData(:,3);
Z        = nodeData(:,4);
nNodes   = numel(nodeID);
fprintf('%d nodes.\n', nNodes);

maxID          = max(nodeID);
id2row         = zeros(maxID, 1, 'int32');
id2row(nodeID) = 1:nNodes;

fprintf('Reading elements ... ');
elemData = readmatrix(elemFile, 'FileType','text');
E        = elemData(:, 2:5);
nElem    = size(E, 1);
fprintf('%d elements.\n', nElem);

fprintf('Extracting outer surface ... ');

faceConn = int32([1 2 3; 1 2 4; 1 3 4; 2 3 4]);
allFaces = zeros(nElem*4, 3, 'int32');
for f = 1:4
    rows = (f-1)*nElem+1 : f*nElem;
    allFaces(rows,:) = int32(E(:, faceConn(f,:)));
end
allFaces = sort(allFaces, 2);

[sortedF, ia] = sortrows(allFaces);
isDup  = all(diff(sortedF,1,1)==0, 2);
dupM   = false(nElem*4, 1);
dupM([find(isDup); find(isDup)+1]) = true;

outerFaces = allFaces(ia(~dupM), :);
fprintf('%d outer faces.\n', size(outerFaces,1));

readSet = @(f) readmatrix(f,'FileType','text');

leftIDs  = readSet(leftFile);
rightIDs = readSet(rightFile);
CMOD1_id = readSet(cmod1File);
CMOD2_id = readSet(cmod2File);

cLeft   = [0.957  0.263  0.212];
cRight  = [0.129  0.588  0.953];
cCMOD1  = [0.612  0.153  0.690];
cCMOD2  = [0.000  0.737  0.831];
cMesh   = [0.780  0.820  0.860];
cEdge   = [0.340  0.380  0.430];

fig = figure('Name','Job-2  |  FE Mesh + BCs + CMOD', ...
             'Color',[0.12 0.13 0.15], ...
             'NumberTitle','off', ...
             'Units','normalized', ...
             'Position',[0.02 0.02 0.94 0.92]);

ax = axes('Parent',fig, ...
          'Color',[0.15 0.16 0.19], ...
          'GridColor',[0.40 0.43 0.50], ...
          'GridAlpha',0.35, ...
          'XColor',[0.80 0.82 0.85], ...
          'YColor',[0.80 0.82 0.85], ...
          'ZColor',[0.80 0.82 0.85], ...
          'FontSize',10, ...
          'FontName','Helvetica Neue');

faceYmean = mean(Y(outerFaces), 2);
faceYnorm = (faceYmean - min(faceYmean)) ./ (max(faceYmean) - min(faceYmean) + eps);

cLow  = [0.30 0.36 0.44];
cHigh = [0.72 0.76 0.82];
faceCols = cLow + faceYnorm .* (cHigh - cLow);

patch(ax, ...
    'Faces',     double(outerFaces), ...
    'Vertices',  [X, Y, Z], ...
    'FaceVertexCData', faceCols, ...
    'FaceColor', 'flat', ...
    'EdgeColor', cEdge, ...
    'FaceAlpha',  0.28, ...
    'EdgeAlpha',  0.18, ...
    'LineWidth',   0.4, ...
    'DisplayName','Mesh surface');

hold(ax,'on');

plotBCnodes(ax, leftIDs,  id2row, X, Y, Z, cLeft,  'Left BC  (XSYMM)',     'o', 220);

plotBCnodes(ax, rightIDs, id2row, X, Y, Z, cRight, 'Right BC (ENCASTRE)',  's', 220);

bbox = [max(X)-min(X), max(Y)-min(Y), max(Z)-min(Z)];
off  = max(bbox) * 0.08;

cmod_info = { CMOD1_id, cCMOD1, 'CMOD1  (+1)', [-off +off +off] ; ...
              CMOD2_id, cCMOD2, 'CMOD2  (−1)', [+off -off -off] };

for c = 1:2
    ids  = cmod_info{c,1};
    col  = cmod_info{c,2};
    name = cmod_info{c,3};
    dxyz = cmod_info{c,4};

    for j = 1:numel(ids)
        r = id2row(ids(j));
        if r == 0, continue; end
        px = X(r);  py = Y(r);  pz = Z(r);

        scatter3(ax, px,py,pz, 680, col, 'd', ...
                 'MarkerFaceAlpha',0.15,'MarkerEdgeAlpha',0.0, ...
                 'HandleVisibility','off');

        sc = scatter3(ax, px,py,pz, 350, col, 'd','filled', ...
                 'MarkerEdgeColor','w','LineWidth',1.8);
        if j == 1
            set(sc,'DisplayName', sprintf('%s   node %d', name, ids(j)));
        else
            set(sc,'HandleVisibility','off');
        end

        lx = px+dxyz(1);  ly = py+dxyz(2);  lz = pz+dxyz(3);
        plot3(ax,[px lx],[py ly],[pz lz],'-','Color',[col 0.7], ...
              'LineWidth',1.2,'HandleVisibility','off');

        text(ax, lx, ly, lz, ...
             sprintf(' %s\n node %d ', name, ids(j)), ...
             'FontSize',9,'FontWeight','bold','Color',col, ...
             'FontName','Helvetica Neue', ...
             'BackgroundColor',[0.10 0.11 0.13], ...
             'EdgeColor',col,'LineWidth',1.5,'Margin',5, ...
             'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

if ~isempty(CMOD1_id) && ~isempty(CMOD2_id)
    r1 = id2row(CMOD1_id(1));
    r2 = id2row(CMOD2_id(1));
    if r1>0 && r2>0

        nSeg = 60;
        tt   = linspace(0,1,nSeg).';
        lx   = X(r1) + tt*(X(r2)-X(r1));
        ly   = Y(r1) + tt*(Y(r2)-Y(r1));
        lz   = Z(r1) + tt*(Z(r2)-Z(r1));
        colM = cCMOD1 + tt .* (cCMOD2-cCMOD1);
        for s = 1:nSeg-1
            plot3(ax,[lx(s) lx(s+1)],[ly(s) ly(s+1)],[lz(s) lz(s+1)], ...
                  '--','Color',[colM(s,:) 0.9],'LineWidth',2.2, ...
                  'HandleVisibility','off');
        end

        plot3(ax,NaN,NaN,NaN,'--w','LineWidth',2.2, ...
              'DisplayName','CMOD gauge line');

        mx   = (X(r1)+X(r2))/2;
        my   = (Y(r1)+Y(r2))/2;
        mz   = (Z(r1)+Z(r2))/2;
        dist = norm([X(r2)-X(r1), Y(r2)-Y(r1), Z(r2)-Z(r1)]);
        text(ax, mx, my, mz+off*0.7, ...
             sprintf(' CMOD gauge \n d₀ = %.3f mm ', dist), ...
             'FontSize',9,'FontWeight','bold', ...
             'FontName','Helvetica Neue', ...
             'Color',[1 1 1], ...
             'BackgroundColor',[0.20 0.18 0.05], ...
             'EdgeColor',[0.95 0.85 0.10],'LineWidth',1.5,'Margin',5, ...
             'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

axis(ax,'equal');  grid(ax,'on');
xlabel(ax,'X  (mm)','FontSize',11,'FontWeight','bold','Color',[0.85 0.87 0.90]);
ylabel(ax,'Y  (mm)','FontSize',11,'FontWeight','bold','Color',[0.85 0.87 0.90]);
zlabel(ax,'Z  (mm)','FontSize',11,'FontWeight','bold','Color',[0.85 0.87 0.90]);

title(ax, ...
    { '\fontsize{13}\bf Job-2  –  FE Mesh  |  Boundary Conditions  |  CMOD', ...
      sprintf('\\fontsize{10}%d nodes    %d C3D10 elements    Left XSYMM  +  Right ENCASTRE  +  CMOD1/CMOD2', ...
              nNodes, nElem) }, ...
    'Color',[0.95 0.96 0.98],'Interpreter','tex');

leg = legend(ax,'Location','bestoutside','FontSize',9, ...
             'NumColumns',1,'TextColor',[0.90 0.92 0.95], ...
             'Color',[0.13 0.14 0.17],'EdgeColor',[0.35 0.38 0.45]);
title(leg,'\bf Legend','Color',[1 1 1]);

view(ax,35,25);
ax.Clipping = 'off';

btnDefs = { 'Isometric',  [35, 25] ; ...
            'Top  (XY)',  [ 0, 90] ; ...
            'Front (XZ)', [ 0,  0] ; ...
            'Side  (YZ)', [-90, 0] };
bW = 0.065;  bH = 0.038;  bX = 0.004;
for k = 1:size(btnDefs,1)
    bY = 0.28 - (k-1)*0.05;
    uicontrol('Style','pushbutton','String',btnDefs{k,1}, ...
        'Units','normalized','Position',[bX bY bW bH], ...
        'BackgroundColor',[0.22 0.24 0.30], ...
        'ForegroundColor',[0.92 0.94 0.97], ...
        'FontSize',8,'FontWeight','bold', ...
        'Callback',@(~,~) view(ax, btnDefs{k,2}(1), btnDefs{k,2}(2)));
end

uicontrol('Style','pushbutton','String','Export PNG', ...
    'Units','normalized','Position',[bX, 0.06, bW, bH], ...
    'BackgroundColor',[0.10 0.45 0.20], ...
    'ForegroundColor',[1 1 1], ...
    'FontSize',8,'FontWeight','bold', ...
    'Callback',@(~,~) exportgraphics(fig, ...
        fullfile(dataDir,[prefix '_mesh_view.png']), ...
        'Resolution',300));

fprintf('\nVisualization ready. Rotate with mouse or use the view buttons.\n');
fprintf('Left  BC nodes : %d\n', numel(leftIDs));
fprintf('Right BC nodes : %d\n', numel(rightIDs));
fprintf('CMOD1 node     : %s\n', mat2str(CMOD1_id(:).'));
fprintf('CMOD2 node     : %s\n', mat2str(CMOD2_id(:).'));

function plotBCnodes(ax, ids, id2row, X, Y, Z, col, label, marker, msz)

    firstDone = false;
    for j = 1:numel(ids)
        r = id2row(ids(j));
        if r == 0, continue; end
        px = X(r);  py = Y(r);  pz = Z(r);

        scatter3(ax, px,py,pz, msz*3.5, col, marker, ...
                 'MarkerFaceAlpha',0.10,'MarkerEdgeAlpha',0.0, ...
                 'HandleVisibility','off');

        sc = scatter3(ax, px,py,pz, msz, col, marker, 'filled', ...
                 'MarkerEdgeColor','w','LineWidth',1.0);

        if ~firstDone
            set(sc,'DisplayName', sprintf('%s', label));
            firstDone = true;
        else
            set(sc,'HandleVisibility','off');
        end

        if numel(ids) <= 50
            text(ax, px, py, pz, sprintf('  %d', ids(j)), ...
                 'FontSize',6,'Color',col, ...
                 'FontName','Helvetica Neue', ...
                 'HorizontalAlignment','left');
        end
    end
end