%% =========================================================
%  extract.m
%  Extracts nodes, elements, two BC node sets (Left / Right),
%  and two CMOD node sets from an Abaqus .inp file.
%
%  Writes SIX files:
%    <prefix>_nodes.txt          – all nodes  [id  x  y  z]
%    <prefix>_elements.txt       – all elements
%    <prefix>_left_nodes.txt     – Left BC  node set  (Set-7 / XSYMM)
%    <prefix>_right_nodes.txt    – Right BC node set  (Set-8 / ENCASTRE)
%    <prefix>_CMOD1_nodes.txt    – single node with +1 coefficient
%    <prefix>_CMOD2_nodes.txt    – single node with -1 coefficient
%
%  CMOD (Crack Mouth Opening Displacement):
%    delta_X = U1(CMOD1) - U1(CMOD2)
% =========================================================

clear; clc;

%% ── SETTINGS ─────────────────────────────────────────────────────────────
inpFile   = 'Job-2.inp';   % <-- your Abaqus input file
prefix    = 'Job-2';       % <-- prefix for all output files

% Node-set names in the .inp that correspond to Left and Right BCs
leftSetName  = 'Set-7';   % XSYMM (left)
rightSetName = 'Set-8';   % ENCASTRE (right)

%% ── READ FILE ────────────────────────────────────────────────────────────
fid = fopen(inpFile,'r');
if fid == -1, error('Cannot open: %s', inpFile); end
lines  = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines  = lines{1};
nLines = numel(lines);

%% ── CONTAINERS ───────────────────────────────────────────────────────────
nodes    = zeros(300000, 4);   nNode = 0;   % [id  x  y  z]
elements = cell(200000, 1);    nElem = 0;

leftIDs  = [];
rightIDs = [];
CMOD1IDs = [];
CMOD2IDs = [];

%% ── HELPERS ──────────────────────────────────────────────────────────────
isKeyword = @(s) ~isempty(s) && s(1) == '*';
toNums    = @(s) sscanf(strrep(s,',',' '),'%f').';
lowerTrim = @(s) lower(strtrim(s));

%% ── PARSE ────────────────────────────────────────────────────────────────
i = 1;
while i <= nLines

    line = strtrim(lines{i});

    % ---- *Node block ------------------------------------------------
    if startsWith(line,'*Node','IgnoreCase',true) && ...
       ~startsWith(line,'*Nset','IgnoreCase',true)
        i = i + 1;
        while i <= nLines && ~isKeyword(strtrim(lines{i}))
            vals = toNums(lines{i});
            if numel(vals) >= 3
                nNode = nNode + 1;
                if numel(vals) >= 4
                    nodes(nNode,:) = vals(1:4);
                else
                    nodes(nNode,:) = [vals(1:3), 0];   % pad Z if 2-D
                end
            end
            i = i + 1;
        end
        continue;
    end

    % ---- *Element block ---------------------------------------------
    if startsWith(line,'*Element','IgnoreCase',true)
        i = i + 1;
        while i <= nLines && ~isKeyword(strtrim(lines{i}))
            vals = toNums(lines{i});
            if ~isempty(vals)
                nElem = nElem + 1;
                elements{nElem} = vals;
            end
            i = i + 1;
        end
        continue;
    end

    % ---- *Nset block ------------------------------------------------
    if startsWith(line,'*Nset','IgnoreCase',true)
        head  = lowerTrim(line);
        parts = regexp(head,'\s*,\s*','split');
        isGen = any(strcmpi(strtrim(parts),'generate'));

        % Extract set name (handles nset= and set=, with or without quotes)
        setName = '';
        for p = 2:numel(parts)
            tok = strtrim(parts{p});
            if startsWith(tok,'nset='), setName = strtrim(tok(6:end)); end
            if startsWith(tok,'set='),  setName = strtrim(tok(5:end)); end
        end
        setName = strrep(strrep(setName,'"',''),'''','');

        % Identify which of the four sets this is
        isLeft  = strcmpi(setName, leftSetName);
        isRight = strcmpi(setName, rightSetName);
        isCM1   = strcmpi(setName, 'CMOD1');
        isCM2   = strcmpi(setName, 'CMOD2');

        if ~(isLeft || isRight || isCM1 || isCM2)
            % Skip all other sets
            i = i + 1;
            while i <= nLines && ~isKeyword(strtrim(lines{i})), i = i+1; end
            continue;
        end

        % Read node-set body
        ids = [];
        i   = i + 1;
        while i <= nLines && ~isKeyword(strtrim(lines{i}))
            v = toNums(lines{i});
            if ~isempty(v)
                if isGen
                    a   = v(1);  b = v(2);
                    stp = 1;
                    if numel(v) >= 3 && ~isnan(v(3)) && v(3) ~= 0
                        stp = v(3);
                    end
                    ids = [ids, a:stp:b]; %#ok<AGROW>
                else
                    ids = [ids, v]; %#ok<AGROW>
                end
            end
            i = i + 1;
        end
        ids = unique(ids(~isnan(ids)));

        if isLeft,  leftIDs  = [leftIDs;  ids(:)]; end %#ok<AGROW>
        if isRight, rightIDs = [rightIDs; ids(:)]; end %#ok<AGROW>
        if isCM1,   CMOD1IDs = [CMOD1IDs; ids(:)]; end %#ok<AGROW>
        if isCM2,   CMOD2IDs = [CMOD2IDs; ids(:)]; end %#ok<AGROW>
        continue;
    end

    i = i + 1;
end

%% ── TRIM & DEDUPLICATE ───────────────────────────────────────────────────
nodes    = nodes(1:nNode, :);
elements = elements(1:nElem);

[~, firstIdx] = unique(nodes(:,1),'stable');
nodes = nodes(sort(firstIdx), :);

leftIDs  = unique(leftIDs);
rightIDs = unique(rightIDs);
CMOD1IDs = unique(CMOD1IDs);
CMOD2IDs = unique(CMOD2IDs);

%% ── VALIDATION: CMOD must be exactly 1 node each ─────────────────────────
if numel(CMOD1IDs) ~= 1
    warning('CMOD1 has %d node(s); expected exactly 1.', numel(CMOD1IDs));
end
if numel(CMOD2IDs) ~= 1
    warning('CMOD2 has %d node(s); expected exactly 1.', numel(CMOD2IDs));
end

%% ── SUMMARY ──────────────────────────────────────────────────────────────
fprintf('Nodes:         %d\n',   size(nodes,1));
fprintf('Elements:      %d\n',   nElem);
fprintf('Left  BC (%s): %d node(s)\n',  leftSetName,  numel(leftIDs));
fprintf('Right BC (%s): %d node(s)\n',  rightSetName, numel(rightIDs));
fprintf('CMOD1 (+1):    %d node(s) -> %s\n', numel(CMOD1IDs), mat2str(CMOD1IDs(:).'));
fprintf('CMOD2 (-1):    %d node(s) -> %s\n', numel(CMOD2IDs), mat2str(CMOD2IDs(:).'));

%% ── WRITE: nodes (id  x  y  z) ──────────────────────────────────────────
writeIDs([prefix '_nodes.txt'], nodes, true);

%% ── WRITE: elements ──────────────────────────────────────────────────────
fid = fopen([prefix '_elements.txt'],'w');
if fid < 0, error('Cannot write %s_elements.txt', prefix); end
for k = 1:nElem
    v = elements{k};
    fprintf(fid,'%d', v(1));
    fprintf(fid,' %d', v(2:end));
    fprintf(fid,'\n');
end
fclose(fid);

%% ── WRITE: four node-set files ───────────────────────────────────────────
writeIDs([prefix '_left_nodes.txt'],  leftIDs);
writeIDs([prefix '_right_nodes.txt'], rightIDs);
writeIDs([prefix '_CMOD1_nodes.txt'], CMOD1IDs);
writeIDs([prefix '_CMOD2_nodes.txt'], CMOD2IDs);

disp('Export complete.');

%% ── HELPER FUNCTION ──────────────────────────────────────────────────────
function writeIDs(filename, data, isNodeMatrix)
% Write a vector of node IDs, or an Nx4 node matrix when isNodeMatrix==true.
    if nargin < 3, isNodeMatrix = false; end
    fid = fopen(filename,'w');
    if fid < 0, error('Cannot write %s', filename); end
    if isNodeMatrix
        fprintf(fid,'%d %.6f %.6f %.6f\n', data.');
    else
        if ~isempty(data)
            fprintf(fid,'%d\n', data(:));
        end
    end
    fclose(fid);
end