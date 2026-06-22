clear; clc;

inpFile   = 'Job-2.inp';
prefix    = 'Job-2';

leftSetName  = 'Set-7';
rightSetName = 'Set-8';

fid = fopen(inpFile,'r');
if fid == -1, error('Cannot open: %s', inpFile); end
lines  = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
lines  = lines{1};
nLines = numel(lines);

nodes    = zeros(300000, 4);   nNode = 0;
elements = cell(200000, 1);    nElem = 0;

leftIDs  = [];
rightIDs = [];
CMOD1IDs = [];
CMOD2IDs = [];

isKeyword = @(s) ~isempty(s) && s(1) == '*';
toNums    = @(s) sscanf(strrep(s,',',' '),'%f').';
lowerTrim = @(s) lower(strtrim(s));

i = 1;
while i <= nLines

    line = strtrim(lines{i});

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
                    nodes(nNode,:) = [vals(1:3), 0];
                end
            end
            i = i + 1;
        end
        continue;
    end

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

    if startsWith(line,'*Nset','IgnoreCase',true)
        head  = lowerTrim(line);
        parts = regexp(head,'\s*,\s*','split');
        isGen = any(strcmpi(strtrim(parts),'generate'));

        setName = '';
        for p = 2:numel(parts)
            tok = strtrim(parts{p});
            if startsWith(tok,'nset='), setName = strtrim(tok(6:end)); end
            if startsWith(tok,'set='),  setName = strtrim(tok(5:end)); end
        end
        setName = strrep(strrep(setName,'"',''),'''','');

        isLeft  = strcmpi(setName, leftSetName);
        isRight = strcmpi(setName, rightSetName);
        isCM1   = strcmpi(setName, 'CMOD1');
        isCM2   = strcmpi(setName, 'CMOD2');

        if ~(isLeft || isRight || isCM1 || isCM2)

            i = i + 1;
            while i <= nLines && ~isKeyword(strtrim(lines{i})), i = i+1; end
            continue;
        end

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
                    ids = [ids, a:stp:b];
                else
                    ids = [ids, v];
                end
            end
            i = i + 1;
        end
        ids = unique(ids(~isnan(ids)));

        if isLeft,  leftIDs  = [leftIDs;  ids(:)]; end
        if isRight, rightIDs = [rightIDs; ids(:)]; end
        if isCM1,   CMOD1IDs = [CMOD1IDs; ids(:)]; end
        if isCM2,   CMOD2IDs = [CMOD2IDs; ids(:)]; end
        continue;
    end

    i = i + 1;
end

nodes    = nodes(1:nNode, :);
elements = elements(1:nElem);

[~, firstIdx] = unique(nodes(:,1),'stable');
nodes = nodes(sort(firstIdx), :);

leftIDs  = unique(leftIDs);
rightIDs = unique(rightIDs);
CMOD1IDs = unique(CMOD1IDs);
CMOD2IDs = unique(CMOD2IDs);

if numel(CMOD1IDs) ~= 1
    warning('CMOD1 has %d node(s); expected exactly 1.', numel(CMOD1IDs));
end
if numel(CMOD2IDs) ~= 1
    warning('CMOD2 has %d node(s); expected exactly 1.', numel(CMOD2IDs));
end

fprintf('Nodes:         %d\n',   size(nodes,1));
fprintf('Elements:      %d\n',   nElem);
fprintf('Left  BC (%s): %d node(s)\n',  leftSetName,  numel(leftIDs));
fprintf('Right BC (%s): %d node(s)\n',  rightSetName, numel(rightIDs));
fprintf('CMOD1 (+1):    %d node(s) -> %s\n', numel(CMOD1IDs), mat2str(CMOD1IDs(:).'));
fprintf('CMOD2 (-1):    %d node(s) -> %s\n', numel(CMOD2IDs), mat2str(CMOD2IDs(:).'));

writeIDs([prefix '_nodes.txt'], nodes, true);

fid = fopen([prefix '_elements.txt'],'w');
if fid < 0, error('Cannot write %s_elements.txt', prefix); end
for k = 1:nElem
    v = elements{k};
    fprintf(fid,'%d', v(1));
    fprintf(fid,' %d', v(2:end));
    fprintf(fid,'\n');
end
fclose(fid);

writeIDs([prefix '_left_nodes.txt'],  leftIDs);
writeIDs([prefix '_right_nodes.txt'], rightIDs);
writeIDs([prefix '_CMOD1_nodes.txt'], CMOD1IDs);
writeIDs([prefix '_CMOD2_nodes.txt'], CMOD2IDs);

disp('Export complete.');

function writeIDs(filename, data, isNodeMatrix)

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