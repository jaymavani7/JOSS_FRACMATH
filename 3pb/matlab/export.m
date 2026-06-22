function export_mesh_data(inpFile, outFolder)

    clc;

    if nargin < 1 || isempty(inpFile)
        hits = dir('Gregoire_3PB/*.inp');
        if isempty(hits)
            hits = dir('*.inp');
        end
        if isempty(hits)
            error('No .inp file found. Run your Python script first or provide path.');
        end
        inpFile = fullfile(hits(1).folder, hits(1).name);
    end

    if nargin < 2 || isempty(outFolder)
        outFolder = 'Gregoire_3PB';
    end

    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    fprintf('Processing: %s\n', inpFile);

    try
        perform_extraction(inpFile, outFolder);
        fprintf('Success! Files written to %s/\n', outFolder);
    catch ME
        fprintf('Error during extraction: %s\n', ME.message);
    end

end

function perform_extraction(inpFile, outFolder)
    fid = fopen(inpFile, 'r');
    if fid == -1, error('File not found.'); end
    C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = C{1};

    nodes = [];
    elems = [];
    sets = struct('Load_Nodes',[],'Support_Left',[],'Support_Right',[],'CMOD1',[],'CMOD2',[]);
    bc_data = {};

    i = 1;
    while i <= numel(lines)
        line = strtrim(lines{i});
        low = lower(line);

        if startsWith(low, '*node') && ~contains(low, 'output')
            i = i + 1;
            while i <= numel(lines) && ~startsWith(lines{i}, '*')
                v = sscanf(strrep(lines{i}, ',', ' '), '%f');
                if numel(v) >= 3
                    nodes = [nodes; v(1:3)'];
                end
                i = i + 1;
            end
            continue;
        end

        if startsWith(low, '*element') && ~contains(low, 'output')
            i = i + 1;
            while i <= numel(lines) && ~startsWith(lines{i}, '*')
                v = sscanf(strrep(lines{i}, ',', ' '), '%f');
                if numel(v) >= 4
                    elems = [elems; v(1:4)'];
                end
                i = i + 1;
            end
            continue;
        end

        if startsWith(low, '*nset')
            nsetName = extract_val(line, 'nset');
            isGen = contains(low, 'generate');
            i = i + 1;
            acc = [];
            while i <= numel(lines) && ~startsWith(lines{i}, '*')
                v = sscanf(strrep(lines{i}, ',', ' '), '%f')';
                if isGen
                    for g = 1:3:numel(v)
                        acc = [acc, v(g):v(g+2):v(g+1)];
                    end
                else
                    acc = [acc, v];
                end
                i = i + 1;
            end

            fn = fieldnames(sets);
            for f = 1:numel(fn)
                if strcmpi(nsetName, fn{f})
                    sets.(fn{f}) = unique([sets.(fn{f}), acc]);
                end
            end
            continue;
        end

        if startsWith(low, '*boundary')
            i = i + 1;
            while i <= numel(lines) && ~startsWith(lines{i}, '*')
                pts = strsplit(strtrim(lines{i}), ',');
                if numel(pts) >= 3
                    val = 0;
                    if numel(pts) >= 4, val = str2double(pts{4}); end
                    bc_data{end+1} = {str2double(pts{2}), str2double(pts{3}), val};
                end
                i = i + 1;
            end
            continue;
        end
        i = i + 1;
    end

    [~, idx] = unique(nodes(:,1), 'first');
    nodes = nodes(idx, :);

    writemtx(fullfile(outFolder, 'nodes.txt'), nodes, '%d %.6f %.6f');
    writemtx(fullfile(outFolder, 'elements.txt'), elems, '%d %d %d %d');

    writeids(fullfile(outFolder, 'top_nodes.txt'), sets.Load_Nodes);
    writeids(fullfile(outFolder, 'left_nodes.txt'), sets.Support_Left);
    writeids(fullfile(outFolder, 'right_nodes.txt'), sets.Support_Right);

    writeids(fullfile(outFolder, 'cmod1.txt'), min(sets.CMOD1));
    writeids(fullfile(outFolder, 'cmod2.txt'), max(sets.CMOD2));

    fid = fopen(fullfile(outFolder, 'boundary_conditions.txt'), 'w');
    for b = 1:numel(bc_data)
        fprintf(fid, '%d %d %.6f\n', bc_data{b}{1}, bc_data{b}{2}, bc_data{b}{3});
    end
    fclose(fid);
end

function v = extract_val(line, param)
    v = '';
    t = regexp(line, [param ' *= *([^, ]+)'], 'tokens');
    if ~isempty(t), v = t{1}{1}; end
end

function writemtx(path, M, fmt)
    fid = fopen(path, 'w');
    if ~isempty(M), fprintf(fid, [fmt '\n'], M'); end
    fclose(fid);
end

function writeids(path, ids)
    fid = fopen(path, 'w');
    if ~isempty(ids), fprintf(fid, '%d\n', ids); end
    fclose(fid);
end
