function [A, info] = psd_adjacency(level, labels, S)
%PSD_ADJACENCY Which channels, or which regions, are neighbours (for the cluster tests).
%   [A, info] = psd_adjacency('sensor', channels, S): A (channels x
%   channels, logical) from the positions of the cap S.montages{1}, as
%   FieldTrip's ft_prepare_neighbours with method 'triangulation' makes
%   them (the cluster tests of FieldTrip: a Delaunay triangulation of the
%   electrodes projected on a plane, also compressed in x and in y, links
%   longer than the mean + 3 SD left out); info.n_neighbours per channel.
%   [A, info] = psd_adjacency('source', regions, S): two regions of the
%   atlas S.atlases{1} are neighbours when they touch on the cortex (two
%   adjacent vertices of the mesh, one in each: the cortex file's VertConn,
%   Brainstorm's tess_vertconn; from the edges of its triangles when a
%   surface has no VertConn), never across the hemispheres (labels ending
%   in L and R, or the sides of the midline).
%   The cortex is that of the shipped ICBM152, or of S.anatomy (one head
%   for every recording); with individual heads the shipped one, whose
%   regions are the same (as atlas_regions).
%   info.pairs: the neighbour pairs as text ('Fz-FCz', ...), for the log;
%   info.hemi: each region's hemisphere (1 left, 2 right; [] for channels).
labels = cellstr(labels(:)');
K = numel(labels);
hemi = [];
switch level
    case 'sensor'
        elec = get_ftElec(S.montages{1});
        [ok, row] = ismember(lower(labels), lower(elec.label));
        if ~all(ok)
            error('psd_adjacency:channels', 'Channels not in the cap %s: %s', S.montages{1}, strjoin(labels(~ok), ', '));
        end
        sub = struct('label', {elec.label(row)}, 'chanpos', elec.chanpos(row, :), 'elecpos', elec.elecpos(row, :), ...
            'unit', elec.unit);
        cfg = struct('method', 'triangulation', 'elec', sub, 'feedback', 'no', 'compress', 'yes');
        [~, nb] = evalc('ft_prepare_neighbours(cfg)');
        A = false(K);
        for k = 1:numel(nb)
            i = find(strcmpi(labels, nb(k).label));
            j = find(ismember(lower(labels), lower(nb(k).neighblabel)));
            A(i, j) = true;
        end
    case 'source'
        An = [];
        if ~isempty(S.anatomy) && ~S.anatomy_individual, An = get_anatomy(S); end
        [Sc, C] = get_atlas(S.atlases{1}, An, {'Vertices', 'Faces', 'VertConn'});
        [ok, row] = ismember(lower(strtrim(labels)), lower(strtrim({Sc.Label})));
        if ~all(ok)
            error('psd_adjacency:regions', 'Regions not in the atlas %s: %s', S.atlases{1}, strjoin(labels(~ok), ', '));
        end
        Sc = Sc(row);
        region = zeros(size(C.Vertices, 1), 1);
        for k = 1:K, region(Sc(k).Vertices) = k; end
        if isfield(C, 'VertConn') && ~isempty(C.VertConn)
            [vi, vj] = find(C.VertConn);              % Brainstorm's vertex adjacency (tess_vertconn)
            pairs = [region(vi), region(vj)];
        else                                          % the same from the faces' edges, as tess_vertconn makes it
            R = region(C.Faces);
            pairs = [R(:, [1 2]); R(:, [2 3]); R(:, [1 3])];
        end
        pairs = pairs(all(pairs > 0, 2) & pairs(:, 1) ~= pairs(:, 2), :);
        A = false(K);
        A(sub2ind([K K], pairs(:, 1), pairs(:, 2))) = true;
        hemi = hemispheres(labels, Sc, C);
        A(hemi ~= hemi') = false;
    otherwise
        error('psd_adjacency:level', 'level: ''sensor'' or ''source''');
end
A = A | A';
A(1:K + 1:end) = false;
[i, j] = find(triu(A));
info = struct('n_neighbours', sum(A, 2), 'pairs', {strcat(labels(i), '-', labels(j))}, 'hemi', hemi);
end

function h = hemispheres(labels, Sc, C)
% 1 left, 2 right: the label's last word (L, R, lh, rh); a region without
% one, the side of its centroid relative to the midline of the others;
% labels without sides, the connected parts of the cortex mesh (the two
% hemispheres are separate surfaces in Brainstorm's cortices)
h = zeros(numel(labels), 1);
last = lower(regexp(labels, '[A-Za-z]+$', 'match', 'once'));
h(ismember(last, {'l', 'lh', 'left'})) = 1;
h(ismember(last, {'r', 'rh', 'right'})) = 2;
if all(h > 0), return; end
if all(h == 0)
    F = double(C.Faces);
    G = graph([F(:, 1); F(:, 2); F(:, 1)], [F(:, 2); F(:, 3); F(:, 3)], [], size(C.Vertices, 1));
    comp = conncomp(simplify(G));
    for k = 1:numel(Sc), h(k) = mode(comp(Sc(k).Vertices)); end
    [~, ~, h] = unique(h);                         % one number per part of the mesh
    return
end
c = cell2mat(arrayfun(@(s) mean(C.Vertices(s.Vertices, :), 1), Sc(:), 'UniformOutput', false));
d = mean(c(h == 1, :), 1) - mean(c(h == 2, :), 1);  % from right to left
m = (mean(c(h == 1, :), 1) + mean(c(h == 2, :), 1)) / 2;
side = (c - m) * d(:);
h(h == 0 & side >= 0) = 1;
h(h == 0 & side < 0) = 2;
end
