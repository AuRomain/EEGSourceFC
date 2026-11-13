function geo = cortex_geometry(S, Lk)
%CORTEX_GEOMETRY What the maps of the PSD step need: the 2-D layout of the channels, or the cortex with its regions.
%   geo = cortex_geometry(S, Lk) for a level Lk of the PSD step (Lk.level
%   'sensor' or 'source', Lk.units its channels or regions): for the
%   channels their 2-D layout (FieldTrip's ft_prepare_layout); for the
%   regions the cortex of the head (the shipped ICBM152, or S.anatomy when
%   it is one head for everyone): vertices V, faces F, normals N, each
%   vertex's region (0: none) and hemisphere (vh), the left-right axis lr,
%   each region's hemisphere (psd_adjacency) and its most outward vertex
%   (point). Drawn by cortex_view. Used by psd_figures, psd_vertex_maps
%   and other maps of the regions on the cortex.
geo = struct('level', Lk.level);
if strcmp(Lk.level, 'sensor')
    elec = get_ftElec(S.montages{1});
    [~, row] = ismember(lower(Lk.units), lower(elec.label));
    sub = struct('label', {elec.label(row)}, 'chanpos', elec.chanpos(row, :), 'elecpos', elec.elecpos(row, :), ...
        'unit', elec.unit);
    [~, lay] = evalc('ft_prepare_layout(struct(''elec'', sub, ''skipscale'', ''yes'', ''skipcomnt'', ''yes''))');
    [~, i] = ismember(lower(Lk.units), lower(lay.label));
    geo.xy = lay.pos(i, :);
    geo.lay = lay;
else
    An = [];
    if ~isempty(S.anatomy) && ~S.anatomy_individual, An = get_anatomy(S); end
    [Sc, C] = get_atlas(S.atlases{1}, An, {'Vertices', 'Faces'});
    [~, row] = ismember(lower(strtrim(Lk.units)), lower(strtrim({Sc.Label})));
    Sc = Sc(row);
    [~, ai] = psd_adjacency('source', Lk.units, S);
    V = C.Vertices;
    region = zeros(size(V, 1), 1);
    for k = 1:numel(Sc), region(Sc(k).Vertices) = k; end
    cen = cell2mat(arrayfun(@(s) mean(V(s.Vertices, :), 1), Sc(:), 'UniformOutput', false));
    lr = mean(cen(ai.hemi == 1, :), 1) - mean(cen(ai.hemi == 2, :), 1);     % right to left
    lr = lr / norm(lr);
    mid = (mean(cen(ai.hemi == 1, :), 1) + mean(cen(ai.hemi == 2, :), 1)) / 2;
    vh = 1 + ((V - mid) * lr(:) < 0);                 % vertex hemisphere: 1 left, 2 right
    vh(region > 0) = ai.hemi(region(region > 0));
    point = zeros(numel(Sc), 1);                      % each region's most outward vertex (from the
    centre = mean(V, 1);                              % brain's centre through the region's centroid)
    for k = 1:numel(Sc)
        u = cen(k, :) - centre;  u = u / norm(u);
        [~, j] = max((V(Sc(k).Vertices, :) - centre) * u(:));
        point(k) = Sc(k).Vertices(j);
    end
    N = C.VertNormals;
    geo = struct('level', 'source', 'V', V, 'F', C.Faces, 'N', N, 'region', region, 'vh', vh, 'lr', lr, ...
        'hemi', ai.hemi, 'point', point);
end
end
