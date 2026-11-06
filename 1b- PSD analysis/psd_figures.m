function psd_figures(S, L, stats, fig_dir)
%PSD_FIGURES Figures of the PSD analysis: group spectra, box plots, scalp and cortex maps, clusters.
%   psd_figures(S, L, stats, fig_dir) (run_psd_analysis, with S.figures).
%   The groups are those of each comparison of S.comparisons ("A - B":
%   names A and B), each name in the same colour in every figure; without
%   comparisons, all recordings together. Per level (sensor,
%   source_<inverse>):
%   spectra_<level>_<comparison>.png   each group's mean spectrum with its
%       95% CI (t over recordings): log10 PSD averaged over the channels or
%       regions, and the same above the 1/f fit (log10 PSD minus the
%       specparam aperiodic fit), the edges of S.psd.bands marked
%   boxplots_<level>_<measure>_<comparison>.png   the measure averaged over
%       the channels or regions, per band (or exponent, offset, alpha
%       frequency): box, each recording, n; Hedges' g in the title is
%       descriptive (the tests are per channel or region, in the maps)
%   map_<level>_<measure>_<comparison>.png   Hedges' g of each channel
%       (<measure> aperiodic_alpha: exponent, offset and alpha frequency)
%       (scalp map, FieldTrip's ft_plot_topo, interpolated inside the head)
%       or region (cortex: lateral and medial views of each hemisphere,
%       and the whole cortex from above and from below, L and R marked),
%       a row or panel per band, one colour scale per figure (blue: A < B,
%       red: A > B). Filled black dots: significant channel or region (p_bb <=
%       S.alpha); open circles: in a significant cluster (cluster_q <=
%       S.alpha). On the cortex the regions neither significant nor in a
%       significant cluster are paler
%   map_source_<inverse>_..._<comparison>.fig   the same cortex maps in 3-D,
%       one brain per band, to turn with the mouse in MATLAB (openfig, or
%       a double click on the file)
%   map_<level>_<measure>_<set>_<scale>.png   the same with partial r, for
%       the correlations with a significant channel, region or cluster
%   cluster_<level>_<n>.png   each significant cluster (clusters.csv, row
%       n): the mean of its channels or regions per group (comparisons), or
%       against the score (correlations), with its members named
if ~isfolder(fig_dir), mkdir(fig_dir); end
for f = [dir(fullfile(fig_dir, '*.png')); dir(fullfile(fig_dir, '*.fig'))]'   % none left from an earlier run
    delete(fullfile(f.folder, f.name));
end
o = S.psd;
T = [];
try
    T = read_subject_table(S);
catch
end
cmap = diverging_map();
for k = 1:numel(L)
    Gr = groups(S, T, L(k));
    geo = geometry(S, L(k));
    for c = 1:numel(Gr)
        spectra_figure(L(k), Gr(c), o, fig_dir);
        for m = intersect({'absolute', 'relative', 'periodic'}, o.measures, 'stable')
            box_figure(L(k), Gr(c), m{1}, o, fig_dir);
        end
        if any(ismember({'aperiodic', 'alpha_frequency'}, o.measures))
            box_figure(L(k), Gr(c), 'other', o, fig_dir);
        end
    end
    if isempty(stats.res), continue; end
    res = stats.res(strcmp({stats.res.level}, L(k).name));
    for m = map_measures(o)
        sel = res(strcmp({res.measure}, m{1}) | (strcmp(m{1}, 'other') & ismember({res.measure}, ...
            {'exponent', 'offset', 'alpha_frequency'})));
        if isempty(sel), continue; end
        for cmp = S.comparisons(:)'                   % group comparisons
            P = panels(sel, 'G', @(X) X.Comparison == string(cmp.name), 'g');
            if isempty(P), continue; end
            map_figure(L(k), geo, P, cmap, sprintf('%s, %s: Hedges'' g (%s)', L(k).name, map_title(m{1}), ...
                cmp.name), 'Hedges'' g', fullfile(fig_dir, sprintf('map_%s_%s_%s.png', L(k).name, map_file(m{1}), ...
                safe(cmp.name))));
        end
        tests = strings(0, 2);                        % correlations with something significant
        for r = sel(:)'
            X = r.C;
            if isempty(X) || height(X) == 0, continue; end
            hit = X.significant;
            if ismember('cluster_significant', X.Properties.VariableNames), hit = hit | X.cluster_significant; end
            tests = [tests; [X.CorrelationSet(hit), X.Scale(hit)]]; %#ok<AGROW>
        end
        if isempty(tests), continue; end
        tests = unique(tests, 'rows');
        for t = 1:size(tests, 1)
            P = panels(sel, 'C', @(X) X.CorrelationSet == tests(t, 1) & X.Scale == tests(t, 2), 'r');
            map_figure(L(k), geo, P, cmap, sprintf('%s, %s: partial r with %s (%s)', L(k).name, map_title(m{1}), ...
                tests(t, 2), tests(t, 1)), 'r', fullfile(fig_dir, sprintf('map_%s_%s_%s_%s.png', L(k).name, ...
                map_file(m{1}), safe(tests(t, 1)), safe(tests(t, 2)))));
        end
    end
    if isfield(stats, 'clusters') && height(stats.clusters)
        cluster_figures(S, T, L(k), Gr, stats.clusters, fig_dir);
    end
end
end


%% ===== groups and geometry =====
function Gr = groups(S, T, Lk)
% per comparison: its two groups over the level's recordings (names,
% members, colours: Okabe-Ito hues, one per group name across comparisons)
pal = [213 94 0; 0 114 178; 170 51 119; 0 158 115; 230 159 0; 86 180 233; 102 102 102] / 255;
n = numel(Lk.ids);
Gr = struct('name', {}, 'names', {}, 'in', {}, 'colours', {});
if ~isempty(T) && ~isempty(S.comparisons)
    [ok, row] = ismember(string(Lk.ids), string(T.(S.id_column)));
    seen = strings(0);
    for cmp = S.comparisons(:)'
        names = split(string(cmp.name), " - ");
        if numel(names) ~= 2, names = string(cmp.name) + [" (A)"; " (B)"]; end
        in = false(n, 2);
        in(ok, 1) = select_subjects(T(row(ok), :), cmp.A);
        in(ok, 2) = select_subjects(T(row(ok), :), cmp.B);
        col = zeros(2, 3);
        for j = 1:2
            c = find(seen == names(j), 1);
            if isempty(c), seen(end + 1) = names(j); c = numel(seen); end %#ok<AGROW>
            col(j, :) = pal(min(c, end), :);
        end
        Gr(end + 1) = struct('name', string(cmp.name), 'names', names(:)', 'in', in, 'colours', col); %#ok<AGROW>
    end
end
if isempty(Gr)
    Gr = struct('name', "all recordings", 'names', "all", 'in', true(n, 1), 'colours', [0.3 0.3 0.3]);
end
end

function geo = geometry(S, Lk)
% what the maps need: the 2-D layout of the channels, or the cortex with
% its regions, hemispheres and the surface point of each region
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

%% ===== spectra and box plots =====
function spectra_figure(Lk, g, o, fig_dir)
f = Lk.freq;
A = squeeze(mean(log10(double(Lk.pow)), 2, 'omitnan'));      % recordings x F
B = squeeze(mean(double(Lk.flat), 2, 'omitnan'));
if isvector(A), A = A(:)'; B = B(:)'; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1100 430]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
unit = 'microvolts^2/Hz';
if strcmp(Lk.level, 'source'), unit = 'a.u.'; end
ylab = {sprintf('log10 PSD (%s), mean of the %s', unit, unit_word(Lk.level)), ...
    'log10 PSD above the 1/f fit (specparam)'};
for p = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    if p == 1, Y = A; else, Y = B; end
    h = gobjects(0);  lab = strings(0);
    for j = 1:size(g.in, 2)
        y = Y(g.in(:, j), :);
        y = y(all(~isnan(y(:, ~all(isnan(y), 1))), 2), :);
        nn = size(y, 1);
        if nn == 0, continue; end
        m = mean(y, 1, 'omitnan');
        if nn > 1
            ci = tinv(0.975, nn - 1) * std(y, 0, 1, 'omitnan') / sqrt(nn);
            ok = ~isnan(m);
            fill(ax, [f(ok), fliplr(f(ok))], [m(ok) - ci(ok), fliplr(m(ok) + ci(ok))], g.colours(j, :), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h(end + 1) = plot(ax, f, m, '-', 'Color', g.colours(j, :), 'LineWidth', 1.5); %#ok<AGROW>
        lab(end + 1) = sprintf('%s (n = %d)', g.names(j), nn); %#ok<AGROW>
    end
    band_marks(ax, o, Lk);
    xlim(ax, [f(1), f(end)]);
    xlabel(ax, 'Frequency (Hz)');  ylabel(ax, ylab{p});
    if p == 2, yline(ax, 0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off'); end
    legend(ax, h, lab, 'Location', 'northeast', 'Interpreter', 'none');
    box(ax, 'on');
end
title(tl, sprintf('%s: mean spectra and 95%% CI (%s)', Lk.name, g.name), 'Interpreter', 'none');
exportgraphics(fig, fullfile(fig_dir, sprintf('spectra_%s_%s.png', Lk.name, safe(g.name))), 'Resolution', 150);
close(fig);
end

function band_marks(ax, o, Lk)
% the fixed edges of the bands (anchored bands: the median of their edges)
% as dotted lines; the names of the bands that no other band contains
e = squeeze(median(Lk.M.edges, 1, 'omitnan'));
if isvector(e), e = reshape(e, [], 2); end
for x = unique(e(~isnan(e)))'
    xline(ax, x, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
end
yl = ylim(ax);
for b = 1:size(e, 1)
    inside = e(:, 1) <= e(b, 1) & e(:, 2) >= e(b, 2) & (1:size(e, 1))' ~= b;
    if any(inside) || any(isnan(e(b, :))), continue; end
    text(ax, mean(e(b, :)), yl(2), o.bands(b).name, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.35 0.35 0.35], 'Interpreter', 'none');
end
end

function box_figure(Lk, g, m, o, fig_dir)
if strcmp(m, 'other')
    names = {};
    if ismember('aperiodic', o.measures), names = {'exponent', 'offset'}; end
    if ismember('alpha_frequency', o.measures), names{end + 1} = 'alpha_frequency'; end
    Y = cellfun(@(x) mean(Lk.M.(x), 2, 'omitnan'), names, 'UniformOutput', false);
    ylab = strrep(names, '_', ' ');
    file = sprintf('boxplots_%s_aperiodic_alpha_%s.png', Lk.name, safe(g.name));
else
    names = {o.bands.name};
    Y = arrayfun(@(b) mean(Lk.M.(m)(:, :, b), 2, 'omitnan'), 1:numel(names), 'UniformOutput', false);
    ylab = repmat({measure_label(m, Lk.level)}, 1, numel(names));
    file = sprintf('boxplots_%s_%s_%s.png', Lk.name, m, safe(g.name));
end
np = numel(names);
nc = min(np, 4);  nr = ceil(np / nc);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 300 * nc, 300 * nr]);
tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:np
    ax = nexttile(tl);  hold(ax, 'on');
    y = Y{p};
    xt = strings(1, size(g.in, 2));
    for j = 1:size(g.in, 2)
        v = y(g.in(:, j));  v = v(~isnan(v));
        xt(j) = sprintf('%s (n = %d)', g.names(j), numel(v));
        if isempty(v), continue; end
        boxchart(ax, repmat(j, numel(v), 1), v, 'BoxFaceColor', g.colours(j, :), 'WhiskerLineColor', ...
            g.colours(j, :), 'MarkerStyle', 'none', 'BoxWidth', 0.5);
        swarmchart(ax, repmat(j, numel(v), 1), v, 14, g.colours(j, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'XJitterWidth', 0.35);
    end
    set(ax, 'XTick', 1:size(g.in, 2), 'XTickLabel', xt, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    xlim(ax, [0.4, size(g.in, 2) + 0.6]);
    ylabel(ax, ylab{p}, 'Interpreter', 'none');
    t = strrep(names{p}, '_', ' ');
    if size(g.in, 2) == 2
        t = sprintf('%s: g = %.2f', t, hedges(y(g.in(:, 1)), y(g.in(:, 2))));   % A - B, as the tests
    end
    title(ax, t, 'Interpreter', 'none', 'FontSize', 9);
    box(ax, 'on');
end
title(tl, sprintf('%s, mean of the %s (%s); g descriptive, tests per %s in the maps', Lk.name, ...
    unit_word(Lk.level), g.name, unit_word(Lk.level, 1)), 'Interpreter', 'none', 'FontSize', 10);
exportgraphics(fig, fullfile(fig_dir, file), 'Resolution', 150);
close(fig);
end

%% ===== maps =====
function ms = map_measures(o)
ms = intersect({'absolute', 'relative', 'periodic'}, o.measures, 'stable');
if any(ismember({'aperiodic', 'alpha_frequency'}, o.measures)), ms{end + 1} = 'other'; end
end

function s = map_title(m)
s = m;
if strcmp(m, 'other'), s = 'exponent, offset, alpha frequency'; end
end

function s = map_file(m)
s = m;
if strcmp(m, 'other'), s = 'aperiodic_alpha'; end
end

function P = panels(sel, kind, pick, col)
% one panel per family of sel: values (col: g or r), significant units,
% members of significant clusters
P = struct('title', {}, 'val', {}, 'sig', {}, 'clu', {});
for r = sel(:)'
    X = r.(kind);
    if isempty(X) || height(X) == 0, continue; end
    X = X(pick(X), :);
    if height(X) == 0, continue; end
    K = max(X.state);
    v = nan(K, 1);  s = false(K, 1);  c = false(K, 1);
    v(X.state) = X.(col);
    s(X.state) = X.significant;
    if ismember('cluster_significant', X.Properties.VariableNames), c(X.state) = X.cluster_significant; end
    t = r.band;
    if isempty(t), t = r.measure; end
    P(end + 1) = struct('title', strrep(t, '_', ' '), 'val', v, 'sig', s, 'clu', c); %#ok<AGROW>
end
end

function map_figure(Lk, geo, P, cmap, ttl, clabel, file)
K = numel(Lk.units);
for p = 1:numel(P)                                    % every unit (one never tested: NaN)
    for f = {'val', 'sig', 'clu'}
        x = P(p).(f{1});
        if numel(x) < K
            if islogical(x), x(end + 1:K) = false; else, x(end + 1:K) = NaN; end
        end
        P(p).(f{1}) = x(:);
    end
end
m = max(abs(vertcat(P.val)), [], 'omitnan');
if isempty(m) || ~(m > 0), m = 1; end
clim = [-m m];
np = numel(P);
if strcmp(geo.level, 'sensor')
    nc = min(np, 4);  nr = ceil(np / nc);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 280 * nc + 80, 290 * nr + 40]);
    tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:np
        ax = nexttile(tl);
        topo(ax, geo, P(p), clim, cmap);
        title(ax, P(p).title, 'Interpreter', 'none', 'FontSize', 9);
    end
else
    views = {'left lateral', 'left medial', 'right medial', 'right lateral', 'dorsal', 'ventral'};
    names = [views(1:4), {'dorsal (L | R)', 'ventral (R | L)'}];   % from below the left is on the right
    nv = numel(views);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 nv * 215 + 80, 190 * np + 50]);
    tl = tiledlayout(fig, np, nv, 'TileSpacing', 'tight', 'Padding', 'compact');
    for p = 1:np
        for v = 1:nv
            ax = nexttile(tl);
            cortex(ax, geo, P(p), clim, cmap, views{v});
            if v == 1
                text(ax, 0, 1.02, P(p).title, 'Units', 'normalized', 'FontSize', 9, 'FontWeight', 'bold', ...
                    'Interpreter', 'none', 'VerticalAlignment', 'bottom');
            end
            if p == 1, title(ax, names{v}, 'FontSize', 8, 'FontWeight', 'normal'); end
        end
    end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = clabel;
title(tl, {strrep(ttl, '_', '\_'), 'filled dot: significant (q_{bb} \leq \alpha); open circle: in a significant cluster'}, ...
    'FontSize', 10, 'Interpreter', 'tex');
exportgraphics(fig, file, 'Resolution', 150);
close(fig);
if ~strcmp(geo.level, 'sensor'), cortex_3d(geo, P, cmap, clim, ttl, clabel, strrep(file, '.png', '.fig')); end
end

function cortex_3d(geo, P, cmap, clim, ttl, clabel, file)
% the cortex maps in 3-D, one brain per band, saved as a MATLAB figure that
% opens visible and can be turned with the mouse
np = numel(P);
nc = min(np, 4);  nr = ceil(np / nc);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 330 * nc + 90, 300 * nr + 70], ...
    'Name', ttl, 'NumberTitle', 'off');
tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:np
    ax = nexttile(tl);
    cortex(ax, geo, P(p), clim, cmap, '3d');
    title(ax, P(p).title, 'Interpreter', 'none', 'FontSize', 10);
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = clabel;
title(tl, {strrep(ttl, '_', '\_'), ['filled dot: significant (q_{bb} \leq \alpha); open circle: in a ' ...
    'significant cluster. Turn a brain with the mouse']}, 'FontSize', 10, 'Interpreter', 'tex');
set(fig, 'CreateFcn', ['set(gcbo, ''Visible'', ''on''); r = rotate3d(gcbo); ' ...   % when opened: visible,
    'r.ActionPostCallback = @(~, e) camlight(findobj(e.Axes, ''Type'', ''light''), ''headlight''); ' ...   % turned
    'r.Enable = ''on''; clear r']);                   % with the mouse, the light following the camera
savefig(fig, file);                                % ('compact': no effect from R2026a, which warns)
close(fig);
end

function topo(ax, geo, P, clim, cmap)
hold(ax, 'on');
lay = geo.lay;
ok = ~isnan(P.val);
if nnz(ok) >= 3
    axes(ax); %#ok<LAXES>                          % ft_plot_topo draws in the current axes
    ft_plot_topo(geo.xy(ok, 1), geo.xy(ok, 2), P.val(ok), 'mask', lay.mask, 'outline', lay.outline, ...
        'interplim', 'mask', 'gridscale', 120, 'clim', clim, 'parent', ax);
end
for i = 1:numel(lay.outline)
    plot(ax, lay.outline{i}(:, 1), lay.outline{i}(:, 2), 'k-', 'LineWidth', 1);
end
plot(ax, geo.xy(:, 1), geo.xy(:, 2), '.', 'Color', [0.25 0.25 0.25], 'MarkerSize', 4);
c = P.clu & ~P.sig;
plot(ax, geo.xy(c, 1), geo.xy(c, 2), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 5, ...
    'LineWidth', 0.9);
plot(ax, geo.xy(P.sig, 1), geo.xy(P.sig, 2), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
colormap(ax, cmap);
caxis(ax, clim); %#ok<CAXIS>
axis(ax, 'equal', 'off');
end

function cortex(ax, geo, P, clim, cmap, view_name)
% one view of the cortex: 'left lateral', 'left medial', 'right medial',
% 'right lateral' (that hemisphere), 'dorsal', 'ventral' (both, the front
% at the top), or '3d' (both, lit from every side, to be turned)
hold(ax, 'on');
both = ismember(view_name, {'dorsal', 'ventral', '3d'});
left = startsWith(view_name, 'left');
lateral = endsWith(view_name, 'lateral');
h = 1 + ~left;
grey = [0.78 0.78 0.78];
col = repmat(grey, size(geo.V, 1), 1);
in = geo.region > 0;
v = P.val(geo.region(in));
c = colour(v, clim, cmap);
c(isnan(v), :) = repmat(grey, nnz(isnan(v)), 1);
pale = ~(P.sig(geo.region(in)) | P.clu(geo.region(in)));
c(pale, :) = 0.6 * c(pale, :) + 0.4 * repmat(grey, nnz(pale), 1);
col(in, :) = c;
if both, F = geo.F; else, F = geo.F(all(geo.vh(geo.F) == h, 2), :); end
patch(ax, 'Faces', F, 'Vertices', geo.V, 'FaceVertexCData', col, 'FaceColor', 'interp', 'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', 'AmbientStrength', 0.45, 'DiffuseStrength', 0.65, 'SpecularStrength', 0.05);
up = [0 0 1] - geo.lr * geo.lr(3);  up = up / norm(up);   % superior
ant = cross(geo.lr, up);                              % anterior (left x superior)
switch view_name
    case 'dorsal',  cam = up;  top = ant;
    case 'ventral', cam = -up;  top = ant;
    case '3d',      cam = geo.lr + 0.35 * up - 0.25 * ant;  top = up;   % from the left, a little above
    otherwise                                         % from the left (+lr) or the right
        if left == lateral, cam = geo.lr; else, cam = -geo.lr; end
        top = up;
end
view(ax, cam);
camup(ax, top);
axis(ax, 'equal', 'off');
camlight(ax, 'headlight');                            % in 3-D it follows the camera (cortex_3d)
% marks: in a fixed view on each region's vertex nearest the camera (a region
% that cannot be seen from there is hidden by the surface in front of it);
% in 3-D on its most outward vertex
if both, k = find(P.sig | P.clu); else, k = find(geo.hemi(:) == h & (P.sig | P.clu)); end
k = k(:);
pt = geo.point(k);
if ~strcmp(view_name, '3d')
    toward = cam(:) / norm(cam);
    for i = 1:numel(k)
        vk = find(geo.region == k(i));
        [~, j] = max(geo.V(vk, :) * toward);
        pt(i) = vk(j);
    end
    x = geo.V(pt, :) + 0.002 * toward';               % 2 mm towards the camera
else
    x = geo.V(pt, :) + 0.002 * geo.N(pt, :);
end
s = P.sig(k);
plot3(ax, x(~s, 1), x(~s, 2), x(~s, 3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot3(ax, x(s, 1), x(s, 2), x(s, 3), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
colormap(ax, cmap);
caxis(ax, clim); %#ok<CAXIS>
end

function c = colour(v, clim, cmap)
n = size(cmap, 1);
i = round((v - clim(1)) / diff(clim) * (n - 1)) + 1;
i = min(max(i, 1), n);
i(isnan(i)) = 1;
c = cmap(i, :);
end

function cmap = diverging_map()
% blue - white - red (ColorBrewer RdBu, 11 classes, interpolated)
k = [5 48 97; 33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; 253 219 199; 244 165 130; ...
    214 96 77; 178 24 43; 103 0 31] / 255;
cmap = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

%% ===== clusters =====
function cluster_figures(S, T, Lk, Gr, Cl, fig_dir)
for i = find(Cl.level == string(Lk.name) & Cl.cluster_significant)'
    units = split(Cl.units(i), " ");
    in = ismember(lower(string(Lk.units)), lower(units));
    b = find(strcmp({S.psd.bands.name}, Cl.band(i)));
    if isempty(b), V = Lk.M.(char(Cl.measure(i))); else, V = Lk.M.(char(Cl.measure(i)))(:, :, b); end
    y = mean(V(:, in), 2, 'omitnan');
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 520 430]);
    ax = axes(fig);  hold(ax, 'on');
    if Cl.kind(i) == "comparison"
        g = Gr(strcmp([Gr.name], Cl.test(i)));
        if isempty(g), close(fig); continue; end
        xt = strings(1, 2);
        for j = 1:2
            v = y(g.in(:, j));  v = v(~isnan(v));
            xt(j) = sprintf('%s (n = %d)', g.names(j), numel(v));
            boxchart(ax, repmat(j, numel(v), 1), v, 'BoxFaceColor', g.colours(j, :), 'WhiskerLineColor', ...
                g.colours(j, :), 'MarkerStyle', 'none', 'BoxWidth', 0.5);
            swarmchart(ax, repmat(j, numel(v), 1), v, 18, g.colours(j, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
                'XJitterWidth', 0.35);
        end
        set(ax, 'XTick', 1:2, 'XTickLabel', xt, 'TickLabelInterpreter', 'none');
        xlim(ax, [0.4 2.6]);
        eff = sprintf('mean g %.2f', Cl.mean_effect(i));
    else
        if isempty(T), close(fig); continue; end
        [ok, row] = ismember(string(Lk.ids), string(T.(S.id_column)));
        x = nan(numel(y), 1);
        x(ok) = double(T.(char(Cl.scale(i)))(row(ok)));
        cs = S.correlations(strcmp({S.correlations.name}, Cl.test(i)));
        if ~isempty(cs)
            tested = false(numel(y), 1);
            tested(ok) = select_subjects(T(row(ok), :), cs.subjects);
            x(~tested) = NaN;
        end
        scatter(ax, x, y, 26, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.7);
        xlabel(ax, Cl.scale(i), 'Interpreter', 'none');
        eff = sprintf('mean partial r %.2f', Cl.mean_effect(i));
    end
    ylabel(ax, sprintf('%s, cluster mean', measure_label(char(Cl.measure(i)), Lk.level)), 'Interpreter', 'none');
    members = strjoin(units, ' ');
    if strlength(members) > 90, members = extractBefore(members, 88) + " ..."; end
    title(ax, {sprintf('%s %s, %s: cluster %d (%d %s), %s, p = %.4f, q = %.4f', Lk.name, ...
        strtrim(Cl.measure(i) + " " + Cl.band(i)), Cl.test(i), Cl.cluster(i), Cl.size(i), unit_word(Lk.level), ...
        eff, Cl.cluster_p(i), Cl.cluster_q(i)), members}, 'Interpreter', 'none', 'FontSize', 8);
    box(ax, 'on');
    exportgraphics(fig, fullfile(fig_dir, sprintf('cluster_%s_%d.png', Lk.name, i)), 'Resolution', 150);
    close(fig);
end
end

%% ===== small helpers =====
function g = hedges(a, b)
% Hedges' g of a - b (descriptive; hedges_g)
g = hedges_g(a, b);
end

function s = measure_label(m, level)
unit = 'microvolts^2';
if strcmp(level, 'source'), unit = 'a.u.'; end
switch m
    case 'absolute', s = sprintf('log10 band power (%s)', unit);
    case 'relative', s = 'relative power (%)';
    case 'periodic', s = 'log10 power above the 1/f fit';
    case 'exponent', s = '1/f exponent';
    case 'offset', s = '1/f offset (log10)';
    case 'alpha_frequency', s = 'alpha frequency (Hz)';
    otherwise, s = m;
end
end

function w = unit_word(level, one)
w = 'channels';
if strcmp(level, 'source'), w = 'regions'; end
if nargin > 1, w = w(1:end - 1); end
end

function s = safe(x)
% a name for a file: 'patient - control' -> 'patient_control'
s = char(regexprep(strtrim(string(x)), '[^A-Za-z0-9]+', '_'));
end
