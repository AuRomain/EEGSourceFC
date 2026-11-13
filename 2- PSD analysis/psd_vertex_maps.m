function V = psd_vertex_maps(S, L, out)
%PSD_VERTEX_MAPS Band power of every source of the cortex: values and maps (descriptive, not tested).
%   V = psd_vertex_maps(S, L, out) (run_psd_analysis, S.psd.vertex_maps)
%   for every inverse method of the regions' levels of L: the PSD of each
%   source (vertex) of each recording (pow_vertex of psd_spectra: w' C(f) w,
%   w the source's inverse filter, C the channels' cross-spectra: each
%   source's power from the cross-spectral matrix, as LORETA's band images,
%   Frei et al. 2001), and from it, per source and band of S.psd.bands
%   (psd_measures, the recording's band edges): 'absolute' (log10 band
%   power) and 'relative' (% of S.psd.total), those of S.psd.measures.
%   Written to <out>/vertex_values_<inverse>.mat (struct V): the
%   recordings (stems, ids), the bands, the values (recordings x sources x
%   bands, single), per comparison of S.comparisons each group's mean and
%   Hedges' g (A - B, no covariates), per scale of S.correlations (its
%   subjects, zero_as_missing) the r with the source's value
%   (S.correlation_type, no covariates). These describe; they are not
%   tested (the tests are on the channels and regions).
%   With S.figures, in <out>/figures:
%     vertex_<measure>_<inverse>_<comparison>.png  one row per band: each
%       group's mean (left, top and right), then g (blue: A < B, red: A > B)
%     vertex_<measure>_<inverse>.png  without comparisons: all recordings'
%       mean
%     vertex_<measure>_<inverse>_<comparison>.fig  the g maps in 3-D, a
%       brain per band, to turn with the mouse
%     vertex_<measure>_<inverse>_r_<set>_<scale>.png  r per band (up to
%       10 scales in all; more: in the .mat only)
%   Needs one cortex for every recording: the shipped head, or S.anatomy
%   without S.anatomy_individual; otherwise, or when the spectra were
%   computed without S.psd.vertex_maps, nothing is written (a message says
%   why).
V = struct([]);
o = S.psd;
if ~o.vertex_maps, return; end
if S.anatomy_individual
    fprintf('Vertex maps: not made, each recording has its own cortex (S.anatomy_individual)\n');
    return
end
meas = intersect({'absolute', 'relative'}, o.measures, 'stable');
if isempty(meas), return; end
T = [];
try
    T = read_subject_table(S);
catch
end
bands = string({o.bands.name});
nb = numel(bands);
for k = find(strcmp({L.level}, 'source'))
    Lk = L(k);
    nr = numel(Lk.stems);
    vals = struct();
    nv = [];
    for i = 1:nr
        c = fullfile(S.results, 'PSD', 'spectra', 'source', [Lk.stems{i} '.mat']);
        w = load(c, 'r');
        x = w.r.inverses(strcmp({w.r.inverses.inverse}, Lk.inverse));
        if ~isfield(x, 'pow_vertex') || isempty(x.pow_vertex)
            fprintf('Vertex maps (%s): not made, %s has no sources'' spectra (S.psd.vertex_maps was false)\n', ...
                Lk.inverse, Lk.stems{i});
            vals = [];
            break
        end
        if isempty(nv)
            nv = size(x.pow_vertex, 1);
            for m = meas, vals.(m{1}) = nan(nr, nv, nb, 'single'); end
        elseif size(x.pow_vertex, 1) ~= nv
            fprintf('Vertex maps (%s): not made, the recordings have different cortices\n', Lk.inverse);
            vals = [];
            break
        end
        M = psd_measures(x.freq, double(x.pow_vertex), [], o, Lk.iaf(i), S.iaf);
        for m = meas, vals.(m{1})(i, :, :) = single(M.(m{1})); end
    end
    if isempty(vals), continue; end
    geo = cortex_geometry(S, Lk);
    if size(geo.V, 1) ~= nv
        fprintf('Vertex maps (%s): not made, %d sources but %d vertices in the cortex of the maps\n', ...
            Lk.inverse, nv, size(geo.V, 1));
        continue
    end
    Vk = struct('inverse', Lk.inverse, 'stems', {Lk.stems}, 'ids', Lk.ids, 'bands', bands, ...
        'edges', Lk.M.edges, 'values', vals, 'groups', struct([]), 'correlations', struct([]));
    % the comparisons: means and g
    Gr = psd_groups(S, T, Lk);
    for c = 1:numel(Gr)
        G = struct('name', Gr(c).name, 'names', Gr(c).names, 'n', sum(Gr(c).in, 1));
        for m = meas
            X = double(vals.(m{1}));
            mu = nan(nv, nb, size(Gr(c).in, 2));
            for j = 1:size(Gr(c).in, 2), mu(:, :, j) = squeeze(mean(X(Gr(c).in(:, j), :, :), 1, 'omitnan')); end
            G.(m{1}) = struct('mean', single(mu), 'g', []);
            if size(Gr(c).in, 2) == 2
                G.(m{1}).g = single(g_columns(X(Gr(c).in(:, 1), :, :), X(Gr(c).in(:, 2), :, :)));
            end
        end
        Vk.groups = [Vk.groups, G];
    end
    % the correlations: r with each scale
    if ~isempty(T)
        [okT, row] = ismember(string(Lk.ids), string(T.(S.id_column)));
        names_T = string(T.Properties.VariableNames);
        for cs = S.correlations(:)'
            scales = names_T(~cellfun(@isempty, regexp(cellstr(names_T), cs.columns, 'once')));
            scales = setdiff(scales, [string(cs.covariates), string(S.id_column), string(S.group_column)], 'stable');
            scales = scales(arrayfun(@(sc) isnumeric(T.(sc)), scales));
            insub = false(nr, 1);
            insub(okT) = select_subjects(T(row(okT), :), cs.subjects);
            for sc = scales
                xs = nan(nr, 1);
                xs(okT) = double(T.(sc)(row(okT)));
                if cs.zero_as_missing, xs(xs == 0) = NaN; end
                use = insub & ~isnan(xs);
                if nnz(use) < 4, continue; end
                C = struct('set', string(cs.name), 'scale', sc, 'type', string(S.correlation_type), 'n', nnz(use));
                for m = meas
                    X = double(vals.(m{1})(use, :, :));
                    r = nan(nv, nb);
                    for b = 1:nb
                        Y = X(:, :, b);
                        okc = all(~isnan(Y), 1);
                        r(okc, b) = corr(xs(use), Y(:, okc), 'Type', type_name(S.correlation_type));
                    end
                    C.(m{1}) = single(r);
                end
                Vk.correlations = [Vk.correlations, C];
            end
        end
    end
    V = [V, Vk]; %#ok<AGROW>
    file = fullfile(out, sprintf('vertex_values_%s.mat', Lk.inverse));
    save(file, 'Vk', '-v7.3');
    fprintf('Vertex maps (%s): %d recordings x %d sources x %d bands (%s), %d comparison(s), %d correlation(s): %s\n', ...
        Lk.inverse, nr, nv, nb, strjoin(meas, ', '), numel(Vk.groups), numel(Vk.correlations), file);
    if S.figures
        try
            figures(Vk, geo, meas, fullfile(out, 'figures'));
        catch err
            fprintf(2, 'Vertex map figures (%s) failed: %s\n', Lk.inverse, err.message);
        end
    end
end
end

function g = g_columns(A, B)
% Hedges' g of every column (sources x bands) of A against B (recordings
% first), as hedges_g: the pooled SD (df = na + nb - 2), Hedges' exact J
na = sum(~isnan(A), 1);  nb = sum(~isnan(B), 1);
df = na + nb - 2;
sp = sqrt(((na - 1) .* var(A, 0, 1, 'omitnan') + (nb - 1) .* var(B, 0, 1, 'omitnan')) ./ df);
J = exp(gammaln(df / 2) - log(sqrt(df / 2)) - gammaln((df - 1) / 2));
g = J .* (mean(A, 1, 'omitnan') - mean(B, 1, 'omitnan')) ./ sp;
g(na < 2 | nb < 2 | ~(sp > 0)) = NaN;
g = reshape(g, size(A, 2), size(A, 3));
end

function t = type_name(x)
t = 'Pearson';
if strcmpi(x, 'spearman'), t = 'Spearman'; end
end

function figures(Vk, geo, meas, fig_dir)
if ~isfolder(fig_dir), mkdir(fig_dir); end
views = {'left lateral', 'dorsal', 'right lateral'};
names = {'left', 'top (L | R)', 'right'};
seq = parula(256);
div = diverging_map();
bands = Vk.bands;
n_r = numel(Vk.correlations);
for m = meas
    for G = Vk.groups                                 % (no comparisons: one group, all recordings)
        mu = double(G.(m{1}).mean);
        P = arrayfun(@(j) mu(:, :, j), 1:size(mu, 3), 'UniformOutput', false);
        heads = cellstr(compose('%s (n = %d)', G.names(:), G.n(:)))';
        g = double(G.(m{1}).g);
        stem = sprintf('vertex_%s_%s', m{1}, Vk.inverse);
        what = 'mean';
        if ~isempty(g)
            stem = [stem '_' safe(G.name)];
            what = sprintf('group means and Hedges'' g (%s)', G.name);
        end
        draw(geo, P, heads, g, bands, views, names, seq, div, unit_label(m{1}), 'Hedges'' g', ...
            sprintf('%s, %s power of every source: %s; descriptive, not tested', Vk.inverse, m{1}, what), ...
            fullfile(fig_dir, [stem '.png']));
        if ~isempty(g), brains_3d(geo, g, bands, div, sprintf('%s, %s power: Hedges'' g (%s; descriptive)', ...
                Vk.inverse, m{1}, G.name), fullfile(fig_dir, [stem '.fig'])); end
    end
    if n_r > 10
        if strcmp(m{1}, meas{1})
            fprintf('Vertex maps (%s): %d correlations, their maps not drawn (more than 10; in the .mat)\n', Vk.inverse, n_r);
        end
        continue
    end
    for C = Vk.correlations
        rname = 'r';
        if strcmpi(C.type, 'spearman'), rname = 'Spearman r'; end
        draw(geo, {}, {}, double(C.(m{1})), bands, views, names, seq, div, '', rname, ...
            sprintf('%s, %s power of every source: %s with %s (%s, n = %d; descriptive, not tested)', Vk.inverse, ...
            m{1}, rname, C.scale, C.set, C.n), fullfile(fig_dir, sprintf('vertex_%s_%s_r_%s_%s.png', m{1}, ...
            Vk.inverse, safe(C.set), safe(C.scale))));
    end
end
end

function draw(geo, P, heads, E, bands, views, names, seq, div, plabel, elabel, ttl, file)
% one row per band, a panel (three views and its colour bar) per map: the
% maps of P (sources x bands each, one sequential scale per band, shared
% by them) and of E (sources x bands, one diverging scale for all bands)
nv = numel(views);  nb = numel(bands);
np = numel(P);  ne = ~isempty(E);
npan = np + ne;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 npan * (nv * 185 + 70) + 60, 175 * nb + 110]);
tl = tiledlayout(fig, nb, npan, 'TileSpacing', 'compact', 'Padding', 'compact');
if ne
    e = prctile(abs(E(:)), 99);
    if ~(e > 0), e = 1; end
end
for b = 1:nb
    if np
        v = cell2mat(cellfun(@(x) x(:, b), P, 'UniformOutput', false));
        lim = prctile(v(:), [2 98]);
        if ~(lim(2) > lim(1)), lim = lim(1) + [-1 1]; end
    end
    for j = 1:npan
        t = tiledlayout(tl, 1, nv, 'TileSpacing', 'none', 'Padding', 'tight');
        t.Layout.Tile = (b - 1) * npan + j;
        for w = 1:nv
            ax = nexttile(t);
            if j <= np
                cortex_view(ax, geo, struct('vertex', P{j}(:, b)), lim, seq, views{w});
            else
                cortex_view(ax, geo, struct('vertex', E(:, b)), [-e e], div, views{w});
            end
            if b == 1, title(ax, names{w}, 'FontSize', 8, 'FontWeight', 'normal'); end
        end
        cb = colorbar(ax);
        cb.Layout.Tile = 'east';
        cb.FontSize = 7;
        if j <= np, cb.Label.String = plabel; else, cb.Label.String = elabel; end
        if b == 1
            h = elabel;
            if j <= np, h = heads{j}; end
            title(t, h, 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
        end
        if j == 1, ylabel(t, bands(b), 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none'); end
    end
end
title(tl, ttl, 'FontSize', 10, 'Interpreter', 'none');
exportgraphics(fig, file, 'Resolution', 150);
close(fig);
end

function brains_3d(geo, E, bands, cmap, ttl, file)
% E in 3-D, one brain per band, saved as a MATLAB figure that opens visible
% and can be turned with the mouse (as psd_figures' cortex maps)
nb = numel(bands);
nc = min(nb, 4);  nr = ceil(nb / nc);
e = prctile(abs(E(:)), 99);
if ~(e > 0), e = 1; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 330 * nc + 90, 300 * nr + 70], ...
    'Name', ttl, 'NumberTitle', 'off');
tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
for b = 1:nb
    ax = nexttile(tl);
    cortex_view(ax, geo, struct('vertex', E(:, b)), [-e e], cmap, '3d');
    title(ax, bands(b), 'Interpreter', 'none', 'FontSize', 10);
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'Hedges'' g';
title(tl, {ttl, 'Turn a brain with the mouse'}, 'FontSize', 10, 'Interpreter', 'none');
set(fig, 'CreateFcn', ['set(gcbo, ''Visible'', ''on''); r = rotate3d(gcbo); ' ...
    'r.ActionPostCallback = @(~, e) camlight(findobj(e.Axes, ''Type'', ''light''), ''headlight''); ' ...
    'r.Enable = ''on''; clear r']);
savefig(fig, file);
close(fig);
end

function s = unit_label(m)
s = 'log10 power';
if strcmp(m, 'relative'), s = '% of total power'; end
end

function cmap = diverging_map()
% blue - white - red (ColorBrewer RdBu, 11 classes, interpolated), as psd_figures
k = [5 48 97; 33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; 253 219 199; 244 165 130; ...
    214 96 77; 178 24 43; 103 0 31] / 255;
cmap = interp1(linspace(0, 1, size(k, 1)), k, linspace(0, 1, 256));
end

function s = safe(x)
% a name for a file: 'patient - control' -> 'patient_control' (as psd_figures)
s = char(regexprep(strtrim(string(x)), '[^A-Za-z0-9]+', '_'));
end
