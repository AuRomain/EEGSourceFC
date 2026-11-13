function C = psd_region_check(S, L, out)
%PSD_REGION_CHECK The regions' power from their signal against the mean of their sources' power (region_power_check.csv).
%   C = psd_region_check(S, L, out) (run_psd_analysis) for every level of
%   the regions (source_<inverse>, or source_<inverse>_vertexmean with
%   S.psd.region_spectra = 'vertex_mean'), from both spectra of each
%   region (psd_spectra): 'signal', the PSD of the region's signal as
%   stage 1 makes it (the sign-flipped mean of its sources with
%   S.scout_function 'mean'), and 'vertex_mean', the mean of its sources'
%   PSDs. In the signal, sources of opposite or shifted phase cancel; in
%   the mean of the sources' power they do not, so the ratio of the two
%   band powers is at most 1 with the 'mean' scout function (Cauchy-Schwarz)
%   and is 1 when the region's sources are coherent. A region whose
%   group difference or correlation depends on which of the two is used
%   holds sources that do not act together.
%   One row per inverse x region x band (x comparison of S.comparisons):
%     ratio_median, ratio_q1, ratio_q3   the band power of the signal over
%       the mean of the sources' band power, over the recordings
%     r_absolute, r_relative, r_periodic   Pearson's r over the recordings
%       between the two (log10 absolute power, relative power, power above
%       the 1/f fit: psd_measures)
%     g_<measure>_signal, g_<measure>_vertex_mean   Hedges' g of the
%       comparison (hedges_g, without covariates) with each of the two:
%       descriptive, not tested
%   Written to <out>/region_power_check.csv; with S.figures, the median
%   ratio on the cortex, figures/region_power_check_<inverse>.png. The log
%   gives each inverse's median ratio per band and its five lowest regions.
C = table();
k_src = find(strcmp({L.level}, 'source') & ~arrayfun(@(x) isempty(x.check), L));
if isempty(k_src), return; end
o = S.psd;
T = [];
try
    T = read_subject_table(S);
catch
end
bands = string({o.bands.name});
rows = {};
for k = k_src
    Lk = L(k);
    sig = Lk.check.signal;  vm = Lk.check.vertex_mean;
    ratio = 10 .^ (sig.absolute - vm.absolute);      % recordings x regions x bands
    Gr = psd_groups(S, T, Lk);
    if isempty(S.comparisons) || isempty(T), Gr = Gr([]); end
    K = numel(Lk.units);
    for b = 1:numel(bands)
        q = quantile(ratio(:, :, b), [0.25 0.5 0.75], 1);   % 3 x regions
        r = nan(K, 3);
        meas = {'absolute', 'relative', 'periodic'};
        for j = 1:3
            for u = 1:K
                r(u, j) = corr_ok(sig.(meas{j})(:, u, b), vm.(meas{j})(:, u, b));
            end
        end
        base = table(repmat(string(Lk.inverse), K, 1), string(Lk.units(:)), repmat(bands(b), K, 1), ...
            q(2, :)', q(1, :)', q(3, :)', r(:, 1), r(:, 2), r(:, 3), 'VariableNames', {'inverse', 'region', ...
            'band', 'ratio_median', 'ratio_q1', 'ratio_q3', 'r_absolute', 'r_relative', 'r_periodic'});
        if isempty(Gr)
            rows{end + 1} = base; %#ok<AGROW>
            continue
        end
        for c = 1:numel(Gr)
            X = base;
            X.comparison = repmat(Gr(c).name, K, 1);
            for j = 1:3
                for kind = ["signal", "vertex_mean"]
                    V = Lk.check.(kind).(meas{j})(:, :, b);
                    g = nan(K, 1);
                    for u = 1:K
                        a = V(Gr(c).in(:, 1), u);  z = V(Gr(c).in(:, 2), u);
                        a = a(~isnan(a));  z = z(~isnan(z));
                        if numel(a) >= 2 && numel(z) >= 2, g(u) = hedges_g(a, z); end
                    end
                    X.("g_" + meas{j} + "_" + kind) = g;
                end
            end
            rows{end + 1} = X; %#ok<AGROW>
        end
    end
    med = squeeze(median(median(ratio, 1, 'omitnan'), 2, 'omitnan'));
    reg = median(median(ratio, 1, 'omitnan'), 3, 'omitnan');
    [~, lo] = sort(reg(:));
    lo = lo(1:min(5, K));
    fprintf(['Region power check (%s): band power of the region signal / mean of its sources'' power, ' ...
        'median over recordings and regions: %s; lowest regions: %s\n'], Lk.inverse, ...
        strjoin(compose('%s %.2f', bands(:), med(:)), ', '), ...
        strjoin(compose('%s %.2f', string(Lk.units(lo))', reg(lo)'), ', '));
    if S.figures
        try
            check_figure(S, Lk, squeeze(median(ratio, 1, 'omitnan')), bands, fullfile(out, 'figures'));
        catch err
            fprintf(2, 'Region power check figure (%s) failed: %s\n', Lk.inverse, err.message);
        end
    end
end
C = vertcat(rows{:});
writetable(C, fullfile(out, 'region_power_check.csv'));
end

function r = corr_ok(x, y)
ok = ~isnan(x) & ~isnan(y);
r = NaN;
if nnz(ok) >= 3 && std(x(ok)) > 0 && std(y(ok)) > 0, r = corr(x(ok), y(ok)); end
end

function check_figure(S, Lk, R, bands, fig_dir)
% the median ratio (regions x bands) on the cortex, one row per band
if ~isfolder(fig_dir), mkdir(fig_dir); end
geo = cortex_geometry(S, Lk);
R = reshape(R, numel(Lk.units), []);
views = {'left lateral', 'left medial', 'right medial', 'right lateral', 'dorsal', 'ventral'};
names = [views(1:4), {'dorsal (L | R)', 'ventral (R | L)'}];
nv = numel(views);  nb = numel(bands);
cmap = parula(256);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 nv * 215 + 80, 190 * nb + 60]);
tl = tiledlayout(fig, nb, nv, 'TileSpacing', 'tight', 'Padding', 'compact');
for b = 1:nb
    vtx = nan(size(geo.V, 1), 1);
    in = geo.region > 0;
    vtx(in) = R(geo.region(in), b);
    for v = 1:nv
        ax = nexttile(tl);
        cortex_view(ax, geo, struct('vertex', vtx), [0 1], cmap, views{v});
        if v == 1
            text(ax, 0, 1.02, bands(b), 'Units', 'normalized', 'FontSize', 9, 'FontWeight', 'bold', ...
                'Interpreter', 'none', 'VerticalAlignment', 'bottom');
        end
        if b == 1, title(ax, names{v}, 'FontSize', 8, 'FontWeight', 'normal'); end
    end
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'signal / mean of the sources'' power';
title(tl, {sprintf('%s: band power of each region''s signal over the mean of its sources'' power (median over recordings)', ...
    Lk.inverse), '1: the region''s sources act together; low: they cancel in the region''s signal'}, ...
    'FontSize', 10, 'Interpreter', 'none');
exportgraphics(fig, fullfile(fig_dir, sprintf('region_power_check_%s.png', Lk.inverse)), 'Resolution', 150);
close(fig);
end
