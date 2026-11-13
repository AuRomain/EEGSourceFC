function group_matrices(S, D, k, T, nodes, file)
%GROUP_MATRICES Figure: the mean connectivity matrix of every group, regions ordered by network.
%   group_matrices(S, D, k, T, nodes, file) for measure k of D (load_static),
%   the groups of S.groups (D.group), one panel each on a common colour
%   scale, and for each comparison of S.comparisons the difference of its
%   group means (A - B). Networks (network_assignment) are outlined.
%   Recordings whose matrix of this measure is NaN (stage 1 could not
%   compute it, e.g. the symmetric orthogonalisation with fewer
%   independent signals than regions, or no window of the measure's
%   length in any piece of the recording) are left out, as in the tests:
%   n in each title counts the recordings drawn. Without any, no figure,
%   and a warning.
n = numel(D.roi_labels);
nr = numel(D.ids);
W = reshape(D.conn(k, :, :, :), n, n, nr);
low = tril(true(n), -1);
Wl = reshape(W, n * n, nr);
has = any(isfinite(Wl(low(:), :)), 1)';                      % recordings with this measure
[lab, nets] = network_assignment(nodes);
order = 1:n;
if ~isempty(nets)
    [~, g] = ismember(lab, [nets; "Other"]);
    [~, order] = sort(g);
end
panels = {};  titles = {};  kinds = [];
for g = 1:numel(S.groups)
    in = D.group(:) == g & has;
    if ~any(in), continue; end
    panels{end + 1} = mean(W(:, :, in), 3, 'omitnan'); %#ok<AGROW>
    titles{end + 1} = sprintf('%s (n = %d)', S.groups{g}, nnz(in)); %#ok<AGROW>
    kinds(end + 1) = 1; %#ok<AGROW>
end
for c = 1:numel(S.comparisons)
    a = select_subjects(T, S.comparisons(c).A);
    b = select_subjects(T, S.comparisons(c).B);
    a = a(:) & has;  b = b(:) & has;
    if nnz(a) < 2 || nnz(b) < 2, continue; end
    panels{end + 1} = mean(W(:, :, a), 3, 'omitnan') - mean(W(:, :, b), 3, 'omitnan'); %#ok<AGROW>
    titles{end + 1} = S.comparisons(c).name; %#ok<AGROW>
    kinds(end + 1) = 2; %#ok<AGROW>
end
v = cellfun(@(P) P(low), panels(kinds == 1), 'UniformOutput', false);
v = vertcat(v{:});
lim1 = [min(v, [], 'omitnan'), max(v, [], 'omitnan')];
if ~any(kinds == 1) || any(~isfinite(lim1))
    warning('group_matrices:no_values', 'Group matrices of %s: no recording of the groups has this measure; no figure', ...
        D.measures{k});
    return
end
if diff(lim1) == 0, lim1 = lim1 + [-1 1]; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 400 * numel(panels) 400]);
tl = tiledlayout(fig, 1, numel(panels), 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:numel(panels)
    ax = nexttile(tl);
    P = panels{p}(order, order);
    P(1:n + 1:end) = NaN;
    if kinds(p) == 1
        h = imagesc(ax, P, lim1);
        colormap(ax, parula(256));
    else
        lim = max(abs(P(:)), [], 'omitnan');
        if ~(lim > 0), lim = 1; end
        h = imagesc(ax, P, [-lim lim]);
        colormap(ax, diverging_map());
    end
    set(h, 'AlphaData', ~isnan(P));
    axis(ax, 'square');
    colorbar(ax);
    title(ax, titles{p}, 'Interpreter', 'none');
    if ~isempty(nets)
        hold(ax, 'on');
        lg = lab(order);
        for a = 1:numel(nets)
            m = find(lg == nets(a));
            if isempty(m), continue; end
            rectangle(ax, 'Position', [m(1) - 0.5, m(1) - 0.5, numel(m), numel(m)], 'EdgeColor', [0.2 0.2 0.2]);
        end
        centre = arrayfun(@(x) mean(find(lg == x)), nets);
        ok = ~isnan(centre);
        set(ax, 'XTick', centre(ok), 'XTickLabel', nets(ok), 'YTick', centre(ok), 'YTickLabel', nets(ok), 'FontSize', 7);
    else
        set(ax, 'XTick', [], 'YTick', []);
    end
end
hz = unique(D.band_hz, 'rows');
band = sprintf('%.4g-%.4g Hz', hz(1, 1), hz(1, 2));
if size(hz, 1) > 1, band = 'individual bands'; end
title(tl, sprintf('%s, %s', D.measures{k}, band), 'Interpreter', 'none');
exportgraphics(fig, file);
close(fig);
end

function m = diverging_map()
% blue - white - red, 256 steps
t = linspace(0, 1, 128)';
blue = [0.13 0.40 0.67];  red = [0.70 0.09 0.17];
m = [blue + t * (1 - blue); 1 - t * (1 - red)];
end
