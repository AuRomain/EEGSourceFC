function [R, E] = compute_NBS_analysis(S, D, k, T, nodes, Gatlas, out, which_static)
%COMPUTE_NBS_ANALYSIS The edges of one static measure, per group comparison: NBS, network enrichment, edge table.
%   [R, E] = compute_NBS_analysis(S, D, k, T, nodes, Gatlas, out, which_static)
%   with D from load_static (D.conn: the whole recording, or the mean of the
%   windows), k the measure, T the subject table rows of D's recordings,
%   nodes and Gatlas the network table rows and atlas_regions of the
%   regions (may be empty), out the output folder (NBS/<tag> of the
%   analysis folder) and which_static 'full' or 'window_mean' (the suffix
%   of the file names). Called by static_analysis. For each comparison of
%   S.comparisons, on the recordings of its two groups with the covariates
%   S.covariates (those usable for them, usable_covariates), without the
%   recordings for which the measure could not be computed (NaN edges, a
%   warning names them):
%   - network-based statistic (Zalesky et al. 2010) with the NBS 1.2
%     toolbox's NBSstats:
%     design [1, A, covariates] (ANCOVA with S.covariates), Welch-type G of
%     every edge (welch_G: each group's residual variance on its own;
%     Welch's t without covariates; NBSglm's pooled-variance t is liberal
%     with unequal groups) as z on its own Aspin-Welch degrees of
%     freedom, in the data and in every permutation (t_to_z, as PALM),
%     with Freedman-Lane permutations (as
%     NBSglm), primary threshold S.nbs.thresh (3.1, on z), component size = number
%     of edges ('extent') or summed statistic above threshold
%     ('intensity'), both contrasts A > B and B > A, S.nbs.perms
%     permutations with a fixed seed; every component reported with
%     its p = (1 + permutations whose largest component is at least as
%     large) / (perms + 1) (NBS 1.2 reports b / perms, which can be 0),
%     and p_both = min(1, 2p): each direction is tested at its
%     own alpha, so the error over the two is up to 2 alpha; p_both is the
%     one to compare with alpha;
%   - network enrichment (network_enrichment): which networks and pairs of
%     networks concentrate the edge statistics, against random and spin
%     (rotated) network labels;
%   - every edge's group means, difference, G and its z.
%   Writes in out: NBS_<which>.csv (components: comparison, contrast,
%   component, n_edges, n_nodes, p, p_both, regions), NBS_<which>.mat,
%   enrichment_<which>.csv and edges_<comparison>_<which>.csv; figures of
%   the edge statistics and of the components with p_both < S.alpha
%   (S.figures: a connectogram and a brain view).
%   R: the NBS table; E: the enrichment table.
if nargin < 8 || isempty(which_static), which_static = 'full'; end
n = numel(D.roi_labels);
nr = numel(D.ids);
W = reshape(D.conn(k, :, :, :), n, n, nr);
ind_upper = find(triu(true(n), 1));                % edge order of NBSstats
y = reshape(W, n * n, nr)';
y = y(:, ind_upper);
not_computed = any(isnan(y), 2);                   % the measure not computed for these recordings: left out
if any(not_computed)
    warning('compute_NBS_analysis:not_computed', ['NBS of %s: %d recording(s) without this measure (NaN edges) ' ...
        'left out: %s'], D.measures{k}, nnz(not_computed), strjoin(string(D.ids(not_computed)), ', '));
end
[Z, znames] = covariate_matrix(T, S.covariates, S.categorical_covariates);
rows = {};
E = table();
nbs = struct('comparison', {}, 'contrast', {}, 'components', {}, 'p', {}, 'n', {}, 'covariates', {});
if ~isfolder(out), mkdir(out); end
for c = 1:numel(S.comparisons)
    cmp = S.comparisons(c);
    inA = select_subjects(T, cmp.A);
    inB = select_subjects(T, cmp.B);
    inAB = inA | inB;
    Zk = nan(nr, 0);
    zused = {};
    if ~isempty(Z)
        [Zc, zused] = usable_covariates(Z(inAB, :), znames, ['NBS ' cmp.name], double(inA(inAB)));
        Zk = nan(nr, size(Zc, 2));  Zk(inAB, :) = Zc;
    end
    sel = inAB & all(~isnan(Zk), 2) & ~not_computed;
    if nnz(inA & sel) < 2 || nnz(inB & sel) < 2
        warning('compute_NBS_analysis:n', 'NBS %s skipped: fewer than 2 recordings in a group with the covariates', cmp.name);
        continue
    end
    X = [ones(nnz(sel), 1), double(inA(sel)), Zk(sel, :)];
    names = split(string(cmp.name), " - ");
    if numel(names) ~= 2, names = ["A"; "B"]; end
    contrasts = {sprintf('%s > %s', names(1), names(2)), 1; sprintf('%s > %s', names(2), names(1)), -1};
    cvec = [0; 1; zeros(size(X, 2) - 2, 1)];
    G = welch_perm_stats(X, y(sel, :), cvec, inA(sel), S.nbs.perms, S.seed);
    for d = 1:2
        nbs_size = 'Extent';
        if strcmpi(S.nbs.size, 'intensity'), nbs_size = 'Intensity'; end   % NBSstats: exact, case-sensitive
        STATS = struct('thresh', S.nbs.thresh, 'alpha', 1, 'N', n, 'test_stat', G * contrasts{d, 2}, 'size', nbs_size);
        [~, n_cnt, con_mat, pval] = evalc('NBSstats(STATS)');   % quiet: NBSstats prints a line per permutation
        pval = (round(pval * S.nbs.perms) + 1) / (S.nbs.perms + 1);   % b / K -> (b + 1) / (K + 1)
        for m = 1:n_cnt
            A = con_mat{m} | con_mat{m}';
            reg = find(any(A, 1));
            rows(end + 1, :) = {string(cmp.name), string(contrasts{d, 1}), m, nnz(triu(A, 1)), numel(reg), ...
                pval(m), min(1, 2 * pval(m)), strjoin(string(D.roi_labels(reg)), '; ')}; %#ok<AGROW>
        end
        nbs(end + 1) = struct('comparison', cmp.name, 'contrast', contrasts{d, 1}, 'components', {con_mat}, ...
            'p', pval, 'n', [nnz(inA & sel), nnz(inB & sel)], 'covariates', {zused}); %#ok<AGROW>
    end
    % where the difference lies: networks against random and spin labels
    [P, stat, Gm] = network_enrichment(W(:, :, sel), X, inA(sel), nodes, Gatlas, ...
        struct('n_perm', S.enrichment.n_perm, 'seed', S.seed));
    if height(P)
        P = addvars(P, repmat(string(cmp.name), height(P), 1), 'Before', 1, 'NewVariableNames', 'comparison');
        E = [E; P]; %#ok<AGROW>
    end
    % every edge
    [ii, jj] = find(tril(true(n), -1));
    lin = sub2ind([n n], ii, jj);
    mA = mean(W(:, :, inA & sel), 3);  mB = mean(W(:, :, inB & sel), 3);
    lbl = string(D.roi_labels(:));
    Te = table(lbl(ii), lbl(jj), mA(lin), mB(lin), mA(lin) - mB(lin), Gm(lin), stat(lin), ...
        'VariableNames', {'region_a', 'region_b', 'mean_A', 'mean_B', 'difference', 'G', 'z'});
    Te = sortrows(Te, 'z', 'descend', 'ComparisonMethod', 'abs');
    writetable(Te, fullfile(out, ['edges_' file_name(cmp.name) '_' which_static '.csv']));
    if S.figures
        edge_figures(S, D, stat, nodes, Gatlas, nbs(end - 1:end), rows, cmp.name, contrasts, out, which_static);
    end
end
vars = {'comparison', 'contrast', 'component', 'n_edges', 'n_nodes', 'p', 'p_both', 'regions'};
if isempty(rows)
    R = cell2table(cell(0, numel(vars)), 'VariableNames', vars);
else
    R = cell2table(rows, 'VariableNames', vars);
end
writetable(R, fullfile(out, ['NBS_' which_static '.csv']));
if height(E), writetable(E, fullfile(out, ['enrichment_' which_static '.csv'])); end
ids = D.ids; roi_labels = D.roi_labels; settings_nbs = S.nbs; covariates = znames; measure = D.measures{k}; %#ok<NASGU>
save(fullfile(out, ['NBS_' which_static '.mat']), 'nbs', 'ids', 'roi_labels', 'settings_nbs', 'covariates', 'measure');
end

function T = welch_perm_stats(X, Y, c, grp, perms, seed)
% (perms + 1) x edges, as NBSglm's test_stat: the Welch-type G of contrast c
% (welch_G) as z on its own Aspin-Welch degrees of freedom (t_to_z, as
% PALM) for the data (first row) and for perms Freedman-Lane
% permutations (residuals of the model without the contrasted column)
X0 = X(:, c == 0);
fit0 = X0 * (X0 \ Y);
r0 = Y - fit0;
rs = RandStream('mrg32k3a', 'Seed', seed);
T = zeros(perms + 1, size(Y, 2));
[g, ~, df] = welch_G(Y, X, c, grp);
T(1, :) = t_to_z(g, df);
for i = 1:perms
    [g, ~, df] = welch_G(fit0 + r0(randperm(rs, size(Y, 1)), :), X, c, grp);
    T(i + 1, :) = t_to_z(g, df);
end
T(isnan(T)) = 0;                                 % as NBSglm: an undefined edge is not suprathreshold
end

function edge_figures(S, D, stat, nodes, Gatlas, nbs, rows, cname, contrasts, out, which_static)
% the edge statistics as a matrix ordered by network; the components with
% p_both < alpha as a connectogram and a brain view
n = numel(D.roi_labels);
[lab, nets] = network_assignment(nodes);
order = 1:n;
if ~isempty(nets)
    [~, g] = ismember(lab, [nets; "Other"]);
    [~, order] = sort(g);
end
stem = fullfile(out, [file_name(cname) '_' which_static]);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 640 560]);
ax = axes(fig);
lim = max(abs(stat(:)));
if lim == 0, lim = 1; end
imagesc(ax, stat(order, order), [-lim lim]);
axis(ax, 'square');
colormap(ax, diverging_map());
cb = colorbar(ax);
cb.Label.String = 'z of the Welch-type G, A - B';
title(ax, sprintf('%s: edge statistics', cname), 'Interpreter', 'none');
if ~isempty(nets)
    hold(ax, 'on');
    lg = lab(order);
    for a = 1:numel(nets)
        m = find(lg == nets(a));
        if isempty(m), continue; end
        rectangle(ax, 'Position', [m(1) - 0.5, m(1) - 0.5, numel(m), numel(m)], 'EdgeColor', [0.2 0.2 0.2]);
    end
    centre = arrayfun(@(a) mean(find(lg == a)), nets);
    ok = ~isnan(centre);
    set(ax, 'XTick', centre(ok), 'XTickLabel', nets(ok), 'YTick', centre(ok), 'YTickLabel', nets(ok));
else
    set(ax, 'XTick', [], 'YTick', []);
end
exportgraphics(fig, [stem '_edges.png']);
close(fig);
if isempty(nodes) || isempty(rows), return; end
R = cell2table(rows, 'VariableNames', {'comparison', 'contrast', 'component', 'n_edges', 'n_nodes', 'p', 'p_both', 'regions'});
for d = 1:numel(nbs)
    hits = find(R.comparison == string(cname) & R.contrast == string(contrasts{d, 1}) & R.p_both < S.alpha);
    for h = hits'
        A = double(nbs(d).components{R.component(h)});
        A = (A | A') .* stat;
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 420]);
        tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
        circular_graph(nexttile(tl), A, nodes, char(network_column(nodes)), 'signed');
        if ~isempty(Gatlas.vertices) && all(ismember({'x_mni', 'y_mni', 'z_mni'}, nodes.Properties.VariableNames))
            glass_brain(nexttile(tl), Gatlas, [nodes.x_mni, nodes.y_mni, nodes.z_mni], A, ...
                network_colours(lab), 'top', 'signed');
        end
        title(tl, sprintf('%s, %s: component %d, %d edges, p = %.4f (both directions %.4f)', cname, ...
            contrasts{d, 1}, R.component(h), R.n_edges(h), R.p(h), R.p_both(h)), 'Interpreter', 'none');
        exportgraphics(fig, sprintf('%s_NBS_%d_component_%d.png', stem, d, R.component(h)));
        close(fig);
    end
end
end

function c = network_column(nodes)
[~, ~, c] = network_assignment(nodes);
end

function m = diverging_map()
% blue - white - red, 256 steps
t = linspace(0, 1, 128)';
blue = [0.13 0.40 0.67];  red = [0.70 0.09 0.17];
m = [blue + t * (1 - blue); 1 - t * (1 - red)];
end

function f = file_name(name)
% a comparison's name in a file name: 'patient - control' -> 'patient_control'
f = regexprep(char(name), '[^A-Za-z0-9]+', '_');
f = regexprep(f, '^_|_$', '');
end
