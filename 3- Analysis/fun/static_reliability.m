function R = static_reliability(D, k, nodes, families, rec_groups)
%STATIC_RELIABILITY Split-half reliability of the static features and edges of one connectivity measure.
%   R = static_reliability(D, k, nodes, families) with D from load_static,
%   k the measure, nodes and families as static_features; rec_groups
%   (optional: each recording's group, in the order of D.ids) adds
%   r_within_groups.
%   Stage 1 computed
%   every whole-recording matrix also on the odd and on the even blocks of
%   each recording (S.reliability_block s, split_halves): two estimates
%   from interleaved halves of the same recording. Across recordings:
%     features  for every feature of static_features: r, the Pearson
%               correlation of the two halves, and r_sb = 2r / (1 + r)
%               (Spearman-Brown), the reliability of the value of the whole
%               recording, the one that is tested; for r <= 0 r_sb is r
%               unchanged (the correction assumes a positive reliability;
%               it would give values below -1)
%     edges     every edge's r_sb: median and quartiles over the edges
%     pattern   within each recording, the correlation of the edges of its
%               two halves (median over recordings: how reproducible the
%               pattern is), and identification: the share of recordings
%               whose odd-block pattern correlates best with their own
%               even-block pattern among all recordings' (Finn et al. 2015;
%               chance 1 / recordings)
%   A group difference or correlation of a feature with low reliability
%   (r_sb below about 0.5) is attenuated: read its results with care.
%   A reliability is computed only from at least max(3, half of the
%   recordings) recordings with both halves; with fewer, r and r_sb are
%   NaN, the note gives n of N, and a warning says so (so that a handful
%   of the recordings cannot stand for all of them). Interleaved blocks of one session give
%   within-session consistency, an upper bound of test-retest reliability,
%   and the groups are pooled (a group difference raises r):
%   r_within_groups is r after removing each group's mean from both
%   halves (groups of at least 2 recordings). A stable
%   bias is reliable too: the window mean of a phase measure, biased by
%   its short windows, can be reliable without any coupling, so a high
%   reliability does not show that a matrix measures coupling, and is not
%   a reason to choose one matrix or to keep a feature (docs/static_mode.md).
%   R: table family, feature, n, r, r_sb, r_within_groups, note.
N = numel(D.ids);
n_min = max(3, ceil(N / 2));                       % recordings with both halves needed
few = sprintf("%%d of %d recordings with both halves, fewer than %d: not computed", N, n_min);
F1 = static_features(D, k, nodes, families, 1);
F2 = static_features(D, k, nodes, families, 2);    % the same rows, in the same order
grp = strings(height(F1), 0);                      % each row's recording's group
if nargin > 4 && ~isempty(rec_groups)
    [~, at] = ismember(string(F1.id), string(D.ids));
    rg = string(rec_groups(:));
    grp = rg(at);
end
[groups, fam, feat] = findgroups(F1.family, F1.feature);
ng = max([0; groups(:)]);                          % 0: no feature family (S.static_families empty)
rows = cell(ng, 7);
order = zeros(ng, 1);
for g = 1:ng
    m = groups == g;
    a = F1.value(m);  b = F2.value(m);
    ok = isfinite(a) & isfinite(b);
    r = NaN;  rw = NaN;  note = "";
    if nnz(ok) >= n_min
        r = corr(a(ok), b(ok));
        if ~isempty(grp)
            gm = grp(m);
            rw = within_r(a(ok), b(ok), gm(ok));
        end
    else
        note = sprintf(few, nnz(ok));
    end
    rows(g, :) = {fam(g), feat(g), nnz(ok), r, sb(r), rw, note};
    order(g) = find(m, 1);
end
[~, s] = sort(order);                              % the order of static_features
rows = rows(s, :);
n = numel(D.roi_labels);
nr = numel(D.ids);
low = find(tril(true(n), -1));
E1 = reshape(D.halves(1, k, :, :, :), n * n, nr);  E1 = E1(low, :);
E2 = reshape(D.halves(2, k, :, :, :), n * n, nr);  E2 = E2(low, :);
okr = all(isfinite(E1), 1) & all(isfinite(E2), 1);  % recordings with both halves (a half can be NaN: too short, measure not computed)
E1 = E1(:, okr);  E2 = E2(:, okr);  nr = nnz(okr);
q = nan(1, 3);  within = NaN;  ident = NaN;
notes = ["median over the edges of each edge's Spearman-Brown r", "", "", ...
    "median over recordings of the correlation of the edges of the two halves", ...
    sprintf("share of recordings identified from their other half (chance %.3f)", 1 / nr)];
if nr >= n_min
    re = zeros(numel(low), 1);
    for e = 1:numel(low), re(e) = corr(E1(e, :)', E2(e, :)'); end
    q = quantile(sb(re), [0.25, 0.5, 0.75]);
    C = corr(E1, E2);                              % recordings x recordings: odd of i with even of j
    [~, best] = max(C, [], 2);
    within = median(diag(C));  ident = mean(best == (1:nr)');
else
    notes(:) = sprintf(few, nr);
    which = 'full';  if isfield(D, 'matrix'), which = D.matrix; end
    warning('static_reliability:few', ['%s, %s matrix: %d of %d recordings have both halves (at least %d ' ...
        'needed): no split-half reliability'], D.measures{k}, which, nr, N, n_min);
end
rows(end + 1, :) = {"edges", "edge r_sb, median", nr, NaN, q(2), NaN, notes(1)};
rows(end + 1, :) = {"edges", "edge r_sb, lower quartile", nr, NaN, q(1), NaN, notes(2)};
rows(end + 1, :) = {"edges", "edge r_sb, upper quartile", nr, NaN, q(3), NaN, notes(3)};
rows(end + 1, :) = {"pattern", "within-recording r of the halves", nr, within, NaN, NaN, notes(4)};
rows(end + 1, :) = {"pattern", "identification", nr, ident, NaN, NaN, notes(5)};
R = cell2table(rows, 'VariableNames', {'family', 'feature', 'n', 'r', 'r_sb', 'r_within_groups', 'note'});
end

function r = within_r(a, b, g)
% r of the halves after removing each group's mean (groups of at least 2
% recordings; recordings without a group left out)
r = NaN;
g = string(g(:));  a = a(:);  b = b(:);
use = ~ismissing(g) & g ~= "";
[~, ~, j] = unique(g(use));
n = accumarray(j, 1);
keep = false(size(g));
keep(use) = n(j) >= 2;
if nnz(keep) < 3, return; end
[~, ~, j] = unique(g(keep));
a = a(keep);  b = b(keep);
ma = accumarray(j, a, [], @mean);  mb = accumarray(j, b, [], @mean);
a = a - ma(j);  b = b - mb(j);
if std(a) > 0 && std(b) > 0, r = corr(a, b); end
end

function s = sb(r)
% Spearman-Brown for r > 0; r <= 0 unchanged
s = r;
s(r > 0) = 2 * r(r > 0) ./ (1 + r(r > 0));
end
