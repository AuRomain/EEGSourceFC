function R = static_reliability(D, k, nodes, families)
%STATIC_RELIABILITY Split-half reliability of the static features and edges of one connectivity measure.
%   R = static_reliability(D, k, nodes, families) with D from load_static,
%   k the measure, nodes and families as static_features. Stage 1 computed
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
%   R: table family, feature, n, r, r_sb, note.
F1 = static_features(D, k, nodes, families, 1);
F2 = static_features(D, k, nodes, families, 2);    % the same rows, in the same order
[groups, fam, feat] = findgroups(F1.family, F1.feature);
rows = cell(max(groups), 6);
order = zeros(max(groups), 1);
for g = 1:max(groups)
    m = groups == g;
    a = F1.value(m);  b = F2.value(m);
    ok = isfinite(a) & isfinite(b);
    r = NaN;
    if nnz(ok) >= 3, r = corr(a(ok), b(ok)); end
    rows(g, :) = {fam(g), feat(g), nnz(ok), r, sb(r), ""};
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
re = zeros(numel(low), 1);
for e = 1:numel(low), re(e) = corr(E1(e, :)', E2(e, :)'); end
q = quantile(sb(re), [0.25, 0.5, 0.75]);
C = corr(E1, E2);                                  % recordings x recordings: odd of i with even of j
[~, best] = max(C, [], 2);
rows(end + 1, :) = {"edges", "edge r_sb, median", nr, NaN, q(2), "median over the edges of each edge's Spearman-Brown r"};
rows(end + 1, :) = {"edges", "edge r_sb, lower quartile", nr, NaN, q(1), ""};
rows(end + 1, :) = {"edges", "edge r_sb, upper quartile", nr, NaN, q(3), ""};
rows(end + 1, :) = {"pattern", "within-recording r of the halves", nr, median(diag(C)), NaN, ...
    "median over recordings of the correlation of the edges of the two halves"};
rows(end + 1, :) = {"pattern", "identification", nr, mean(best == (1:nr)'), NaN, ...
    sprintf("share of recordings identified from their other half (chance %.3f)", 1 / nr)};
R = cell2table(rows, 'VariableNames', {'family', 'feature', 'n', 'r', 'r_sb', 'note'});
end

function s = sb(r)
% Spearman-Brown for r > 0; r <= 0 unchanged
s = r;
s(r > 0) = 2 * r(r > 0) ./ (1 + r(r > 0));
end
