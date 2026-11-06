function [tabs, fams] = family_adjust(tabs, fam, alpha)
%FAMILY_ADJUST Tests adjusted for the selection of families (Benjamini-Bogomolov), max-statistic within families.
%   [tabs, fams] = family_adjust(tabs, fam, alpha)
%   tabs   cell of result tables, each with p_unc, q_fdr, and (as
%          permu_ttest_permu_corr_BNS writes them) p_maxT_family and family_p
%   fam    cell of string arrays, one per table: the family of every row
%          (names unique across tables: band and inverse, feature family,
%          comparison or correlation set)
%   alpha  level
%   Each family is tested as a whole: family_p is the max-statistic
%   (Westfall-Young) p of its tests, computed by permu_ttest_permu_corr_BNS
%   on the permutations all the family's tests share: the share of
%   permutations whose largest |z| over the family reaches the data's
%   largest. The F families with a family_p are selected by Benjamini-
%   Hochberg at alpha over their family_p (family_q, family_selected; R
%   selected). In a selected family a test is significant when its
%   family-wise p (p_maxT_family: its |z| against the family's largest in
%   each permutation) is at most alpha * R / F:
%   p_bb = min(1, p_maxT_family * F / R), NaN in the families not selected,
%   significant = p_bb <= alpha. This keeps the family-wise error averaged
%   over the selected families at alpha (Benjamini & Bogomolov 2014, JRSS-B
%   76:303, with family-wise error within the families), whereas the tests
%   of each family read alone, across the many families of a run, give
%   false discoveries in most null runs (audit V3-17, M4). Family and test
%   agree: a selected family has at least one significant test. Both
%   p-values are at least 1 / (n_perm + 1) whatever the size of the family,
%   so one strong test can make its family selected and itself significant
%   when 1 / (n_perm + 1) <= alpha / F (resolution_note). q_fdr, the
%   Benjamini-Hochberg q within the family, is kept as a description; it
%   decides nothing.
%   A family whose tables lack p_maxT_family (tables made otherwise) is
%   tested as before 5 November 2026: Simes p of its p_unc, and q_bb =
%   min(1, q_fdr * F / R), significant = q_bb <= alpha (FDR within the
%   family). A family put together from tables of different family_p
%   (not done by the pipeline) gets Bonferroni over its parts: family_p and
%   p_bb times their number.
%   fams   one row per family: family, n_tests, family_p, family_test
%          ("max-statistic" or "Simes"), family_q, family_selected, n_q_fdr
%          (q_fdr <= alpha), n_significant, F, R
names = strings(0, 1);  p = zeros(0, 1);  pf = zeros(0, 1);  pw = zeros(0, 1);
for i = 1:numel(tabs)
    names = [names; string(fam{i}(:))]; %#ok<AGROW>
    p = [p; tabs{i}.p_unc(:)]; %#ok<AGROW>
    if all(ismember({'family_p', 'p_maxT_family'}, tabs{i}.Properties.VariableNames))
        pf = [pf; tabs{i}.family_p(:)]; %#ok<AGROW>
        pw = [pw; tabs{i}.p_maxT_family(:)]; %#ok<AGROW>
    else
        pf = [pf; -ones(height(tabs{i}), 1)]; %#ok<AGROW>   % -1: not tested by the maximum
        pw = [pw; nan(height(tabs{i}), 1)]; %#ok<AGROW>
    end
end
[u, ~, idx] = unique(names, 'stable');
fp = nan(numel(u), 1);  how = strings(numel(u), 1);  parts = ones(numel(u), 1);
for f = 1:numel(u)
    [fp(f), how(f), parts(f)] = family_test(p(idx == f), pf(idx == f));
end
fq = fdr_bh(fp);
F = nnz(~isnan(fp));
sel = fq <= alpha;
R = nnz(sel);
n_tests = accumarray(idx, ~isnan(p), [numel(u), 1]);
n_q = zeros(numel(u), 1);  n_sig = n_q;
for i = 1:numel(tabs)
    [~, fi] = ismember(string(fam{i}(:)), u);
    T = tabs{i};
    T.family_p = fp(fi);
    T.family_q = fq(fi);
    T.family_selected = sel(fi);
    s = sel(fi);
    if all(how(fi) ~= "Simes")
        T.p_bb = nan(height(T), 1);
        T.p_bb(s) = min(1, T.p_maxT_family(s) .* parts(fi(s)) * F / R);
        T.significant = T.p_bb <= alpha;
    else
        T.q_bb = nan(height(T), 1);
        T.q_bb(s) = min(1, T.q_fdr(s) * F / R);
        T.significant = T.q_bb <= alpha;
    end
    n_q = n_q + accumarray(fi, T.q_fdr <= alpha, [numel(u), 1]);
    n_sig = n_sig + accumarray(fi, T.significant, [numel(u), 1]);
    tabs{i} = T;
end
fams = table(u, n_tests, fp, how, fq, sel, n_q, n_sig, repmat(F, numel(u), 1), repmat(R, numel(u), 1), ...
    'VariableNames', {'family', 'n_tests', 'family_p', 'family_test', 'family_q', 'family_selected', 'n_q_fdr', ...
                      'n_significant', 'F', 'R'});
end

function [fp, how, parts] = family_test(p, pf)
% the p of a family as a whole: the max-statistic p of its tests (column
% family_p, the same in every row of a table) when every row has it, with
% Bonferroni over the parts when a family has several; otherwise Simes
parts = 1;
if all(pf >= 0 | isnan(pf))
    v = unique(pf(~isnan(pf)));
    how = "max-statistic";
    if isempty(v), fp = NaN; return; end
    parts = numel(v);
    fp = min(1, parts * min(v));
    return
end
fp = simes(p);
how = "Simes";
end

function ps = simes(p)
% Simes combination: min over i of m p_(i) / i (NaN when no p-value)
p = sort(p(~isnan(p)));
m = numel(p);
if m == 0, ps = NaN; return; end
ps = min(1, min(m * p(:) ./ (1:m)'));
end
