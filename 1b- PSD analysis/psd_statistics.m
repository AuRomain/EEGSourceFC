function stats = psd_statistics(S, L, out)
%PSD_STATISTICS Group comparisons and correlations of the PSD measures, per channel or region, with clusters.
%   stats = psd_statistics(S, L, out) (run_psd_analysis; L its levels).
%   Families: one per level x measure x band (psd_<level>_<measure>_<band>;
%   exponent, offset and alpha_frequency have no band), whose states are
%   the channels or regions. Each family is tested by stage 3's
%   permu_ttest_permu_corr_BNS, unchanged: every comparison of
%   S.comparisons (Welch t, or ANCOVA with S.covariates, Freedman-Lane
%   permutations; Hedges' g with a bootstrap CI) and every correlation set
%   (S.correlations, S.covariate_sets; partial r), S.psd.n_perm
%   permutations (S.n_perm when empty), p_maxT (family-wise over the
%   channels or regions of one measure), p_maxT_family and family_p
%   (family-wise over the family: comparison or correlation set) and q_fdr
%   (Benjamini-Hochberg within the family, a description). The recordings are those with a row in the subject
%   table and a group of S.groups (others named in a warning); a subject
%   with two recordings is an error.
%   Over the families of the PSD analysis only (not those of stage 3):
%   family_adjust (max-statistic family_p, families selected by Benjamini-
%   Hochberg, Benjamini-Bogomolov p_bb = p_maxT_family x F / R;
%   significant = p_bb <= S.alpha). The console lines count q_fdr < alpha
%   within each family; the last line the tests significant after the
%   adjustment over families. Before the tests: the sample sizes with the
%   covariates each comparison can use (check_sample_sizes), the
%   permutation resolution (resolution_note: a family can be selected only
%   when 1 / (n_perm + 1) <= alpha / families; the PSD has many families,
%   hence S.psd.n_perm), and, when a covariate is a
%   measure of the spectral profile (e.g. alpha_pw_post), the profile itself
%   (prepare_spectral_covariates: stage 3, which computes it otherwise,
%   runs later).
%   Clusters (S.psd.cluster.run): on the same permutations, clusters of
%   neighbouring channels or regions (psd_adjacency) whose statistic
%   passes the two-sided p < S.psd.cluster.alpha with the same sign, tested
%   by their mass against the largest cluster of the same sign in each
%   permutation, p doubled for the two tails (cluster_perm_test, FieldTrip's
%   cluster test: cluster_p, family-wise over the map); then
%   Benjamini-Hochberg over every cluster of the PSD analysis (cluster_q;
%   cluster_significant = cluster_q <= S.alpha).
%   Writes in out: <level>/<family>_group_comparison.csv, _correlations.csv
%   (.xlsx too; columns Level, Measure, Band, Unit added), summary.csv (one
%   row per family), families.csv (family_adjust), clusters.csv (one row
%   per cluster: its channels or regions, mass, p, q, mean g or r) and
%   neighbours_<level>.csv.
%   stats: res (one entry per family: G, C and what it is), families,
%   summary, clusters.
o = S.psd;
if ~isfolder(out), mkdir(out); end
prepare_spectral_covariates(S);                     % covariates from the spectral profile: computed first
T = read_subject_table(S);
Sx = S;
Sx.figures = false;                                 % the PSD figures are psd_figures'
if isfield(o, 'n_perm') && ~isempty(o.n_perm), Sx.n_perm = o.n_perm; end    % its own number of permutations (many families)
cov = strjoin(S.covariates, ', ');
if isempty(S.covariates), cov = 'none'; end
fprintf('PSD statistics, %d permutations, seed %d; covariates of the group comparisons: %s\n', Sx.n_perm, ...
    S.seed, cov);
list = measure_list(o);
try                                                  % with the covariates each comparison can use
    check_sample_sizes(S, unique(vertcat(L.ids)));
catch err
    warning('psd:sample_size', 'Sample sizes not checked: %s', err.message);
end
nF = family_sizes(S, T, L, size(list, 1));
fprintf('%s\n', resolution_note(Sx.n_perm, S.alpha, nF, "at most: the families are known only once tested"));
res = struct('family', {}, 'level', {}, 'measure', {}, 'band', {}, 'dir', {}, 'G', {}, 'C', {});
for k = 1:numel(L)
    keep = in_groups(S, T, L(k).ids, L(k).name);
    ids = L(k).ids(keep);
    [u, ~, j] = unique(ids);
    dup = u(accumarray(j, 1) > 1);
    if ~isempty(dup)
        error('psd:subjects', '%s: subjects with more than one recording (one per subject): %s', L(k).name, ...
            strjoin(dup, ', '));
    end
    copts = struct();
    if o.cluster.run
        [A, ai] = psd_adjacency(L(k).level, L(k).units, S);
        copts = struct('adjacency', A, 'cluster_alpha', o.cluster.alpha);
        nb = arrayfun(@(i) strjoin(L(k).units(A(i, :)), ' '), (1:numel(L(k).units))', 'UniformOutput', false);
        writetable(table(string(L(k).units(:)), ai.n_neighbours, string(nb), 'VariableNames', ...
            {'unit', 'n_neighbours', 'neighbours'}), fullfile(out, ['neighbours_' L(k).name '.csv']));
        fprintf('%s: %d recordings in the tests; neighbours per %s: median %g (%d-%d)\n', L(k).name, nnz(keep), ...
            unit_word(L(k).level), median(ai.n_neighbours), min(ai.n_neighbours), max(ai.n_neighbours));
        if any(ai.n_neighbours == 0)
            warning('psd:neighbours', '%s: without neighbours (clusters of one): %s', L(k).name, ...
                strjoin(L(k).units(ai.n_neighbours == 0), ', '));
        end
    end
    d = fullfile(out, L(k).name);
    for q = 1:size(list, 1)
        [m, b] = list{q, :};
        V = values(L(k), m, b);
        V = V(keep, :);
        feat = m;
        if b > 0, feat = [m '_' o.bands(b).name]; end
        family = ['psd_' L(k).name '_' feat];
        [nr, K] = size(V);
        F = table(repmat(ids(:), K, 1), repmat(string(feat), nr * K, 1), repelem((1:K)', nr), V(:), ...
            'VariableNames', {'id', 'feature', 'state', 'value'});
        F = F(~isnan(F.value), :);
        if isempty(F)
            fprintf('%-44s no values (band outside the spectra or the fit range)\n', family);
            continue
        end
        t0 = tic;
        [G, C] = permu_ttest_permu_corr_BNS(Sx, F, family, d, copts);
        G = annotate(G, L(k), m, b, o);
        C = annotate(C, L(k), m, b, o);
        res(end + 1) = struct('family', family, 'level', L(k).name, 'measure', m, 'band', band_name(o, b), ...
            'dir', d, 'G', G, 'C', C); %#ok<AGROW>
        fprintf('%-44s group tests %4d (within-family q < %.2f: %d), correlations %5d (within-family q < %.2f: %d)  [%.0f s]\n', family, ...
            n_tests(G), S.alpha, count(G, 'q_fdr', S.alpha), n_tests(C), S.alpha, count(C, 'q_fdr', S.alpha), toc(t0));
    end
end
[res, families] = adjust(res, S.alpha);
[res, clusters] = cluster_fdr(res, S.alpha, o.cluster.run);
summary = table();
for r = 1:numel(res)
    x = res(r);
    row = table(string(x.family), string(x.level), string(x.measure), string(x.band), ...
        n_tests(x.G), count(x.G, 'q_fdr', S.alpha), count(x.G, 'significant'), ...
        n_tests(x.C), count(x.C, 'q_fdr', S.alpha), count(x.C, 'significant'), ...
        'VariableNames', {'family', 'level', 'measure', 'band', 'group_tests', 'group_q_below_alpha', ...
        'group_significant', 'correlations', 'corr_q_below_alpha', 'corr_significant'});
    if o.cluster.run
        in = clusters.family == string(x.family);
        row.clusters = nnz(in);
        row.clusters_significant = nnz(in & clusters.cluster_significant);
    end
    summary = [summary; row]; %#ok<AGROW>
    write_both(x.G, fullfile(x.dir, [x.family '_group_comparison']));
    write_both(x.C, fullfile(x.dir, [x.family '_correlations']));
end
writetable(summary, fullfile(out, 'summary.csv'));
writetable(families, fullfile(out, 'families.csv'));
if o.cluster.run, writetable(clusters, fullfile(out, 'clusters.csv')); end
if height(families)
    fprintf(['PSD: %d families with tests, %d selected (FDR %.2f over the PSD families); %d significant tests ' ...
        'after the selection adjustment'], families.F(1), families.R(1), S.alpha, sum(families.n_significant));
    if o.cluster.run
        fprintf('; %d clusters, %d significant (q <= %.2f)', height(clusters), nnz(clusters.cluster_significant), S.alpha);
    end
    fprintf('\n');
end
stats = struct('res', res, 'families', families, 'summary', summary, 'clusters', clusters);
end


function F = family_sizes(S, T, L, n_measures)
% the number of families, at most: level x measure x band x (comparison or
% correlation set with at least one score: a numeric column of the table
% matching it)
sets = {};
for k = 1:numel(S.correlations), sets(end + 1, :) = {S.correlations(k).columns, S.correlations(k).covariates}; end %#ok<AGROW>
for k = 1:numel(S.covariate_sets)
    cs = S.covariate_sets(k);
    b = S.correlations(strcmp({S.correlations.name}, cs.correlations));
    if isempty(b), continue; end
    col = b.columns;
    if isfield(cs, 'columns') && ~isempty(cs.columns), col = cs.columns; end
    sets(end + 1, :) = {col, cs.covariates}; %#ok<AGROW>
end
names_T = string(T.Properties.VariableNames);
ns = zeros(1, size(sets, 1));
for k = 1:size(sets, 1)
    sc = names_T(~cellfun(@isempty, regexp(cellstr(names_T), sets{k, 1}, 'once')));
    sc = setdiff(sc, [string(sets{k, 2}), string(S.id_column), string(S.group_column)]);
    ns(k) = nnz(arrayfun(@(v) isnumeric(T.(v)), sc));
end
F = numel(L) * n_measures * (numel(S.comparisons) + nnz(ns));
end

function list = measure_list(o)
% {measure, band index (0: none)} of every family, in the order of S.psd.measures
list = cell(0, 2);
for m = o.measures(:)'
    switch m{1}
        case {'absolute', 'relative', 'periodic'}
            for b = 1:numel(o.bands), list(end + 1, :) = {m{1}, b}; end %#ok<AGROW>
        case 'aperiodic'
            list(end + 1, :) = {'exponent', 0}; %#ok<AGROW>
            list(end + 1, :) = {'offset', 0}; %#ok<AGROW>
        case 'alpha_frequency'
            list(end + 1, :) = {'alpha_frequency', 0}; %#ok<AGROW>
    end
end
end

function V = values(Lk, m, b)
% recordings x channels (regions) of one measure (and band)
if b > 0, V = Lk.M.(m)(:, :, b); else, V = Lk.M.(m); end
end

function s = band_name(o, b)
s = '';
if b > 0, s = o.bands(b).name; end
end

function X = annotate(X, Lk, m, b, o)
% the level, measure, band and channel or region of every row
if isempty(X) || height(X) == 0, return; end
X = addvars(X, repmat(string(Lk.name), height(X), 1), repmat(string(m), height(X), 1), ...
    repmat(string(band_name(o, b)), height(X), 1), string(Lk.units(X.state))', ...
    'After', 'Parameter', 'NewVariableNames', {'Level', 'Measure', 'Band', 'Unit'});
end

function keep = in_groups(S, T, ids, what)
% recordings with a row in the subject table and a group of S.groups (as
% run_statistics); the others named in a warning
keep = ismember(string(ids), string(T.(S.id_column)));
if ~isempty(S.group_column)
    [~, row] = ismember(string(ids), string(T.(S.id_column)));
    g = strings(numel(ids), 1);
    g(keep) = string(T.(S.group_column)(row(keep)));
    keep = keep & ismember(lower(strtrim(g)), lower(strtrim(string(S.groups))));
end
if any(~keep)
    warning('psd:no_group', '%s: %s left out (no row in the subject table, or in none of S.groups)', ...
        what, strjoin(unique(string(ids(~keep))), ', '));
end
end

function [res, fams] = adjust(res, alpha)
% Benjamini-Bogomolov over the PSD families (family_adjust), as stage 3
% does over its own: a family is one table's comparison or correlation set
tabs = {};  fam = {};  where = zeros(0, 2);
for r = 1:numel(res)
    for kind = 1:2
        if kind == 1, X = res(r).G; else, X = res(r).C; end
        if isempty(X) || height(X) == 0, continue; end
        if kind == 1, unit = X.Comparison; else, unit = X.CorrelationSet + " (correlations)"; end
        tabs{end + 1} = X; %#ok<AGROW>
        fam{end + 1} = X.Family + " | " + unit; %#ok<AGROW>
        where(end + 1, :) = [r, kind]; %#ok<AGROW>
    end
end
fams = table();
if isempty(tabs), return; end
[tabs, fams] = family_adjust(tabs, fam, alpha);
for i = 1:numel(tabs)
    if where(i, 2) == 1, res(where(i, 1)).G = tabs{i}; else, res(where(i, 1)).C = tabs{i}; end
end
end

function [res, clusters] = cluster_fdr(res, alpha, run)
% one row per cluster of every test; Benjamini-Hochberg over all of them
clusters = table();
if ~run, return; end
rows = {};
for r = 1:numel(res)
    for kind = 1:2
        if kind == 1, X = res(r).G; else, X = res(r).C; end
        if isempty(X) || height(X) == 0 || ~ismember('cluster', X.Properties.VariableNames), continue; end
        if kind == 1
            test = X.Comparison;  scale = strings(height(X), 1);  eff = X.g;  what = "comparison";
        else
            test = X.CorrelationSet;  scale = X.Scale;  eff = X.r;  what = "correlation";
        end
        key = test + "|" + scale;
        for t = unique(key(X.cluster > 0))'
            in0 = key == t;
            for c = unique(X.cluster(in0 & X.cluster > 0))'
                in = in0 & X.cluster == c;
                i = find(in, 1);
                rows(end + 1, :) = {string(res(r).family), string(res(r).level), string(res(r).measure), ...
                    string(res(r).band), what, test(i), scale(i), c, X.cluster_sign(i), nnz(in), ...
                    X.cluster_mass(i), X.cluster_p(i), mean(eff(in), 'omitnan'), strjoin(X.Unit(in), ' '), r, kind}; %#ok<AGROW>
            end
        end
    end
end
names = {'family', 'level', 'measure', 'band', 'kind', 'test', 'scale', 'cluster', 'sign', 'size', ...
    'mass', 'cluster_p', 'mean_effect', 'units', 'res', 'table'};
if isempty(rows)                                     % typed columns, as with clusters
    rows = [repmat({strings(0, 1)}, 1, 7), repmat({zeros(0, 1)}, 1, 6), {strings(0, 1)}, repmat({zeros(0, 1)}, 1, 2)];
    clusters = table(rows{:}, 'VariableNames', names);
    clusters.cluster_q = zeros(0, 1);
else
    clusters = cell2table(rows, 'VariableNames', names);
    clusters.cluster_q = fdr_bh(clusters.cluster_p);
end
clusters.cluster_significant = clusters.cluster_q <= alpha;
for r = 1:numel(res)                                 % the q of each row's cluster in the tables
    for kind = 1:2
        if kind == 1, X = res(r).G; else, X = res(r).C; end
        if isempty(X) || height(X) == 0 || ~ismember('cluster', X.Properties.VariableNames), continue; end
        X.cluster_q = nan(height(X), 1);
        X.cluster_significant = false(height(X), 1);
        if kind == 1, key = X.Comparison + "|"; else, key = X.CorrelationSet + "|" + X.Scale; end
        for i = find(clusters.res == r & clusters.table == kind)'
            in = key == clusters.test(i) + "|" + clusters.scale(i) & X.cluster == clusters.cluster(i);
            X.cluster_q(in) = clusters.cluster_q(i);
            X.cluster_significant(in) = clusters.cluster_significant(i);
        end
        if kind == 1, res(r).G = X; else, res(r).C = X; end
    end
end
clusters = removevars(clusters, {'res', 'table'});
clusters.mean_effect_is = repmat("g", height(clusters), 1);
clusters.mean_effect_is(clusters.kind == "correlation") = "r";
end

function n = n_tests(X)
n = 0;
if ~isempty(X) && height(X), n = nnz(~isnan(X.p_unc)); end
end

function n = count(X, col, alpha)
n = 0;
if isempty(X) || height(X) == 0 || ~ismember(col, X.Properties.VariableNames), return; end
if nargin > 2, n = nnz(X.(col) < alpha); else, n = nnz(X.(col)); end
end

function w = unit_word(level)
w = 'channel';
if strcmp(level, 'source'), w = 'region'; end
end

function write_both(T, stem)
if isempty(T) || height(T) == 0, return; end
writetable(T, [stem '.csv']);
try
    writetable(T, [stem '.xlsx']);
catch
end
end
