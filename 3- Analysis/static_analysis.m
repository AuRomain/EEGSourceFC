function res = static_analysis(S, inv, band, adir)
%STATIC_ANALYSIS The static connectivity of one inverse method and band: features, tests, NBS, reliability, figures.
%   res = static_analysis(S, inv, band, adir) with adir the analysis folder
%   (analysis_dir: Analysis/ of the brain-network-state run, or
%   results/Statistics/<run> in a static-only run). For every measure of
%   S.static_measures (tag <inv>_<band>_<measure>) and every matrix of
%   S.static_matrices ('full': the measure over the whole recording;
%   'window_mean': the mean of its windows, the original "static" matrix):
%   - the feature families of S.static_families (static_features):
%       global    Strength, ClustCoeff, CPL (graph metrics; family
%                 static_GT_<matrix>, the original static_GT_metrics_stats_permu)
%       networks  mean connectivity inside each network and between each
%                 pair (S.network_table; family static_networks_<matrix>)
%       regions   node strength of every region (static_regions_<matrix>)
%       power     log10 band power, overall and per network (static_power;
%                 once per band and inverse method, in the folder of the
%                 first measure and matrix: it depends on neither)
%     each with every group comparison (S.comparisons, covariates
%     S.covariates: Welch t, or ANCOVA with the Welch-type G and
%     Freedman-Lane permutations) and every correlation set
%     (S.correlations, S.covariate_sets), by permu_ttest_permu_corr_BNS, in
%     Statistical analysis/<tag>/ (with features_<matrix>.csv);
%   - the edges (compute_NBS_analysis): the network-based statistic of every
%     comparison, the network enrichment and every edge's statistics, in
%     NBS/<tag>/;
%   - the split-half reliability of every feature and of the edges
%     (static_reliability, S.reliability), for each matrix:
%     Statistical analysis/<tag>/static_reliability.csv (whole recording)
%     and static_reliability_window_mean.csv (the window mean on each half;
%     reliability.csv is that of the states, state_reliability);
%   - figures (S.figures): the group mean matrices
%     (group_matrices_<matrix>.png).
%   A warning when the groups of a comparison differ in the seconds
%   analysed (connectivity estimates depend on the amount of data).
%   res: one element per family tested (tag, family, dir, G, C, seconds;
%   G and C empty and failed true when the family failed), for the
%   adjustment over families of run_statistics.
res = struct('tag', {}, 'family', {}, 'dir', {}, 'G', {}, 'C', {}, 'seconds', {}, 'failed', {});
Ttab = read_subject_table(S);
for st = 1:numel(S.static_matrices)
    which = S.static_matrices{st};
    D = load_static(S, inv, band, which);
    [ok, row] = ismember(D.ids, Ttab.(S.id_column));
    if ~all(ok), error('static_analysis:table', 'No subject table row for: %s', strjoin(D.ids(~ok), ', ')); end
    T = Ttab(row, :);
    nodes = load_network_table(S, D.roi_labels);
    Gatlas = atlas_regions(S, D.roi_labels, nodes);
    if st == 1
        fprintf('%s %s: %d recordings (%s), %.0f-%.0f s analysed\n', inv, band, numel(D.ids), ...
            strjoin(arrayfun(@(g) sprintf('%s %d', S.groups{g}, nnz(D.group == g)), 1:numel(S.groups), ...
            'UniformOutput', false), ', '), min(D.analysed_s), max(D.analysed_s));
        duration_check(S, D, T);
    end
    for m = 1:numel(S.static_measures)
        k = find(strcmp(D.measures, S.static_measures{m}));
        tag = sprintf('%s_%s_%s', inv, band, D.measures{k});
        fprintf('\n--- %s, %s ---\n', tag, which);
        tag_dir = fullfile(adir, 'Statistical analysis', tag);
        if ~isfolder(tag_dir), mkdir(tag_dir); end
        fams = S.static_families;
        if st > 1 || m > 1, fams = setdiff(fams, {'power'}, 'stable'); end   % power: once (one family in the adjustment)
        F = static_features(D, k, nodes, fams);
        for f = fams
            Ff = F(F.family == f{1}, {'id', 'feature', 'state', 'value'});
            name = family_name(f{1}, which);
            if isempty(Ff), fprintf('%s: no features (no network table?)\n', name); continue; end
            t0 = tic;
            try
                [G, C] = permu_ttest_permu_corr_BNS(S, Ff, name, tag_dir);
                res(end + 1) = struct('tag', string(tag), 'family', string(name), 'dir', tag_dir, 'G', G, 'C', C, ...
                    'seconds', toc(t0), 'failed', false); %#ok<AGROW>
            catch err                             % one family failing does not stop the others
                fprintf(2, '%s failed: %s\n', name, err.message);
                res(end + 1) = struct('tag', string(tag), 'family', string(name), 'dir', tag_dir, 'G', [], 'C', [], ...
                    'seconds', toc(t0), 'failed', true); %#ok<AGROW>
            end
        end
        writetable(F, fullfile(tag_dir, ['features_' which '.csv']));
        if S.reliability && ~any(isfinite(D.halves(:, k, :, :, :)), 'all')
            fprintf('Split-half reliability of the %s matrix: not available (stage 1 outputs made before it existed)\n', which);
        elseif S.reliability
            R = static_reliability(D, k, nodes, fams);
            rel_file = 'static_reliability.csv';
            if ~strcmp(which, 'full'), rel_file = ['static_reliability_' which '.csv']; end
            writetable(R, fullfile(tag_dir, rel_file));
            fprintf('Split-half reliability, %s (blocks of %g s): edges median r_sb %.2f, identification %.2f\n', ...
                which, S.reliability_block, R.r_sb(R.feature == "edge r_sb, median"), R.r(R.feature == "identification"));
        end
        try
            Rn = compute_NBS_analysis(S, D, k, T, nodes, Gatlas, fullfile(adir, 'NBS', tag), which);
            for i = find(Rn.p_both < S.alpha)'
                fprintf('NBS %s, %s: component of %d edges, %d regions, p = %.4f (both directions %.4f)\n', ...
                    Rn.comparison(i), Rn.contrast(i), Rn.n_edges(i), Rn.n_nodes(i), Rn.p(i), Rn.p_both(i));
            end
            if ~any(Rn.p_both < S.alpha), fprintf('NBS: no component with p (both directions) < %.2f\n', S.alpha); end
        catch err
            fprintf(2, 'NBS (%s, %s) failed: %s\n', tag, which, err.message);
        end
        if S.figures
            try
                group_matrices(S, D, k, T, nodes, fullfile(tag_dir, ['group_matrices_' which '.png']));
            catch err
                fprintf(2, 'group matrices figure failed: %s\n', err.message);
            end
        end
    end
end
end

function name = family_name(family, which)
% static_GT_full (the original static_GT_metrics_stats_permu), static_networks_full, ...; static_power
switch family
    case 'global', name = ['static_GT_' which];
    case 'power',  name = 'static_power';
    otherwise,     name = ['static_' family '_' which];
end
end

function duration_check(S, D, T)
% connectivity estimates depend on the amount of data (wPLI's bias, the
% spread of every estimate): a warning when the groups of a comparison
% differ in the seconds analysed
for c = 1:numel(S.comparisons)
    a = D.analysed_s(select_subjects(T, S.comparisons(c).A));
    b = D.analysed_s(select_subjects(T, S.comparisons(c).B));
    if numel(a) < 2 || numel(b) < 2, continue; end
    p = ranksum(a, b);
    if p < 0.05 || abs(median(a) - median(b)) > 0.1 * median([a; b])
        warning('static_analysis:duration', ['%s: the groups differ in the duration analysed (median %.0f vs ' ...
            '%.0f s, rank-sum p = %.3f). Connectivity estimates depend on it: set S.core_duration to the ' ...
            'shortest recording, or add the duration as a covariate'], S.comparisons(c).name, median(a), median(b), p);
    end
end
end
