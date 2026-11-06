function run_statistics(S)
%RUN_STATISTICS  Stage 3: statistics and network figures.
%   run_statistics() or run(...) uses pipeline_settings(); run_statistics(S)
%   uses the settings S. A variable S in the workspace is not used (V3-5).
% The original folder had no main script: 21 scripts, one per feature
% family and covariate choice, each with its own copy of the code and
% hard-coded paths (O13). This main script runs, for every band and
% inverse method of the settings:
%   the five original analyses (functions of the same names):
%     permu_ttest_permu_corr_microfeatures_multcomp_BNS  temporal features of the states
%     permu_ttest_permu_corr_GT_multcomp_BNS             graph metrics of the states
%     permu_ttest_permu_corr_language_GT_multcomp_BNS    language networks of the states
%       (each for the raw and the smoothed labels, S.label_sets)
%     static_GT_metrics_stats_permu                      graph metrics of the static matrix
%       (now the family 'global' of static_analysis, below)
%     static_language_GT_metrics_stats_permu             language networks of the static matrix
%       (each for the whole-recording and the window-mean matrix, S.static_matrices)
%   the state-free features (state_free_features, S.state_free): how the
%     connectivity of the windows varies without states (variability,
%     autocorrelation, gradient scores of node strength, components)
%   the window measures (window_measures, S.window_measures): band power,
%     frequency, muscle and eye indices of every window, per recording, and
%     how much the states follow them
%   the variability against surrogates (surrogate_variability,
%     S.surrogate_variability): the window-to-window SD and CV of mean
%     connectivity against stationary surrogates with the same spectra
%   each with every group comparison (S.comparisons, covariates
%   S.covariates) and every correlation set (S.correlations,
%   S.covariate_sets: these replace the resid_* copies);
%   the static analyses (static_analysis): for every measure of
%     S.static_measures and matrix of S.static_matrices, the feature
%     families of S.static_families (global: Strength, ClustCoeff, CPL;
%     networks; regions; power), the network-based statistic with network
%     enrichment and edge tables (compute_NBS_analysis), the split-half
%     reliability (S.reliability) and the group mean matrices;
%   then the network scripts: compute_export_nodes_edges,
%   node_strength_RSN, state_network_contrast (the networks each state
%   strengthens or weakens, tested over recordings with a spin test),
%   plot_state_networks (connectograms and brain views of the states),
%   compute_language_network_strength; the split-half reliability of the
%   recording-level measures (state_reliability, S.reliability); and checks
%   of the duration differences by ranks, without outliers and beyond
%   coverage (duration_checks, S.duration_checks; not a test family);
%   with individual bands (S.individual_bands), also the individual alpha
%   frequency itself (group comparisons, correlations; Statistical analysis/IAF);
%   first of all, the spectral profile of every recording
%   (spectral_profile_features, S.spectral_profile: alpha frequency and
%   peak, 1/f background, muscle and eye indices at the scalp), joined to
%   the subject table, so that its measures can be covariates of every
%   family, and tested as the family 'spectral' (Statistical
%   analysis/spectral_profile; a measure that is a covariate is tested
%   adjusted for the other covariates, not for itself);
%   before the families, the number of recordings of every comparison and
%   correlation set against what its test needs, with the covariates each
%   can use (check_sample_sizes: a warning names those that cannot be
%   tested, whose rows would be NaN or absent; the counts of the summary
%   are the tests with a p-value), and what the number of permutations
%   allows (resolution_note).
% Multiple testing: a family is one feature family x one comparison or
% correlation set, for one band and inverse method. Every test of a family
% runs on the same permutations, so the family's largest statistic can be
% compared with the largest of each permutation (max-statistic, Westfall-
% Young): family_p, the family as a whole, and p_maxT_family, each test
% against it. The families of the whole run are selected by Benjamini-
% Hochberg over their family_p, and a test of a selected family is
% significant when p_bb = p_maxT_family x F / R <= alpha (Benjamini-
% Bogomolov: family_q, family_selected, p_bb, significant; family_adjust,
% audit M4). q_fdr, the Benjamini-Hochberg q within the family (the console
% lines "within-family q < alpha" count these), describes; it does not
% decide. The tables are rewritten with these columns at the end of the
% stage; families.csv lists every family with F and R; the figures of tests
% that are not significant after this adjustment are removed (those drawn
% with q_fdr < alpha), and the last lines give the number of significant
% tests and the permutation resolution (resolution_note).
% Outputs in results/Connectivity/BNS analysis/<run>/Analysis/; the log,
% the settings and a summary of every family in Statistical analysis/.
% Static only (S.dynamic = false: no windows, no stage 2): the spectral
% profile, the IAF and the static analyses, in results/Statistics/<run>/
% (analysis_dir).
% Settings from pipeline_settings.m, or the argument S.

stage_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(stage_dir));
if nargin < 1, S = []; end
S = setup_toolboxes(stage_settings(S, 'run_statistics'));
adir = analysis_dir(S);
out = fullfile(adir, 'Statistical analysis');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out, 'analysis_log.txt'));
stage_name = 'Stage 3';
if ~S.dynamic, stage_name = 'Statistics (static only)'; end
fprintf('%s, %s, %d permutations, seed %d\nSubject table: %s\nCovariates of the group comparisons: %s\n', ...
    stage_name, char(datetime('now')), S.n_perm, S.seed, S.subs_info, strjoin([{''}, S.covariates], ' '));
t_start = tic;
fprintf('%s\n', resolution_note(S.n_perm, S.alpha));
do_language = ~isempty(S.language_file) && isfile(S.language_file);
summary = {};
done = struct('tag', {}, 'dir', {}, 'G', {}, 'C', {});   % every family's tables, for the adjustment over families
% the spectral profile (S.spectral_profile): before every other family,
% which can use its measures as covariates (read_subject_table)
P = [];
if S.spectral_profile
    fprintf('\n=== spectral profile ===\n');
    try
        [P, tested] = spectral_profile_features(S);
    catch err
        fprintf(2, 'Spectral profile failed: %s\n', err.message);
    end
end
try                                             % with the spectral measures, which can be covariates
    check_sample_sizes(S);
catch err
    warning('run_statistics:sample_size', 'Sample sizes not checked: %s', err.message);
end
if ~isempty(P)                                  % a measure that is a covariate: tested adjusted for the
    try                                         % others, not for itself (permu_ttest_permu_corr_BNS)
        fprintf('\n=== spectral profile: tests ===\n');
        F = table(repmat(P.ID, numel(tested), 1), repelem(string(tested(:)), height(P), 1), ...
            zeros(height(P) * numel(tested), 1), reshape(P{:, tested}, [], 1), ...
            'VariableNames', {'id', 'feature', 'state', 'value'});
        F = in_groups(S, F(~isnan(F.value), :), 'spectral profile');
        [G, C] = permu_ttest_permu_corr_BNS(S, F, 'spectral', fullfile(out, 'spectral_profile'));
        done(end + 1) = struct('tag', "spectral", 'dir', fullfile(out, 'spectral_profile'), 'G', G, 'C', C);
        nG = 0; nC = 0;
        if height(G), nG = nnz(G.q_fdr < S.alpha); end
        if height(C), nC = nnz(C.q_fdr < S.alpha); end
        summary(end + 1, :) = {"spectral", "spectral", n_tests(G), nG, n_tests(C), nC};
        fprintf('%-26s group tests %4d (within-family q < %.2f: %d), correlations %5d (within-family q < %.2f: %d)\n', ...
            'spectral', n_tests(G), S.alpha, nG, n_tests(C), S.alpha, nC);
    catch err
        fprintf(2, 'Spectral profile tests failed: %s\n', err.message);
    end
end
% the individual alpha frequency (only with individual bands): group
% comparisons and correlations, as a feature without states
iaf_file = fullfile(S.results, 'Connectivity', 'iaf.csv');
if S.iaf_used && isfile(iaf_file)
    o = detectImportOptions(iaf_file);
    o = setvartype(o, intersect({'recording', 'ID', 'IAF_source', 'IAF_flag'}, o.VariableNames), 'string');
    I = readtable(iaf_file, o);
    I = I(~isnan(I.IAF), :);
    F = table(I.ID, repmat("IAF", height(I), 1), zeros(height(I), 1), I.IAF, ...
        'VariableNames', {'id', 'feature', 'state', 'value'});
    F = in_groups(S, F, 'IAF');
    fprintf('\n=== individual alpha frequency (%d recordings) ===\n', height(I));
    try
        [G, C] = permu_ttest_permu_corr_BNS(S, F, 'IAF', fullfile(out, 'IAF'));
        done(end + 1) = struct('tag', "IAF", 'dir', fullfile(out, 'IAF'), 'G', G, 'C', C);
        nG = 0; nC = 0;
        if height(G), nG = nnz(G.q_fdr < S.alpha); end
        if height(C), nC = nnz(C.q_fdr < S.alpha); end
        summary(end + 1, :) = {"IAF", "IAF", n_tests(G), nG, n_tests(C), nC};
        fprintf('%-26s group tests %4d (within-family q < %.2f: %d), correlations %5d (within-family q < %.2f: %d)\n', ...
            'IAF', n_tests(G), S.alpha, nG, n_tests(C), S.alpha, nC);
    catch err
        fprintf(2, 'IAF statistics failed: %s\n', err.message);
    end
end
for b = 1:numel(S.bands)
    for iv = 1:numel(S.inv)
        band = S.bands(b).name;
        inv = S.inv{iv};
        tag = sprintf('%s_%s_%s', inv, band, S.conn{S.co});
        if S.dynamic, fprintf('\n=== %s ===\n', tag); else, fprintf('\n=== %s %s ===\n', inv, band); end
        runs = {};
        if S.dynamic                                % the states and the windows (stage 2)
            lang_files = dir(fullfile(bns_run_dir(S), 'Results', ['language_GT_metrics_*_' tag '.mat']));
            for ls = S.label_sets
                runs(end + 1, :) = {['microfeatures_' ls{1}], @() permu_ttest_permu_corr_microfeatures_multcomp_BNS(S, tag, ls{1})}; %#ok<AGROW>
                runs(end + 1, :) = {['GT_' ls{1}], @() permu_ttest_permu_corr_GT_multcomp_BNS(S, tag, ls{1})}; %#ok<AGROW>
                if do_language && ~isempty(lang_files)          % stage 2 made language features for this atlas
                    runs(end + 1, :) = {['language_' ls{1}], @() permu_ttest_permu_corr_language_GT_multcomp_BNS(S, tag, ls{1})}; %#ok<AGROW>
                end
            end
            for st = S.static_matrices
                if do_language && ~isempty(lang_files)
                    runs(end + 1, :) = {['static_language_' st{1}], @() static_language_GT_metrics_stats_permu(S, inv, band, st{1})}; %#ok<AGROW>
                end
            end
            if S.state_free                         % window-to-window variation without states
                runs(end + 1, :) = {'state_free', @() state_free_features(S, inv, band, tag)}; %#ok<AGROW>
            end
            if S.window_measures                    % power, frequency and artefacts of the windows
                runs(end + 1, :) = {'window_measures', @() window_measures(S, inv, band, tag)}; %#ok<AGROW>
            end
            if S.surrogate_variability              % variability beyond what the spectra imply
                runs(end + 1, :) = {'surrogate_variability', @() surrogate_variability(S, inv, band, tag)}; %#ok<AGROW>
            end
        end
        for k = 1:size(runs, 1)
            t0 = tic;
            try
                [G, C] = runs{k, 2}();
            catch err                              % one family failing does not stop the others
                fprintf(2, '%s failed: %s\n', runs{k, 1}, err.message);
                summary(end + 1, :) = {string(tag), string(runs{k, 1}), NaN, NaN, NaN, NaN}; %#ok<AGROW>
                continue
            end
            done(end + 1) = struct('tag', string(tag), 'dir', fullfile(out, tag), 'G', G, 'C', C); %#ok<AGROW>
            nG = 0; nC = 0; nGt = n_tests(G); nCt = n_tests(C);
            if height(G), nG = nnz(G.q_fdr < S.alpha); end
            if height(C), nC = nnz(C.q_fdr < S.alpha); end
            summary(end + 1, :) = {string(tag), string(runs{k, 1}), nGt, nG, nCt, nC}; %#ok<AGROW>
            fprintf('%-26s group tests %4d (within-family q < %.2f: %d), correlations %5d (within-family q < %.2f: %d)  [%.0f s]\n', runs{k, 1}, ...
                nGt, S.alpha, nG, nCt, S.alpha, nC, toc(t0));
        end
        % the static analyses: features, NBS, reliability, figures of every
        % measure of S.static_measures
        t0 = tic;
        try
            res = static_analysis(S, inv, band, adir);
        catch err
            fprintf(2, 'static analyses failed: %s\n', err.message);
            res = struct('tag', {}, 'family', {}, 'dir', {}, 'G', {}, 'C', {}, 'seconds', {}, 'failed', {});
            summary(end + 1, :) = {string(tag), "static", NaN, NaN, NaN, NaN}; %#ok<AGROW>
        end
        for r = res
            if r.failed
                summary(end + 1, :) = {r.tag, r.family, NaN, NaN, NaN, NaN}; %#ok<AGROW>
                continue
            end
            done(end + 1) = struct('tag', r.tag, 'dir', r.dir, 'G', r.G, 'C', r.C); %#ok<AGROW>
            nG = 0; nC = 0; nGt = n_tests(r.G); nCt = n_tests(r.C);
            if height(r.G), nG = nnz(r.G.q_fdr < S.alpha); end
            if height(r.C), nC = nnz(r.C.q_fdr < S.alpha); end
            summary(end + 1, :) = {r.tag, r.family, nGt, nG, nCt, nC}; %#ok<AGROW>
            fprintf('%-26s %-24s group tests %4d (within-family q < %.2f: %d), correlations %5d (within-family q < %.2f: %d)\n', r.tag, r.family, ...
                nGt, S.alpha, nG, nCt, S.alpha, nC);
        end
        fprintf('static analyses: %.0f s\n', toc(t0));
        if ~S.dynamic, continue; end
        if S.duration_checks                        % robustness of the duration differences (not a family)
            for ls = S.label_sets
                try
                    duration_checks(S, tag, ls{1});
                catch err
                    fprintf(2, 'duration checks (%s) failed: %s\n', ls{1}, err.message);
                end
            end
        end
        steps = {@() compute_export_nodes_edges(S, tag), @() node_strength_RSN(S, tag), ...
                 @() state_network_contrast(S, inv, band, tag), @() plot_state_networks(S, tag), ...
                 @() compute_language_network_strength(S, tag)};
        if S.reliability, steps{end + 1} = @() state_reliability(S, inv, band, tag); end %#ok<AGROW>
        for k = 1:numel(steps)
            try
                steps{k}();
            catch err
                fprintf(2, '%s failed: %s\n', func2str(steps{k}), err.message);
            end
        end
    end
end
summary = cell2table(summary, 'VariableNames', {'tag', 'family', 'group_tests', 'group_q_below_alpha', ...
    'correlations', 'corr_q_below_alpha'});
[summary, families] = adjust_over_families(done, summary, S.alpha, S.figures);
writetable(summary, fullfile(out, 'summary.csv'));
writetable(families, fullfile(out, 'families.csv'));
if height(families)
    fprintf(['\n%d families with tests, %d selected (FDR %.2f over families); %d significant tests after the ' ...
        'selection adjustment (group comparisons %d, correlations %d; within-family FDR q <= %.2f, not a decision: %d)\n'], ...
        families.F(1), families.R(1), S.alpha, sum(families.n_significant), sum(summary.group_significant), ...
        sum(summary.corr_significant), S.alpha, sum(families.n_q_fdr));
    fprintf('%s\n', resolution_note(S.n_perm, S.alpha, families.F(1)));
end
provenance = struct('stage', 3, 'date', char(datetime('now')), 'settings', S, 'summary', summary, 'families', families);
save(fullfile(out, 'provenance_stage3.mat'), 'provenance');
fprintf('%s finished in %.1f min: %s\n', stage_name, toc(t_start) / 60, fileparts(out));
diary off
end

function [summary, families] = adjust_over_families(done, summary, alpha, figures)
% Benjamini-Bogomolov over all families of the run (family_adjust); every
% table rewritten with the new columns; the summary gets the number of
% significant tests of each feature family; the figures of the tests that
% are not significant after the adjustment are removed (drawn for
% p_maxT_family <= alpha)
tabs = {};  fam = {};  where = zeros(0, 3);
for d = 1:numel(done)
    for kind = 1:2
        if kind == 1, T = done(d).G; else, T = done(d).C; end
        if isempty(T) || height(T) == 0, continue; end
        if kind == 1, unit = T.Comparison; else, unit = T.CorrelationSet + " (correlations)"; end
        tabs{end + 1} = T; %#ok<AGROW>
        fam{end + 1} = done(d).tag + " | " + T.Family + " | " + unit; %#ok<AGROW>
        where(end + 1, :) = [d, kind, 0]; %#ok<AGROW>
    end
end
summary.group_significant = zeros(height(summary), 1);
summary.corr_significant = zeros(height(summary), 1);
if isempty(tabs)
    families = table();
    return
end
[tabs, families] = family_adjust(tabs, fam, alpha);
kept = 0;  removed = 0;
for i = 1:numel(tabs)
    d = done(where(i, 1));
    T = tabs{i};
    if where(i, 2) == 1, suffix = '_group_comparison'; col = 'group_significant'; else, suffix = '_correlations'; col = 'corr_significant'; end
    write_both(T, fullfile(d.dir, [char(T.Family(1)) suffix]));
    row = summary.tag == d.tag & summary.family == T.Family(1);
    summary.(col)(row) = nnz(T.significant);
    if figures                                     % the figures drawn for q_fdr < alpha: the significant kept
        [k, r] = prune_family_figures(T, d.dir, alpha);
        kept = kept + k;  removed = removed + r;
    end
end
if figures && kept + removed > 0
    fprintf(['Figures: %d of the tests significant after the adjustment over families kept, %d removed ' ...
        '(family-wise p <= %.2f only)\n'], kept, removed, alpha);
end
end

function write_both(T, stem)
writetable(T, [stem '.csv']);
try
    writetable(T, [stem '.xlsx']);
catch
end
end

function F = in_groups(S, F, what)
% the recordings of S.data with a row in the subject table and a group of
% S.groups, as in the other families (stage1_outputs); the others are left
% out with a warning (before, one recording without a row stopped the family)
if isempty(S.subs_info) || isempty(S.group_column), return, end
T = read_subject_table(S);
[ok, row] = ismember(string(F.id), string(T.(S.id_column)));
g = strings(height(F), 1);
g(ok) = string(T.(S.group_column)(row(ok)));
keep = ok & ismember(lower(strtrim(g)), lower(strtrim(string(S.groups))));
if any(~keep)
    warning('run_statistics:no_group', '%s: %s left out (no row in the subject table, or in none of S.groups)', ...
        what, strjoin(unique(string(F.id(~keep))), ', '));
end
F = F(keep, :);
end

function n = n_tests(X)
% tests of a family's table that have a p-value (a test its subjects could
% not support is a row of NaN, never a p-value: perm_ancova_maxT)
n = height(X);
if n && ismember('p_unc', X.Properties.VariableNames), n = nnz(~isnan(X.p_unc)); end
end
