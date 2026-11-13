function [G, C] = permu_ttest_permu_corr_BNS(S, F, family, out_dir, opts)
%PERMU_TTEST_PERMU_CORR_BNS Group comparisons and clinical correlations of one feature family.
%   [G, C] = permu_ttest_permu_corr_BNS(S, F, family, out_dir)
%   [G, C] = permu_ttest_permu_corr_BNS(S, F, family, out_dir, opts) with
%           opts.adjacency (states x states, true for neighbours: channels
%           or regions as states, the PSD analysis) and opts.cluster_alpha
%           also tests clusters of neighbouring states (cluster_perm_test,
%           on the permutations of each test): G and C get the columns
%           cluster, cluster_sign, cluster_size, cluster_mass, cluster_p
%   F       long table (id, feature, state, value) built by one of the five
%           analysis functions (permu_ttest_permu_corr_*_multcomp_BNS,
%           static_*_stats_permu); state 0 = a feature without states.
%           Optional columns covariate and covariate_name: a covariate of
%           that feature only (e.g. its value in stationary surrogates,
%           surrogate_variability), the same for every state of a subject,
%           added to the covariates of its group comparisons and
%           correlations; covariate_name "" = none
%   family  name used in the output files, e.g. 'microfeatures_smoothed'
%   The statistics:
%   - subjects joined to the table by ID; groups and subgroups
%     defined on the table (S.comparisons), not by fixed indices;
%   - group comparison of every feature and state (perm_ancova_maxT):
%     Welch t without covariates, ANCOVA with the Welch-type G statistic
%     (each group's residual variance on its own) and
%     Freedman-Lane permutations with covariates (S.covariates);
%     each statistic converted to z with its own Welch-Satterthwaite or
%     Aspin-Welch degrees of freedom in the data and in every permutation,
%     as PALM does (t_to_z), and the max-statistic across the states of a
%     feature. The covariates are in the model, not regressed out within
%     each group, which would remove the group difference itself.
%     PERMUTOOLS (permuttest2, permucorr, booteffectsize) is not used: with
%     missing values ('rows','all') its permuted statistics keep the counts
%     and means of the unpermuted data (in a null simulation half of the
%     correlations had p < 0.05), and a constant feature gets the smallest
%     possible p; with 'rows','complete' a subject missing one state is
%     dropped from every state. The tests here drop missing values per
%     state and give NaN when a statistic is undefined;
%   - covariates: text columns (e.g. sex), and numeric columns named in
%     S.categorical_covariates (e.g. a site coded 1, 2, 3), coded as
%     indicators; within each comparison or correlation set, a covariate
%     that is missing for every subject, constant, collinear, or (in a
%     comparison) of a single value in one group that the other group
%     (almost) never has is left out, with a warning (usable_covariates),
%     and the covariates actually used are written in every row; a
%     covariate is not used for its own measure (a feature named as the
%     covariate, e.g. alpha_pw_post of the family 'spectral' when it is in
%     S.covariates; IAF, IAF_cov and the IAF_PAF or IAF_CoG that the IAF
%     is, S.iaf.measure, count as one measure), with a message: adjusted
%     for itself, the model fits perfectly and the test means nothing;
%   - each group's mean and SD of the raw values (mean_A, sd_A, mean_B,
%     sd_B, on the subjects of the test; diff = mean_A - mean_B without
%     covariates, the adjusted difference with them);
%   - the effect size of the difference tested, with a bootstrap CI, on
%     the subjects of the test: g, Hedges' g without covariates, the
%     adjusted g with them (the adjusted difference over the pooled SD of
%     the raw values, Hedges' J; perm_ancova_maxT); g_unadjusted, Hedges'
%     g of the raw values, and its CI;
%   - correlations with clinical scores (S.correlations: columns chosen by
%     regular expression, subjects by condition, e.g. patients only):
%     Pearson or Spearman correlation (S.correlation_type), or partial
%     correlation (covariates removed from feature and score), tested with
%     the slope over its HC3 standard error and permutations (with
%     covariates, of the score residualised on them: Smith;
%     perm_partial_corr_maxT), or a wild bootstrap (S.correlation_test),
%     StatisticalTest naming the type and the test; zeros treated as
%     missing only where asked; S.covariate_sets re-run a correlation set
%     with more covariates;
%   - features constant across all subjects are not tested; every other
%     test is written, with its N; the family (this feature
%     family and one comparison or correlation set) tested as a whole on
%     the permutations its tests share: each test's |z| (the group
%     comparison's or the correlation's statistic as the |z| of the same
%     tail probability on its degrees of freedom: t_to_z, PALM's palm_gtoz),
%     and the largest over the family compared with the largest of each
%     permutation (max-statistic, Westfall-Young): p_maxT_family (each
%     test, family-wise over the family) and family_p (the smallest, the
%     family as a whole); Benjamini-Hochberg q within the family
%     (q_fdr, a description); run_statistics then selects families and
%     decides (family_adjust: p_bb, significant); seeds fixed;
%     figures closed after saving, old figures and tables of the
%     family removed first (a family that fails, or no longer has a
%     comparison or correlation set, leaves no table of an earlier run);
%     a figure that fails is a warning, and the family keeps its tables and
%     tests (S.stop_on_error: the error).
%   Writes <family>_group_comparison.csv/.xlsx and <family>_correlations.csv/.xlsx
%   in out_dir; figures of the results with p_maxT_family <= alpha (the
%   tests that can be significant) in out_dir/figures (family_figure_name;
%   the p in their titles is the family-wise p, and run_statistics removes,
%   after the adjustment over families, those of the tests that are not
%   significant):
%   box plots, each group in its own colour (the same colour for a group
%   name in every figure) with its subjects' values and n; scatter plots
%   with one panel per comparison of S.comparisons (e.g. patients and
%   controls, or two subgroups of patients), each group's points in
%   its colour and marker with its least-squares line, n and raw r (of
%   S.correlation_type) in the legend, the subjects the correlation was
%   tested in marked (a dashed line when they are the panel's two groups
%   together), the same axes in every panel, and the tested r (partial with
%   covariates; Spearman's when S.correlation_type is 'spearman'), the
%   family-wise p, tested subjects and covariates in the title; a
%   group not tested whose values on the scale are all the same (controls
%   scored 0 on a symptom scale) is not drawn, a line under the axis label
%   says so. The panels show raw values;
%   the r tested is partial (covariates removed from feature and score).
T = read_subject_table(S);
ids = unique(string(F.id), 'stable');
[ok, row] = ismember(ids, T.(S.id_column));
if ~all(ok)
    error('permu:join', 'Subjects with EEG features but no row in the subject table: %s', strjoin(ids(~ok), ', '));
end
T = T(row, :);
feats = unique(F.feature, 'stable');
Y_all = containers.Map();
for f = feats'
    Y = feature_matrix(F, ids, f);
    v = Y(~isnan(Y));
    if isempty(v) || all(v == v(1)), continue; end   % constant for everyone (e.g. OrgTM diagonal)
    Y_all(char(f)) = Y;
end
feats = string(Y_all.keys)';
feats = feats(ismember(feats, unique(F.feature, 'stable')));
[~, order] = ismember(feats, unique(F.feature, 'stable'));
[~, s] = sort(order);  feats = feats(s);            % keep the order of F
Zf_all = feature_covariates(F, ids, feats);         % covariates of one feature (optional)
if nargin < 5, opts = struct(); end
clusters = isfield(opts, 'adjacency') && ~isempty(opts.adjacency);
if ~isfolder(out_dir), mkdir(out_dir); end
for stem = string(family) + ["_group_comparison", "_correlations"]   % none left from an earlier run: a family
    for ext = [".csv", ".xlsx"]                                     % that fails, or no longer has comparisons
        if isfile(fullfile(out_dir, stem + ext)), delete(fullfile(out_dir, stem + ext)); end   % or sets, has no
    end                                                             % tables, not an old one
end
fig_dir = fullfile(out_dir, 'figures');
if S.figures
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    for f = dir(fullfile(fig_dir, [family '_*.png']))', delete(fullfile(f.folder, f.name)); end   % none left from an earlier run
end

%% Group comparisons
G = {};
[Zg, zg_names] = covariate_matrix(T, S.covariates, categorical_list(S));
own_notes = strings(0);
for k = 1:numel(S.comparisons)
    cmp = S.comparisons(k);
    inA = select_subjects(T, cmp.A);
    inB = select_subjects(T, cmp.B);
    if any(inA & inB), error('permu:comparison', '%s: subjects in both groups', cmp.name); end
    if nnz(inA) < 2 || nnz(inB) < 2
        warning('permu:comparison', '%s: %d and %d subjects, skipped', cmp.name, nnz(inA), nnz(inB));
        continue
    end
    sel = inA | inB;
    [Z, names] = usable_covariates(Zg(sel, :), zg_names, cmp.name, double(inA(sel)));
    first = numel(G) + 1;  fmax = -inf(S.n_perm, 1);   % the family's largest |z| in each permutation
    for f = feats'
        Y = Y_all(char(f));
        own = own_covariates(f, names, S);           % not adjusted for itself
        [Zf, names_f] = deal(Z(:, ~own), names(~own));
        own_notes = [own_notes, f + " (" + string(names(own)) + ")"]; %#ok<AGROW>
        if isKey(Zf_all, char(f))                   % the feature's own covariate
            zf = Zf_all(char(f));
            [Zf, names_f] = usable_covariates([Zf, zf.z(sel)], [names_f, {zf.name}], cmp.name + ": " + f, double(inA(sel)));
        end
        if isempty(Zf)
            test = "Welch t, permutations";
        else
            test = "ANCOVA (" + strjoin(string(names_f), ", ") + "), Welch-type G, Freedman-Lane permutations";
        end
        po = struct('n_perm', S.n_perm, 'seed', S.seed, 'n_boot', 2000);
        [R, D] = perm_ancova_maxT(Y(sel, :), inA(sel), Zf, po);
        [fmax, z_obs] = add_to_family_max(fmax, D);
        R = R(:, {'state', 'diff', 't', 'df', 'p_unc', 'p_maxT', 'g', 'g_ci_low', 'g_ci_high', 'g_unadjusted', ...
            'g_unadjusted_ci_low', 'g_unadjusted_ci_high', 'nB', 'nA', ...
            'meanB', 'sdB', 'meanA', 'sdA'});
        R.z_family = z_obs(:);
        R = renamevars(R, {'nB', 'nA', 'meanB', 'sdB', 'meanA', 'sdA'}, ...
            {'N_A', 'N_B', 'mean_A', 'sd_A', 'mean_B', 'sd_B'});   % perm_ancova_maxT's B is comparison A
        if clusters, R = [R, cluster_columns(D, F, f, opts)]; end %#ok<AGROW>
        R = addvars(R, repmat(string(cmp.name), height(R), 1), repmat(f, height(R), 1), repmat(test, height(R), 1), ...
            'Before', 1, 'NewVariableNames', {'Comparison', 'Parameter', 'StatisticalTest'});
        R.state(:) = state_numbers(F, f, height(R));
        G{end + 1} = R; %#ok<AGROW>
    end
    G(first:end) = family_wise(G(first:end), fmax);
end
if ~isempty(G)
    G = vertcat(G{:});
    G.q_fdr = nan(height(G), 1);
    for c = unique(G.Comparison)'
        m = G.Comparison == c;
        G.q_fdr(m) = fdr_bh(G.p_unc(m));
    end
    G = addvars(G, repmat(string(family), height(G), 1), 'Before', 1, 'NewVariableNames', 'Family');
    write_both(G, fullfile(out_dir, [family '_group_comparison']));
else
    G = table();
end

%% Correlations with clinical scores
sets = struct('name', {}, 'columns', {}, 'subjects', {}, 'covariates', {}, 'zero_as_missing', {});
for k = 1:numel(S.correlations), sets(end + 1) = S.correlations(k); end %#ok<AGROW>
for k = 1:numel(S.covariate_sets)             % correlation sets with more covariates
    cs = S.covariate_sets(k);
    base = S.correlations(strcmp({S.correlations.name}, cs.correlations));
    if isempty(base), warning('permu:covset', 'Covariate set %s: no correlation set %s', cs.name, cs.correlations); continue; end
    base.name = cs.name;
    base.covariates = cs.covariates;
    if isfield(cs, 'columns') && ~isempty(cs.columns), base.columns = cs.columns; end
    absent = setdiff(base.covariates, T.Properties.VariableNames);
    if ~isempty(absent)
        warning('permu:covset', 'Covariate set %s skipped: no column %s in the subject table', cs.name, strjoin(absent, ', '));
        continue
    end
    sets(end + 1) = base; %#ok<AGROW>
end
C = {};
names_T = string(T.Properties.VariableNames);
for k = 1:numel(sets)
    cs = sets(k);
    scales = names_T(~cellfun(@isempty, regexp(cellstr(names_T), cs.columns, 'once')));
    scales = setdiff(scales, [string(cs.covariates), string(S.id_column), string(S.group_column)], 'stable');
    text = scales(~arrayfun(@(c) isnumeric(T.(c)), scales));
    if ~isempty(text)                              % e.g. scores read as text ("12,5", "n/a")
        warning('permu:scale_text', ['Correlation set %s: column(s) %s match its columns but are not numeric ' ...
            '(text in the subject table): not tested. Write them as numbers, missing values empty'], ...
            cs.name, strjoin(text, ', '));
    end
    scales = setdiff(scales, text, 'stable');
    if isempty(scales), continue; end
    insub = select_subjects(T, cs.subjects);
    [Zc, zc_names] = covariate_matrix(T(insub, :), cs.covariates, categorical_list(S));
    [Zc, zc_names] = usable_covariates(Zc, zc_names, cs.name, []);
    Ck = {};  fmax = -inf(S.n_perm, 1);         % the family's largest |z| in each permutation
    for f = feats'
        Y = Y_all(char(f));
        Y = Y(insub, :);
        own = own_covariates(f, zc_names, S);        % not adjusted for itself
        [Zf, names_f] = deal(Zc(:, ~own), zc_names(~own));
        own_notes = [own_notes, f + " (" + string(zc_names(own)) + ")"]; %#ok<AGROW>
        if isKey(Zf_all, char(f))                   % the feature's own covariate
            zf = Zf_all(char(f));
            [Zf, names_f] = usable_covariates([Zf, zf.z(insub)], [names_f, {zf.name}], cs.name + ": " + f, []);
        end
        for sc = scales
            x = double(T.(sc)(insub));
            if cs.zero_as_missing, x(x == 0) = NaN; end
            if nnz(~isnan(x)) < 4 + size(Zf, 2), continue; end
            po = struct('n_perm', S.n_perm, 'seed', S.seed, 'type', S.correlation_type, 'method', S.correlation_test);
            [R, D] = perm_partial_corr_maxT(Y, x, Zf, po);
            [fmax, z_obs] = add_to_family_max(fmax, D);
            R.z_family = z_obs(:);
            if clusters, R = [R, cluster_columns(D, F, f, opts)]; end %#ok<AGROW>
            R = addvars(R, repmat(string(list_or_none(names_f)), height(R), 1), repmat(f, height(R), 1), ...
                repmat(sc, height(R), 1), repmat(correlation_test(S, ~isempty(Zf)), height(R), 1), ...
                'Before', 1, 'NewVariableNames', {'Covariates', 'Parameter', 'Scale', 'StatisticalTest'});
            R.state(:) = state_numbers(F, f, height(R));
            Ck{end + 1} = R; %#ok<AGROW>
        end
    end
    if isempty(Ck), continue; end
    Ck = family_wise(Ck, fmax);
    Ck = vertcat(Ck{:});
    Ck.q_fdr = fdr_bh(Ck.p_unc);
    Ck = addvars(Ck, repmat(string(cs.name), height(Ck), 1), 'Before', 1, 'NewVariableNames', 'CorrelationSet');
    C{end + 1} = Ck; %#ok<AGROW>
end
if ~isempty(C)
    C = vertcat(C{:});
    C = addvars(C, repmat(string(family), height(C), 1), 'Before', 1, 'NewVariableNames', 'Family');
    write_both(C, fullfile(out_dir, [family '_correlations']));
else
    C = table();
end

if ~isempty(own_notes)
    fprintf('%s: a covariate is not used for its own measure: %s\n', family, strjoin(unique(own_notes, 'stable'), ', '));
end

%% Figures of the results with a family-wise p <= alpha (closed after saving)
if ~S.figures, return; end
try                                     % a failing figure keeps the family's tables and tests
    P = figure_groups(S, T);
    if ~isempty(G)
        for i = find(G.p_maxT_family <= S.alpha)'
            gk = P([P.comparison] == find(strcmp({S.comparisons.name}, G.Comparison(i)), 1));
            Y = Y_all(char(G.Parameter(i)));
            col = find(state_numbers(F, G.Parameter(i), size(Y, 2)) == G.state(i));
            fig = figure('Visible', 'off', 'Color', 'w');
            ax = axes(fig);  hold(ax, 'on');
            xt = strings(1, 2);
            for j = 1:2
                y = Y(gk(j).in, col);  y = y(~isnan(y));  n = numel(y);
                xt(j) = sprintf('%s (n = %d)', gk(j).name, n);
                boxchart(ax, repmat(j, n, 1), y, 'BoxFaceColor', gk(j).colour, 'WhiskerLineColor', gk(j).colour, ...
                    'MarkerStyle', 'none', 'BoxWidth', 0.5);
                swarmchart(ax, repmat(j, n, 1), y, 22, gk(j).colour, gk(j).marker, 'filled', ...
                    'MarkerFaceAlpha', 0.75, 'XJitterWidth', 0.35);
            end
            set(ax, 'XTick', 1:2, 'XTickLabel', xt, 'TickLabelInterpreter', 'none');
            xlim(ax, [0.4 2.6]);
            ylabel(ax, G.Parameter(i), 'Interpreter', 'none');
            title(ax, sprintf('%s%s: %s, g = %.2f, family-wise p = %.4f', G.Parameter(i), state_text(G.state(i)), ...
                G.Comparison(i), G.g(i), G.p_maxT_family(i)), 'Interpreter', 'none');
            exportgraphics(fig, fullfile(fig_dir, family_figure_name(G, i)));
            close(fig);
        end
    end
    if ~isempty(C)
        for i = find(C.p_maxT_family <= S.alpha)'
            cs = sets(strcmp({sets.name}, C.CorrelationSet(i)));
            Y = Y_all(char(C.Parameter(i)));
            col = find(state_numbers(F, C.Parameter(i), size(Y, 2)) == C.state(i));
            x = double(T.(C.Scale(i)));
            if cs.zero_as_missing, x(x == 0) = NaN; end
            y = Y(:, col);
            has = ~isnan(x) & ~isnan(y);
            tested = select_subjects(T, cs.subjects) & has;
            nc = max(1, numel(S.comparisons));
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 470 * nc 430]);
            tl = tiledlayout(fig, 1, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
            axs = gobjects(1, nc);
            for k = 1:nc
                ax = nexttile(tl);  hold(ax, 'on');  axs(k) = ax;
                gk = P([P.comparison] == k);
                if isempty(gk)                                  % no comparisons: the tested subjects alone
                    gk = struct('name', "tested", 'spec', {cs.subjects}, 'in', tested, 'colour', [0.3 0.3 0.3], ...
                        'marker', 'o', 'comparison', 1);
                end
                h = gobjects(0);  lab = strings(0);  is_tested = false;  notes = strings(0);
                drawn = false(numel(gk), 1);
                for j = 1:numel(gk)
                    m = gk(j).in & has;
                    same = isequal(m, tested);
                    if ~same && nnz(m) > 1 && all(x(m) == x(find(m, 1)))   % no spread on the scale: not drawn
                        notes(end + 1) = sprintf('%s: all %d at %g, not shown', gk(j).name, nnz(m), x(find(m, 1))); %#ok<AGROW>
                        continue
                    end
                    drawn(j) = true;
                    h(end + 1) = scatter(ax, x(m), y(m), 30, gk(j).colour, gk(j).marker, 'filled', 'MarkerFaceAlpha', 0.75); %#ok<AGROW>
                    is_tested = is_tested || same;
                    lab(end + 1) = group_label(gk(j).name, x(m), y(m), same, S.correlation_type); %#ok<AGROW>
                    fit_line(ax, x(m), y(m), gk(j).colour, '-');
                end
                shown = false(size(has));
                for j = find(drawn)', shown = shown | (gk(j).in & has); end
                if ~is_tested && isequal(shown, tested)         % the tested subjects are the panel's groups together
                    h(end + 1) = fit_line(ax, x(tested), y(tested), [0 0 0], '--'); %#ok<AGROW>
                    lab(end + 1) = group_label("all tested", x(tested), y(tested), true, S.correlation_type); %#ok<AGROW>
                end
                legend(ax, h, lab, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 8);
                if k <= numel(S.comparisons), title(ax, S.comparisons(k).name, 'Interpreter', 'none'); end
                xlabel(ax, [string(C.Scale(i)), notes], 'Interpreter', 'none');   % notes: groups not drawn
                ylabel(ax, C.Parameter(i), 'Interpreter', 'none');
            end
            linkaxes(axs, 'xy');
            rname = r_name(S.correlation_type);
            if C.Covariates(i) ~= "none", rname = "partial " + rname; end
            title(tl, sprintf('%s%s: %s = %.2f, family-wise p = %.4f (tested in %s; covariates: %s)', ...
                C.CorrelationSet(i), state_text(C.state(i)), rname, C.r(i), C.p_maxT_family(i), describe_subjects(cs.subjects, P), ...
                C.Covariates(i)), 'Interpreter', 'none', 'FontSize', 10);
            exportgraphics(fig, fullfile(fig_dir, family_figure_name(C, i)));
            close(fig);
        end
    end
catch err
    if isfield(S, 'stop_on_error') && S.stop_on_error, rethrow(err); end
    warning('permu_ttest_permu_corr_BNS:figures', '%s: figures not made (tables and tests kept): %s', ...
        family, err.message);
end
end

function [fmax, z_obs] = add_to_family_max(fmax, D)
% the tests of one feature (D of perm_ancova_maxT or perm_partial_corr_maxT,
% on the permutations every test of the family shares) added to the
% family's largest |z| in each permutation; z_obs: their |z| in the data
% (t_to_z; perm_ancova_maxT gives z already, with df Inf)
z_obs = abs(t_to_z(D.t_obs, D.df));
fmax = max(fmax, max(abs(t_to_z(D.T, repmat(D.df, size(D.T, 1), 1))), [], 2));   % max omits NaN (tests not done)
end

function Rs = family_wise(Rs, fmax)
% p_maxT_family of every test of a family (the result tables Rs of its
% features, each with z_family), from the family's largest |z| in each
% permutation: the share of permutations whose largest |z| is at least the
% test's |z| (one plus, over n_perm + 1); family_p, the same for every row:
% the smallest of them, the p of the family as a whole (its largest |z|
% against the permutations' largest). The single-step max-statistic
% (Westfall & Young 1993) as PERMUTOOLS' permuttest2 computes it with
% 'correct' true (Crosse et al. 2024; permuttest2.m lines 208-221: the
% largest |t| over the variables in each permutation, p = (1 + #{max >=
% |t|}) / (n_perm + 1)), and FieldTrip's ft_statistics_montecarlo with
% correctm 'max' (lines 347-352); validated against permuttest2 in
% docs/validation.md ("Multiple testing: the max-statistic family test")
z = cell2mat(cellfun(@(R) R.z_family, Rs(:), 'UniformOutput', false));
n_perm = numel(fmax);
p = nan(size(z));
for i = find(~isnan(z))'
    p(i) = (1 + sum(fmax >= z(i) * (1 - 1e-10))) / (n_perm + 1);   % a permutation equal to the data counts despite rounding
end
fp = min(p, [], 'omitnan');
if isempty(fp), fp = NaN; end
i0 = 0;
for i = 1:numel(Rs)
    m = height(Rs{i});
    Rs{i}.p_maxT_family = p(i0 + (1:m));
    Rs{i}.family_p = repmat(fp, m, 1);
    Rs{i}.z_family = [];
    i0 = i0 + m;
end
end

function Cl = cluster_columns(D, F, f, opts)
% the cluster test of one feature's states (cluster_perm_test), with the
% adjacency of its states (the state numbers index opts.adjacency)
st = unique(F.state(F.feature == f));
st = st(1:min(numel(D.t_obs), end));
Cl = cluster_perm_test(D.t_obs, D.T, D.df, opts.adjacency(st, st), opts.cluster_alpha);
end

function P = figure_groups(S, T)
% the two groups of every comparison, named as in "A - B": their subjects, a
% colour per group name (the same in every figure; Okabe-Ito hues, distinct
% under colour blindness) and a marker per side (A circle, B triangle)
pal = [213 94 0; 0 114 178; 170 51 119; 0 158 115; 230 159 0; 86 180 233; 102 102 102] / 255;
P = struct('name', {}, 'spec', {}, 'in', {}, 'colour', {}, 'marker', {}, 'comparison', {});
seen = strings(0);
for k = 1:numel(S.comparisons)
    cmp = S.comparisons(k);
    names = split(string(cmp.name), " - ");
    if numel(names) ~= 2, names = string(cmp.name) + [" (A)"; " (B)"]; end
    specs = {cmp.A, cmp.B};  marker = {'o', '^'};
    for j = 1:2
        c = find(seen == names(j), 1);
        if isempty(c), seen(end + 1) = names(j); c = numel(seen); end %#ok<AGROW>
        P(end + 1) = struct('name', names(j), 'spec', {specs{j}}, 'in', select_subjects(T, specs{j}), ...
            'colour', pal(min(c, end), :), 'marker', marker{j}, 'comparison', k); %#ok<AGROW>
    end
end
end

function s = state_text(k)
% ', state k' in a figure title; nothing for a feature without states (0)
s = '';
if k > 0, s = sprintf(', state %d', k); end
end

function h = fit_line(ax, x, y, colour, style)
% least-squares line over the range of x (an empty line when it is undefined)
h = plot(ax, NaN, NaN, style, 'Color', colour, 'LineWidth', 1.5);
if numel(x) >= 3 && std(x) > 0
    b = [ones(numel(x), 1), x(:)] \ y(:);
    set(h, 'XData', [min(x), max(x)], 'YData', b(1) + b(2) * [min(x), max(x)]);
end
end

function s = group_label(name, x, y, tested, type)
if isempty(x), s = name + ": no values"; return; end
r = NaN;
if numel(x) >= 3 && std(x) > 0 && std(y) > 0, r = corr(x(:), y(:), 'Type', type); end
s = sprintf('%s (n = %d): raw %s = %.2f', name, numel(x), r_name(type), r);
if isnan(r), s = sprintf('%s (n = %d): %s undefined', name, numel(x), r_name(type)); end
if tested, s = s + ", tested"; end
end

function s = r_name(type)
% 'r' for Pearson, 'Spearman r' for Spearman
s = "r";
if strcmpi(type, 'spearman'), s = "Spearman r"; end
end

function s = describe_subjects(spec, P)
% a group's name when the condition is one of the comparisons' groups
for g = P
    if isequal(g.spec, spec), s = g.name; return; end
end
if isempty(spec), s = "all subjects"; return; end
s = strjoin(compose('%s %s %s', string(spec(:, 1)), string(spec(:, 2)), string(spec(:, 3))), ' and ');
end

function Y = feature_matrix(F, ids, f)
% subjects x states (one column for a feature without states)
sel = F(F.feature == f, :);
st = unique(sel.state);
[~, col] = ismember(sel.state, st);
[~, row] = ismember(string(sel.id), ids);
Y = nan(numel(ids), numel(st));
Y(sub2ind(size(Y), row, col)) = sel.value;
end

function Z = feature_covariates(F, ids, feats)
% the covariate of each feature that has one (F columns covariate and
% covariate_name): name and subjects x 1 values, the same for every state
Z = containers.Map();
if ~all(ismember({'covariate', 'covariate_name'}, F.Properties.VariableNames)), return; end
for f = feats'
    sel = F(F.feature == f, :);
    nm = unique(string(sel.covariate_name));
    nm = nm(~ismissing(nm) & nm ~= "");
    if isempty(nm), continue; end
    if numel(nm) > 1, error('permu:covariate', '%s: more than one covariate of its own (%s)', f, strjoin(nm, ', ')); end
    sel.value = sel.covariate;
    z = feature_matrix(sel, ids, f);
    if any(max(z, [], 2) > min(z, [], 2))
        error('permu:covariate', '%s: its covariate %s differs between states of a subject', f, nm);
    end
    Z(char(f)) = struct('name', char(nm), 'z', min(z, [], 2));
end
end

function st = state_numbers(F, f, n)
st = unique(F.state(F.feature == f));
st = st(1:min(n, end));
end

function own = own_covariates(f, names, S)
% the covariate columns that are the feature F itself (a text covariate's
% indicators 'col=level' by their column): the same name, or both among the
% names of the IAF (IAF; IAF_cov, the IAF or the CoG without a peak; the
% IAF_PAF or IAF_CoG that the IAF is, S.iaf.measure)
base = regexprep(string(names), '=.*$', '');
own = base == string(f);
iaf = ["IAF", "IAF_cov"];
if isfield(S, 'iaf') && isfield(S.iaf, 'measure')
    if strcmpi(S.iaf.measure, 'paf'), iaf(end + 1) = "IAF_PAF"; else, iaf(end + 1) = "IAF_CoG"; end
end
if ismember(string(f), iaf), own = own | ismember(base, iaf); end
end

function c = categorical_list(S)
% numeric covariates to code as indicators (S.categorical_covariates)
c = {};
if isfield(S, 'categorical_covariates'), c = cellstr(S.categorical_covariates); end
end

function s = list_or_none(c)
if isempty(c), s = 'none'; else, s = strjoin(c, ', '); end
end

function s = correlation_test(S, partial)
% the StatisticalTest of a correlation: its type and its null distribution
type = S.correlation_type;
s = string(upper(type(1))) + type(2:end) + " r, slope / HC3 SE";
if partial, s = "partial " + s; end
if strcmp(S.correlation_test, 'wild_bootstrap')
    s = s + ", wild bootstrap (Rademacher, null imposed)";
elseif partial
    s = s + ", Smith permutations";               % the score's residuals
else
    s = s + ", permutations";
end
end

function write_both(T, stem)
writetable(T, [stem '.csv']);
try
    writetable(T, [stem '.xlsx']);
catch
end
end
