function report = check_sample_sizes(S, ids)
%CHECK_SAMPLE_SIZES Recordings of every comparison and correlation set against what their tests need.
%   report = check_sample_sizes(S), called by run_statistics before the
%   families (after the spectral profile, whose measures can be
%   covariates): the recordings of every comparison and correlation set,
%   against what their tests need. A comparison needs at least 2
%   recordings per group and the covariates + 4 in all (perm_ancova_maxT:
%   2 residual degrees of freedom); a correlation the covariates + 4
%   recordings with the score (perm_partial_corr_maxT). Those that cannot
%   be tested give rows of NaN (comparisons) or no rows (correlations) in
%   every family: a warning names them here, once, with their numbers.
%   The covariates are counted as the tests use them: per comparison or
%   correlation set, those usable_covariates keeps (a covariate missing in
%   a group, constant, collinear or of one value in one group is left out
%   there, and named), a text covariate one column per level but the first,
%   a numeric one in S.categorical_covariates one per value but the first.
%   Recordings: those of stage 1 with outputs (stage1_summary.csv, status
%   ok, or skipped (up to date) when run_all reused them) with a row
%   in the subject table and every covariate used.
%   report = check_sample_sizes(S, ids) counts the recordings of subject
%   IDs ids instead (check_subject_table, before stage 1; the PSD
%   statistics). A covariate the pipeline makes later (the spectral
%   profile, the IAF of stage 1) and not yet in the table is not counted,
%   with a note.
%   Also:
%   - a note when one group of a comparison is at least 3 times larger
%     than the other: with unequal groups and variances the groups are not
%     exchangeable under the null; with the z of each permutation on its
%     own degrees of freedom (as PALM, since 5 November 2026) the p-values
%     stayed at or below their level in a simulation of 8 against 30
%     recordings with 3 times the SD (the raw statistics gave up to about
%     twice too many below 0.01; docs/validation.md, "Unequal groups");
%   - a warning naming the numeric covariates not in
%     S.categorical_covariates with 3 to 5 distinct whole values (a site
%     coded 1, 2, 3 enters as one straight line, not as categories);
%   - the IDs of the table matched to recordings ignoring case or leading
%     zeros (read_subject_table).
%   report: one row per comparison or correlation set (kind, name, n,
%   needed, testable); empty without recordings or a subject table.
report = table('Size', [0 5], 'VariableTypes', {'string', 'string', 'string', 'double', 'logical'}, ...
    'VariableNames', {'kind', 'name', 'n', 'needed', 'testable'});
if isempty(S.subs_info) || ~isfile(S.subs_info), return; end
if nargin < 2
    f = fullfile(S.results, 'Connectivity', 'stage1_summary.csv');
    if ~isfile(f), return; end
    o = detectImportOptions(f, 'FileType', 'text', 'Delimiter', ',');   % run_all writes commas
    o = setvartype(o, intersect({'file', 'status'}, o.VariableNames), 'string');
    R = readtable(f, o);
    if ~all(ismember({'file', 'status'}, R.Properties.VariableNames)), return; end
    done = R.status == "ok" | R.status == "skipped (up to date)";      % outputs made now or reused
    [~, stems] = cellfun(@fileparts, cellstr(R.file(done)), 'UniformOutput', false);
    ids = subject_ids(S, stems);
end
[T, matched] = read_subject_table(S);
if height(matched)
    fprintf('Subject IDs matched to the recordings ignoring case or leading zeros: %s\n', ...
        strjoin(matched.table_id + " -> " + matched.recording_id, ', '));
end
T = T(ismember(T.(S.id_column), string(ids)), :);
catg = {};
if isfield(S, 'categorical_covariates'), catg = cellstr(S.categorical_covariates); end
fprintf('Sample sizes (recordings in the subject table, with the covariates each test can use):\n');
categorical_warning(T, S, catg);
covs = available(T, S.covariates, 'S.covariates');
Zg = [];  zg_names = {};
if ~isempty(T), [Zg, zg_names] = covariate_matrix(T, covs, catg); end
for k = 1:numel(S.comparisons)
    c = S.comparisons(k);
    try
        inA = select_subjects(T, c.A);
        inB = select_subjects(T, c.B);
    catch err
        warning('run_statistics:sample_size', 'Comparison %s: %s', c.name, err.message);
        continue
    end
    sel = inA | inB;
    note = '';  Z = zeros(nnz(sel), 0);  names = {};
    if nnz(inA) >= 2 && nnz(inB) >= 2
        [Z, names, note] = usable_covariates(Zg(sel, :), zg_names, c.name, double(inA(sel)));
    end
    okr = all(~isnan(Z), 2);
    a = inA(sel);  nA = nnz(a & okr);  nB = nnz(~a & okr);
    need = size(Z, 2) + 4;
    ok = nA >= 2 && nB >= 2 && nA + nB >= need;
    fprintf('  comparison %s: %d and %d, %d covariate column(s) (%s)%s\n', c.name, nA, nB, size(Z, 2), ...
        list_or_none(names), verdict(ok));
    if ~isempty(note), fprintf('    %s\n', note); end
    if ok && max(nA, nB) >= 3 * min(nA, nB)
        fprintf(['    groups of %d and %d (one at least 3 times the other): with unequal variances too, the groups ' ...
            'are not exchangeable; the p-values held their level in a simulation of such groups, read them ' ...
            'with the group SDs (docs/validation.md, "Unequal groups")\n'], nA, nB);
    end
    report(end + 1, :) = {"comparison", string(c.name), sprintf('%d + %d', nA, nB), need, ok}; %#ok<AGROW>
    if ~ok
        warning('run_statistics:sample_size', ['Comparison %s: %d and %d recordings with %d covariate column(s): ' ...
            'not testable (needs at least 2 per group and %d in all, the covariates + 4); its rows will be NaN'], ...
            c.name, nA, nB, size(Z, 2), need);
    end
end
sets = S.correlations;                       % with the covariate sets, as permu_ttest_permu_corr_BNS
for k = 1:numel(S.covariate_sets)
    cs = S.covariate_sets(k);
    base = S.correlations(strcmp({S.correlations.name}, cs.correlations));
    if isempty(base), continue; end              % warned there
    base.name = cs.name;
    base.covariates = cs.covariates;
    if isfield(cs, 'columns') && ~isempty(cs.columns), base.columns = cs.columns; end
    sets(end + 1) = base; %#ok<AGROW>
end
names_T = string(T.Properties.VariableNames);
for k = 1:numel(sets)
    cs = sets(k);
    sc = names_T(~cellfun(@isempty, regexp(cellstr(names_T), cs.columns, 'once')));
    sc = setdiff(sc, [string(cs.covariates), string(S.id_column), string(S.group_column)], 'stable');
    sc = sc(arrayfun(@(v) isnumeric(T.(v)), sc));
    if isempty(sc)
        fprintf('  correlations %s: no numeric column matches %s%s\n', cs.name, cs.columns, verdict(false));
        report(end + 1, :) = {"correlations", string(cs.name), "0", NaN, false}; %#ok<AGROW>
        warning('run_statistics:sample_size', 'Correlation set %s: no numeric column of the subject table matches %s: not tested', ...
            cs.name, cs.columns);
        continue
    end
    try
        insub = select_subjects(T, cs.subjects);
        [Zc, zc_names] = covariate_matrix(T(insub, :), available(T, cs.covariates, cs.name), catg);
    catch err
        warning('run_statistics:sample_size', 'Correlation set %s: %s', cs.name, err.message);
        continue
    end
    [Zc, zc_names, note] = usable_covariates(Zc, zc_names, cs.name, []);
    okc = all(~isnan(Zc), 2);
    n = zeros(1, numel(sc));
    for j = 1:numel(sc)
        x = double(T.(sc(j))(insub));
        if cs.zero_as_missing, x(x == 0) = NaN; end
        n(j) = nnz(~isnan(x) & okc);
    end
    needc = size(Zc, 2) + 4;
    ok = max(n) >= needc;
    fprintf('  correlations %s: %d score(s), %d to %d recordings, %d covariate column(s) (%s)%s\n', cs.name, numel(sc), ...
        min(n), max(n), size(Zc, 2), list_or_none(zc_names), verdict(ok));
    if ~isempty(note), fprintf('    %s\n', note); end
    report(end + 1, :) = {"correlations", string(cs.name), string(max(n)), needc, ok}; %#ok<AGROW>
    if ~ok
        warning('run_statistics:sample_size', ['Correlation set %s: at most %d recordings with a score and %d covariate ' ...
            'column(s): not testable (needs the covariates + 4 = %d); it will have no rows'], cs.name, max(n), size(Zc, 2), needc);
    elseif min(n) < needc
        fprintf('    %d of its %d scores have fewer than %d recordings: not tested\n', nnz(n < needc), numel(n), needc);
    end
end
end

function covs = available(T, covs, what)
% the covariates in the table; those the pipeline makes later (spectral
% profile, IAF of stage 1) and not yet there are not counted (a note), any
% other absent one is named in a warning (its tests will fail)
covs = cellstr(covs);
absent = setdiff(covs, T.Properties.VariableNames, 'stable');
if isempty(absent), return; end
later = absent(ismember(absent, [spectral_profile_measures(), {'IAF', 'IAF_PAF', 'IAF_CoG', 'IAF_IRASA', 'IAF_specparam'}]));
if ~isempty(later)
    fprintf('  %s: covariate(s) %s made later by the pipeline (spectral profile or IAF): not counted here\n', ...
        what, strjoin(later, ', '));
end
other = setdiff(absent, later, 'stable');
if ~isempty(other)
    warning('run_statistics:sample_size', '%s: covariate(s) %s not in the subject table', what, strjoin(other, ', '));
end
covs = setdiff(covs, absent, 'stable');
end

function categorical_warning(T, S, catg)
% numeric covariates with 3 to 5 distinct whole values, not declared
% categorical: probably codes (site, scanner) entered as one straight line
c = cellstr(S.covariates);
for f = {'correlations', 'covariate_sets'}
    for k = 1:numel(S.(f{1})), c = [c, cellstr(S.(f{1})(k).covariates)]; end %#ok<AGROW>
end
c = setdiff(intersect(unique(c), T.Properties.VariableNames), catg);
odd = strings(0);
for v = reshape(string(c), 1, [])
    x = T.(v);
    if ~isnumeric(x), continue; end
    x = x(~isnan(x));
    u = unique(x);
    if numel(u) >= 3 && numel(u) <= 5 && all(u == round(u))
        odd(end + 1) = sprintf('%s (%s)', v, strjoin(string(u'), ', ')); %#ok<AGROW>
    end
end
if ~isempty(odd)
    warning('check_sample_sizes:categorical', ['Numeric covariate(s) with 3 to 5 distinct whole values: %s. ' ...
        'Each enters the model as one straight line; if it codes categories (a site, a scanner), name it in ' ...
        'S.categorical_covariates (one indicator per value but the first)'], strjoin(odd, '; '));
end
end

function s = verdict(ok)
if ok, s = ''; else, s = ': NOT TESTABLE'; end
end

function s = list_or_none(c)
if isempty(c), s = 'none'; else, s = strjoin(c, ', '); end
end
