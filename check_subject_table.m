function report = check_subject_table(S)
%CHECK_SUBJECT_TABLE The subject table against the recordings, before anything is computed.
%   report = check_subject_table(S), with S checked (setup_toolboxes), is
%   called by run_full_pipeline before stage 1 when a stage that uses the
%   table will run (1b, 2 or 3), so that a table that does not fit the
%   recordings stops the run at once instead of after hours of stage 1. It
%   reads the table (read_subject_table: IDs matched to the recordings
%   exactly, or one to one ignoring case and leading zeros, printed) and the
%   recordings of S.data (get_subs, subject_ids), and
%   - lists the recordings without a row and the rows without a recording;
%   - checks that S.group_column is a column, and counts the recordings of
%     each group of S.groups and those in none (left out of stages 2 and 3);
%   - checks that every column named by S.comparisons and by the
%     covariates of S.covariates and S.correlations is in the table (a
%     covariate made by the pipeline, a spectral-profile or IAF measure,
%     excepted), and that every correlation set matches a numeric column;
%   - prints the sample size of every comparison and correlation set, with
%     the covariates each can use (check_sample_sizes).
%   Errors (check_subject_table:*): the table not found; no recording with a
%   row (with examples of both kinds of ID, to fix S.id_regexp or
%   S.id_column); S.group_column not a column; no recording in S.groups;
%   a column of a comparison condition or of a covariate (S.covariates,
%   S.correlations) not in the table. Otherwise a short summary.
%   Without a subject table (S.subs_info empty) only a note: the groups then
%   come from the start of the IDs (S.group_column = ''). Without recordings
%   in S.data, the recordings of stage1_summary.csv are checked; without
%   either, nothing (stage 1 reports the missing recordings).
%   report: struct with the recordings (ids), those without a row
%   (no_row), the rows without a recording (no_recording), the IDs matched
%   loosely (matched), the recordings per group (n_group) and the sample
%   sizes (sizes, check_sample_sizes).
report = struct('ids', strings(0, 1), 'no_row', strings(0, 1), 'no_recording', strings(0, 1), ...
    'matched', table(), 'n_group', [], 'sizes', table());
fprintf('Subject table: ');
if isempty(S.subs_info)
    fprintf('none (S.subs_info is empty): groups from the start of the IDs; stage 3 needs a table\n');
    return
end
if ~isfile(S.subs_info)
    error('check_subject_table:file', 'Subject table not found: %s (S.subs_info)', S.subs_info);
end
fprintf('%s\n', S.subs_info);
try
    subs = get_subs(S.data, S.file_pattern);
catch err                                         % no recordings there: those of stage 1, if any
    f = fullfile(S.results, 'Connectivity', 'stage1_summary.csv');
    if ~isfile(f)                                 % stage 1 will say so
        fprintf('  not checked: %s\n', err.message);
        return
    end
    Q = readtable(f, 'Delimiter', ',', 'TextType', 'string');
    subs = cellstr(Q.file);
    fprintf('  (%s: the recordings of stage1_summary.csv)\n', err.message);
end
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
ids = subject_ids(S, stems);
[T, matched] = read_subject_table(S);
tid = T.(S.id_column);
has = ismember(ids, tid);
report.ids = ids;  report.matched = matched;
report.no_row = ids(~has);
report.no_recording = tid(~ismember(tid, ids));
if height(matched)
    fprintf('  IDs matched ignoring case or leading zeros: %s\n', strjoin(matched.table_id + " -> " + matched.recording_id, ', '));
end
if ~any(has)
    error('check_subject_table:none', ['No recording has a row in the subject table: recording IDs (S.id_regexp ' ...
        'on the file names) such as %s, table IDs (column %s) such as %s. Fix S.id_regexp or S.id_column'], ...
        strjoin(ids(1:min(3, end)), ', '), S.id_column, strjoin(tid(1:min(3, end)), ', '));
end
fprintf('  %d recordings, %d with a row', numel(ids), nnz(has));
if any(~has), fprintf('; without a row (left out of stages 2 and 3): %s', some(ids(~has))); end
fprintf('\n');
if ~isempty(report.no_recording)
    fprintf('  %d row(s) without a recording (not used): %s\n', numel(report.no_recording), some(report.no_recording));
end
% groups
Tr = T(ismember(tid, ids), :);
if ~isempty(S.group_column)
    if ~ismember(S.group_column, T.Properties.VariableNames)
        error('check_subject_table:group_column', 'S.group_column: %s is not a column of %s; columns: %s', ...
            S.group_column, S.subs_info, strjoin(T.Properties.VariableNames, ', '));
    end
    g = lower(strtrim(string(Tr.(S.group_column))));
    [~, gi] = ismember(g, lower(strtrim(string(S.groups))));
else
    gi = zeros(height(Tr), 1);
    for k = 1:numel(S.groups)
        gi(startsWith(Tr.(S.id_column), S.groups{k}, 'IgnoreCase', true) & gi == 0) = k;
    end
end
report.n_group = accumarray(gi(gi > 0), 1, [numel(S.groups), 1])';
if ~any(gi)
    src = 'the start of the IDs';
    if ~isempty(S.group_column), src = ['column ' S.group_column]; end
    error('check_subject_table:groups', 'No recording is in a group of S.groups (%s) by %s', ...
        strjoin(S.groups, ', '), src);
end
fprintf('  groups: %s', strjoin(compose('%s %d', string(S.groups(:)), report.n_group(:)), ', '));
if any(gi == 0), fprintf('; in none (left out): %s', some(Tr.(S.id_column)(gi == 0))); end
fprintf('\n');
% columns named by the comparisons and the covariates
made = [spectral_profile_measures(), {'IAF', 'IAF_PAF', 'IAF_CoG', 'IAF_IRASA', 'IAF_specparam'}];
cols = T.Properties.VariableNames;
for k = 1:numel(S.comparisons)
    c = S.comparisons(k);
    spec = [c.A; c.B];
    if isempty(spec), continue; end
    bad = setdiff(spec(:, 1), cols);
    if ~isempty(bad)
        error('check_subject_table:column', 'Comparison %s: column(s) %s not in the subject table', ...
            c.name, strjoin(bad, ', '));
    end
end
need = S.covariates;
for k = 1:numel(S.correlations), need = [need, S.correlations(k).covariates]; end %#ok<AGROW>
bad = setdiff(setdiff(need, cols), made);
if ~isempty(bad)
    error('check_subject_table:column', 'Covariate(s) %s not in the subject table (S.covariates, S.correlations)', ...
        strjoin(bad, ', '));
end
for k = 1:numel(S.covariate_sets)                 % a set with an absent column is skipped (with a message)
    bad = setdiff(setdiff(S.covariate_sets(k).covariates, cols), made);
    if ~isempty(bad)
        warning('check_subject_table:covariate_set', 'Covariate set %s: column(s) %s not in the subject table: it will be skipped', ...
            S.covariate_sets(k).name, strjoin(bad, ', '));
    end
end
report.sizes = check_sample_sizes(S, ids(has));
end

function s = some(x)
% up to 10 of x, and how many more
x = string(x(:));
s = strjoin(x(1:min(10, end)), ', ');
if numel(x) > 10, s = sprintf('%s and %d more', s, numel(x) - 10); end
end
