function [T, matched] = read_subject_table(S)
%READ_SUBJECT_TABLE The subject table (phenotypes, scores, covariates) as a MATLAB table.
%   T = read_subject_table(S) reads S.subs_info: .csv (any delimiter), .tsv
%   (e.g. BIDS participants.tsv) or .xlsx/.xls. Column names are kept as
%   written ('VariableNamingRule', 'preserve'); 'n/a', 'NA', 'NaN' and
%   empty cells are missing values. The ID column (S.id_column) becomes a
%   string column, without surrounding spaces; IDs must be unique. Text
%   columns holding only numbers - e.g. with decimal commas ("286,5") from
%   Excel with European settings - become numeric. A covariate (of
%   S.covariates, S.correlations, S.covariate_sets; not in
%   S.categorical_covariates) whose values are mostly numbers but not all
%   is an error naming the others ("12 years", "<5"): read as text it would
%   enter as one indicator per value, read as numbers they would be missing
%   values. Numeric columns holding a usual
%   missing-value code (-9999, -999, -99, -9, 99, 999, 9999) far from their
%   other values (more than 10 robust SD from their median) give a warning,
%   once per file and MATLAB session: write them as empty or NA.
%   IDs are matched to the recordings (the subject IDs of the files of
%   S.data, subject_ids) exactly; an ID of the table without an exact match
%   that matches one recording, and only one, ignoring case (sub_001 for
%   SUB_001), or then also ignoring leading zeros of its numbers (sub-1 or
%   1 for sub-01 or 001: Excel drops them), is replaced by the recording's
%   ID, so that every stage finds it. A match that is not one to one (two
%   rows for one recording, two recordings for one row) is not made. Without
%   S.data (not found) IDs are matched exactly only.
%   [T, matched] = read_subject_table(S) also returns those replacements
%   (table: table_id, recording_id), which check_subject_table and
%   check_sample_sizes print.
%   Every stage reads the table through this function, and joins it to the
%   EEG features by ID (never by row order). With individual bands
%   (S.individual_bands), the individual alpha frequencies of stage 1
%   (results/Connectivity/iaf.csv) are joined too: columns IAF, IAF_PAF,
%   IAF_CoG, IAF_IRASA, IAF_specparam, ... (usable as covariates); with
%   the spectral profile of the statistics (S.spectral_profile), its measures
%   (results/Connectivity/spectral_profile.csv: alpha_pw_post, exponent_all,
%   ...; spectral_profile_features), once the statistics (or the PSD statistics,
%   for their covariates) has computed them. A covariate of S that needs
%   one of these files while it does not exist yet gives a warning, once
%   per file and MATLAB session (and covariate_matrix an error naming the
%   file if it is used).
%   A joined column never replaces a column of the subject table.
%   With the groups taken from the start of the IDs (S.group_column = ''),
%   a column group holds each subject's group (S.groups), for the
%   conditions of S.comparisons and S.correlations.
if isempty(S.subs_info) || ~isfile(S.subs_info)
    error('read_subject_table:file', 'Subject table not found: %s (S.subs_info)', S.subs_info);
end
[~, ~, ext] = fileparts(S.subs_info);
% the ID column is read as text, so IDs such as 001 keep their zeros
switch lower(ext)
    case {'.xlsx', '.xls'}
        opts = detectImportOptions(S.subs_info, 'VariableNamingRule', 'preserve', 'TextType', 'string');
        if ismember(S.id_column, opts.VariableNames), opts = setvartype(opts, S.id_column, 'string'); end
        T = readtable(S.subs_info, opts);
    otherwise
        if strcmpi(ext, '.tsv')
            delim = {'\t'};
        else
            o = detectImportOptions(S.subs_info, 'FileType', 'text');
            delim = cellstr(o.Delimiter);
        end
        opts = detectImportOptions(S.subs_info, 'FileType', 'text', 'Delimiter', delim, ...
            'DecimalSeparator', decimal_separator(S.subs_info, delim), 'VariableNamingRule', 'preserve', ...
            'TextType', 'string');
        opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'n/a', 'NA', 'NaN', ''});
        if ismember(S.id_column, opts.VariableNames), opts = setvartype(opts, S.id_column, 'string'); end
        T = readtable(S.subs_info, opts);
end
T = numeric_text_columns(T, S.id_column);
text_covariates(T, S, opts);
missing_codes(T, S);
if ~ismember(S.id_column, T.Properties.VariableNames)
    error('read_subject_table:id', 'Column %s (S.id_column) not in %s; columns: %s', S.id_column, ...
        S.subs_info, strjoin(T.Properties.VariableNames, ', '));
end
T.(S.id_column) = strtrim(string(T.(S.id_column)));
d = unique(T.(S.id_column)(duplicated_ids(T.(S.id_column))));
if ~isempty(d)
    error('read_subject_table:duplicates', 'Duplicate IDs in %s: %s', S.subs_info, strjoin(d, ', '));
end
[T, matched] = match_recording_ids(T, S);
% groups from the start of the IDs (S.group_column = ''): a column group
% holding them, as get_subject_groups assigns them, so that S.comparisons
% and the correlation sets can name the groups ({'group', '==', 'P'})
if isfield(S, 'group_column') && isempty(S.group_column) && isfield(S, 'groups') && ~isempty(S.groups)
    if ismember('group', T.Properties.VariableNames)
        warning('read_subject_table:group', ['S.group_column is empty (groups from the start of the IDs), but %s ' ...
            'has a column group: the comparisons use that column'], S.subs_info);
    else
        g = strings(height(T), 1);  g(:) = missing;
        for k = 1:numel(S.groups)
            g(startsWith(T.(S.id_column), S.groups{k}, 'IgnoreCase', true) & ismissing(g)) = string(S.groups{k});
        end
        T.group = g;
    end
end
% individual alpha frequencies of stage 1 (only with bands anchored on
% them: IAF, IAF_PAF, IAF_CoG, ...) and the spectral profile of the
% statistics (alpha_pw_post, exponent_all, ...): joined by ID, so that the statistics can
% test them and use them as covariates
if isfield(S, 'results')
    if isfield(S, 'iaf_used') && S.iaf_used
        T = join_measures(T, S, fullfile(S.results, 'Connectivity', 'iaf.csv'), ...
            {'IAF', 'IAF_PAF', 'IAF_CoG', 'IAF_IRASA', 'IAF_specparam'}, 'stage 1 (run_all)');
    end
    if isfield(S, 'spectral_profile') && S.spectral_profile
        m = {};
        if ~isempty(which('spectral_profile_measures')), m = spectral_profile_measures(); end   % Analysis/fun
        T = join_measures(T, S, fullfile(S.results, 'Connectivity', 'spectral_profile.csv'), m, ...
            'the statistics, or the PSD statistics when a covariate needs it');
    end
end
end

function T = join_measures(T, S, f, measures, maker)
% the columns of a per-recording table of the pipeline (recording, ID,
% measures), joined by ID; missing for an ID it lacks. A missing file that
% covariates of S need: a warning (once per file and session)
persistent warned
if isempty(warned), warned = strings(0); end
if ~isfile(f)
    need = setdiff(intersect(all_covariates(S), measures), T.Properties.VariableNames);
    if ~isempty(need) && ~ismember(string(f), warned)
        warning('read_subject_table:not_joined', ['%s not found (yet): covariate(s) %s, which it holds, are not ' ...
            'in the subject table until %s writes it'], f, strjoin(need, ', '), maker);
        warned(end + 1) = string(f);
    end
    return
end
o = detectImportOptions(f);
o = setvartype(o, intersect({'recording', 'ID', 'IAF_source', 'IAF_flag'}, o.VariableNames), 'string');
I = readtable(f, o);
cols = setdiff(I.Properties.VariableNames, [{'recording', 'ID'}, T.Properties.VariableNames], 'stable');
[ok, row] = ismember(T.(S.id_column), I.ID);
for c = cols
    if isnumeric(I.(c{1})), v = nan(height(T), 1); else, v = strings(height(T), 1); end
    v(ok) = I.(c{1})(row(ok));
    T.(c{1}) = v;
end
end

function [T, matched] = match_recording_ids(T, S)
% IDs of the table without an exact recording, matched one to one ignoring
% case (pass 1), then also leading zeros of numbers (pass 2), replaced by
% the recording's ID
matched = table(strings(0, 1), strings(0, 1), 'VariableNames', {'table_id', 'recording_id'});
if ~isfield(S, 'data') || isempty(S.data) || ~isfolder(S.data) || ~isfield(S, 'file_pattern') || ...
        isempty(which('get_subs'))
    return
end
try
    subs = get_subs(S.data, S.file_pattern);
catch
    return                                          % no recordings there: exact matching only
end
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
rec = unique(subject_ids(S, stems));
tid = T.(S.id_column);
free_r = rec(~ismember(rec, tid));                 % recordings without an exact row
free_t = find(~ismember(tid, rec) & ~ismissing(tid) & tid ~= "");
keys = {@(x) lower(strtrim(x)), @(x) regexprep(lower(strtrim(x)), '(?<![0-9])0+(?=[0-9])', '')};
for pass = 1:2
    kr = keys{pass}(free_r);  kt = keys{pass}(tid(free_t));
    done_t = false(size(free_t));  done_r = false(size(free_r));
    for i = 1:numel(free_t)
        h = find(kr == kt(i));
        if isscalar(h) && nnz(kt == kr(h)) == 1    % one to one
            matched(end + 1, :) = {tid(free_t(i)), free_r(h)}; %#ok<AGROW>
            tid(free_t(i)) = free_r(h);
            done_t(i) = true;  done_r(h) = true;
        end
    end
    free_t = free_t(~done_t);  free_r = free_r(~done_r);
end
T.(S.id_column) = tid;
end

function c = all_covariates(S)
% the covariates named anywhere in S (comparisons, correlation sets, covariate sets)
c = {};
if isfield(S, 'covariates'), c = cellstr(S.covariates); end
for f = {'correlations', 'covariate_sets'}
    if isfield(S, f{1}) && isstruct(S.(f{1})) && isfield(S.(f{1}), 'covariates')
        for k = 1:numel(S.(f{1})), c = [c, cellstr(S.(f{1})(k).covariates)]; end %#ok<AGROW>
    end
end
c = unique(c);
end

function d = decimal_separator(file, delimiter)
% ',' when the file is not comma-separated and holds numbers with a decimal
% comma (European Excel export): readtable then reads "286,5" correctly;
% with '.', a column holding "400" and "286,5" gets NaN for "286,5", silently
d = '.';
delimiter = cellstr(delimiter);
if any(strcmp(delimiter, ',')), return; end
txt = fileread(file);
sep = regexptranslate('escape', strjoin(delimiter, ''));
if ~isempty(regexp(txt, ['(^|[' sep '])\s*[-+]?\d+,\d+\s*(?=[' sep ']|\r?$)'], 'once', 'lineanchors'))
    d = ',';
end
end

function T = numeric_text_columns(T, id_column)
% Columns read as text whose values are all numbers become numeric. This
% covers decimal commas (tables exported by Excel with Italian or other
% European settings: "286,5"), which readtable leaves as text or, with
% str2double, would read as 2865 (comma taken as a thousands separator).
for v = string(T.Properties.VariableNames)
    if v == id_column || isnumeric(T.(v)) || islogical(T.(v)), continue; end
    x = strtrim(string(T.(v)));
    miss = ismissing(x) | x == "" | ismember(lower(x), ["n/a", "na", "nan"]);
    if all(miss), continue; end
    ok = ~cellfun(@isempty, regexp(cellstr(x(~miss)), '^[-+]?(\d+([.,]\d*)?|[.,]\d+)([eE][-+]?\d+)?$', 'once'));
    if all(ok)
        y = nan(size(x));
        y(~miss) = str2double(strrep(x(~miss), ",", "."));
        T.(v) = y;
    end
end
end

function text_covariates(T, S, opts)
% a covariate mostly numeric with values that are not numbers: an error naming them. Read as text
% (readtable's choice when enough values are not numbers), the covariate would enter as one
% indicator per value; read as numbers (its choice otherwise), the others would be missing,
% silently: the file is read again with these columns as text to see them. In a column read as
% numbers, a value without digits ("-", "?") stays a missing value
catg = {};
if isfield(S, 'categorical_covariates'), catg = cellstr(S.categorical_covariates); end
cols = setdiff(intersect(all_covariates(S), intersect(T.Properties.VariableNames, opts.VariableNames)), ...
    [catg, {S.id_column}]);
if isempty(cols), return; end
o = setvartype(opts, cols, 'string');
o.SelectedVariableNames = cols;
R = readtable(S.subs_info, o);
for i = 1:numel(cols)
    c = cols(i);
    v = T.(c{1});
    x = strtrim(string(R.(c{1})));
    use = ~(ismissing(x) | x == "" | ismember(lower(x), ["n/a", "na", "nan"]));
    if isnumeric(v) || islogical(v)
        num = ~isnan(double(v(:))) | ~contains(x, digitsPattern);
    else
        num = ~isnan(str2double(strrep(x, ",", ".")));
    end
    num = num(use);  x = x(use);
    if ~isempty(x) && mean(num) >= 0.5 && ~all(num)
        bad = unique(x(~num), 'stable');
        error('read_subject_table:text_covariate', ['Covariate %s in %s: %d of %d values are numbers, the others ' ...
            'are not (%s). Read as text, it would enter as one indicator per value: correct them (or leave them ' ...
            'empty), or name it in S.categorical_covariates if it codes categories'], c{1}, S.subs_info, ...
            nnz(num), numel(num), strjoin('"' + bad(1:min(5, end)) + '"', ', '));
    end
end
end

function missing_codes(T, S)
% numeric columns with a usual missing-value code far from their other values: a warning, once per file and session
persistent warned
if isempty(warned), warned = strings(0, 1); end
if ismember(string(S.subs_info), warned), return; end
codes = [-9999, -999, -99, -9, 99, 999, 9999];
found = strings(0, 1);
for v = string(T.Properties.VariableNames)
    if v == string(S.id_column) || ~isnumeric(T.(v)), continue; end
    x = double(T.(v)(:));
    x = x(~isnan(x));
    hit = ismember(x, codes);
    if ~any(hit) || nnz(~hit) < 3, continue; end
    o = x(~hit);
    s = 1.4826 * median(abs(o - median(o)));
    if s == 0, s = std(o); end
    if s == 0, continue; end
    far = hit & abs(x - median(o)) > 10 * s;
    for k = unique(x(far))'
        found(end + 1) = sprintf('%s: %g (%d times)', v, k, nnz(x == k)); %#ok<AGROW>
    end
end
if ~isempty(found)
    warned(end + 1) = string(S.subs_info);
    warning('read_subject_table:missing_codes', ['%s: values that look like missing-value codes, far from the ' ...
        'other values of their column: %s. They are used as numbers; write missing values as empty or NA'], ...
        S.subs_info, strjoin(found, ', '));
end
end

function tf = duplicated_ids(ids)
[~, first] = unique(ids, 'stable');
tf = true(size(ids));
tf(first) = false;
end
