function [ids, group, T] = get_subject_groups(S, names)
%GET_SUBJECT_GROUPS Subject ID and group of each recording.
%   [ids, group, T] = get_subject_groups(S, names) with names the stage 1
%   folders (file names without extension, e.g. SUB_001_eo_ready).
%   ids    subject ID: first token of S.id_regexp on the name, or the
%          whole name when it does not match
%   group  index into S.groups (0 = in no group: left out of the group analyses)
%   T      the subject table rows of ids (empty when there is no table)
%   The group comes from the subject table (column S.group_column) when one
%   is given, otherwise from the start of the ID.
names = string(names(:));
ids = subject_ids(S, names);
if numel(unique(ids)) < numel(ids)
    error('get_subject_groups:ids', ['Several recordings give the same subject ID (S.id_regexp = ''%s''): %s. ' ...
        'BIDS data with several sessions: S.id_regexp = ''^(sub-[^_]+_ses-[^_]+)'' (one subject-table row ' ...
        'per session), or one session per run'], S.id_regexp, strjoin(unique(ids(duplicated(ids))), ', '));
end
group = zeros(numel(ids), 1);
T = table();
if ~isempty(S.group_column) && ~isempty(S.subs_info)
    T = read_subject_table(S);
    [ok, row] = ismember(ids, T.(S.id_column));
    if ~all(ok)
        warning('get_subject_groups:table', 'No row in the subject table for: %s (left out)', strjoin(ids(~ok), ', '));
    end
    g = strings(numel(ids), 1);
    g(ok) = string(T.(S.group_column)(row(ok)));
    [~, group] = ismember(lower(strtrim(g)), lower(strtrim(string(S.groups))));   % case ignored
    T = T(row(ok), :);
else
    for k = 1:numel(S.groups)
        group(startsWith(ids, S.groups{k}, 'IgnoreCase', true) & group == 0) = k;
    end
end
if any(group == 0)
    warning('get_subject_groups:none', 'Recordings in none of the groups %s (left out): %s', ...
        strjoin(S.groups, ', '), strjoin(names(group == 0), ', '));
end
end

function tf = duplicated(x)
[~, first] = unique(x, 'stable');
tf = true(size(x));
tf(first) = false;
end
