function [subs, ids, group, files] = stage1_outputs(S, inv, band)
%STAGE1_OUTPUTS The stage 1 outputs that stages 2 and 3 use, checked.
%   [subs, ids, group, files] = stage1_outputs(S, inv, band)
%   subs   folder of each recording (file name without extension)
%   ids    subject IDs; group: index into S.groups
%   files  cmats_<inv>_<montage>_<atlas>_<band>.mat of each recording
%          (static_cmats_*.mat in a static-only run, S.dynamic = false)
%   Recordings are those of S.data (S.file_pattern), not every folder of
%   results/Connectivity/cmats: a folder left by a recording that failed, or
%   by a file since removed, is not used. Recordings without outputs (stage
%   1 failed) are left out with a warning; recordings in none of the groups
%   too (get_subject_groups). Stops when an output was made with other
%   stage 1 settings, code or inputs than now, or from a recording that has
%   changed since (settings_key), or when the measure S.state_conn
%   is not at the position S.co of its cmats: the windows would not be the
%   ones the settings describe. Recordings with fewer windows of
%   S.state_conn than max(S.NCs), or for which S.state_conn could not be
%   computed (NaN, info.not_computed: e.g. a symmetric orthogonalisation
%   of more regions than the data rank), are left out with a warning (5 Nov 2026).
names = get_subs(S.data, S.file_pattern);
[~, stems] = cellfun(@fileparts, names, 'UniformOutput', false);
subs = unique(stems);
tag_ma = [file_tag(S.montages{1}) '_' file_tag(S.atlases{1})];
prefix = 'cmats';
if isfield(S, 'dynamic') && ~S.dynamic, prefix = 'static_cmats'; end    % static only: no windows
files = fullfile(S.results, 'Connectivity', 'cmats', subs, sprintf('%s_%s_%s_%s.mat', prefix, inv, tag_ma, band));
have = cellfun(@isfile, files);
if any(~have)
    warning('stage1_outputs:missing', ['No stage 1 output (failed, not run, or no individual alpha ' ...
        'frequency for a band anchored on it) for %d recording(s), left out: %s'], nnz(~have), strjoin(subs(~have), ', '));
end
subs = subs(have);  files = files(have);
[ids, group] = get_subject_groups(S, subs);
keep = group > 0;
subs = subs(keep);  ids = ids(keep);  group = group(keep);  files = files(keep);
if isempty(subs)
    error('stage1_outputs:none', 'No recording with stage 1 outputs in the groups %s', strjoin(S.groups, ', '));
end
S1 = setfield(setfield(S, 'montages', S.montages(1)), 'atlases', S.atlases(1)); %#ok<SFLD>
[~, ~, accepted] = settings_key(S1, '', inv, band);   % the current key and those of equivalent earlier code
n_windows = inf(size(files));
not_done = false(size(files));
for i = 1:numel(files)
    w = load(files{i}, 'info');
    if ~isfield(w.info, 'settings_key') || ~any(strcmp(w.info.settings_key, accepted))
        error('stage1_outputs:settings', ['%s was made with other stage 1 settings, code or inputs than the ' ...
            'current ones: run stage 1 again'], files{i});
    end
    file = fullfile(S.data, names{find(strcmp(stems, subs{i}), 1)});
    [~, input_key] = settings_key(S1, file);
    if ~isfield(w.info, 'input_key') || ~strcmp(w.info.input_key, input_key)
        error('stage1_outputs:input', '%s changed after its stage 1 outputs were made: run stage 1 again', file);
    end
    if numel(w.info.conn) < S.co || ~strcmpi(w.info.conn{S.co}, S.state_conn)
        error('stage1_outputs:conn', '%s: cmats{%d} is %s, not %s', files{i}, S.co, w.info.conn{min(S.co, end)}, S.state_conn);
    end
    if prefix == "cmats" && isfield(w.info, 'windows') && numel(w.info.windows) >= S.co
        n_windows(i) = numel(w.info.windows{S.co}.first_sample);
    end
    if isfield(w.info, 'not_computed') && isfield(w.info.not_computed, lower(S.state_conn))
        not_done(i) = true;                              % NaN for this recording (e.g. data rank below the regions)
    end
end
if any(not_done)
    warning('stage1_outputs:not_computed', ['%s could not be computed for %d recording(s) (NaN; the reason is in ' ...
        'their info.not_computed and stage1_summary.csv), left out of stages 2 and 3: %s'], S.state_conn, ...
        nnz(not_done), strjoin(subs(not_done), ', '));
    subs = subs(~not_done);  ids = ids(~not_done);  group = group(~not_done);  files = files(~not_done);
    n_windows = n_windows(~not_done);
    if isempty(subs)
        error('stage1_outputs:none', 'No recording for which %s could be computed', S.state_conn);
    end
end
% a recording with fewer windows of the clustered measure than the largest
% number of states (none when every piece is shorter than a window, e.g. the
% 6 s AEC windows of 4-6 s pieces) has no state features: left out of
% stages 2 and 3, with a warning, instead of stopping stage 2 late
few = n_windows < max(S.NCs);
if any(few)
    warning('stage1_outputs:few_windows', ['%d recording(s) with fewer than %d windows of %s (the largest ' ...
        'number of states), left out of stages 2 and 3: %s'], nnz(few), max(S.NCs), S.state_conn, ...
        strjoin(compose('%s (%d)', string(subs(few)), n_windows(few)), ', '));
    subs = subs(~few);  ids = ids(~few);  group = group(~few);  files = files(~few);
    if isempty(subs)
        error('stage1_outputs:none', 'No recording with at least %d windows of %s', max(S.NCs), S.state_conn);
    end
end
for g = 1:numel(S.groups)
    if ~any(group == g)
        warning('stage1_outputs:group', 'No recording in group %s', S.groups{g});
    end
end
end
