function V = stage2_settings(S)
%STAGE2_SETTINGS The settings that shape the results of stage 2.
%   V = stage2_settings(S): the fields of S used by stage 2 (states and their
%   features), for the name of the run folder (bns_run_dir) and the run
%   fingerprint stored with every stage 2 output (V3-1, V3-6).
f = {'state_conn', 'NCs', 'k_selection', 'stability', 'reps', 'max_iterations', 'clustering', ...
     'clustering_engine', 'null_test', 'minTime', 'minTime_s', 'energy_entropy', 'seed', ...
     'network_table', 'language_file', 'groups', 'group_column', 'group_files'};
V = struct();
for i = 1:numel(f)
    if isfield(S, f{i}), V.(f{i}) = S.(f{i}); end
end
if isfield(S, 'k_fixed') && ~isempty(S.k_fixed)      % only when set, so earlier runs keep their folders
    V.k_fixed = S.k_fixed;
end
if isfield(V, 'group_files') && isempty(V.group_files) && isfield(S, 'groups')
    V.group_files = default_group_files(S.groups);  % as check_settings fills it, so that the settings
end                                                  % before and after check_settings name the same run folder
for g = {'network_table', 'language_file'}           % tables named by path: by content
    if isfield(V, g{1}) && ~isempty(V.(g{1}))
        V.(g{1}) = file_fingerprint(V.(g{1}));
    end
end
end
