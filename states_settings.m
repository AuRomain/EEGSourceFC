function V = states_settings(S)
%STATES_SETTINGS The settings that shape the results of the brain network states.
%   V = states_settings(S): the fields of S used by the states (stage 3; their
%   features), for the name of the run folder (bns_run_dir) and the run
%   fingerprint stored with every output of the states.
f = {'state_conn', 'NCs', 'k_selection', 'stability', 'reps', 'max_iterations', 'clustering', ...
     'clustering_engine', 'null_test', 'minTime', 'minTime_s', 'energy_entropy', 'seed', ...
     'network_table', 'language_file', 'groups', 'group_column', 'group_files'};
V = struct();
for i = 1:numel(f)
    if isfield(S, f{i}), V.(f{i}) = S.(f{i}); end
end
if isfield(S, 'state_conn') && strcmpi(S.state_conn, 'dwpli')   % dwPLI states: their own smoothing
    V = rmfield(V, intersect(fieldnames(V), {'minTime', 'minTime_s'}));   % only for them, so that
    V.minTime_dwpli = 0;                                              % the run folders of the
    if isfield(S, 'minTime_dwpli'), V.minTime_dwpli = S.minTime_dwpli; end   % other measures do not
end                                                                   % depend on it
if isfield(S, 'k_fixed') && ~isempty(S.k_fixed)      % only when set (unset: the folder name does not change)
    V.k_fixed = S.k_fixed;
end
seq = {'sequence_features', true; 'sequence_windows', 'shortest'; 'sequence_surrogates', 200};
for i = 1:size(seq, 1)                               % the settings of the sequence features: only when
    f = seq{i, 1};                                   % not the default, so that the default does not
    if isfield(S, f) && ~isequal(S.(f), seq{i, 2})   % change the folder name
        V.(f) = S.(f);
    end
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
