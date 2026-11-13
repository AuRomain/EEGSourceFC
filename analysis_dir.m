function d = analysis_dir(S)
%ANALYSIS_DIR Folder of the outputs of the statistics.
%   d = analysis_dir(S): <bns_run_dir(S)>/Analysis, beside the states of
%   the brain-network-state run that the statistics test; in a static-only run
%   (S.dynamic = false, no states), results/Statistics/<S.run_name>
%   ('main' when S.run_name is empty).
if ~isfield(S, 'dynamic') || S.dynamic
    d = fullfile(bns_run_dir(S), 'Analysis');
else
    name = 'main';
    if isfield(S, 'run_name') && ~isempty(S.run_name), name = char(S.run_name); end
    d = fullfile(S.results, 'Statistics', name);
end
end
