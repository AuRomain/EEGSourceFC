function Q = nbs_run_q(nbs_dir, alpha)
%NBS_RUN_Q Benjamini-Hochberg q over every component of every NBS analysis of a run.
%   Q = nbs_run_q(nbs_dir, alpha), nbs_dir the run's NBS folder (one
%   folder per band and inverse, each with NBS_<matrix>.csv of
%   compute_NBS_analysis). Each component's p_both (its family-wise p over
%   both contrasts of its comparison) enters one Benjamini-Hochberg
%   adjustment over all the components of the run (every band, inverse,
%   measure, matrix and comparison): q_run, added to every NBS_<matrix>.csv,
%   and Q, one row per component (folder, file and the file's columns).
%   The NBS stays exploratory, outside the adjustment over families of the
%   other analyses: p_both controls the error of one analysis, q_run the
%   false discoveries among the components of all of them.
Q = table();
if ~isfolder(nbs_dir), return; end
d = dir(fullfile(nbs_dir, '**', 'NBS_*.csv'));
T = cell(numel(d), 1);
for k = 1:numel(d)
    T{k} = readtable(fullfile(d(k).folder, d(k).name), 'TextType', 'string');
end
p = cell2mat(cellfun(@(t) pboth(t), T, 'UniformOutput', false));
if isempty(p)
    return
end
q = fdr_bh(p);
at = 0;
for k = 1:numel(d)
    m = height(T{k});
    if ~ismember('p_both', T{k}.Properties.VariableNames), continue; end
    T{k}.q_run = q(at + 1:at + m);
    at = at + m;
    writetable(T{k}, fullfile(d(k).folder, d(k).name));
    [~, folder] = fileparts(d(k).folder);
    Q = [Q; [table(repmat(string(folder), m, 1), repmat(string(d(k).name), m, 1), ...
        'VariableNames', {'folder', 'file'}), T{k}]]; %#ok<AGROW>
end
fprintf(['NBS over the run: %d component(s) in %d analysis file(s); %d with q_run <= %.2f (Benjamini-Hochberg ' ...
    'over all of them, exploratory)\n'], numel(p), numel(d), nnz(q <= alpha), alpha);
end

function p = pboth(t)
p = zeros(0, 1);
if ismember('p_both', t.Properties.VariableNames), p = double(t.p_both(:)); end
end
