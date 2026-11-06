function run_all(S)
%RUN_ALL  Stage 1: source connectivity of every recording.
%   run_all() or run('1- EEG source connectivity/run_all.m') uses
%   pipeline_settings(); run_all(S) uses the settings S (another dataset:
%   run_all(settings_mystudy())). A variable S in the workspace is not used.
% Same flow as the original run_all.m: list the subjects, then run
% run_pipeline on each. Changes (see CHANGES.md):
% - settings from pipeline_settings.m or the argument S (a function since the
%   audit of Oct 2026: a script read any S left in the workspace, V3-5),
%   paths relative to the pipeline (N15);
% - toolboxes added by setup_toolboxes.m (FieldTrip, Brainstorm, ...);
% - recordings run in parallel when the Parallel Computing Toolbox is present;
% - a recording whose outputs already exist is skipped (delete them to redo);
% - a recording that fails is reported and the others go on; the summary
%   of every recording (rate, duration, channels, data rank, cuts, windows)
%   is written to results/Connectivity/stage1_summary.csv (the range of the
%   data ranks is printed), the settings and versions
%   to results/logs/provenance_stage1.mat;
% - the head model (S.anatomy, get_anatomy): the shipped ICBM152, another
%   head for every recording, or one head per subject (S.anatomy_individual);
%   recordings without a head stop the run before anything is computed,
%   unless S.anatomy_fallback (they then use the shipped ICBM152, named in
%   the anatomy column of stage1_summary.csv);
% - optional individual bands (edges anchored on each recording's individual
%   alpha frequency, S.bands and S.iaf): the IAF of every recording in
%   results/Connectivity/iaf.csv (also joined to the subject table as the
%   columns IAF, IAF_PAF, ...: a feature and a possible covariate in stage 3);
% - every recording needs its own subject ID (S.id_regexp on its file name):
%   two recordings with the same ID (a BIDS dataset with several sessions
%   and the default '^([^_]+)', which gives sub-01 for each session) stop
%   the run before anything is computed (stage 2 stopped on them, after
%   stage 1), with the fix: an S.id_regexp that keeps the session, or one
%   session per run.

stage_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(stage_dir));
if nargin < 1, S = []; end
S = setup_toolboxes(stage_settings(S, 'run_all'));

data_path = S.data;                       % 'path' would shadow MATLAB's path() (N15)
path_to_save = fullfile(S.results, 'Connectivity');
if ~isfolder(path_to_save), mkdir(path_to_save); end
if ~isfolder(fullfile(S.results, 'logs')), mkdir(fullfile(S.results, 'logs')); end

% get all subjects in directory
subs = get_subs(data_path, S.file_pattern);
nb_subs = length(subs);
fprintf('Stage 1: %d recordings in %s\n', nb_subs, data_path);
check_ids(S, subs);

% head models: every recording's, before anything runs (get_anatomy)
if ~S.anatomy_individual
    heads = {get_anatomy(S)};
else
    heads = cell(1, nb_subs);
    problems = strings(0, 1);
    for s = 1:nb_subs
        [heads{s}, problems(end + 1)] = get_anatomy(S, subs{s}); %#ok<AGROW>
    end
    problems = problems(strlength(problems) > 0);
    if ~isempty(problems)
        error('run_all:anatomy', ['%d of %d recordings without a head model (S.anatomy = %s):\n  %s\n' ...
            'Make them with tools/make_anatomy.m, or set S.anatomy_fallback = true to use the shipped ICBM152 ' ...
            'for these subjects'], numel(problems), nb_subs, S.anatomy, strjoin(problems, '\n  '));
    end
    fallback = cellfun(@(A) A.fallback, heads);
    if any(fallback)
        warning('run_all:anatomy_fallback', '%d of %d recordings use the shipped ICBM152 (no folder in %s): %s', ...
            nnz(fallback), nb_subs, S.anatomy, strjoin(subs(fallback), ', '));
    end
    [~, first] = unique(cellfun(@(A) A.dir, heads, 'UniformOutput', false), 'stable');
    heads = heads(first);
end
if ~heads{1}.shipped || numel(heads) > 1
    fprintf('Head models: %d (%s)\n', numel(heads), strjoin(cellfun(@(A) A.name, heads, 'UniformOutput', false), ', '));
end

% leadfields: loaded (or computed, for a new cap or head) once here, not
% by every parallel worker at the same time
for h = 1:numel(heads)
    for m = 1:numel(S.montages)
        for a = 1:numel(S.atlases)
            get_ftHeadmodel_and_ftLeadfield(S.montages{m}, S.atlases{a}, S, heads{h});
        end
    end
end

% loop over all subjects in directory
t0 = tic;
qc = cell(nb_subs, 1);
use_pool = S.parallel && nb_subs > 12 && ~isempty(ver('parallel'));   % a pool takes about a minute to start
if use_pool
    if isempty(gcp('nocreate')), parpool('Processes', min(S.max_workers, nb_subs)); end
    parfor s = 1:nb_subs
        qc{s} = run_pipeline(subs{s}, data_path, path_to_save, S);
    end
else
    for s = 1:nb_subs
        qc{s} = run_pipeline(subs{s}, data_path, path_to_save, S);
    end
end
qc = struct2table(vertcat(qc{:}), 'AsArray', true);
writetable(qc, fullfile(path_to_save, 'stage1_summary.csv'));
ok = ~isnan(qc.data_rank);
if any(ok)              % independent signals in the channels used (run_pipeline: data_rank)
    fprintf('Data rank %d-%d (median %g) of %d-%d channels; below S.min_channels (%d): %d (stage1_summary.csv)\n', ...
        min(qc.data_rank(ok)), max(qc.data_rank(ok)), median(qc.data_rank(ok)), min(qc.n_channels(ok)), ...
        max(qc.n_channels(ok)), S.min_channels, nnz(qc.data_rank(ok) < S.min_channels));
end
if S.iaf_used          % individual alpha frequencies (bands anchored on them): iaf.csv, by subject ID
    [~, stems] = cellfun(@fileparts, cellstr(qc.file), 'UniformOutput', false);
    I = table(string(stems), subject_ids(S, stems), qc.iaf, qc.iaf_paf, qc.iaf_cog, qc.iaf_irasa, ...
        qc.iaf_specparam, qc.iaf_channels, string(qc.iaf_source), string(qc.iaf_flag), ...
        'VariableNames', {'recording', 'ID', 'IAF', 'IAF_PAF', 'IAF_CoG', 'IAF_IRASA', 'IAF_specparam', ...
        'IAF_channels', 'IAF_source', 'IAF_flag'});
    writetable(I, fullfile(path_to_save, 'iaf.csv'));
    fprintf('Individual alpha frequency: %d of %d recordings (%.2f-%.2f Hz); iaf.csv\n', nnz(~isnan(I.IAF)), ...
        height(I), min(I.IAF), max(I.IAF));
end
provenance = struct('stage', 1, 'date', char(datetime('now')), 'settings', S, 'summary', qc);
save(fullfile(S.results, 'logs', 'provenance_stage1.mat'), 'provenance');
failed = qc.file(startsWith(qc.status, "failed"));
fprintf('Stage 1 finished in %.1f min; outputs in %s\n', toc(t0) / 60, path_to_save);
if ~isempty(failed)
    warning('run_all:failed', '%d recording(s) failed: %s (see stage1_summary.csv)', numel(failed), strjoin(failed, ', '));
end

function check_ids(S, subs)
% every recording its own subject ID: stages 2 and 3, iaf.csv and the
% subject table identify the recordings by it (get_subject_groups)
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
ids = subject_ids(S, stems);
[u, ~, j] = unique(ids);
twice = find(accumarray(j(:), 1) > 1);
if isempty(twice), return; end
list = arrayfun(@(k) sprintf('%s (%s)', u(k), strjoin(string(stems(j == k)), ', ')), twice, 'UniformOutput', false);
error('run_all:ids', ['%d subject ID(s) given to several recordings by S.id_regexp = ''%s'': %s.\n' ...
    'Each recording needs its own ID. BIDS data with several sessions: keep the session in the ID, ' ...
    'S.id_regexp = ''^(sub-[^_]+_ses-[^_]+)'' (sub-01_ses-1, ...; the subject table then needs one row per ' ...
    'session), or analyse one session per run (S.file_pattern = ''sub-*/ses-1/eeg/*_eeg.set'')'], ...
    numel(twice), S.id_regexp, strjoin(list, '; '));
