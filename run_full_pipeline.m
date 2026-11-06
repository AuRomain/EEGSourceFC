function run_full_pipeline(S, varargin)
%RUN_FULL_PIPELINE  The whole pipeline, stage after stage.
%   run_full_pipeline() uses pipeline_settings(); run_full_pipeline(S) uses
%   the settings S of a study, e.g. run_full_pipeline(settings_mystudy()).
%   A variable S left in the workspace is not used.
%   The stages, each the same function as when it is called alone:
%     1  run_all                             connectivity of every recording
%     1b run_psd_analysis                    power spectra of the channels and
%                                            regions, statistics and maps (only
%                                            with S.psd.run = true, the default;
%                                            it reads the recordings itself)
%     2  compute_modkmeans_states_meta_crit  brain network states (only with
%                                            S.dynamic = true; a static-only
%                                            run skips it)
%     3  run_statistics                      statistics and figures
%   What is up to date is reused (the recordings stage 1 has already done,
%   stage 2's cached clustering), so the same command, run again after an
%   interruption or a change of settings, redoes only what is missing or
%   changed.
%   run_full_pipeline(S, 'from', 2) starts at stage 2 (stage 1 done and
%   checked); run_full_pipeline(S, 'to', 1) stops after stage 1;
%   'from' and 'to' take 1, '1b', 2 or 3.
%   The settings are checked and the toolboxes found before anything runs;
%   when a stage that uses the subject table will run (1b, 2 or 3), the
%   table is checked against the recordings too (check_subject_table:
%   recordings without a row, rows without a recording, the group and
%   comparison columns, the sample sizes), and a table that fits no
%   recording stops the run before stage 1.
%   An error stops the run, since each stage needs the outputs of the one
%   before; the message names the stage and how to resume from it. Stages
%   2 and 3 do not start without a recording with stage 1 outputs. A
%   recording that fails in stage 1 is left out of the next stages, as when
%   run_all is called alone, and named at the end. At the end: the time of
%   each stage and the files to check before reading the results (README,
%   section 4).

root = fileparts(mfilename('fullpath'));
addpath(root);
if nargin < 1, S = []; end
S = stage_settings(S, 'run_full_pipeline');
p = inputParser;
p.FunctionName = 'run_full_pipeline';
names = {'1', '1b', '2', '3'};             % the stages in their order
p.addParameter('from', 1, @(x) ~isempty(stage_index(x)));
p.addParameter('to', 3, @(x) ~isempty(stage_index(x)));
p.parse(varargin{:});
from = stage_index(p.Results.from);
to = stage_index(p.Results.to);
if from > to
    error('run_full_pipeline:stages', '''from'' (stage %s) comes after ''to'' (stage %s)', names{from}, names{to});
end
Sc = setup_toolboxes(S);                % the settings checked, the stages on the path; each stage gets S as given
if to >= 2                              % a stage that uses the subject table will run: checked first
    fprintf('\n========== Subject table ==========\n');
    check_subject_table(Sc);
end
stages = {'run_all', 'connectivity of every recording'
    'run_psd_analysis', 'PSD analysis'
    'compute_modkmeans_states_meta_crit', 'brain network states'
    'run_statistics', 'statistics'};
if ~Sc.dynamic, stages{4, 2} = 'statistics (static connectivity only)'; end
minutes = nan(1, 4);
t0 = tic;
for k = from:to
    if k == 2 && ~Sc.psd.run
        fprintf('\n========== Stage 1b (PSD analysis): skipped, S.psd.run = false ==========\n');
        continue
    end
    if k == 3 && ~Sc.dynamic
        fprintf('\n========== Stage 2 of 3 (brain network states): skipped, S.dynamic = false ==========\n');
        continue
    end
    if k >= 3 && stage1_status(Sc) == 0             % stage 3 alone would only log its families as failed
        error('run_full_pipeline:stage1', ['No recording has stage 1 outputs in %s (no stage1_summary.csv, ' ...
            'or every recording failed: see it): stage %s needs them. Stage 1: run_full_pipeline(S, ''to'', 1)'], ...
            fullfile(Sc.results, 'Connectivity'), names{k});
    end
    if k == 2
        fprintf('\n========== Stage 1b: %s (%s) ==========\n', stages{k, 2}, stages{k, 1});
    else
        fprintf('\n========== Stage %s of 3: %s (%s) ==========\n', names{k}, stages{k, 2}, stages{k, 1});
    end
    t = tic;
    try
        feval(stages{k, 1}, S);
    catch err
        resume = names{k};
        if k == 2, resume = ['''' resume '''']; end
        e = MException('run_full_pipeline:stage', ['Stage %s (%s) stopped after %.1f min: %s\n' ...
            'After fixing it, resume with run_full_pipeline(S, ''from'', %s)'], names{k}, stages{k, 1}, ...
            toc(t) / 60, err.message, resume);
        throw(addCause(e, err));
    end
    minutes(k) = toc(t) / 60;
end

fprintf('\n========== Done in %.1f min ==========\n', toc(t0) / 60);
for k = find(~isnan(minutes))
    fprintf('Stage %s (%s): %.1f min\n', names{k}, stages{k, 2}, minutes(k));
end
fprintf('Check before reading the results:\n');
[done, n, failed, summary1] = stage1_status(Sc);
if ~isnan(n)
    fprintf('  %s: %d recordings', summary1, n);
    if isempty(failed)
        fprintf(', all made or reused\n');
    else
        fprintf(', %d made or reused; %d failed, left out: %s\n', done, numel(failed), strjoin(failed, ', '));
    end
end
if ~isnan(minutes(2))
    run = 'main';
    if isfield(Sc, 'run_name') && ~isempty(Sc.run_name), run = char(Sc.run_name); end
    fprintf(['  %s: psd_log.txt (channels left out, failures), recordings.csv (fit quality, IAF), ' ...
        'Statistical analysis/summary.csv and clusters.csv\n'], fullfile(Sc.results, 'PSD', run));
end
if Sc.dynamic && to >= 3
    fprintf(['  %s: k_criteria_*.csv (stability of the number of states), null_test_*.csv (states against ' ...
        'surrogates), continuum_test_*.csv (clusters or a continuum)\n'], fullfile(bns_run_dir(Sc), 'Results'));
end
if to == 4
    fprintf('  %s: analysis_log.txt (recordings per comparison), summary.csv, families.csv\n', ...
        fullfile(analysis_dir(Sc), 'Statistical analysis'));
end
end


function k = stage_index(x)
% the position of a stage (1, '1b', 2, 3) in the order 1, 1b, 2, 3; [] if none
k = [];
if isnumeric(x) && isscalar(x) && ismember(x, 1:3)
    k = x + (x > 1);
elseif (ischar(x) || (isstring(x) && isscalar(x))) && any(strcmpi(strtrim(char(x)), {'1', '1b', '2', '3'}))
    k = find(strcmpi(strtrim(char(x)), {'1', '1b', '2', '3'}));
end
end

function [done, n, failed, file] = stage1_status(S)
% the recordings of stage1_summary.csv with outputs (made now or reused), all
% of them (NaN without the file) and those that failed
file = fullfile(S.results, 'Connectivity', 'stage1_summary.csv');
done = 0;  n = NaN;  failed = strings(0, 1);
if ~isfile(file), return; end
Q = readtable(file, 'Delimiter', ',', 'TextType', 'string');
if ~all(ismember({'file', 'status'}, Q.Properties.VariableNames)), return; end
n = height(Q);
done = nnz(Q.status == "ok" | Q.status == "skipped (up to date)");
failed = Q.file(startsWith(Q.status, "failed"));
end
