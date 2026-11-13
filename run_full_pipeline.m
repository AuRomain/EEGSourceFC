function run_full_pipeline(S, varargin)
%RUN_FULL_PIPELINE  The whole pipeline, stage after stage.
%   run_full_pipeline() uses pipeline_settings(); run_full_pipeline(S) uses
%   the settings S of a study, e.g. run_full_pipeline(settings_mystudy()).
%   A variable S left in the workspace is not used.
%   The stages, each the same function as when it is called alone, numbered
%   as the pipeline's folders (without the brain network states, the
%   statistics are stage 3):
%     1  run_all                             connectivity of every recording
%     2  run_psd_analysis                    power spectra of the channels and
%                                            regions, statistics and maps (only
%                                            with S.psd.run = true, the default;
%                                            it reads the recordings itself)
%     3  compute_modkmeans_states_meta_crit  brain network states (only with
%                                            S.dynamic = true; a static-only
%                                            run skips it)
%     4  run_statistics                      statistics and figures
%   What is up to date is reused (the recordings stage 1 has already done,
%   the cached clustering of the states), so the same command, run again
%   after an interruption or a change of settings, redoes only what is
%   missing or changed.
%   run_full_pipeline(S, 'from', 2) starts at stage 2 (stage 1 done and
%   checked); run_full_pipeline(S, 'to', 1) stops after stage 1; 'from'
%   and 'to' take the number of a stage; '1b', the PSD analysis's number
%   before version 1.2.0, is taken as 2, with a warning.
%   The settings are checked and the toolboxes found before anything runs;
%   when a stage that uses the subject table will run (any after stage 1),
%   the table is checked against the recordings too (check_subject_table:
%   recordings without a row, rows without a recording, the group and
%   comparison columns, the sample sizes), and a table that fits no
%   recording stops the run before stage 1.
%   An error stops the run, since each stage needs the outputs of the one
%   before; the message names the stage and how to resume from it. Inside
%   the statistics a failing family, NBS, reliability or figure does not
%   stop the others (unless S.stop_on_error = true), nor inside the PSD
%   analysis a failing statistics, figure, region check or vertex map step:
%   the run then ends "Done with errors" and lists them (failures.csv of
%   the statistics and of the PSD analysis). When only figures failed, it
%   ends "every result written, but N figure step(s) failed" instead
%   (failure_summary).
%   The states and the statistics do not start without a recording with
%   stage 1 outputs. A recording that fails in stage 1 is left out of the
%   next stages, as when run_all is called alone, and named at the end. At
%   the end: the time of each stage and the files to check before reading
%   the results (README, section 4).

root = fileparts(mfilename('fullpath'));
addpath(root);
if nargin < 1, S = []; end
S = stage_settings(S, 'run_full_pipeline');
% the stages of this copy of the pipeline: function, what it does, folder
stages = {'run_all', 'connectivity of every recording', 'EEG source connectivity'
    'run_psd_analysis', 'PSD analysis', 'PSD analysis'
    'compute_modkmeans_states_meta_crit', 'brain network states', 'Dynamic connectivity'
    'run_statistics', 'statistics', 'Analysis'};
num = nan(1, 4);                           % the number of each stage (NaN: not in this copy)
for k = 1:4
    [~, n] = stage_folder(root, stages{k, 3});
    if ~isempty(n), num(k) = n; end
end
last = max(num);
p = inputParser;
p.FunctionName = 'run_full_pipeline';
p.addParameter('from', 1, @(x) ~isempty(stage_index(x, num)));
p.addParameter('to', last, @(x) ~isempty(stage_index(x, num)));
p.parse(varargin{:});
from = stage_index(p.Results.from, num);
to = stage_index(p.Results.to, num);
if any(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(strtrim(char(x)), '1b'), {p.Results.from, p.Results.to}))
    warning('run_full_pipeline:stage_1b', ['Stage ''1b'' (the PSD analysis) is stage 2 since 1.2.0, ' ...
        'and the stages after it are numbered on: taken as 2']);
end
if from > to
    error('run_full_pipeline:stages', '''from'' (stage %d) comes after ''to'' (stage %d)', num(from), num(to));
end
Sc = setup_toolboxes(S);                % the settings checked, the stages on the path; each stage gets S as given
if to >= 2                              % a stage that uses the subject table will run: checked first
    fprintf('\n========== Subject table ==========\n');
    check_subject_table(Sc);
end
if ~Sc.dynamic, stages{4, 2} = 'statistics (static connectivity only)'; end
minutes = nan(1, 4);
failures = strings(0, 2);                % the steps of the PSD analysis and the statistics that failed
fail_files = {};                         % their failures.csv
t0 = tic;
for k = from:to
    if isnan(num(k)), continue; end                 % not in this pipeline (no states)
    if k == 2 && ~Sc.psd.run
        fprintf('\n========== Stage %d of %d (PSD analysis): skipped, S.psd.run = false ==========\n', num(k), last);
        continue
    end
    if k == 3 && ~Sc.dynamic
        fprintf('\n========== Stage %d of %d (brain network states): skipped, S.dynamic = false ==========\n', ...
            num(k), last);
        continue
    end
    if k >= 3 && stage1_status(Sc) == 0             % the statistics alone would only log their families as failed
        error('run_full_pipeline:stage1', ['No recording has stage 1 outputs in %s (no stage1_summary.csv, ' ...
            'or every recording failed: see it): stage %d needs them. Stage 1: run_full_pipeline(S, ''to'', 1)'], ...
            fullfile(Sc.results, 'Connectivity'), num(k));
    end
    fprintf('\n========== Stage %d of %d: %s (%s) ==========\n', num(k), last, stages{k, 2}, stages{k, 1});
    t = tic;
    try
        if k == 2
            f = run_psd_analysis(S);
            if ~isempty(f)
                failures = [failures; "PSD: " + f(:, 1), f(:, 2)]; %#ok<AGROW>
                fail_files{end + 1} = fullfile(Sc.results, 'PSD', psd_run(Sc), 'failures.csv'); %#ok<AGROW>
            end
        elseif k == 4
            f = run_statistics(S);
            if ~isempty(f)
                failures = [failures; f]; %#ok<AGROW>
                fail_files{end + 1} = fullfile(analysis_dir(Sc), 'Statistical analysis', 'failures.csv'); %#ok<AGROW>
            end
        else
            feval(stages{k, 1}, S);
        end
    catch err
        e = MException('run_full_pipeline:stage', ['Stage %d (%s) stopped after %.1f min: %s\n' ...
            'After fixing it, resume with run_full_pipeline(S, ''from'', %d)'], num(k), stages{k, 1}, ...
            toc(t) / 60, err.message, num(k));
        throw(addCause(e, err));
    end
    minutes(k) = toc(t) / 60;
end

if isempty(failures)
    fprintf('\n========== Done in %.1f min ==========\n', toc(t0) / 60);
else
    [head, lines] = failure_summary(failures, toc(t0) / 60);   % figure-only failures apart
    fprintf(2, '\n========== %s ==========\n', head);
    fprintf(2, '  %s\n', lines);
    fprintf(2, '  (%s)\n', fail_files{:});
end
for k = find(~isnan(minutes))
    fprintf('Stage %d (%s): %.1f min\n', num(k), stages{k, 2}, minutes(k));
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
    fprintf(['  %s: psd_log.txt (channels left out, failures), recordings.csv (fit quality, IAF), ' ...
        'Statistical analysis/summary.csv and clusters.csv\n'], fullfile(Sc.results, 'PSD', psd_run(Sc)));
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


function run = psd_run(S)
% the PSD analysis's run folder name (results/PSD/<run>)
run = 'main';
if isfield(S, 'run_name') && ~isempty(S.run_name), run = char(S.run_name); end
end

function k = stage_index(x, num)
% the position (1 connectivity, 2 PSD, 3 states, 4 statistics) of the stage
% numbered x (a number, or its text) in this copy of the pipeline, whose
% numbers are num (NaN: not in this copy); '1b', the PSD analysis until
% 1.2.0, is 2; [] if none
k = [];
if (ischar(x) || (isstring(x) && isscalar(x)))
    x = strtrim(char(x));
    if strcmpi(x, '1b'), k = 2; return; end
    x = str2double(x);
end
if isnumeric(x) && isscalar(x) && ~isnan(x)
    k = find(num == x, 1);
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
