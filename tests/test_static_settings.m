function tests = test_static_settings
%TEST_STATIC_SETTINGS The settings of the static-only mode (S.dynamic = false) and their checks.
%   From the pipeline folder:  results = runtests('tests/test_static_settings.m')
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
tc.TestData.root = root;
end

function test_defaults_of_a_static_run(tc)
S = check_settings(static_test_settings());
verifyFalse(tc, S.dynamic);
verifyEqual(tc, S.static_measures, {'wpli', 'aec_orth'});      % every measure of S.conn
verifyEqual(tc, S.static_matrices, {'full'});                  % as static_test_settings sets it
S = static_test_settings();
S.static_matrices = {};                                        % the default: the window mean (6 Nov 2026)
verifyEqual(tc, check_settings_quietly(S).static_matrices, {'window_mean'});
verifyEqual(tc, S.static_families, {'global', 'networks', 'regions', 'power'});
S = static_test_settings();
S.state_conn = 'plv';                                          % not in S.conn: no states anyway
S = check_settings(S);
verifyEqual(tc, S.state_conn, 'wpli');
verifyEqual(tc, S.co, 1);
S = static_test_settings();
S.dynamic = 0;                                                 % 0 / 1 accepted
verifyFalse(tc, check_settings(S).dynamic);
end

function test_the_default_follows_the_folder(tc)
% S.dynamic left out: true where the folder of stage 2 is, false in a
% static-only copy (v_static); true there is refused
S = rmfield(static_test_settings(), 'dynamic');
has_stage2 = isfolder(fullfile(tc.TestData.root, '2- Dynamic connectivity'));
if has_stage2
    S.save_scouts = true;
    verifyTrue(tc, check_settings(S).dynamic);
else
    verifyFalse(tc, check_settings(S).dynamic);
    S.dynamic = true;
    verifyError(tc, @() check_settings(S), 'check_settings:dynamic');
end
end

function test_errors_name_the_field(tc)
S = static_test_settings();
S.static_matrices = {'full', 'windows'};
verifyError(tc, @() check_settings(S), 'check_settings:static_matrices');
S.static_matrices = {'full', 'full'};
verifyError(tc, @() check_settings(S), 'check_settings:static_matrices');
S = static_test_settings();
S.static_measures = {'plv'};
verifyError(tc, @() check_settings(S), 'check_settings:static_measures');
S = static_test_settings();
S.static_families = {'global', 'modularity'};
verifyError(tc, @() check_settings(S), 'check_settings:static_families');
S = static_test_settings();
S.dynamic = 'no';
verifyError(tc, @() check_settings(S), 'check_settings:dynamic');
S = static_test_settings();
S.run_name = '../elsewhere';
verifyError(tc, @() check_settings(S), 'check_settings:run_name');
S = static_test_settings();
S.static_familes = {'global'};                                 % misspelled: stops, naming the field
verifyError(tc, @() check_settings(S), 'check_settings:unknown');
end

function test_settings_of_the_windows_are_reported(tc)
% static only: a setting of the windows or the states changed from its
% default is not used, and the warning says which; the template's values,
% and the window settings that the filter's 'window' rule reads, are not reported
verifyWarningFree(tc, @() check_settings(static_test_settings()));
S = pipeline_settings();
S.dynamic = false;
verifyWarningFree(tc, @() check_settings_quietly(S));        % (the window mean's windows: test_window_mean_of_a_static_run)
S = static_test_settings();
S.NCs = 4:8;
S.window_overlap = 0.75;
verifyWarning(tc, @() check_settings(S), 'check_settings:not_used');
txt = evalc('check_settings(S);');
verifySubstring(tc, txt, 'S.window_overlap, S.NCs');
S = static_test_settings();
S.filter_transition = 'window';
S.window_cycles = 5;
verifyWarningFree(tc, @() check_settings(S));
S.filter_transition = 'brainstorm';
verifyWarning(tc, @() check_settings(S), 'check_settings:not_used');
end

function test_window_mean_of_a_static_run(tc)
% static only, S.static_matrices with 'window_mean': accepted; the window
% settings are used then (not reported); a warning when the windows are
% longer than the shortest piece kept (S.min_segment)
S = static_test_settings();
S.static_matrices = {'full', 'window_mean'};
S.min_segment = 6;                                             % the AEC windows (6 s) fit
S.window_overlap = 0.75;  S.aec_step = 1;  S.window_cycles = 4;
S = verifyWarningFree(tc, @() check_settings(S));
verifyEqual(tc, S.static_matrices, {'full', 'window_mean'});
S = static_test_settings();
S.static_matrices = 'window_mean';                             % one name: a cell
verifyEqual(tc, check_settings_quietly(S).static_matrices, {'window_mean'});
S.min_segment = 4;                                             % the default: shorter than the AEC windows
verifyWarning(tc, @() check_settings(S), 'check_settings:window_vs_segment');
txt = evalc('check_settings(S);');
verifySubstring(tc, txt, 'The windows of aec_orth (up to 6 s)');
S.conn = {'wpli'};  S.state_conn = 'wpli';                     % phase only, 6 cycles of 10 Hz: 0.6 s
verifyWarningFree(tc, @() check_settings(S));
S.window_length = 5;                                           % fixed windows of 5 s
verifyWarning(tc, @() check_settings(S), 'check_settings:window_vs_segment');
txt = evalc('check_settings(S);');
verifySubstring(tc, txt, 'The windows of wpli (up to 5 s)');
S.static_matrices = {'full'};                                  % no window mean: the window settings unused
verifyWarning(tc, @() check_settings(S), 'check_settings:not_used');
end

function test_windows_of_one_measure(tc)
% S.measure_windows (6 Nov 2026): measures of S.conn only, the settings of
% their kind, values in range; the window warning reads each measure's
% windows; not used without windows; another stage 1 key only when windows
% are computed, whatever the order of the fields
S = static_test_settings();
S.static_matrices = {'window_mean'};
S.min_segment = 10;
S.measure_windows = struct('wpli', struct('window_length', 2), 'aec_orth', struct('aec_window', 8, 'aec_step', 1));
C = verifyWarningFree(tc, @() check_settings(S));
verifyEqual(tc, C.measure_windows, S.measure_windows);
bad = {struct('plv', struct('window_length', 2)), ...          % not in S.conn
    struct('wpli', struct('aec_window', 2)), ...               % an AEC setting for a phase measure
    struct('aec_orth', struct('window_length', 2)), ...        % and the other way round
    struct('wpli', struct('window_overlap', 1)), ...
    struct('wpli', struct('window_cycles', 0)), ...
    struct('aec_orth', struct('aec_step', -1)), ...
    struct('wpli', struct()), struct('wpli', 2), 'wpli'};
for i = 1:numel(bad)
    S2 = S;  S2.measure_windows = bad{i};
    verifyError(tc, @() check_settings(S2), 'check_settings:measure_windows', sprintf('case %d', i));
end
S2 = S;  S2.measure_windows = [];                               % [] = none
verifyEqual(tc, check_settings(S2).measure_windows, struct());
S2 = S;  S2.min_segment = 4;  S2.measure_windows = struct('aec_orth', struct('aec_window', 3));
verifyWarningFree(tc, @() check_settings(S2));                  % AEC windows of 3 s, wPLI of 0.6 s
S2.measure_windows.wpli = struct('window_length', 5);
verifyWarning(tc, @() check_settings(S2), 'check_settings:window_vs_segment');
txt = evalc('check_settings(S2);');
verifySubstring(tc, txt, 'The windows of wpli (up to 5 s)');
S2 = S;  S2.static_matrices = {'full'};                         % no windows: not used
verifyWarning(tc, @() check_settings(S2), 'check_settings:not_used');
txt = evalc('check_settings(S2);');
verifySubstring(tc, txt, 'S.measure_windows');
key = @(X) settings_key(check_settings_quietly(X), '', 'eloreta', 'alpha');
k0 = key(setfield(S, 'measure_windows', struct()));
km = key(S);
verifyNotEqual(tc, km, k0);
S3 = S;  S3.measure_windows = struct('aec_orth', struct('aec_step', 1, 'aec_window', 8), 'wpli', struct('window_length', 2));
verifyEqual(tc, key(S3), km);                                   % the same settings in another order
w = warning('off', 'check_settings:not_used');
c = onCleanup(@() warning(w));
verifyEqual(tc, key(S2), key(setfield(S2, 'measure_windows', struct())));   % no windows: no effect on the key
if isfolder(fullfile(tc.TestData.root, '2- Dynamic connectivity'))   % with the states: the state measure's windows
    Sd = S;  Sd.dynamic = true;  Sd.save_scouts = true;  Sd.min_segment = 4;
    Sd.measure_windows = struct('wpli', struct('window_length', 5));
    verifyWarning(tc, @() check_settings(Sd), 'check_settings:window_vs_segment');
    verifyNotEqual(tc, key(Sd), key(setfield(Sd, 'measure_windows', struct())));
end
end

function S = check_settings_quietly(S)
w = warning('off', 'check_settings:window_vs_segment');
c = onCleanup(@() warning(w));
S = check_settings(S);
end

function test_simulation_needs_the_whole_recording(tc)
% run_static_simulation (and run_static_demo) report on the whole-recording
% results: with the window mean only (the default since 6 November 2026) it
% stops at once, before simulating, instead of after stage 1 and the
% statistics on a missing file (the demo did until 6 November, late)
addpath(fullfile(tc.TestData.root, 'simulation'));
S = static_test_settings();
for m = {{'window_mean'}, {}}
    S.static_matrices = m{1};
    verifyError(tc, @() run_static_simulation(tempname, struct(), S), 'run_static_simulation:full');
end
verifyNotEmpty(tc, regexp(fileread(fullfile(tc.TestData.root, 'run_static_demo.m')), ...
    'S\.static_matrices = \{''full''\};', 'once'));          % the demo asks for the whole recording
end

function test_outputs_folder(tc)
S = check_settings(static_test_settings());
verifyEqual(tc, analysis_dir(S), fullfile(S.results, 'Statistics', 'test'));
S.run_name = '';
verifyEqual(tc, analysis_dir(S), fullfile(S.results, 'Statistics', 'main'));
end

function test_keys_follow_the_settings(tc)
% stage 1 outputs depend on S.dynamic and S.reliability_block (what is
% computed and saved), not on the statistics
S = check_settings(static_test_settings());
k = settings_key(S, '', 'eloreta', 'alpha');
S2 = S;  S2.reliability_block = 5;
verifyNotEqual(tc, settings_key(S2, '', 'eloreta', 'alpha'), k);
S2 = S;  S2.n_perm = 10;  S2.static_families = {'global'};  S2.nbs.perms = 10;
verifyEqual(tc, settings_key(S2, '', 'eloreta', 'alpha'), k);
if isfolder(fullfile(tc.TestData.root, '2- Dynamic connectivity'))
    S2 = S;  S2.dynamic = true;
    verifyNotEqual(tc, settings_key(S2, '', 'eloreta', 'alpha'), k);
end
% the window mean of a static run: another key; its outputs hold the
% whole-recording matrices too, so a run without it accepts them (not the
% other way round)
Sw = S;  Sw.static_matrices = {'full', 'window_mean'};
[kw, ~, acc_w] = settings_key(Sw, '', 'eloreta', 'alpha');
[~, ~, acc] = settings_key(S, '', 'eloreta', 'alpha');
verifyNotEqual(tc, kw, k);
verifyEqual(tc, acc{1}, k);
verifyTrue(tc, ismember(kw, acc));
verifyFalse(tc, ismember(k, acc_w));
Sw.static_matrices = {'window_mean'};                          % the same outputs: the same key
verifyEqual(tc, settings_key(Sw, '', 'eloreta', 'alpha'), kw);
end

function test_scout_function(tc)
% S.scout_function: 'mean' by default, Brainstorm's five functions only;
% the free orientation needs 'pca'; power and RMS with phase measures give
% a warning; the mean keeps the key outputs had before the setting existed,
% any other function changes it
S = check_settings(static_test_settings());
verifyEqual(tc, S.scout_function, 'mean');
k = settings_key(S, '', 'eloreta', 'alpha');
verifyEqual(tc, settings_key(rmfield(S, 'scout_function'), '', 'eloreta', 'alpha'), k);
keys = strings(1, 4);
fun = {'pca', 'max', 'power', 'rms'};
for i = 1:4
    S2 = static_test_settings();  S2.scout_function = upper(fun{i});   % any case
    S2.conn = {'aec', 'aec_orth'};                                       % amplitude only: no warning
    S2 = verifyWarningFree(tc, @() check_settings(S2));
    verifyEqual(tc, S2.scout_function, fun{i});
    S2.conn = S.conn;
    keys(i) = string(settings_key(S2, '', 'eloreta', 'alpha'));
end
verifyEqual(tc, numel(unique([keys, string(k)])), 5);
S2 = static_test_settings();  S2.scout_function = 'median';
verifyError(tc, @() check_settings(S2), 'check_settings:scout_function');
S2 = static_test_settings();  S2.orientation = 'free';
verifyError(tc, @() check_settings(S2), 'check_settings:scout_function');
S2.scout_function = 'pca';
verifyWarningFree(tc, @() check_settings(S2));
for f = {'power', 'rms'}
    S2 = static_test_settings();  S2.scout_function = f{1};             % with wPLI
    verifyWarning(tc, @() check_settings(S2), 'check_settings:scout_function');
end
S2 = static_test_settings();  S2.scout_function = 'max';
verifyWarningFree(tc, @() check_settings(S2));
end

