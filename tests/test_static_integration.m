function tests = test_static_integration
%TEST_STATIC_INTEGRATION The static-only pipeline (S.dynamic = false) on a BIDS dataset: every output, and its checks.
%   From the pipeline folder:  results = runtests('tests/test_static_integration.m')
%   The simulated recordings of static_test_recordings (BIDS layout, three
%   file formats, participants.tsv), stage 1 in the alpha band with wPLI
%   and orthogonalised AEC, the statistics with few permutations (both by
%   run_full_pipeline, which also checks its options and errors), the PSD
%   step (stage 2) in two bands; the stage 1 outputs are kept in a temporary folder between runs (they are
%   recomputed when the code or the settings change).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'));
S = static_test_settings();
S.parallel = false;
S.spectral_profile = true;
S.figures = true;
S = setup_toolboxes(S);
S.data = static_test_recordings(S);
S.results = fullfile(tempdir, 'static_test_results');
S.file_pattern = {'sub-*/eeg/*_eeg.set', 'sub-*/eeg/*_eeg.edf', 'sub-*/eeg/*_eeg.vhdr'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.subs_info = fullfile(S.data, 'participants.tsv');
S.id_column = 'participant_id';
S.covariates = {'age'};
S.correlations = struct('name', 'Score', 'columns', '^score$', 'subjects', {{'group', '==', 'patient'}}, ...
    'covariates', {{}}, 'zero_as_missing', false);
S.n_perm = 200;
S.nbs = struct('thresh', 3.1, 'perms', 200, 'size', 'extent');
S.enrichment = struct('n_perm', 200);
S.psd.bands = struct('name', {'alpha', 'beta'}, 'fmin', {7, 13}, 'fmax', {13, 30});
tc.TestData.S = S;
tc.TestData.out = fullfile(S.results, 'Statistics', S.run_name);
tc.TestData.log = evalc('run_full_pipeline(S);');           % stage 1, the PSD, no states, the statistics
end

function teardownOnce(tc)
if ~isfield(tc.TestData, 'out'), return; end                % setupOnce stopped before the run
if isfolder(tc.TestData.out), rmdir(tc.TestData.out, 's'); end
d = fullfile(tc.TestData.S.results, 'PSD', tc.TestData.S.run_name);
if isfolder(d), rmdir(d, 's'); end
end

function test_stage1_outputs_have_no_windows(tc)
S = tc.TestData.S;
d = fullfile(S.results, 'Connectivity', 'cmats');
verifyEmpty(tc, dir(fullfile(d, '*', 'cmats_*.mat')));
f = dir(fullfile(d, '*', 'static_cmats_eloreta_ASA_10-05_343_Desikan_Killiany_68_alpha.mat'));
verifyEqual(tc, numel(f), 8);
w = load(fullfile(f(1).folder, f(1).name));
verifyEmpty(tc, w.static_cmats);
verifyEqual(tc, size(w.static_full_cmats), [2, 68, 68]);
verifyEqual(tc, size(w.static_halves), [2, 2, 68, 68]);
verifyEqual(tc, size(w.static_power), [68, 1]);
verifyFalse(tc, isfield(w.info, 'window_s'));
verifyEmpty(tc, w.info.windows);
verifyFalse(tc, isfolder(fullfile(S.results, 'Connectivity', 'BNS analysis')));
end

function test_statistics_outputs(tc)
out = tc.TestData.out;
verifyFalse(tc, contains(tc.TestData.log, 'failed'), tc.TestData.log);
sa = fullfile(out, 'Statistical analysis');
prov = sprintf('provenance_stage%d.mat', last_stage(tc.TestData.S));     % 4, or 3 without the states
for f = {'summary.csv', 'families.csv', 'analysis_log.txt', prov}
    verifyTrue(tc, isfile(fullfile(sa, f{1})), f{1});
end
% provenance: the code of the statistics and the subject table
S = tc.TestData.S;
[~, code] = code_fingerprint(check_settings(S), 'statistics');
P = load(fullfile(sa, prov)).provenance;
verifyEqual(tc, P.stage, last_stage(S));
verifyEqual(tc, P.code, code);
verifyEqual(tc, P.subject_table, file_fingerprint({S.subs_info}));
Q = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, unique(Q.code_fingerprint), string(code));
verifyEqual(tc, unique(Q.subject_table_fingerprint), string(P.subject_table));
verifySubstring(tc, tc.TestData.log, ['Code of the statistics: ' code]);
for tag = {'eloreta_alpha_wpli', 'eloreta_alpha_aec_orth'}
    d = fullfile(sa, tag{1});
    files = {'static_GT_full_group_comparison.csv', 'static_networks_full_group_comparison.csv', ...
        'static_regions_full_group_comparison.csv', 'static_networks_full_correlations.csv', ...
        'features_full.csv', 'static_reliability.csv', 'group_matrices_full.png'};
    if strcmp(tag{1}, 'eloreta_alpha_wpli'), files{end + 1} = 'static_power_group_comparison.csv'; end   % first measure only
    for f = files
        verifyTrue(tc, isfile(fullfile(d, f{1})), fullfile(tag{1}, f{1}));
    end
    d = fullfile(out, 'NBS', tag{1});
    for f = {'NBS_full.csv', 'NBS_full.mat', 'enrichment_full.csv', 'edges_patient_control_full.csv', ...
            'patient_control_full_edges.png'}
        verifyTrue(tc, isfile(fullfile(d, f{1})), fullfile('NBS', tag{1}, f{1}));
    end
end
F = readtable(fullfile(sa, 'eloreta_alpha_wpli', 'features_full.csv'), 'TextType', 'string');
verifyEqual(tc, height(F), 8 * (3 + 28 + 68 + 8));
verifyEqual(tc, numel(unique(F.id)), 8);
verifyEqual(tc, sort(unique(F.feature(F.family == "global"))), ["CPL"; "ClustCoeff"; "Strength"]);
G = readtable(fullfile(sa, 'eloreta_alpha_wpli', 'static_networks_full_group_comparison.csv'), 'TextType', 'string');
verifyTrue(tc, all(ismember({'q_fdr', 'p_maxT_family', 'p_bb', 'significant', 'g'}, G.Properties.VariableNames)));
verifyEqual(tc, height(G), 28);
verifyTrue(tc, all(G.q_fdr >= G.p_unc - 1e-12));
% each group's mean and SD of the raw values, on the subjects of the test
% (diff: adjusted for age here), as the README describes
P = readtable(fullfile(tc.TestData.S.data, 'participants.tsv'), 'FileType', 'text', 'TextType', 'string');
x = F(F.feature == "within_DMN", :);
[~, r] = ismember(x.id, P.participant_id);
a = x.value(P.group(r) == "patient");  b = x.value(P.group(r) == "control");
g = G(G.Parameter == "within_DMN", :);
verifyEqual(tc, [g.mean_A, g.sd_A, g.mean_B, g.sd_B, g.N_A, g.N_B], [mean(a), std(a), mean(b), std(b), numel(a), numel(b)], ...
    'AbsTol', 1e-10);
Su = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(Su), 2 * 3 + 1 + 1);                   % 2 measures x 3 families + power (once) + spectral profile
verifyTrue(tc, isfolder(fullfile(sa, 'spectral_profile')));
Fm = readtable(fullfile(sa, 'families.csv'), 'TextType', 'string');
verifyEqual(tc, height(Fm), (2 * 3 + 1) * 2 + 2);             % x (comparison, correlations) + spectral
verifyEqual(tc, nnz(contains(Fm.family, "static_power")), 2);  % power does not depend on the measure
verifyFalse(tc, isfile(fullfile(sa, 'eloreta_alpha_aec_orth', 'static_power_group_comparison.csv')));
end

function test_edges_only_and_failures(tc)
% S.static_families empty: the edges only (NBS, enrichment, edge tables)
% and their reliability, with a warning; the later steps of the band
% still run. Then a step that fails (a network table without its
% columns): the other steps run, the run ends "Done with errors" and
% failures.csv names it; with S.stop_on_error the run stops
S = tc.TestData.S;
S.static_families = {};
S.psd.run = false;  S.figures = false;  S.spectral_profile = false;  S.correlations = S.correlations([]);
S.run_name = 'edges_only';
out = fullfile(S.results, 'Statistics', 'edges_only');
c = onCleanup(@() rmdir_if(out));
verifyWarning(tc, @() check_settings(S), 'check_settings:static_families');
log = evalc(sprintf('run_full_pipeline(S, ''from'', %d);', last_stage(S)));
verifyFalse(tc, contains(log, 'failed') || contains(log, 'Done with errors'), log);
verifySubstring(tc, log, 'Done in');
sa = fullfile(out, 'Statistical analysis');
for tag = {'eloreta_alpha_wpli', 'eloreta_alpha_aec_orth'}
    for f = {'NBS_full.csv', 'enrichment_full.csv', 'edges_patient_control_full.csv'}
        verifyTrue(tc, isfile(fullfile(out, 'NBS', tag{1}, f{1})), fullfile(tag{1}, f{1}));
    end
    R = readtable(fullfile(sa, tag{1}, 'static_reliability.csv'), 'TextType', 'string');
    verifyEqual(tc, unique(R.family), ["edges"; "pattern"]);
end
verifyFalse(tc, isfile(fullfile(sa, 'failures.csv')));
Su = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(Su), 0);
% a step that fails
bad = fullfile(tempdir, 'network_table_without_labels.csv');
writetable(table(["a"; "b"], [1; 2], 'VariableNames', {'name', 'x'}), bad);
c2 = onCleanup(@() delete(bad));
S.network_table = bad;  S.static_families = {'global'};
log = evalc(sprintf('run_full_pipeline(S, ''from'', %d);', last_stage(S)));
verifySubstring(tc, log, 'Done with errors');
verifySubstring(tc, log, 'WITH 1 FAILURE(S)');
Fl = readtable(fullfile(sa, 'failures.csv'), 'TextType', 'string', 'Delimiter', ',');
verifyEqual(tc, height(Fl), 1);
verifySubstring(tc, Fl.what(1), 'eloreta alpha, full matrices');
verifyFalse(tc, isfile(fullfile(sa, 'eloreta_alpha_wpli', 'static_GT_full_group_comparison.csv')));
S.stop_on_error = true;
verifyError(tc, @() run_quietly(S, 'from', last_stage(S)), 'run_full_pipeline:stage');
end

function test_window_mean_statistics(tc)
% static only with both matrices (S.static_matrices): stage 1 again, with
% the windows (its own key), then each matrix analysed on its own:
% features, split-half reliability, NBS; the window mean of wPLI above its
% whole-recording value (the bias of short windows)
S = tc.TestData.S;
S.static_matrices = {'full', 'window_mean'};
S.static_measures = {};
S.static_families = {'global'};
S.psd.run = false;  S.figures = false;  S.spectral_profile = false;
S.run_name = 'window_mean';
S.results = fullfile(tempdir, 'static_test_results_window_mean_halves');   % stage 1 kept between runs
c = onCleanup(@() rmdir_if(fullfile(S.results, 'Statistics')));
w = warning('off', 'check_settings:window_vs_segment');                 % 6 s AEC windows, S.min_segment 4 s
c2 = onCleanup(@() warning(w));
log = evalc('run_full_pipeline(S);');
verifyFalse(tc, contains(log, 'failed'), log);
out = fullfile(S.results, 'Statistics', 'window_mean');
sa = fullfile(out, 'Statistical analysis');
for tag = {'eloreta_alpha_wpli', 'eloreta_alpha_aec_orth'}
    for f = {'features_full.csv', 'features_window_mean.csv', 'static_reliability.csv', ...
            'static_reliability_window_mean.csv', 'static_GT_full_group_comparison.csv', ...
            'static_GT_window_mean_group_comparison.csv'}
        verifyTrue(tc, isfile(fullfile(sa, tag{1}, f{1})), fullfile(tag{1}, f{1}));
    end
    for f = {'NBS_full.csv', 'NBS_window_mean.csv'}
        verifyTrue(tc, isfile(fullfile(out, 'NBS', tag{1}, f{1})), fullfile('NBS', tag{1}, f{1}));
    end
end
d = fullfile(sa, 'eloreta_alpha_wpli');
Ff = readtable(fullfile(d, 'features_full.csv'), 'TextType', 'string');
Fw = readtable(fullfile(d, 'features_window_mean.csv'), 'TextType', 'string');
verifyEqual(tc, height(Fw), height(Ff));
verifyGreaterThan(tc, mean(Fw.value(Fw.feature == "Strength")), mean(Ff.value(Ff.feature == "Strength")));
Rf = readtable(fullfile(d, 'static_reliability.csv'), 'TextType', 'string');
Rw = readtable(fullfile(d, 'static_reliability_window_mean.csv'), 'TextType', 'string');
verifyEqual(tc, Rw.feature, Rf.feature);
verifyTrue(tc, all(isfinite(Rw.r_sb(Rw.family == "edges"))));
verifyNotEqual(tc, Rw.r_sb, Rf.r_sb);
Fm = readtable(fullfile(sa, 'families.csv'), 'TextType', 'string');
verifyEqual(tc, nnz(contains(Fm.family, "window_mean")), nnz(contains(Fm.family, "_full")));
% the halves of each measure in its own blocks: wPLI's 0.6 s windows in
% 10 s blocks, AEC's 6 s windows in 18 s blocks; the log gives them and
% how many recordings have both halves
verifyTrue(tc, ~isempty(regexp(log, ['Split-half reliability, window_mean \(blocks of 10 s\), [1-8] of 8 ' ...
    'recordings with both halves'], 'once')), log);
verifyTrue(tc, ~isempty(regexp(log, ['Split-half reliability, window_mean \(blocks of 18 s\), [1-8] of 8 ' ...
    'recordings with both halves'], 'once')), log);
verifyTrue(tc, ~isempty(regexp(log, 'Split-half reliability, full \(blocks of 10 s\), 8 of 8', 'once')), log);
% stage 1 outputs of earlier code (every measure's windows in 10 s blocks,
% no static_halves_block): the AEC halves of the window mean are left out,
% with the reason, wPLI's kept
old = fullfile(tempname, 'results');
copyfile(fullfile(S.results, 'Connectivity', 'cmats'), fullfile(old, 'Connectivity', 'cmats'));
c3 = onCleanup(@() rmdir(fileparts(old), 's'));
f = dir(fullfile(old, 'Connectivity', 'cmats', '*', 'static_cmats_*.mat'));
for i = 1:numel(f)
    w = rmfield(load(fullfile(f(i).folder, f(i).name)), 'static_halves_block');
    save(fullfile(f(i).folder, f(i).name), '-struct', 'w');
end
So = S;  So.results = old;
D = load_static(So, 'eloreta', 'alpha', 'window_mean');
k = find(strcmp(D.measures, 'aec_orth'));  kw = find(strcmp(D.measures, 'wpli'));
verifyTrue(tc, all(isnan(D.halves(:, k, :, :, :)), 'all'));
verifyTrue(tc, all(contains(D.halves_why(k, :), "stage 1 outputs of earlier code")));
verifyTrue(tc, all(D.halves_why(kw, :) == ""));
verifyTrue(tc, all(isfinite(D.halves(:, kw, 1, 2, :)), 'all'));
verifyEqual(tc, D.halves_block, [10; 10]);
D = load_static(S, 'eloreta', 'alpha', 'window_mean');
verifyEqual(tc, D.halves_block([kw, k]), [10; 18]);
% the window mean of a phase measure: one warning in the run (the log
% file, a diary, stays empty under evalc)
verifyEqual(tc, numel(strfind(log, 'Static analyses of the window mean of wpli')), 1);
% a recording without the state measure (S.state_conn) is left out of
% the states, not of the static analyses: each measure keeps its own
% recordings
f = dir(fullfile(old, 'Connectivity', 'cmats', '*', 'static_cmats_*.mat'));
w = load(fullfile(f(1).folder, f(1).name));
w.info.not_computed.(So.conn{1}) = 'test';
save(fullfile(f(1).folder, f(1).name), '-struct', 'w');
[~, rec] = fileparts(f(1).folder);
subs = verifyWarning(tc, @() stage1_outputs(check_settings(So), 'eloreta', 'alpha'), 'stage1_outputs:not_computed');
verifyFalse(tc, ismember(rec, subs));
verifyWarningFree(tc, @() stage1_outputs(check_settings(So), 'eloreta', 'alpha', 'static'));
D = load_static(check_settings(So), 'eloreta', 'alpha', 'window_mean');
verifyEqual(tc, numel(D.ids), 8);
% the stage 1 outputs of the main run (without the window mean) have no
% window mean: refused
verifyError(tc, @() load_static(setfield(tc.TestData.S, 'static_matrices', {'window_mean'}), ...
    'eloreta', 'alpha', 'window_mean'), 'stage1_outputs:settings');
end

function test_statistics_refuse_other_settings(tc)
% outputs made with other stage 1 settings, or from a changed recording,
% are not combined; no window mean without windows
S = tc.TestData.S;
S2 = S;  S2.edge_trim = 2;
verifyError(tc, @() load_static(S2, 'eloreta', 'alpha'), 'stage1_outputs:settings');
verifyError(tc, @() load_static(S, 'eloreta', 'alpha', 'window_mean'), 'load_static:which');
D = load_static(S, 'eloreta', 'alpha');
verifyEqual(tc, numel(D.ids), 8);
verifyEqual(tc, D.measures, {'wpli', 'aec_orth'});
verifyEqual(tc, size(D.conn), [2, 68, 68, 8]);
S3 = S;  S3.groups = {'control'};                          % recordings of another group left out
D3 = load_static(S3, 'eloreta', 'alpha');
verifyEqual(tc, numel(D3.ids), 4);
end

function test_stale_outputs_and_fresh_folders(tc)
% stage 1 outputs made with other settings stop the statistics, also
% without S.stop_on_error (not one failed family per matrix, which would
% leave fewer families to adjust over); a recording whose
% stage 1 failed after some bands were saved (no info.qc) is left out; each
% run of the statistics starts from empty output folders (the cached state
% contrasts kept), with a log of that run only
S = tc.TestData.S;
S.psd.run = false;  S.figures = false;  S.spectral_profile = false;  S.correlations = S.correlations([]);
S.static_families = {'global'};  S.reliability = false;
S.run_name = 'fresh_folders';
out = fullfile(S.results, 'Statistics', S.run_name);
c = onCleanup(@() rmdir_if(out));
S2 = S;  S2.edge_trim = 2;                                 % other stage 1 settings than the outputs'
try
    evalc(sprintf('run_full_pipeline(S2, ''from'', %d);', last_stage(S2)));
    verifyFail(tc, 'the statistics ran on stage 1 outputs of other settings');
catch err
    verifyEqual(tc, err.identifier, 'run_full_pipeline:stage');
    verifyEqual(tc, err.cause{1}.identifier, 'stage1_outputs:settings');
end
% stale outputs of an earlier run under the same name
stale = {fullfile(out, 'Statistical analysis', 'old_family', 'old_group_comparison.csv'), ...
    fullfile(out, 'NBS', 'old_tag', 'component_1.png'), ...
    fullfile(out, 'Network visualization', 'old_tag', 'old_figure.png')};
cache = fullfile(out, 'Network visualization', 'old_tag', 'state_contrasts_old_tag.mat');
for f = [stale, {cache}]
    if ~isfolder(fileparts(f{1})), mkdir(fileparts(f{1})); end
    fclose(fopen(f{1}, 'w'));
end
evalc(sprintf('run_full_pipeline(S, ''from'', %d);', last_stage(S)));
for f = stale, verifyFalse(tc, isfile(f{1}), f{1}); end
verifyTrue(tc, isfile(cache));
logf = fullfile(out, 'Statistical analysis', 'analysis_log.txt');
fid = fopen(logf, 'a');  fprintf(fid, 'LINE OF AN EARLIER RUN\n');  fclose(fid);
evalc(sprintf('run_full_pipeline(S, ''from'', %d);', last_stage(S)));
verifyTrue(tc, isfile(logf));
verifyFalse(tc, contains(fileread(logf), 'LINE OF AN EARLIER RUN'));   % a new log (not appended to)
% a recording without info.qc (stage 1 failed after this band was saved)
S3 = S;
S3.results = fullfile(tempdir, 'static_test_partial');
c3 = onCleanup(@() rmdir_if(S3.results));
copyfile(fullfile(S.results, 'Connectivity', 'cmats'), fullfile(S3.results, 'Connectivity', 'cmats'));
subs = stage1_outputs(S3, 'eloreta', 'alpha', 'static');
d = dir(fullfile(S3.results, 'Connectivity', 'cmats', subs{1}, 'static_cmats_*_alpha.mat'));
w = load(fullfile(d(1).folder, d(1).name));
w.info = rmfield(w.info, 'qc');
save(fullfile(d(1).folder, d(1).name), '-struct', 'w');
subs3 = verifyWarning(tc, @() stage1_outputs(S3, 'eloreta', 'alpha', 'static'), 'stage1_outputs:partial');
verifyEqual(tc, subs3, subs(2:end));
end

function test_full_pipeline(tc)
% the stages in order, numbered as the folders (1-4, or 1-3 without the
% states), the states skipped without windows; 'from' and
% 'to', '1b' (the PSD's number until 1.2.0) taken as 2; a stage that stops
% names itself and how to resume
log = tc.TestData.log;
S = tc.TestData.S;
n = last_stage(S);
[~, nd] = stage_folder(S.root, 'Dynamic connectivity');
i1 = strfind(log, sprintf('Stage 1 of %d: connectivity of every recording (run_all)', n));
ib = strfind(log, sprintf('Stage 2 of %d: PSD analysis (run_psd_analysis)', n));
i3 = strfind(log, sprintf('Stage %d of %d: statistics (static connectivity only) (run_statistics)', n, n));
if isempty(nd)                                                      % no states
    verifyEqual(tc, n, 3);
    i2 = ib + 1;
else
    verifyEqual(tc, [nd n], [3 4]);
    i2 = strfind(log, sprintf('Stage 3 of 4 (brain network states): skipped, S.dynamic = false'));
end
verifyTrue(tc, isscalar(i1) && isscalar(ib) && isscalar(i2) && isscalar(i3) && i1 < ib && ib < i2 && i2 < i3, log);
verifySubstring(tc, log, 'stage1_summary.csv: 8 recordings, all made or reused');
verifySubstring(tc, log, fullfile(tc.TestData.out, 'Statistical analysis'));
log1 = evalc('run_full_pipeline(S, ''to'', 1);');                 % stage 1 only: everything reused
verifySubstring(tc, log1, sprintf('Stage 1 of %d', n));
verifyFalse(tc, contains(log1, 'Stage 2 of') || contains(log1, 'Stage 3 of') || contains(log1, 'Stage 4 of'), log1);
S4 = S;  S4.psd.run = false;                                        % the PSD step switched off
log4 = evalc('run_full_pipeline(S4, ''to'', 2);');
verifySubstring(tc, log4, sprintf('Stage 2 of %d (PSD analysis): skipped, S.psd.run = false', n));
verifyFalse(tc, contains(log4, sprintf('Stage %d of', n)), log4);
log4 = evalc('run_full_pipeline(S4, ''to'', ''1b'');');             % '1b': 2, with a warning
verifySubstring(tc, log4, 'is stage 2 since 1.2.0');
verifySubstring(tc, log4, sprintf('Stage 2 of %d (PSD analysis): skipped', n));
verifyFalse(tc, contains(log4, sprintf('Stage %d of', n)), log4);
verifyError(tc, @() run_full_pipeline(S, 'from', 3, 'to', 1), 'run_full_pipeline:stages');
verifyError(tc, @() run_full_pipeline(S, 'from', n + 1), 'MATLAB:InputParser:ArgumentFailedValidation');
verifyError(tc, @() run_full_pipeline(S, 'from', '2b'), 'MATLAB:InputParser:ArgumentFailedValidation');
S2 = S;  S2.results = fullfile(tempdir, 'static_test_no_stage1');           % no stage 1 outputs there
c = onCleanup(@() rmdir_if(S2.results));
verifyError(tc, @() run_quietly(S2, 'from', n), 'run_full_pipeline:stage1');
S2.data = fullfile(tempdir, 'static_test_no_recordings');                     % stage 1 stops: no recording
mkdir(S2.data);
c2 = onCleanup(@() rmdir_if(S2.data));
try
    evalc('run_full_pipeline(S2);');
    verifyFail(tc, 'stage 1 without recordings did not stop');
catch err
    verifyEqual(tc, err.identifier, 'run_full_pipeline:stage');
    verifySubstring(tc, err.message, 'Stage 1 (run_all) stopped');
    verifySubstring(tc, err.message, 'run_full_pipeline(S, ''from'', 1)');
    verifyEqual(tc, err.cause{1}.identifier, 'get_subs:none');
end
end

function test_psd_step(tc)
% stage 2 (the PSD) on the BIDS recordings (three formats): every recording at both
% levels, the tests of both bands, the figures
S = tc.TestData.S;
out = fullfile(S.results, 'PSD', S.run_name);
verifySubstring(tc, tc.TestData.log, 'PSD analysis finished in');
Q = readtable(fullfile(out, 'recordings.csv'), 'TextType', 'string');
verifyEqual(tc, height(Q), 8);
verifyEqual(tc, Q.status, repmat("ok", 8, 1));
verifyGreaterThan(tc, min(Q.fit_r2_sensor), 0.8);
V = readtable(fullfile(out, 'values_source_eloreta_absolute.csv'), 'TextType', 'string', 'VariableNamingRule', 'preserve');
verifyEqual(tc, height(V), 8);
verifyEqual(tc, nnz(startsWith(V.Properties.VariableNames, 'alpha__')), 68);
sa = fullfile(out, 'Statistical analysis');
Su = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(Su), 2 * (2 * 3 + 3));                 % levels x (bands x 3 measures + exponent, offset, alpha frequency)
verifyEqual(tc, Su.group_tests(Su.family == "psd_sensor_absolute_alpha"), 64);
verifyEqual(tc, Su.correlations(Su.family == "psd_source_eloreta_relative_beta"), 68);
G = readtable(fullfile(sa, 'sensor', 'psd_sensor_periodic_alpha_group_comparison.csv'), 'TextType', 'string');
verifyTrue(tc, all(ismember({'Level', 'Measure', 'Band', 'Unit', 'g', 'q_fdr', 'p_maxT_family', 'p_bb', 'cluster', 'cluster_q'}, ...
    G.Properties.VariableNames)));
verifyTrue(tc, isfile(fullfile(sa, 'clusters.csv')));
for f = {'spectra_sensor_patient_control.png', 'map_source_eloreta_absolute_patient_control.png', ...
        'map_sensor_relative_patient_control.png', 'boxplots_sensor_absolute_patient_control.png'}
    verifyTrue(tc, isfile(fullfile(out, 'figures', f{1})), f{1});
end
end

function run_quietly(S, varargin) %#ok<INUSD> (used by evalc)
evalc('run_full_pipeline(S, varargin{:})');
end

function rmdir_if(d)
if isfolder(d), rmdir(d, 's'); end
end

function test_states_refuse_a_static_run(tc)
assumeFalse(tc, isempty(stage_folder(tc.TestData.S.root, 'Dynamic connectivity')), 'no states in this pipeline');
verifyError(tc, @() compute_modkmeans_states_meta_crit(tc.TestData.S), 'compute_modkmeans_states_meta_crit:static');
end

function n = last_stage(S)
% the number of the statistics, the last stage: 4, or 3 without the states
[~, n] = stage_folder(S.root, 'Analysis');
end
