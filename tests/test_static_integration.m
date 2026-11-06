function tests = test_static_integration
%TEST_STATIC_INTEGRATION The static-only pipeline (S.dynamic = false) on a BIDS dataset: every output, and its checks.
%   From the pipeline folder:  results = runtests('tests/test_static_integration.m')
%   The simulated recordings of static_test_recordings (BIDS layout, three
%   file formats, participants.tsv), stage 1 in the alpha band with wPLI
%   and orthogonalised AEC, the statistics with few permutations (both by
%   run_full_pipeline, which also checks its options and errors), the PSD
%   step (stage 1b) in two bands; the stage 1 outputs are kept in a temporary folder between runs (they are
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
tc.TestData.log = evalc('run_full_pipeline(S);');           % stage 1, no stage 2, the statistics
tc.TestData.S = S;
tc.TestData.out = fullfile(S.results, 'Statistics', S.run_name);
end

function teardownOnce(tc)
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
for f = {'summary.csv', 'families.csv', 'analysis_log.txt', 'provenance_stage3.mat'}
    verifyTrue(tc, isfile(fullfile(sa, f{1})), f{1});
end
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
Su = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(Su), 2 * 3 + 1 + 1);                   % 2 measures x 3 families + power (once) + spectral profile
verifyTrue(tc, isfolder(fullfile(sa, 'spectral_profile')));
Fm = readtable(fullfile(sa, 'families.csv'), 'TextType', 'string');
verifyEqual(tc, height(Fm), (2 * 3 + 1) * 2 + 2);             % x (comparison, correlations) + spectral
verifyEqual(tc, nnz(contains(Fm.family, "static_power")), 2);  % power does not depend on the measure
verifyFalse(tc, isfile(fullfile(sa, 'eloreta_alpha_aec_orth', 'static_power_group_comparison.csv')));
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
S.results = fullfile(tempdir, 'static_test_results_window_mean');   % stage 1 kept between runs
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

function test_full_pipeline(tc)
% the stages in order, stage 2 skipped without states; 'from' and 'to';
% a stage that stops names itself and how to resume
log = tc.TestData.log;
i1 = strfind(log, 'Stage 1 of 3: connectivity of every recording (run_all)');
ib = strfind(log, 'Stage 1b: PSD analysis (run_psd_analysis)');
i2 = strfind(log, 'Stage 2 of 3 (brain network states): skipped, S.dynamic = false');
i3 = strfind(log, 'Stage 3 of 3: statistics (static connectivity only) (run_statistics)');
verifyTrue(tc, isscalar(i1) && isscalar(ib) && isscalar(i2) && isscalar(i3) && i1 < ib && ib < i2 && i2 < i3, log);
verifySubstring(tc, log, 'stage1_summary.csv: 8 recordings, all made or reused');
verifySubstring(tc, log, fullfile(tc.TestData.out, 'Statistical analysis'));
S = tc.TestData.S;
log1 = evalc('run_full_pipeline(S, ''to'', 1);');                 % stage 1 only: everything reused
verifySubstring(tc, log1, 'Stage 1 of 3');
verifyFalse(tc, contains(log1, 'Stage 1b') || contains(log1, 'Stage 2 of 3') || contains(log1, 'Stage 3 of 3'), log1);
S4 = S;  S4.psd.run = false;                                        % the PSD step switched off
log4 = evalc('run_full_pipeline(S4, ''to'', ''1b'');');
verifySubstring(tc, log4, 'Stage 1b (PSD analysis): skipped, S.psd.run = false');
verifyFalse(tc, contains(log4, 'Stage 3 of 3'), log4);
verifyError(tc, @() run_full_pipeline(S, 'from', 3, 'to', 1), 'run_full_pipeline:stages');
verifyError(tc, @() run_full_pipeline(S, 'from', 4), 'MATLAB:InputParser:ArgumentFailedValidation');
S2 = S;  S2.results = fullfile(tempdir, 'static_test_no_stage1');           % no stage 1 outputs there
c = onCleanup(@() rmdir_if(S2.results));
verifyError(tc, @() run_quietly(S2, 'from', 3), 'run_full_pipeline:stage1');
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
% stage 1b on the BIDS recordings (three formats): every recording at both
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

function test_stage2_refuses_a_static_run(tc)
assumeTrue(tc, isfolder(fullfile(tc.TestData.S.root, '2- Dynamic connectivity')), 'no stage 2 in this copy (v_static)');
verifyError(tc, @() compute_modkmeans_states_meta_crit(tc.TestData.S), 'compute_modkmeans_states_meta_crit:static');
end
