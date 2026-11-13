function tests = test_static_stage1
%TEST_STATIC_STAGE1 Stage 1 without windows (S.dynamic = false) on simulated recordings: outputs, cache, formats, channels, cuts.
%   From the pipeline folder:  results = runtests('tests/test_static_stage1.m')
%   Eight simulated recordings in a BIDS layout (static_test_recordings:
%   EEGLAB, EDF and BrainVision files, 64 channels, 60 s), the alpha band,
%   three measures; outputs in a temporary folder.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'));
S = static_test_settings();
S.parallel = false;
S.conn = {'wpli', 'plv', 'aec'};
S.results = tempname;
S = setup_toolboxes(S);
[S.data, truth] = static_test_recordings(S);
S.file_pattern = {'sub-*/eeg/*_eeg.set', 'sub-*/eeg/*_eeg.edf', 'sub-*/eeg/*_eeg.vhdr'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
tc.TestData.S = S;
tc.TestData.truth = truth;
tc.TestData.tmp = tempname;
mkdir(tc.TestData.tmp);
evalc('run_all(S);');
end

function teardownOnce(tc)
for f = {tc.TestData.S.results, tc.TestData.tmp}
    if isfolder(f{1}), rmdir(f{1}, 's'); end
end
end

function f = output(results, rec)
f = fullfile(results, 'Connectivity', 'cmats', rec, 'static_cmats_eloreta_ASA_10-05_343_Desikan_Killiany_68_alpha.mat');
end

%% Outputs
function test_outputs(tc)
S = tc.TestData.S;
Q = readtable(fullfile(S.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(Q), 8);
verifyTrue(tc, all(Q.status == "ok"), strjoin(string(Q.message), ' | '));
verifyEqual(tc, Q.n_channels, 64 * ones(8, 1));
edf = endsWith(Q.file, ".edf");                                  % average-referenced simulations: rank 63,
verifyEqual(tc, nnz(edf), 3);                                     % also in EDF, whose 16-bit rounding MNE's
verifyEqual(tc, Q.data_rank(~edf), 63 * ones(5, 1));              % rule does not count (connectivity.signal_rank;
verifyEqual(tc, Q.data_rank(edf), 63 * ones(3, 1));               % 64 with a 1e-12 relative threshold)
verifyEqual(tc, Q.analysed_s, 58 * ones(8, 1), 'AbsTol', 0.01);   % 60 s less 1 s at each end
f = output(S.results, 'sub-control01_task-rest_eeg');
verifyTrue(tc, isfile(f), f);
verifyEmpty(tc, dir(fullfile(fileparts(f), 'cmats_*.mat')));      % no windows
w = load(f);
verifyEmpty(tc, w.static_cmats);
verifyEqual(tc, size(w.static_full_cmats), [3, 68, 68]);
verifyEqual(tc, size(w.static_halves), [2, 3, 68, 68]);
verifyEqual(tc, size(w.static_power), [68, 1]);
verifyEqual(tc, size(w.static_power_halves), [68, 2]);
verifyEqual(tc, w.info.conn, {'wpli', 'plv', 'aec'});
verifyEqual(tc, w.info.band_hz, [7, 13]);
verifyEqual(tc, sum(w.info.halves_s), w.info.analysed_s, 'AbsTol', 1e-9);
verifyEqual(tc, w.info.halves_s, [29, 29], 'AbsTol', 0.01);       % 1-10, 20-30, 40-50 s and 10-20, 30-40, 50-59 s
verifyEqual(tc, w.info.settings_key, settings_key(S, '', 'eloreta', 'alpha'));
verifyEqual(tc, w.info.qc.status, "ok");
for k = 1:3
    A = squeeze(w.static_full_cmats(k, :, :));
    verifyEqual(tc, A, A.', 'AbsTol', 1e-12);
    verifyTrue(tc, all(isfinite(A(:))));
end
% the halves agree with the whole recording more than with another recording
v = load(output(S.results, 'sub-patient01_task-rest_eeg'));
low = tril(true(68), -1);
a = squeeze(w.static_full_cmats(1, :, :));  h = squeeze(w.static_halves(1, 1, :, :));
o = squeeze(v.static_full_cmats(1, :, :));
verifyGreaterThan(tc, corr(a(low), h(low)), corr(o(low), h(low)));
end

function test_rerun_skips_and_redoes(tc)
% the same files and settings: skipped; another setting of stage 1: redone
S = tc.TestData.S;
evalc('run_all(S);');
Q = readtable(fullfile(S.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyTrue(tc, all(Q.status == "skipped (up to date)"));
S2 = S;
S2.reliability_block = 5;
S2.data = fullfile(tc.TestData.tmp, 'one');
copy_recording(S.data, 'sub-control01', S2.data);
S2.results = fullfile(tc.TestData.tmp, 'rerun');
evalc('run_all(S2);');
f = output(S2.results, 'sub-control01_task-rest_eeg');
w0 = load(f, 'info');
evalc('run_all(S2);');
Q = readtable(fullfile(S2.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyEqual(tc, Q.status, "skipped (up to date)");
S2.conn = {'wpli'};
S2.static_measures = {};                           % S2 was checked with three measures
evalc('run_all(S2);');
Q = readtable(fullfile(S2.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyEqual(tc, Q.status, "ok");
w1 = load(f);
verifyEqual(tc, size(w1.static_full_cmats), [1, 68, 68]);
verifyNotEqual(tc, w1.info.settings_key, w0.info.settings_key);
if isempty(stage_folder(S.root, 'Dynamic connectivity')), return, end   % no windows in this pipeline
S2.dynamic = true;                                 % windows wanted: redone, with cmats_*.mat
S2.null_test.n_surrogates = 0;                     % (no region time series kept)
S2 = check_settings(S2);
evalc('run_all(S2);');
verifyNotEmpty(tc, dir(fullfile(fileparts(f), 'cmats_*.mat')));
w2 = load(f);
verifyEqual(tc, w2.static_full_cmats, w1.static_full_cmats, 'AbsTol', 1e-12);   % whole recording: the same
verifyNotEmpty(tc, w2.static_cmats);
end

function test_individual_bands_without_alpha(tc)
% the default individual bands on a recording without an alpha rhythm (no
% IAF): the four anchored bands are left out, the fixed ones computed. Both
% edges are checked: delta_iaf is 1 Hz to IAF - 6, and a NaN upper edge
% made the whole recording fail (found by a usability test)
S = tc.TestData.S;
p = sim.static_params(S);
p.out_dir = fullfile(tc.TestData.tmp, 'no_alpha');
p.n_per_group = [1, 1];  p.duration_s = 40;  p.osc_rms = 0;  p.alpha_rms = 0;  p.seed = 8;
evalc('sim.simulate_static_dataset(S, p);');
S.individual_bands = true;
S.results = fullfile(tc.TestData.tmp, 'no_alpha_results');
S = setup_toolboxes(S);                            % restingIAF
rec = 'sub-control01_task-rest_eeg';
evalc('qc = run_pipeline([rec ''.set''], fullfile(p.out_dir, ''sub-control01'', ''eeg''), fullfile(S.results, ''Connectivity''), S);');
verifyEqual(tc, qc.status, "ok", qc.message);
verifyTrue(tc, isnan(qc.iaf));
verifyEqual(tc, qc.bands_skipped, "delta_iaf theta_iaf alpha_iaf beta1_iaf");
d = fullfile(S.results, 'Connectivity', 'cmats', rec);
for b = {'beta2', 'gamma'}
    verifyTrue(tc, isfile(fullfile(d, ['static_cmats_eloreta_ASA_10-05_343_Desikan_Killiany_68_' b{1} '.mat'])), b{1});
end
verifyEmpty(tc, dir(fullfile(d, 'static_cmats_*_iaf.mat')));
end

function test_window_mean_through_stage_1(tc)
% static only with the window mean (S.static_matrices): static_cmats is the
% mean of the windows (6 cycles, 50% overlap; AEC 6 s every 0.5 s) and
% static_halves_window its split halves; the windows are not saved (no
% cmats_*.mat); the whole recording is the same as without the window mean;
% these outputs serve a run without the window mean, not the other way
% round; fixed windows (S.window_length) for every band
S = tc.TestData.S;
rec = 'sub-control01_task-rest_eeg';
S2 = S;
S2.data = fullfile(tc.TestData.tmp, 'one_window_mean');
copy_recording(S.data, 'sub-control01', S2.data);
S2.results = fullfile(tc.TestData.tmp, 'window_mean');
S2.static_matrices = {'full', 'window_mean'};
S2.static_measures = {};
w0 = warning('off', 'check_settings:window_vs_segment');       % 6 s AEC windows, S.min_segment 4 s
c = onCleanup(@() warning(w0));
evalc('run_all(S2);');
f = output(S2.results, rec);
w = load(f);
v = load(output(S.results, rec));
verifyEmpty(tc, dir(fullfile(fileparts(f), 'cmats_*.mat')));
verifyEqual(tc, size(w.static_cmats), [3, 68, 68]);
verifyTrue(tc, all(isfinite(w.static_cmats), 'all'));
verifyEqual(tc, w.static_full_cmats, v.static_full_cmats, 'AbsTol', 1e-12);
verifyEqual(tc, w.static_halves, v.static_halves, 'AbsTol', 1e-12);
verifyEqual(tc, w.static_power, v.static_power, 'AbsTol', 1e-12);
verifyEqual(tc, size(w.static_halves_window), [2, 3, 68, 68]);
verifyTrue(tc, all(isfinite(w.static_halves_window), 'all'));
verifyEqual(tc, w.static_halves_block, [10, 10, 18], 'AbsTol', 1e-12);   % AEC: 3 x its 6 s windows
verifyEqual(tc, [w.info.window_s, w.info.step_s, w.info.aec_window_s, w.info.aec_step_s], [0.6, 0.3, 6, 0.5], 'AbsTol', 1e-12);
verifyEqual(tc, w.info.no_windows, struct());
verifyEqual(tc, w.info.qc.n_windows, numel(w.info.windows{1}.first_sample));
verifyFalse(tc, isfield(v, 'static_halves_window'));
low = tril(true(68), -1);
wm = squeeze(w.static_cmats(1, :, :));  wf = squeeze(w.static_full_cmats(1, :, :));
verifyGreaterThan(tc, mean(wm(low)), mean(wf(low)) + 0.2);       % wPLI of 0.6 s windows: the short-window bias
h = squeeze(mean(w.static_halves_window(:, 1, :, :), 1));
verifyGreaterThan(tc, corr(h(low), wm(low)), 0.9);              % the halves: the windows inside their blocks
S3 = S2;  S3.static_matrices = {'full'};                        % without the window mean: reused
evalc('run_all(S3);');
Q = readtable(fullfile(S3.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyEqual(tc, Q.status, "skipped (up to date)");
S4 = S3;  S4.results = fullfile(tc.TestData.tmp, 'full_then_window_mean');
evalc('run_all(S4);');                                          % without the window mean
S4.static_matrices = {'full', 'window_mean'};  S4.static_measures = {};
evalc('run_all(S4);');                                          % with it: not reused, redone
Q = readtable(fullfile(S4.results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
verifyEqual(tc, Q.status, "ok");
verifyEqual(tc, load(output(S4.results, rec), 'static_cmats').static_cmats, w.static_cmats);
S5 = S2;  S5.window_length = 2;                                 % fixed 2 s windows, 50% overlap
S5.results = fullfile(tc.TestData.tmp, 'window_mean_2s');
evalc('run_all(S5);');
x = load(output(S5.results, rec));
verifyEqual(tc, [x.info.window_s, x.info.step_s], [2, 1], 'AbsTol', 1e-12);
verifyEqual(tc, x.static_full_cmats, v.static_full_cmats, 'AbsTol', 1e-12);
xm = squeeze(x.static_cmats(1, :, :));
verifyLessThan(tc, mean(xm(low)), mean(wm(low)));               % longer windows: less bias
verifyGreaterThan(tc, mean(xm(low)), mean(wf(low)));
verifyEqual(tc, x.static_cmats(3, :, :), w.static_cmats(3, :, :), 'AbsTol', 1e-6);   % AEC: its own windows
S6 = S2;  S6.measure_windows = struct('wpli', struct('window_length', 2));   % 2 s windows for wPLI only
S6.results = fullfile(tc.TestData.tmp, 'window_mean_wpli_2s');
evalc('run_all(S6);');
y = load(output(S6.results, rec));
verifyEqual(tc, y.static_cmats(1, :, :), x.static_cmats(1, :, :), 'AbsTol', 1e-12);          % wPLI: as S5
verifyEqual(tc, y.static_halves_window(:, 1, :, :), x.static_halves_window(:, 1, :, :), 'AbsTol', 1e-12);
verifyEqual(tc, y.static_cmats(2:3, :, :), w.static_cmats(2:3, :, :), 'AbsTol', 1e-12);      % PLV, AEC: as S2
verifyEqual(tc, y.static_halves_window(:, 2:3, :, :), w.static_halves_window(:, 2:3, :, :), 'AbsTol', 1e-12);
verifyEqual(tc, y.info.windows{1}.samples, x.info.windows{1}.samples);
verifyEqual(tc, y.info.measure_windows, S6.measure_windows);
verifyEqual(tc, y.static_full_cmats, v.static_full_cmats, 'AbsTol', 1e-12);
end

function test_lcmv_kappa(tc)
% LCMV's covariance (S.lcmv_kappa): FieldTrip's filters equal the
% beamformer (L' iC L)^-1 L' iC of each source, iC the first kappa
% components of the covariance C, each 1 / (s + lambda), lambda 5% of the
% mean eigenvalue (ft_inv, 'lavrentiev'): kappa FieldTrip's own by default
% (components above 10 x channels x eps of the largest), the data rank with
% 'rank' (connectivity.signal_rank); data of 40 signals on 64 channels
% with a floor between the two limits (cleaned data saved with rounding)
S = tc.TestData.S;
[~, lf] = get_ftHeadmodel_and_ftLeadfield(S.montages{1}, S.atlases{1}, S);
el = get_ftElec(S.montages{1});
root = fileparts(fileparts(mfilename('fullpath')));
C = load(fullfile(root, '1- EEG source connectivity', 'inputs', 'tess_cortex_pial_low.mat'), 'VertNormals');
idx = round(linspace(1, numel(el.label), 64));
el.label = el.label(idx);  el.chanpos = el.chanpos(idx, :);  el.elecpos = el.elecpos(idx, :);  el.tra = eye(64);
lf.label = lf.label(idx);
n = numel(lf.leadfield);
G = zeros(64, n);
for i = 1:n
    g = lf.leadfield{i}(idx, :);
    lf.leadfield{i} = g - mean(g, 1);
    if size(g, 2) == 3, G(:, i) = lf.leadfield{i} * C.VertNormals(i, :)'; else, G(:, i) = lf.leadfield{i}; end
end
rs = RandStream('mrg32k3a', 'Seed', 7);
X = randn(rs, 64, 40) * randn(rs, 40, 20000);
X = X + 5e-6 * std(X(:)) * randn(rs, size(X));                  % eigenvalues ~5e-12 of the largest
X = X - mean(X, 1);  X = X - mean(X, 2);                          % average reference, centred
Cx = cov(X');
s = svd(Cx);
r = connectivity.signal_rank(X);
kappa_ft = nnz(s / s(1) > 10 * 64 * eps);
verifyEqual(tc, r, 40);
assumeGreaterThan(tc, kappa_ft, r, 'FieldTrip''s own limit would keep the same components');
[U, s] = svd(Cx);  s = diag(s);
lambda = 0.05 * trace(Cx) / 64;
v = randperm(rs, n, 50);
F = cell(1, 2);
kappas = {[], 'rank'};
tol = [1e-4, 1e-6];        % FieldTrip's own kappa keeps components of 5e-12 of the largest next to the
for j = 1:2                % average reference's zero (4e-16): their boundary is resolved to about 1e-5
    S.lcmv_kappa = kappas{j};
    F{j} = squeeze(get_inverse_solution(X, 512, eye(64), [], lf, el, C.VertNormals, {'lcmv'}, S));
    k = kappa_ft;  if j == 2, k = r; end
    iC = U(:, 1:k) * diag(1 ./ (s(1:k) + lambda)) * U(:, 1:k)';
    L = G(:, v);
    W = (L' * iC) ./ sum(L .* (iC * L), 1)';                       % (L' iC L)^-1 L' iC, source by source
    verifyLessThan(tc, max(abs(F{j}(v, :) - W), [], 'all') / max(abs(W(:))), tol(j), sprintf('kappa %d', k));
end
verifyGreaterThan(tc, max(abs(F{1}(v, :) - F{2}(v, :)), [], 'all') / max(abs(F{2}(v, :)), [], 'all'), 0.1);
end

%% Formats, channels, cuts
function test_formats_give_the_same_connectivity(tc)
% one recording written as EEGLAB, BrainVision (both float32), EDF (16 bit)
% and BDF (24 bit; read without FieldTrip's mex file where it has none for
% the platform, read_24bit.m)
S = tc.TestData.S;
E = load(fullfile(S.data, 'sub-control01', 'eeg', 'sub-control01_task-rest_eeg.set'), '-mat');
X = double(E.EEG.data);  labels = {E.EEG.chanlocs.labels};
d = fullfile(tc.TestData.tmp, 'formats');
fmts = {'set', 'vhdr', 'edf', 'bdf'};
for i = 1:numel(fmts)
    sim.write_eeg(fullfile(d, ['rec_' fmts{i}]), X, E.EEG.srate, labels, fmts{i});
end
S.results = fullfile(tc.TestData.tmp, 'formats_results');
conn = cell(1, numel(fmts));
for i = 1:numel(fmts)
    file = dir(fullfile(d, ['rec_' fmts{i} '.' fmts{i}]));
    qc = run_pipeline(file.name, d, fullfile(S.results, 'Connectivity'), S);
    verifyEqual(tc, qc.status, "ok", qc.message);
    w = load(output(S.results, ['rec_' fmts{i}]), 'static_full_cmats');
    conn{i} = w.static_full_cmats;
end
verifyEqual(tc, conn{2}, conn{1}, 'AbsTol', 1e-6);
verifyGreaterThan(tc, corr(conn{3}(:), conn{1}(:)), 0.999);
verifyGreaterThan(tc, corr(conn{4}(:), conn{1}(:)), 0.9999);
end

function test_flat_channel_and_channel_subset(tc)
% a flat channel is left out and reported; a 19-channel subset is used
% with a warning (fewer than 32 channels)
S = tc.TestData.S;
E = load(fullfile(S.data, 'sub-control01', 'eeg', 'sub-control01_task-rest_eeg.set'), '-mat');
X = double(E.EEG.data);  labels = {E.EEG.chanlocs.labels};
X(strcmp(labels, 'Cz'), :) = 0;
d = fullfile(tc.TestData.tmp, 'channels');
sim.write_eeg(fullfile(d, 'flat'), X, E.EEG.srate, labels, 'set');
S.results = fullfile(tc.TestData.tmp, 'channels_results');
out = fullfile(S.results, 'Connectivity');
qc = run_pipeline('flat.set', d, out, S);
verifyEqual(tc, qc.status, "ok", qc.message);
verifyEqual(tc, qc.flat_channels, "Cz");
verifyEqual(tc, qc.n_channels, 63);
S.channel_subset = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'C4', 'T8', 'P7', ...
    'P3', 'Pz', 'P4', 'P8', 'O1', 'O2', 'Cz'};
S.drop_flat_channels = false;
qc = verifyWarning(tc, @() run_pipeline('flat.set', d, out, S), 'run_pipeline:few_channels');
verifyEqual(tc, qc.n_channels, 19);
S.min_channels = 20;
S.stop_on_error = false;
qc = run_pipeline('flat.set', d, fullfile(tc.TestData.tmp, 'channels_min'), S);
verifyEqual(tc, qc.status, "failed");
verifySubstring(tc, char(qc.message), 'minimum 20');
end

function test_data_rank(tc)
% the rank of the channels used (stage1_summary.csv: data_rank): the
% recording with 20 of its 64 dimensions removed (as ICA removes
% components) has rank 44; with 50 removed, rank 14 and a warning (below
% S.min_channels, 16). Outputs whose summary row has no rank get it when
% they are reused: the recording is read, nothing is recomputed
S = tc.TestData.S;
E = load(fullfile(S.data, 'sub-control01', 'eeg', 'sub-control01_task-rest_eeg.set'), '-mat');
X = double(E.EEG.data);  labels = {E.EEG.chanlocs.labels};
rs = RandStream('mrg32k3a', 'Seed', 3);
[U, ~] = qr(randn(rs, size(X, 1)));                % a random orthonormal basis
d = fullfile(tc.TestData.tmp, 'rank');
for k = [20, 50]
    P = U(:, k + 1:end) * U(:, k + 1:end)';         % k dimensions removed
    sim.write_eeg(fullfile(d, sprintf('removed%d', k)), P * X, E.EEG.srate, labels, 'set');
end
S.results = fullfile(tc.TestData.tmp, 'rank_results');
out = fullfile(S.results, 'Connectivity');
qc = run_pipeline('removed20.set', d, out, S);
verifyEqual(tc, qc.status, "ok", qc.message);
verifyEqual(tc, qc.n_channels, 64);
verifyEqual(tc, qc.data_rank, 44);
verifyEqual(tc, S.min_channels, 16);
qc = verifyWarning(tc, @() run_pipeline('removed50.set', d, out, S), 'run_pipeline:low_rank');
verifyEqual(tc, qc.data_rank, 14);
f = output(S.results, 'removed20');
s = load(f, 'info');
verifyEqual(tc, s.info.qc.data_rank, 44);
s.info.qc = rmfield(s.info.qc, 'data_rank');       % an output without the rank
save(f, '-struct', 's', '-append');
txt = evalc('qc = run_pipeline(''removed20.set'', d, out, S);');
verifyEqual(tc, qc.status, "skipped (up to date)");
verifyEqual(tc, qc.data_rank, 44);
verifySubstring(tc, txt, 'data rank 44 of 64 channels added');
w = load(f, 'info');
verifyEqual(tc, w.info.qc.data_rank, 44);           % kept with the output
verifyEqual(tc, w.info.settings_key, s.info.settings_key);
txt = evalc('qc = run_pipeline(''removed20.set'', d, out, S);');
verifyEqual(tc, qc.data_rank, 44);
verifyFalse(tc, contains(txt, 'added'));            % read once only
end

function test_boundary_cuts_the_recording(tc)
% an EEGLAB boundary event at 30 s: two pieces, each trimmed at both ends
S = tc.TestData.S;
E = load(fullfile(S.data, 'sub-control01', 'eeg', 'sub-control01_task-rest_eeg.set'), '-mat');
EEG = E.EEG;
EEG.event = struct('type', 'boundary', 'latency', 30 * EEG.srate + 0.5, 'duration', 0);
d = fullfile(tc.TestData.tmp, 'boundary');
mkdir(d);
EEG.filepath = d;  EEG.filename = 'cut.set';
save(fullfile(d, 'cut.set'), 'EEG', '-v7');
S.results = fullfile(tc.TestData.tmp, 'boundary_results');
qc = run_pipeline('cut.set', d, fullfile(S.results, 'Connectivity'), S);
verifyEqual(tc, qc.status, "ok", qc.message);
verifyEqual(tc, qc.n_boundaries, 1);
verifyEqual(tc, qc.n_pieces, 2);
verifyEqual(tc, qc.analysed_s, 56, 'AbsTol', 0.02);
w = load(output(S.results, 'cut'), 'info');
verifyEqual(tc, w.info.piece_onset_s(:)', [1, 31], 'AbsTol', 0.01);
end

function test_measures_through_stage_1(tc)
% seven measures through run_pipeline (the symmetric orthogonalisation
% needs fewer regions than the 63 independent channels:
% test_static_estimators): the stored matrices, halves and power equal
% static_connectivity of the stored region time series after
% bandpass_pieces, so stage 1 passes every measure the signal, band, pieces
% and dwPLI epochs, and keeps their order (test_static_estimators checks the
% measures, tools/compare_with_brainstorm.m compares them with Brainstorm)
S = tc.TestData.S;
S.conn = {'plv', 'wpli', 'pli', 'ciplv', 'dwpli', 'aec', 'aec_orth'};
S.save_scouts = true;
S.data = fullfile(tc.TestData.tmp, 'measures');
copy_recording(tc.TestData.S.data, 'sub-control01', S.data);
S.results = fullfile(tc.TestData.tmp, 'measures_results');
S = check_settings(S);
evalc('run_all(S);');
rec = 'sub-control01_task-rest_eeg';
w = load(output(S.results, rec));
sc = load(fullfile(S.results, 'Connectivity', 'Scouts time-series', rec, ...
    'Scouts_eloreta_ASA_10-05_343_Desikan_Killiany_68.mat'));
verifyEqual(tc, w.info.conn, S.conn);
[Z, seg] = bandpass_pieces(sc.src_data, sc.srate, 7, 13, sc.segments, sc.onset_s, S);
verifyEqual(tc, seg, w.info.segments);
[M, P] = static_connectivity(Z, seg, S.conn, sc.srate, 7, 13);
verifyEqual(tc, w.static_full_cmats, M, 'AbsTol', 1e-12);
verifyEqual(tc, w.static_power, P, 'RelTol', 1e-12);
for h = 1:2
    verifyEqual(tc, squeeze(w.static_halves(h, :, :, :)), static_connectivity(Z, w.info.halves{h}, S.conn, sc.srate, 7, 13), ...
        'AbsTol', 1e-12);
end
low = tril(true(68), -1);
ci = squeeze(w.static_full_cmats(4, :, :));
verifyTrue(tc, all(ci(low) >= 0 & ci(low) <= 1));                  % ciPLV in [0, 1]
dw = squeeze(w.static_full_cmats(5, :, :));
verifyTrue(tc, all(dw(low) <= 1) && any(dw(low) < 0));            % dwPLI: at most 1, can be negative
end

function test_scout_functions_equal_their_definitions(tc)
% get_scouts_timeseries against the definitions, on random filters and
% data: each region's sources, sign-flipped to the region's dominant
% orientation (Brainstorm's rule), combined by their mean, first principal
% component (of their covariance, sign: weights summing to a positive
% number, scaled by 1/sqrt(number of sources)), signed maximum, power and
% RMS; the linear ones (mean, PCA) as an operator on the filters too
S = tc.TestData.S;
in = fullfile(S.root, '1- EEG source connectivity', 'inputs');
C = load(fullfile(in, 'tess_cortex_pial_low.mat'), 'VertNormals');
A = load(fullfile(in, 'scout_Desikan-Killiany_68.mat'));
rs = RandStream('mrg32k3a', 'Seed', 7);
nc = 12;  T = 300;
F = randn(rs, size(C.VertNormals, 1), nc);
X = randn(rs, nc, T) + 3;                          % an offset: the PCA uses the (demeaned) covariance
Xc = X - mean(X, 2);  Cx = Xc * Xc' / T;
fun = {'mean', 'pca', 'max', 'power', 'rms'};
Y = struct();
for k = 1:numel(fun)
    [Y.(fun{k}), labels] = get_scouts_timeseries(F, 'Desikan_Killiany_68', Cx, fun{k}, X);
end
verifyEqual(tc, labels, {A.Scouts.Label});
for r = [1, 17, 40, 68]
    v = A.Scouts(r).Vertices;
    [U, ~] = svd(C.VertNormals(v, :), 'econ');
    m = sign(U(:, 1));
    if nnz(m > 0) < nnz(m < 0), m = -m; end
    m(m == 0) = 1;
    Ys = m .* (F(v, :) * X);                       % the sign-flipped sources
    verifyEqual(tc, Y.mean(r, :), mean(Ys, 1), 'AbsTol', 1e-10);
    verifyEqual(tc, Y.power(r, :), mean(Ys .^ 2, 1), 'AbsTol', 1e-9);
    verifyEqual(tc, Y.rms(r, :), sqrt(mean(Ys .^ 2, 1)), 'AbsTol', 1e-10);
    [~, j] = max(abs(Ys), [], 1);
    verifyEqual(tc, Y.max(r, :), Ys(sub2ind(size(Ys), j, 1:T)), 'AbsTol', 1e-10);
    Yc = Ys - mean(Ys, 2);
    [V, D] = eig(Yc * Yc' / T, 'vector');
    [~, j] = max(D);  w = V(:, j);
    if sum(w) < 0, w = -w; end
    verifyEqual(tc, Y.pca(r, :), w' * Ys / sqrt(numel(v)), 'AbsTol', 1e-8 * max(abs(Y.pca(r, :))));
end
for f = {'mean', 'pca'}                            % operators: the same on the filters, then the data
    verifyEqual(tc, get_scouts_timeseries(F, 'Desikan_Killiany_68', Cx, f{1}) * X, Y.(f{1}), 'AbsTol', 1e-9);
end
verifyEqual(tc, get_scouts_timeseries(F, 'Desikan_Killiany_68'), ...
    get_scouts_timeseries(F, 'Desikan_Killiany_68', [], 'mean'));    % the default: the mean
verifyError(tc, @() get_scouts_timeseries(F, 'Desikan_Killiany_68', Cx, 'max'), 'get_scouts_timeseries:data');
verifyError(tc, @() get_scouts_timeseries(F, 'Desikan_Killiany_68', [], 'pca'), 'get_scouts_timeseries:pca');
verifyError(tc, @() get_scouts_timeseries(F, 'Desikan_Killiany_68', Cx, 'median', X), 'get_scouts_timeseries:function');
end

function test_scout_functions_through_stage_1(tc)
% S.scout_function through run_pipeline on one recording: each function
% has its own settings key, recorded in info; the stored region time
% series keep the definitions' relations (power = RMS^2; |max| >= |mean|
% at every sample; the PCA's variance at least the mean's, region by
% region) and the matrices are computed from them. With windows, the
% same whole-recording matrices and windowed ones from the same signals
S = tc.TestData.S;
S.save_scouts = true;
S.data = fullfile(tc.TestData.tmp, 'scouts');
copy_recording(tc.TestData.S.data, 'sub-control01', S.data);
rec = 'sub-control01_task-rest_eeg';
scouts = @(R) load(fullfile(R, 'Connectivity', 'Scouts time-series', rec, ...
    'Scouts_eloreta_ASA_10-05_343_Desikan_Killiany_68.mat'));
fun = {'mean', 'pca', 'max', 'power', 'rms'};
src = struct();  keys = strings(1, numel(fun));  full = struct();
for k = 1:numel(fun)
    S2 = S;
    S2.scout_function = fun{k};
    S2.results = fullfile(tc.TestData.tmp, ['scouts_' fun{k}]);
    evalc('S2 = check_settings(S2); run_all(S2);');   % power, RMS: a warning (phase measures)
    w = load(output(S2.results, rec));
    verifyEqual(tc, w.info.scout_function, fun{k});
    keys(k) = string(w.info.settings_key);
    sc = scouts(S2.results);
    src.(fun{k}) = sc.src_data;
    full.(fun{k}) = w.static_full_cmats;
    [Z, seg] = bandpass_pieces(sc.src_data, sc.srate, 7, 13, sc.segments, sc.onset_s, S2);
    verifyEqual(tc, w.static_full_cmats, static_connectivity(Z, seg, S2.conn, sc.srate, 7, 13), 'AbsTol', 1e-12);
end
verifyEqual(tc, numel(unique(keys)), numel(fun));
verifyEqual(tc, src.power, src.rms .^ 2, 'AbsTol', 1e-9 * max(src.power(:)));
verifyTrue(tc, all(src.rms(:) >= 0));
verifyTrue(tc, all(abs(src.max(:)) >= abs(src.mean(:)) - 1e-9 * max(abs(src.max(:)))));
verifyTrue(tc, all(var(src.pca, 1, 2) >= var(src.mean, 1, 2) * (1 - 1e-9)));
verifyGreaterThan(tc, mean(var(src.pca, 1, 2) ./ var(src.mean, 1, 2)), 1);
verifyNotEqual(tc, full.pca, full.mean);
if isempty(stage_folder(S.root, 'Dynamic connectivity')), return, end   % no windows in this pipeline
S2 = S;
S2.scout_function = 'max';
S2.dynamic = true;
S2.results = fullfile(tc.TestData.tmp, 'scouts_max_windows');
evalc('S2 = check_settings(S2); run_all(S2);');
w = load(output(S2.results, rec));
verifyEqual(tc, w.static_full_cmats, full.max, 'AbsTol', 1e-12);    % whole recording: the same
c = load(fullfile(S2.results, 'Connectivity', 'cmats', rec, 'cmats_eloreta_ASA_10-05_343_Desikan_Killiany_68_alpha.mat'));
verifyEqual(tc, c.info.scout_function, 'max');
verifyTrue(tc, all(cellfun(@(m) all(isfinite(m(:))), c.cmats)));
verifyEqual(tc, scouts(S2.results).src_data, src.max, 'AbsTol', 1e-12);
end

function copy_recording(from, id, to)
src = fullfile(from, id);
mkdir(fullfile(to, id));
copyfile(src, fullfile(to, id));
end
