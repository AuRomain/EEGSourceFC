function tests = test_reference_connectivity
%TEST_REFERENCE_CONNECTIVITY Stage 1's connectivity measures against verified toolboxes.
%   From the pipeline folder:  results = runtests('tests/test_reference_connectivity.m');  table(results)
%   Simulated region signals (10 regions, 90 s at 250 Hz, alpha 8-12 Hz
%   with lagged, zero-lag and envelope coupling, a broadband background).
%   This pipeline's measures (get_all_connectivity_mats: the windows, cmats,
%   and the whole recording, static_full_cmats = static_connectivity) are
%   computed from the same analytic signal as the reference's, and compared
%   with:
%   - Brainstorm (the copy found by setup_toolboxes, its whole tree):
%     PLV, ciPLV, wPLI, PLI from its NxN phase-synchronisation process
%     (process_plv1n -> bst_connectivity, Hilbert, magnitude) over the whole
%     recording and in windows; AEC from bst_henv ('penv'), whole recording
%     and windows. Brainstorm band-passes the broadband signals itself
%     (process_bandpass 'bst-hfilter-2019', no mirror, then hilbert); this
%     pipeline is given hilbert of the same filter's output
%     (bst_bandpass_hfilter). Brainstorm runs headless with a temporary user
%     folder and database (as tools/compare_with_brainstorm.m), deleted with
%     the path restored at the end of setupOnce: ~/.brainstorm and the
%     user's own databases are not used.
%   - FieldTrip: AEC against ft_connectivityanalysis 'amplcorr', each
%     sample of the analytic signal taken as one repetition of a
%     one-frequency spectrum.
%   - ROInets (OHBA, Colclough et al. 2015; as used in var-EEG-FC, Allouch
%     et al. 2023, in ../Toolboxes/GitHub resources): pairwise
%     orthogonalised AEC against do_pairwise_calculation (Brookes et al.
%     2012), symmetric orthogonalisation against remove_source_leakage
%     'closest', plain PLV, PLI and AEC against its plv, pli and aec.
%   - FieldTrip: dwPLI (the debiased wPLI of Vinck et al. 2011, epochs of
%     6 cycles as trials) against ft_freqanalysis and
%     ft_connectivityanalysis 'wpli_debiased' on the same epochs, over the
%     whole recording and in each 60-cycle window; and how FieldTrip's dwPLI of
%     the band-passed signals (stage 1's) compares with that of the
%     broadband ones.
%   - Brainstorm's orthogonalised envelope correlation ('oenv', Hipp et al.
%     2012) is another estimator: the documented relation between the two
%     is tested (test_brainstorm_oenv_relation).
%   Tolerances: 1e-12 in double precision (the whole recording); the
%   windows are stored in single precision (cmats): 6e-8, one unit in
%   the last place of a single in [0.5, 1).
%   A reference toolbox that is not installed skips its tests (Incomplete,
%   with the reason); so does a Brainstorm that cannot start (e.g. a git
%   copy without java/brainstorm.jar, which Brainstorm downloads at its
%   first start: tools/fetch_toolboxes.sh gets it), without stopping the
%   FieldTrip and ROInets tests.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
S = static_test_settings();
tc.TestData.S = setup_toolboxes(S);
[X, fs, labels] = simulate();
band = [8 12];
Xb = bst_bandpass_hfilter(X, fs, band(1), band(2), 0, 0, [], [], 'bst-hfilter-2019');   % Brainstorm's connectivity filter
Z = hilbert(Xb.').';
tc.TestData.X = X;  tc.TestData.fs = fs;  tc.TestData.labels = labels;  tc.TestData.band = band;
tc.TestData.Z = Z;
% this pipeline, as stage 1 calls it with the default windows: 6 cycles of
% the band centre (0.6 s) every 0.3 s for the phase measures, 6 s every
% 0.5 s for AEC, 60 cycles (6 s) every 30 cycles for dwPLI (epochs of 6 cycles)
conn = {'plv', 'ciplv', 'wpli', 'pli', 'aec', 'aec_orth', 'aec_symorth', 'plv_symorth', 'dwpli'};
opts = struct('conn', {conn}, 'window_cycles', 6, 'window_length', [], 'window_overlap', 0.5, ...
    'aec_window', 6, 'aec_step', 0.5);
[P.cmats, ~, P.full, P.info] = get_all_connectivity_mats(Z, fs, band(1), band(2), [1, size(Z, 2)], opts);
P.conn = conn;
tc.TestData.P = P;
tc.TestData.off = ~eye(size(X, 1));
tc.TestData.roinets = fullfile(tc.TestData.S.toolboxes, 'GitHub resources', 'var-EEG-FC-main', ...
    'var-EEG-FC-main', '+ROInets');
tc.TestData.B = brainstorm_reference(X, fs, band, labels, P.info);
end

%% Brainstorm: phase measures
function test_phase_measures_equal_brainstorm(tc)
% PLV, ciPLV, wPLI, PLI over the whole recording (static_full_cmats) and in
% the 299 windows of 0.6 s (cmats) against process_plv1n -> bst_connectivity
B = brainstorm_or_skip(tc);
P = tc.TestData.P;  off = tc.TestData.off;
for m = {'plv', 'ciplv', 'wpli', 'pli'}
    k = find(strcmp(P.conn, m{1}));
    A = squeeze(P.full(k, :, :));
    R = B.(m{1}).none;
    report(m{1}, 'whole recording', A(off), R(off));
    verifyEqual(tc, A(off), R(off), 'AbsTol', 1e-12, [m{1} ', whole recording']);
    W = permute(double(P.cmats{k}), [2 3 1]);           % nROI x nROI x nWindows
    R = B.(m{1}).windows;
    verifyEqual(tc, size(W), size(R), [m{1} ': number of windows']);
    w = P.info.windows{k};
    verifyEqual(tc, (w.centre_sample - 1) / tc.TestData.fs, B.(m{1}).time, 'AbsTol', 1e-9, ...
        [m{1} ': window centres']);
    offw = repmat(off, 1, 1, size(W, 3));
    report(m{1}, sprintf('%d windows', size(W, 3)), W(offw), R(offw));
    verifyEqual(tc, W(offw), R(offw), 'AbsTol', 6e-8, [m{1} ', windows (single precision)']);
end
end

%% AEC: Brainstorm, FieldTrip
function test_aec_equals_brainstorm_penv(tc)
% plain AEC (static_full_cmats and the 169 windows of 6 s) against bst_henv
% 'penv', which filters the broadband signals itself (bst_timefreq,
% Hilbert: the same filter and hilbert as this pipeline's input here)
B = brainstorm_or_skip(tc);
P = tc.TestData.P;  off = tc.TestData.off;
k = find(strcmp(P.conn, 'aec'));
A = squeeze(P.full(k, :, :));
report('aec', 'whole recording (bst_henv penv)', A(off), B.penv.none(off));
verifyEqual(tc, A(off), B.penv.none(off), 'AbsTol', 1e-12);
W = permute(double(P.cmats{k}), [2 3 1]);
verifyEqual(tc, size(W), size(B.penv.windows));
offw = repmat(off, 1, 1, size(W, 3));
report('aec', sprintf('%d windows (bst_henv penv)', size(W, 3)), W(offw), B.penv.windows(offw));
verifyEqual(tc, W(offw), B.penv.windows(offw), 'AbsTol', 6e-8);
end

function test_aec_equals_fieldtrip_amplcorr(tc)
% plain AEC against FieldTrip's amplitude correlation: the analytic
% signal's samples as repetitions of a one-frequency Fourier spectrum
% (ft_connectivityanalysis 'amplcorr': |F| demeaned over repetitions,
% covariance, normalised)
assumeNotEmpty(tc, which('ft_connectivityanalysis'), 'FieldTrip not on the path');
Z = tc.TestData.Z;  [n, N] = size(Z);  off = tc.TestData.off;
freq = struct('label', {tc.TestData.labels}, 'freq', 10, 'dimord', 'rpttap_chan_freq', ...
    'fourierspctrm', reshape(Z.', N, n, 1), 'cumtapcnt', ones(N, 1), 'cumsumcnt', ones(N, 1));
cfg = struct('method', 'amplcorr', 'removemean', 'yes');
C = ft_connectivityanalysis(cfg, freq);
R = C.amplcorrspctrm;
k = find(strcmp(tc.TestData.P.conn, 'aec'));
A = squeeze(tc.TestData.P.full(k, :, :));
report('aec', 'whole recording (FieldTrip amplcorr)', A(off), R(off));
verifyEqual(tc, A(off), R(off), 'AbsTol', 1e-12);
end

%% FieldTrip: debiased wPLI
function test_dwpli_equals_fieldtrip_wpli_debiased(tc)
% dwPLI is FieldTrip's debiased wPLI (Vinck et al. 2011) with epochs as
% trials, here computed through FieldTrip's own entry points on the same
% epochs of the band-passed signals: ft_freqanalysis (mtmfft, Hanning, the
% Fourier frequencies of the epochs inside the band: 8.3, 10, 11.7 Hz) and
% ft_connectivityanalysis 'wpli_debiased' (this pipeline calls
% ft_connectivity_wpli with 'debias'), then the mean over the frequencies.
% Over the whole recording (static_full_cmats: epochs of 6 cycles, 0.6 s,
% one after the other) and in each window of 60 cycles every 30 (6 s
% every 3 s, 10 epochs each: cmats, single precision)
assumeNotEmpty(tc, which('ft_connectivityanalysis'), 'FieldTrip not on the path');
Z = tc.TestData.Z;  fs = tc.TestData.fs;  band = tc.TestData.band;  N = size(Z, 2);
off = tc.TestData.off;  labels = tc.TestData.labels;
P = tc.TestData.P;  k = find(strcmp(P.conn, 'dwpli'));
L = round(6 / mean(band) * fs);
R = fieldtrip_dwpli(real(Z), 1:L:N - L + 1, L, band, fs, labels);
A = squeeze(P.full(k, :, :));
report('dwpli', 'whole recording (FieldTrip wpli_debiased)', A(off), R(off));
verifyEqual(tc, A(off), R(off), 'AbsTol', 1e-12);
w = P.info.windows{k};
verifyEqual(tc, [w.samples, w.step, w.epochs, w.epoch_samples], [max(ceil(6 * fs), 10 * L), ceil(3 * fs), 10, L]);
verifyEqual(tc, numel(w.first_sample), floor((N - w.samples) / w.step) + 1);
d = zeros(numel(w.first_sample), 1);
for t = 1:numel(w.first_sample)
    R = fieldtrip_dwpli(real(Z), w.first_sample(t) + (0:w.epochs - 1) * L, L, band, fs, labels);
    W = squeeze(double(P.cmats{k}(t, :, :)));
    d(t) = max(abs(W(off) - R(off)));
    verifyEqual(tc, W(off), R(off), 'AbsTol', 6e-8, sprintf('window %d (single precision)', t));
end
fprintf('dwpli        %d windows (FieldTrip wpli_debiased)          max |difference| %.2e\n', numel(d), max(d));
end

function test_dwpli_band_passed_and_broadband(tc)
% FieldTrip's usual dwPLI takes the broadband signals; stage 1 gives it
% the band-passed ones. At a frequency inside the band the filter
% (zero-phase, the same gain in every epoch) changes the Fourier
% coefficients of both signals by the same factor, which the debiased wPLI
% ignores, up to the Hanning taper's leakage from neighbouring frequencies
% (the filter's edges are 1.7 Hz from the epochs' outer frequencies): the
% same lagged pairs stand out, near 0 elsewhere, the edges agree. Epochs of
% 6 cycles over 10-80 s of the simulation (middle_stretch)
assumeNotEmpty(tc, which('ft_connectivityanalysis'), 'FieldTrip not on the path');
X = tc.TestData.X;  Z = tc.TestData.Z;  fs = tc.TestData.fs;  band = tc.TestData.band;  n = size(X, 1);
[ii, jj] = find(tril(true(n), -1));
lin = sub2ind([n n], ii, jj);
lagged = ismember([ii, jj], [2 1; 3 2; 8 7], 'rows');  % 2 follows 1 (and 3, which holds 1); 8 follows 7
seg = middle_stretch(fs);
L = round(6 / mean(band) * fs);
first = seg(1):L:seg(2) - L + 1;
M = fieldtrip_dwpli(real(Z), first, L, band, fs, tc.TestData.labels);   % = this pipeline's (test above)
R = fieldtrip_dwpli(X, first, L, band, fs, tc.TestData.labels);
M = M(lin);  R = R(lin);
fprintf(['dwpli        band-passed vs broadband (FieldTrip, %d epochs): r %.3f, max |difference| %.3f; ' ...
    'lagged pairs %.2f-%.2f vs %.2f-%.2f; others |max| %.3f vs %.3f\n'], numel(first), corr(M, R), ...
    max(abs(M - R)), min(M(lagged)), max(M(lagged)), min(R(lagged)), max(R(lagged)), max(abs(M(~lagged))), ...
    max(abs(R(~lagged))));
verifyGreaterThan(tc, corr(M, R), 0.99);                % r 1.000, max |difference| 0.027
verifyLessThan(tc, max(abs(M - R)), 0.05);
verifyGreaterThan(tc, min(M(lagged)), 3 * max(abs(M(~lagged))));
verifyGreaterThan(tc, min(R(lagged)), 3 * max(abs(R(~lagged))));
end

%% ROInets: pairwise and symmetric orthogonalisation
function test_aec_orth_equals_roinets_pairwise(tc)
% 'aec_orth' (static_connectivity and the whole-recording output of
% aec_sliding_window_corrected_pairwise) against ROInets'
% do_pairwise_calculation (Brookes et al. 2012): for every seed x, every
% other signal y orthogonalised to x (Householder QR, = y - (x'y / x'x) x),
% envelopes correlated: envCorrelation(x, y) = corr(env x, env(y orth x)).
% This pipeline averages the two directions: (E + E') / 2.
% The same analytic signal: ROInets' envelope_data demeans, takes the
% Hilbert transform zero-padded to a power of 2 and drops the first sample
% ("1st sample bias"); its moving-average downsampling is switched off
% (windowLength 0). So the signals are given to ROInets with a zero sample
% in front (it adds nothing to the regressions and is the sample dropped),
% and this pipeline gets the analytic signal ROInets takes, without that
% sample.
c = use_reference_roinets(tc); %#ok<NASGU>
Xr = real(tc.TestData.Z);
Xr = Xr - mean(Xr, 2);
[n, N] = size(Xr);
fs = tc.TestData.fs;  off = tc.TestData.off;
Xe = [zeros(n, 1), Xr];
nfft = 2 ^ nextpow2(N + 2);                        % as envelope_data: nextpow2(nSamples + 1)
Ze = hilbert(Xe.', nfft).';
Zp = Ze(:, 2:N + 1);
env = struct('windowLength', 0, 'overlap', 0, 'useHanningWindow', false, 'useFilter', false);
t = (0:N) / fs;
[~, mats] = evalc('ROInets.do_pairwise_calculation(Xe, t, env)');   % its progress lines not shown
E = mats.envCorrelation;
R = (E + E.') / 2;
% for the windows: ROInets' orthogonalisation (over the whole recording, as
% this pipeline's) and envelopes of every ordered pair
envx = ROInets.envelope_data(Xe, t, env);
envo = zeros(n, n, N);                              % envo(i, j, :): envelope of j orthogonalised to i
for i = 1:n
    for j = [1:i - 1, i + 1:n]
        o = ROInets.householder_orthogonalise([Xe(i, :)', Xe(j, :)']);
        envo(i, j, :) = ROInets.envelope_data(o(:, 2)', t, env);
    end
end
clear c                                             % this pipeline's own path again
M = static_connectivity(Zp, [], {'aec_orth'}, fs, 8, 12);
A = squeeze(M);
[Wp, ~, F, win] = connectivity.aec_sliding_window_corrected_pairwise(Zp, fs, 6, 0.5, []);
report('aec_orth', 'whole recording (ROInets pairwise)', A(off), R(off));
verifyEqual(tc, A(off), R(off), 'AbsTol', 1e-12, 'static_connectivity');
verifyEqual(tc, F(off), R(off), 'AbsTol', 1e-12, 'aec_sliding_window_corrected_pairwise');
Et = E.';
verifyGreaterThan(tc, max(abs(E(off) - Et(off))), 1e-3);   % ROInets' directions differ: the average matters
nw = numel(win.first_sample);
Rw = zeros(n, n, nw);
for w = 1:nw
    idx = win.first_sample(w) + (0:win.samples - 1);
    for i = 1:n
        for j = [1:i - 1, i + 1:n]
            Rw(i, j, w) = corr(envx(i, idx)', squeeze(envo(i, j, idx)));
        end
    end
end
Rw = (Rw + permute(Rw, [2 1 3])) / 2;
W = permute(double(Wp), [2 3 1]);
offw = repmat(off, 1, 1, nw);
report('aec_orth', sprintf('%d windows (ROInets envelopes)', nw), W(offw), Rw(offw));
verifyEqual(tc, W(offw), Rw(offw), 'AbsTol', 6e-8, 'windows');
end

function test_symmetric_orthogonalisation_equals_roinets(tc)
% 'aec_symorth' and 'plv_symorth' (whole recording and windows) against
% ROInets: remove_source_leakage(signals, 'closest') (closest orthogonal
% matrix, Colclough et al. 2015), hilbert of the orthogonalised signals,
% ROInets' aec and plv (whole recording and the same windows). Plain AEC,
% PLV and PLI against ROInets' aec, plv and pli. And the pipeline's copy
% of the ROInets functions it uses is the toolbox's, unchanged.
Z = tc.TestData.Z;  P = tc.TestData.P;  off = tc.TestData.off;
src = tc.TestData.roinets;
assumeTrue(tc, isfolder(src), sprintf('ROInets (var-EEG-FC repository) not found: %s', src));
mine = fullfile(tc.TestData.S.root, '1- EEG source connectivity', '+ROInets');
for f = dir(fullfile(mine, '*.m'))'
    verifyEqual(tc, fileread(fullfile(mine, f.name)), fileread(fullfile(src, f.name)), ...
        ['the pipeline''s +ROInets/' f.name ' differs from the toolbox''s']);
end
c = use_reference_roinets(tc); %#ok<NASGU>
L = ROInets.remove_source_leakage(real(Z), 'closest');
Zo = hilbert(L.').';
ref.aec_symorth = ROInets.aec(abs(Zo));
ref.plv_symorth = ROInets.plv(angle(Zo));
ref.aec = ROInets.aec(abs(Z));
ref.plv = ROInets.plv(angle(Z));
ref.pli = ROInets.pli(angle(Z));
% the windows of the pipeline (info.windows), on the orthogonalised signals
for m = {'aec_symorth', 'plv_symorth'}
    w = P.info.windows{strcmp(P.conn, m{1})};
    R = zeros(size(Z, 1), size(Z, 1), numel(w.first_sample));
    for t = 1:numel(w.first_sample)
        idx = w.first_sample(t) + (0:w.samples - 1);
        if strcmp(m{1}, 'aec_symorth')
            R(:, :, t) = ROInets.aec(abs(Zo(:, idx)));
        else
            R(:, :, t) = ROInets.plv(angle(Zo(:, idx)));
        end
    end
    ref.([m{1} '_windows']) = R;
end
clear c
G = L * L.';
G = G ./ sqrt(diag(G) * diag(G).');
verifyLessThan(tc, max(abs(G(off))), 1e-10);       % orthogonal signals
for m = {'aec_symorth', 'plv_symorth', 'aec', 'plv', 'pli'}
    k = strcmp(P.conn, m{1});
    A = squeeze(P.full(k, :, :));
    report(m{1}, 'whole recording (ROInets)', A(off), ref.(m{1})(off));
    verifyEqual(tc, A(off), ref.(m{1})(off), 'AbsTol', 1e-12, m{1});
end
for m = {'aec_symorth', 'plv_symorth'}
    W = permute(double(P.cmats{strcmp(P.conn, m{1})}), [2 3 1]);
    R = ref.([m{1} '_windows']);
    offw = repmat(off, 1, 1, size(W, 3));
    report(m{1}, sprintf('%d windows (ROInets)', size(W, 3)), W(offw), R(offw));
    verifyEqual(tc, W(offw), R(offw), 'AbsTol', 6e-8, [m{1} ', windows']);
end
end

%% Brainstorm's orthogonalised envelope correlation: another estimator
function test_brainstorm_oenv_relation(tc)
% Brainstorm's 'oenv' (bst_henv, HOrthCorr; Hipp et al. 2012) and this
% pipeline's 'aec_orth' (Brookes et al. 2012) both remove the zero-lag
% part of y that x explains and correlate the envelopes, both directions
% averaged; they differ in two ways, by design:
%  - the coefficient: Hipp's orthogonal part |Im(z_y conj(z_x))| / |z_x|
%    is |z_y - beta(t) z_x| with beta(t) = Re(z_y conj(z_x)) / |z_x|^2,
%    a coefficient of the complex analytic signals at every sample; this
%    pipeline uses one least-squares coefficient of the real signals over
%    the recording, beta = x'y / x'x;
%  - Brainstorm averages the absolute values of the two directions, this
%    pipeline the signed values (AEC keeps its sign).
% The same formula with each choice gives each toolbox (1e-12); the two
% agree where it matters: a zero-lag mixture is removed by both (plain AEC
% high, both orthogonalised values near 0), envelope coupling with
% independent phases is kept by both.
B = brainstorm_or_skip(tc);
Z = tc.TestData.Z;  P = tc.TestData.P;  off = tc.TestData.off;
A = squeeze(P.full(strcmp(P.conn, 'aec_orth'), :, :));
R = orthogonalised_aec(Z, 'regression', 'signed');
verifyEqual(tc, A(off), R(off), 'AbsTol', 1e-12, 'aec_orth = the formula with one regression coefficient');
H = orthogonalised_aec(Z, 'instantaneous', 'abs');
report('oenv', 'whole recording (formula, Brainstorm)', H(off), B.oenv.none(off));
verifyEqual(tc, B.oenv.none(off), H(off), 'AbsTol', 1e-12, 'Brainstorm oenv = the formula with beta(t), abs');
D = squeeze(P.full(strcmp(P.conn, 'aec'), :, :));
fprintf(['aec_orth vs Brainstorm oenv: r %.3f over the edges, max |difference| %.3f; ' ...
    'pair 1-3 (zero-lag mixture) aec %.3f, aec_orth %.3f, oenv %.3f; ' ...
    'pair 4-5 (envelope coupling) aec %.3f, aec_orth %.3f, oenv %.3f\n'], ...
    corr(A(off), B.oenv.none(off)), max(abs(A(off) - B.oenv.none(off))), ...
    D(1, 3), A(1, 3), B.oenv.none(1, 3), D(4, 5), A(4, 5), B.oenv.none(4, 5));
verifyGreaterThan(tc, D(1, 3), 0.5);                       % leakage-like coupling seen by plain AEC
verifyLessThan(tc, abs(A(1, 3)), 0.1);                     % removed by both
verifyLessThan(tc, B.oenv.none(1, 3), 0.1);
verifyGreaterThan(tc, A(4, 5), 0.3);                       % envelope coupling kept by both
verifyGreaterThan(tc, B.oenv.none(4, 5), 0.3);
end

%% Simulation and formulas
function [X, fs, labels] = simulate()
% 10 regions, 90 s at 250 Hz: alpha sources (8-12 Hz) coupled as
%   2 follows 1 with a lag of 16 ms (58 degrees at 10 Hz): phase coupling
%   3 = 1 + 0.5 x its own source: a zero-lag mixture (as source leakage)
%   4, 5 share a slow amplitude modulation, phases independent: envelope coupling
%   6 = 0.7 x region 4 + its own source: leakage of an envelope-coupled region
%   8 follows 7 with a lag of 24 ms; 9, 10 independent
% plus a broadband background (AR(1)), which the band-pass removes
rs = RandStream('mrg32k3a', 'Seed', 20261106);
fs = 250;  n = 10;  N = 90 * fs;
[b, a] = butter(4, [8 12] / (fs / 2));
s = filtfilt(b, a, randn(rs, N, n)).';
s = s ./ std(s, 0, 2);
[bl, al] = butter(2, 0.3 / (fs / 2));
m = filtfilt(bl, al, randn(rs, N, 1)).';
m = exp(0.6 * m / std(m));                          % slow positive modulation
lag = @(x, d) [zeros(1, d), x(1:end - d)];
X = zeros(n, N);
X(1, :) = s(1, :);
X(2, :) = 0.8 * lag(s(1, :), 4) + 0.6 * s(2, :);
X(3, :) = s(1, :) + 0.5 * s(3, :);
X(4, :) = m .* s(4, :);
X(5, :) = m .* s(5, :);
X(6, :) = 0.7 * X(4, :) + s(6, :);
X(7, :) = s(7, :);
X(8, :) = 0.7 * lag(s(7, :), 6) + 0.7 * s(8, :);
X(9, :) = s(9, :);
X(10, :) = s(10, :);
X = X + 0.2 * filter(1, [1 -0.95], randn(rs, N, n)).';
labels = arrayfun(@(k) sprintf('R%02d', k), (1:n)', 'UniformOutput', false);
end

function R = orthogonalised_aec(Z, coefficient, directions)
% F(i, j) = corr(|z_i|, |z_j - beta_ij z_i|), averaged over the two
% directions; beta_ij one least-squares coefficient of the real signals
% ('regression') or Re(z_j conj(z_i)) / |z_i|^2 at every sample
% ('instantaneous'); the directions signed or their absolute values
n = size(Z, 1);
env = abs(Z);
X = real(Z);
F = zeros(n);
for i = 1:n
    if strcmp(coefficient, 'regression')
        beta = (X * X(i, :)') / (X(i, :) * X(i, :)');
    else
        beta = real(Z .* conj(Z(i, :))) ./ env(i, :) .^ 2;
    end
    F(i, :) = corr(env(i, :)', abs(Z - beta .* Z(i, :))');
end
if strcmp(directions, 'abs'), F = abs(F); end
R = (F + F.') / 2;
R(1:n + 1:end) = 0;
end

function seg = middle_stretch(fs)
% 10-80 s of the 90 s simulation: the slow amplitude modulation of regions
% 4-6 starts with a large transient (the edge of its low-pass filter, then
% exp: region 4 reaches over 1000 times its median amplitude in the first
% second), which holds almost all of the products of pairs 4-5 and 5-6 over
% the whole simulation; stage 1 trims filter edges in the same way
seg = [10 * fs + 1, 80 * fs];
end

function R = fieldtrip_dwpli(X, first, L, band, fs, labels)
% FieldTrip's debiased wPLI of the epochs X(:, first(e) + (0:L-1)) as
% trials: ft_freqanalysis (mtmfft, Hanning) at the Fourier frequencies
% inside the band, ft_connectivityanalysis 'wpli_debiased', the mean over
% the frequencies (n x n; FieldTrip's own console lines not shown)
f = (ceil(band(1) * L / fs - 1e-9):floor(band(2) * L / fs + 1e-9)) * fs / L;
data = struct('label', {labels(:)}, 'fsample', fs, ...
    'trial', {arrayfun(@(s) X(:, s + (0:L - 1)), first, 'UniformOutput', false)}, ...
    'time', {repmat({(0:L - 1) / fs}, 1, numel(first))}, 'sampleinfo', [first(:), first(:) + L - 1]);
cfg = struct('method', 'mtmfft', 'taper', 'hanning', 'foi', f, 'output', 'fourier', 'keeptrials', 'yes', ...
    'feedback', 'no');
evalc('F = ft_freqanalysis(cfg, data);');
assert(isequal(size(F.freq), size(f)) && max(abs(F.freq - f)) < 1e-9, 'FieldTrip computed other frequencies');
evalc('C = ft_connectivityanalysis(struct(''method'', ''wpli_debiased''), F);');
R = mean(C.wpli_debiasedspctrm, 3);
end

function report(measure, what, a, b)
fprintf('%-12s %-40s max |difference| %.2e (values %.3f to %.3f)\n', measure, what, ...
    max(abs(double(a(:)) - double(b(:)))), min(b(:)), max(b(:)));
end

%% ROInets: the toolbox's copy first on the path while it computes
function c = use_reference_roinets(tc)
src = tc.TestData.roinets;
assumeTrue(tc, isfolder(src), sprintf('ROInets (var-EEG-FC repository) not found: %s', src));
tmp = tempname;
mkdir(tmp);
copyfile(src, fullfile(tmp, '+ROInets'));          % the package alone (the repository's other
addpath(tmp);                                       % folders would shadow this pipeline's functions)
c = onCleanup(@() remove_folder(tmp));
assertTrue(tc, startsWith(which('ROInets.closest_orthogonal_matrix'), tmp), 'ROInets is not the toolbox''s copy');
end

function remove_folder(tmp)
rmpath(tmp);
rmdir(tmp, 's');
end

%% Brainstorm, headless, in setupOnce
function B = brainstorm_or_skip(tc)
B = tc.TestData.B;
assumeEmpty(tc, B.skip, B.skip);                   % not installed, or cannot start here
assertEmpty(tc, B.error, ['Brainstorm failed: ' B.error]);
end

function B = brainstorm_reference(X, fs, band, labels, info)
% Brainstorm's measures of X, computed once; Brainstorm stopped, its
% temporary user folder deleted and the path restored afterwards. Never an
% error (an error would stop every test of this file, also FieldTrip's and
% ROInets'): B.skip says why Brainstorm's tests cannot run (not installed,
% no java/brainstorm.jar, already running, cannot start), B.error what
% failed once it had started
B = struct('skip', '', 'error', '');
bst = fileparts(fileparts(fileparts(which('bst_bandpass_hfilter'))));    % brainstorm3/toolbox/math
if isempty(bst) || ~isfile(fullfile(bst, 'brainstorm.m'))
    B.skip = 'Brainstorm''s whole tree (brainstorm.m) not found next to bst_bandpass_hfilter';
    return
end
if ~isfile(fullfile(bst, 'java', 'brainstorm.jar'))
    B.skip = sprintf(['%s has no java/brainstorm.jar (a git copy of brainstorm3 lacks it; Brainstorm ' ...
        'downloads it at its first start, from github.com/brainstorm-tools/bst-java): run ' ...
        'tools/fetch_toolboxes.sh, or start Brainstorm once (brainstorm nogui, then brainstorm stop)'], bst);
    return
end
p0 = path;
addpath(bst);
try
    running = brainstorm('status');
catch err
    path(p0);
    B.skip = sprintf('Brainstorm cannot start here: %s (%s)', err.message, err.identifier);
    return
end
if running                                         % not stopped by this test: the user's session
    path(p0);
    B.skip = 'Brainstorm is running in this MATLAB session: quit it first (brainstorm stop)';
    return
end
home = char(java.lang.System.getProperty('user.home'));
tmp = tempname;
mkdir(tmp);
c = onCleanup(@() stop_brainstorm(home, tmp, p0));
started = false;
try
    t0 = tic;
    java.lang.System.setProperty('user.home', tmp);    % Brainstorm's user folder: <tmp>/.brainstorm
    brainstorm('setpath');
    bst_startup(bst, -1, 'local', '');                 % server mode, database in <tmp>/.brainstorm/local_db
    gui_brainstorm('CreateProtocol', 'TestReferenceConnectivity', 0, 0);
    t_start = toc(t0);
    started = true;
    db_add_subject('S1');
    iStudy = db_add_condition('S1', 'simulated');
    sMat = db_template('matrixmat');
    sMat.Value = X;
    sMat.Time = (0:size(X, 2) - 1) / fs;
    sMat.Description = labels(:);
    sMat.Comment = 'simulated';
    file = db_add(iStudy, sMat);
    n = size(X, 1);
    k = find(strcmp(info.conn, 'plv'), 1);
    win_s = info.windows{k}.samples / fs;             % 0.6 s, 50% overlap
    for m = {'plv', 'ciplv', 'wpli', 'pli'}
        B.(m{1}).none = plv1n(file, m{1}, 'none', band, win_s, n);
        [B.(m{1}).windows, B.(m{1}).time] = plv1n(file, m{1}, 'windowed', band, win_s, n);
    end
    k = find(strcmp(info.conn, 'aec'), 1);
    aw = info.windows{k};
    opt = struct('SampleRate', fs, 'tfMeasure', 'hilbert', 'Freqs', {{'band', sprintf('%g, %g', band), 'mean'}}, ...
        'tfSplit', 1, 'isParallel', 0, 'WinLen', aw.samples / fs, 'WinOverlap', 1 - aw.step / aw.samples);
    t = sMat.Time;
    opt.CohMeasure = 'penv';  opt.TimeRes = 'none';
    B.penv.none = henv(X, t, opt);
    opt.TimeRes = 'windowed';
    B.penv.windows = henv(X, t, opt);
    opt.CohMeasure = 'oenv';  opt.TimeRes = 'none';
    B.oenv.none = henv(X, t, opt);
    fprintf('Brainstorm: started in %.0f s, measures computed in %.0f s\n', t_start, toc(t0) - t_start);
catch err
    if started
        B.error = sprintf('%s (%s)', err.message, err.identifier);
    else
        B.skip = sprintf('Brainstorm cannot start here: %s (%s)', err.message, err.identifier);
    end
end
end

function [R, time] = plv1n(file, method, timeres, band, win_s, n)
% process_plv1n -> bst_connectivity; R: n x n x time, magnitude, zero diagonal
args = {'plvmethod', method, 'plvmeasure', 2, 'tfmeasure', 'hilbert', 'timeres', timeres, ...
    'tfedit', struct('Freqs', {{'band', sprintf('%g, %g', band), 'mean'}}), 'outputmode', 'input'};
if strcmp(timeres, 'windowed'), args = [args, {'avgwinlength', win_s, 'avgwinoverlap', 50}]; end
out = bst_process('CallProcess', 'process_plv1n', {file}, [], args{:});
if isempty(out), error('Brainstorm returned no result (%s, %s)', method, timeres); end
F = in_bst_timefreq(out(1).FileName, 0);
TF = F.TF;
if size(TF, 1) == n * (n + 1) / 2, TF = process_compress_sym('Expand', TF, n); end
R = double(reshape(TF, n, n, size(TF, 2)));
for k = 1:n, R(k, k, :) = 0; end
time = F.Time;
end

function R = henv(X, t, opt) %#ok<INUSD> (used in evalc)
% bst_henv on the broadband signals (it band-passes them: bst_timefreq, Hilbert);
% its progress line not shown
[~, R] = evalc('bst_henv(X, X, t, opt)');
R = double(reshape(R, size(R, 1), size(R, 2), size(R, 3)));
end

function stop_brainstorm(home, tmp, p0)
try
    brainstorm('stop');
catch
end
java.lang.System.setProperty('user.home', home);
path(p0);                                          % Brainstorm's folders off the path again
if isfolder(tmp), rmdir(tmp, 's'); end
end
