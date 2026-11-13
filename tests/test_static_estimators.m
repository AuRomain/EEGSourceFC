function tests = test_static_estimators
%TEST_STATIC_ESTIMATORS The whole-recording measures against their formulas, and the split halves.
%   From the pipeline folder:  results = runtests('tests');  table(results)
%   Synthetic analytic signals only (no data): every measure of
%   static_connectivity.m computed again with a direct, slow formula, on
%   one continuous piece and on several pieces; known couplings; the
%   halves of split_halves.m; the whole-recording values of the windowed
%   functions of +connectivity and of get_all_connectivity_mats with and
%   without windows (S.dynamic). dwPLI is computed by FieldTrip (put on the
%   path by setup_toolboxes) and checked here against Vinck et al.'s
%   formula written directly.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, '1- EEG source connectivity'));      % static_connectivity, +connectivity
if isempty(which('ft_freqanalysis'))                              % dwPLI is computed by FieldTrip
    evalc('setup_toolboxes(static_test_settings());');
end
rs = RandStream('mrg32k3a', 'Seed', 3);
fs = 128;  n = 5;  N = 40 * fs;
[b, a] = butter(4, [8 12] / (fs / 2));
X = filtfilt(b, a, randn(rs, N, n)).';
X(2, :) = X(2, :) + 0.8 * [zeros(1, 3), X(1, 1:end - 3)];   % region 2 follows region 1 with a lag
X(3, :) = X(3, :) + 0.8 * X(1, :);                            % region 3: zero-lag copy (leakage)
tc.TestData.X = X;
tc.TestData.Z = hilbert(X.').';
tc.TestData.fs = fs;
tc.TestData.seg = [1, 900; 1201, 3000; 3301, N];             % three pieces
tc.TestData.all = {'plv', 'wpli', 'pli', 'ciplv', 'dwpli', 'aec', 'aec_orth', 'aec_symorth', 'plv_symorth'};
end

%% Measures against direct formulas
function test_measures_match_formulas(tc)
Z = tc.TestData.Z;  n = size(Z, 1);
M = static_connectivity(Z, [], {'plv', 'wpli', 'pli', 'ciplv', 'aec', 'aec_orth'}, tc.TestData.fs, 8, 12);
R = direct(Z);
for k = 1:6
    verifyEqual(tc, squeeze(M(k, :, :)), R{k}, 'AbsTol', 1e-10, sprintf('measure %d', k));
end
verifyEqual(tc, size(M), [6, n, n]);
end

function test_pieces_use_only_their_samples(tc)
% PLV, wPLI, PLI, AEC on pieces = the formulas on the samples of the pieces
Z = tc.TestData.Z;  seg = tc.TestData.seg;
in = false(1, size(Z, 2));
for k = 1:size(seg, 1), in(seg(k, 1):seg(k, 2)) = true; end
[M, power] = static_connectivity(Z, seg, {'plv', 'wpli', 'pli', 'aec'}, tc.TestData.fs, 8, 12);
R = direct(Z(:, in));
for k = 1:4
    verifyEqual(tc, squeeze(M(k, :, :)), R{k + (k == 4)}, 'AbsTol', 1e-10);   % PLV, wPLI, PLI, AEC
end
verifyEqual(tc, power, mean(abs(Z(:, in)) .^ 2, 2), 'RelTol', 1e-12);
Z2 = Z;  Z2(:, ~in) = 1e6 * randn(size(Z2(:, ~in)));          % samples outside the pieces are ignored
verifyEqual(tc, static_connectivity(Z2, seg, {'plv', 'wpli', 'aec', 'dwpli'}, tc.TestData.fs, 8, 12), ...
    static_connectivity(Z, seg, {'plv', 'wpli', 'aec', 'dwpli'}, tc.TestData.fs, 8, 12), 'AbsTol', 1e-12);
end

function test_real_input_is_hilbert_of_each_piece(tc)
X = tc.TestData.X;  seg = tc.TestData.seg;
Zp = complex(zeros(size(X)));
for k = 1:size(seg, 1)
    Zp(:, seg(k, 1):seg(k, 2)) = hilbert(X(:, seg(k, 1):seg(k, 2)).').';
end
m = {'plv', 'wpli', 'aec_orth'};
verifyEqual(tc, static_connectivity(X, seg, m, tc.TestData.fs, 8, 12), ...
    static_connectivity(Zp, seg, m, tc.TestData.fs, 8, 12), 'AbsTol', 1e-12);
end

function test_dwpli_matches_vinck_formula(tc)
% debiased wPLI (Vinck et al. 2011: the sum over pairs of different trials)
% written directly, without FieldTrip: epochs of 6 cycles of the band
% centre (0.6 s at 10 Hz, 77 samples) one after the other from the start
% of each piece, each demeaned, Hanning taper, the Fourier transform at
% the frequencies k fs / L inside the band (8.3, 10.0, 11.6 Hz); at each
% frequency X_e = Im(F_i conj(F_j)) of epoch e and dwPLI = ((sum X)^2 -
% sum X^2) / ((sum |X|)^2 - sum X^2); the mean over the frequencies.
% static_connectivity computes it with ft_freqanalysis and
% ft_connectivity_wpli (test_reference_connectivity compares it with
% ft_connectivityanalysis 'wpli_debiased')
Z = tc.TestData.Z;  fs = tc.TestData.fs;  seg = tc.TestData.seg;
L = round(6 / 10 * fs);
f = (ceil(8 * L / fs):floor(12 * L / fs)) * fs / L;
verifyEqual(tc, numel(f), 3);
first = [];
for k = 1:size(seg, 1)
    first = [first, seg(k, 1):L:seg(k, 2) - L + 1]; %#ok<AGROW>
end
M = static_connectivity(Z, seg, {'dwpli'}, fs, 8, 12);
verifyEqual(tc, squeeze(M), vinck(real(Z), first, L, f, fs), 'AbsTol', 1e-10);
% epochs in seconds (S.dwpli_epoch_length): 2 s, frequencies every 0.5 Hz
L = 2 * fs;
M = static_connectivity(Z, [], {'dwpli'}, fs, 8, 12, struct('dwpli_epoch_length', 2));
verifyEqual(tc, squeeze(M), vinck(real(Z), 1:L:size(Z, 2) - L + 1, L, 8:0.5:12, fs), 'AbsTol', 1e-10);
% 10 cycles: 1 s at 10 Hz
M = static_connectivity(Z, [], {'dwpli'}, fs, 8, 12, struct('dwpli_epoch_cycles', 10));
verifyEqual(tc, squeeze(M), vinck(real(Z), 1:fs:size(Z, 2) - fs + 1, fs, 8:12, fs), 'AbsTol', 1e-10);
end

function test_dwpli_needs_ten_epochs(tc)
% fewer than 10 epochs (a short recording or split half), or a band
% without a Fourier frequency of the epochs: dwPLI NaN with a note, not an
% error (which would stop the recording); the other measures are computed
Z = tc.TestData.Z;  fs = tc.TestData.fs;  n = size(Z, 1);
L = round(6 / 10 * fs);
[M, ~, notes] = static_connectivity(Z(:, 1:10 * L), [], {'plv', 'dwpli'}, fs, 8, 12);   % 10 epochs: computed
verifyTrue(tc, all(isfinite(M), 'all'));
verifyEmpty(tc, notes{2});
[M, ~, notes] = static_connectivity(Z(:, 1:10 * L - 1), [], {'plv', 'dwpli'}, fs, 8, 12);   % 9
verifyTrue(tc, all(isfinite(M(1, :, :)), 'all'));
D = squeeze(M(2, :, :));
verifyTrue(tc, all(isnan(D(~eye(n)))) && all(diag(D) == 0));
verifyEmpty(tc, notes{1});
verifySubstring(tc, notes{2}, 'at least 10 epochs');
% epochs of 0.25 s: Fourier frequencies every 4 Hz, none in 9-11 Hz
[M, ~, notes] = static_connectivity(Z, [], {'plv', 'dwpli'}, fs, 9, 11, struct('dwpli_epoch_length', 0.25));
D = squeeze(M(2, :, :));
verifyTrue(tc, all(isnan(D(~eye(n)))) && all(diag(D) == 0) && all(isfinite(M(1, :, :)), 'all'));
verifySubstring(tc, notes{2}, 'no Fourier frequency');
end

function test_symmetric_orthogonalisation(tc)
% the orthogonalised signals are uncorrelated; their analytic signals are
% those of the orthogonalised real signals
Z = tc.TestData.Z;
M = static_connectivity(Z, [], {'aec_symorth', 'plv_symorth'}, tc.TestData.fs, 8, 12);
X = real(Z);
L = ROInets.closest_orthogonal_matrix(X')';
G = L * L.';
verifyLessThan(tc, max(abs(triu(G ./ sqrt(diag(G) * diag(G).'), 1)), [], 'all'), 1e-8);
Zo = hilbert(L.').';
R = direct(Zo);
verifyEqual(tc, squeeze(M(1, :, :)), R{5}, 'AbsTol', 1e-8);     % the Hilbert transform is linear
verifyEqual(tc, squeeze(M(2, :, :)), R{1}, 'AbsTol', 1e-8);
% more regions than independent signals: that measure NaN with a note
% (not an error, which would lose every measure of the recording)
[M, ~, notes] = static_connectivity([Z; Z(1, :) + Z(2, :)], [], {'aec', 'aec_symorth', 'plv_symorth'}, ...
    tc.TestData.fs, 8, 12);
verifyTrue(tc, all(isfinite(M(1, :, :)), 'all'));
for k = 2:3
    A = squeeze(M(k, :, :));
    verifyTrue(tc, all(isnan(A(~eye(6)))) && all(diag(A) == 0));
    verifySubstring(tc, notes{k}, 'needs 6 independent signals; the region signals have rank 5');
end
end

function test_symmetric_orthogonalisation_rank_of_float32_data(tc)
% data of rank 5 stored as 32-bit floats: the rounding (1e-7 of the signal)
% is not an independent signal. MATLAB's rank() counts it (tolerance at
% machine precision) and the orthogonalisation would amplify it; the rank is
% counted as data_rank counts the channels' (connectivity.signal_rank:
% MNE-Python's rule for 32-bit data, singular values above n eps('single')
% times the largest)
X = real(tc.TestData.Z);
X6 = double(single([X; X(1, :) + X(2, :)]));       % rank 5 up to float32 rounding
verifyEqual(tc, rank(X6), 6);                      % what rank() says
Z6 = hilbert(X6.').';
[Zo, note, r] = connectivity.symmetric_orthogonalisation(Z6, '''aec_orth''');
verifyEqual(tc, r, 5);
verifyEmpty(tc, Zo);
verifySubstring(tc, note, 'Use an atlas with fewer regions than that, or ''aec_orth''');
[Zo, note, r] = connectivity.symmetric_orthogonalisation(tc.TestData.Z, '''aec_orth''');   % full rank: done
verifyEqual(tc, r, 5);  verifyEmpty(tc, note);  verifySize(tc, Zo, size(tc.TestData.Z));
% the windowed functions: NaN edges, zero diagonal, the windows kept, and
% get_all_connectivity_mats computes the other measures and says why
fs = tc.TestData.fs;  seg = [1, size(Z6, 2)];
[a, sa, fa, wa, na] = connectivity.aec_sliding_window_corrected(Z6, fs, 6, 0.5, seg);
[p, sp, fp, wp, np] = connectivity.plv_sliding_window_corrected(Z6, fs, 0.6, 0.3, seg);
[a0, ~, ~, wa0] = connectivity.aec_sliding_window(Z6, fs, 6, 0.5, seg);
verifySize(tc, a, size(a0));  verifyEqual(tc, wa, wa0);  verifyNotEmpty(tc, wp.first_sample);
off = repmat(reshape(~eye(6), [1, 6, 6]), size(a, 1), 1, 1);
verifyTrue(tc, all(isnan(a(off))) && all(a(~off) == 0));
verifyTrue(tc, all(isnan(p(repmat(reshape(~eye(6), [1, 6, 6]), size(p, 1), 1, 1)))));
verifyTrue(tc, all(isnan(sa(~eye(6)))) && all(isnan(fa(~eye(6)))) && all(isnan(sp(~eye(6)))) && all(isnan(fp(~eye(6)))));
verifyNotEmpty(tc, na);  verifyNotEmpty(tc, np);
opts = struct('conn', {{'wpli', 'aec', 'aec_symorth', 'plv_symorth'}}, 'window_cycles', 6, 'window_length', [], ...
    'window_overlap', 0.5, 'aec_window', 6, 'aec_step', 0.5);
[c, s, f, info] = get_all_connectivity_mats(Z6, fs, 8, 12, seg, opts);
verifyEqual(tc, sort(fieldnames(info.not_computed))', {'aec_symorth', 'plv_symorth'});
verifyTrue(tc, all(isfinite(c{1}(:))) && all(isfinite(c{2}(:))) && all(isfinite(f(1:2, :, :)), 'all'));
verifyTrue(tc, all(isnan(s(3, ~eye(6)))) && all(isnan(f(4, ~eye(6)))));
end

function test_signal_rank_rule(tc)
% connectivity.signal_rank, MNE-Python's estimate_rank with tol 'float32'
% (mne/rank.py, _estimate_rank_from_s): singular values of the centred
% signals above n x eps('single') x the largest. A simulated case: 129
% channels of rank 117 saved as 32-bit floats after cleaning, with one
% more eigenvalue 1.4e-12 of the largest (rounding) and the smallest real
% one 6e-4. A rule counting the eigenvalues above 1e-12 of the largest
% gives 118
rs = RandStream('mrg32k3a', 'Seed', 11);
N = 20000;
lam = [logspace(0, log10(6e-4), 117), 1.4e-12];   % covariance eigenvalues relative to the largest
[Q, ~] = qr(randn(rs, 129));
S = randn(rs, 118, N);  S = (S - mean(S, 2)) ./ std(S, 1, 2);
[U, ~, V] = svd(S', 'econ');  S = (U * V')' * sqrt(N);   % exactly white: the eigenvalues are lam
X = Q(:, 1:118) * diag(sqrt(lam)) * S;
e = sort(eig(cov(X', 1)), 'descend');
verifyEqual(tc, e(118) / e(1), 1.4e-12, 'RelTol', 1e-3);
verifyEqual(tc, nnz(e > 1e-12 * e(1)), 118);         % eigenvalues above 1e-12 of the largest
verifyEqual(tc, connectivity.signal_rank(X), 117);
verifyEqual(tc, connectivity.signal_rank(X(1:117, :)), 117);
% MNE's tolerance written out: s > n eps32 max(s), s the singular values
s = svd(X - mean(X, 2));
verifyEqual(tc, connectivity.signal_rank(X), nnz(s > 129 * double(eps('single')) * s(1)));
end

function test_long_recordings_in_chunks(tc)
% PLI, ciPLV and pairwise-orthogonalised AEC take their edges (signals) in
% chunks of about 200 MB, as wPLI (128 edges; all signals at once for
% aec_orth would need 7 GB for an hour of 68 regions at 1000 Hz): 68 regions
% of 190,000 samples, several chunks, equal the direct formulas
rs = RandStream('mrg32k3a', 'Seed', 4);
fs = 100;  n = 68;  N = 190000;
[b, a] = butter(4, [8 12] / (fs / 2));
X = filtfilt(b, a, randn(rs, N, n)).';
X(68, :) = X(68, :) + 0.5 * [zeros(1, 2), X(1, 1:end - 2)];
Z = hilbert(X.').';
clear X
verifyGreaterThan(tc, n * (n - 1) / 2, floor(2e8 / (8 * N)));        % more than one chunk of edges
verifyLessThan(tc, floor(2e8 / (16 * N)), n);                        % and of signals
M = static_connectivity(Z, [], {'pli', 'ciplv', 'aec_orth'}, fs, 8, 12);
seg = [1, N];
[pw, ~, pf] = connectivity.pli_sliding_window(Z, fs, 0.6, 0.3, seg);
[cw, ~, cf] = connectivity.ciplv_sliding_window(Z, fs, 0.6, 0.3, seg);
[ow, ~, of] = connectivity.aec_sliding_window_corrected_pairwise(Z, fs, 300, 150, seg);
verifyEqual(tc, double(pf), squeeze(M(1, :, :)), 'AbsTol', 1e-10);
verifyEqual(tc, double(cf), squeeze(M(2, :, :)), 'AbsTol', 1e-10);
verifyEqual(tc, double(of), squeeze(M(3, :, :)), 'AbsTol', 1e-10);
for ij = [1 68; 2 67; 40 41]'                      % pairs in the first and last chunks
    i = ij(1);  j = ij(2);
    s = Z(i, :) .* conj(Z(j, :));
    verifyEqual(tc, M(1, i, j), abs(mean(sign(imag(s)))), 'AbsTol', 1e-12);
    q = mean(exp(1i * angle(s)));
    verifyEqual(tc, M(2, i, j), abs(imag(q)) / sqrt(1 - real(q) ^ 2), 'AbsTol', 1e-10);
    bij = sum(real(Z(j, :)) .* real(Z(i, :))) / sum(real(Z(i, :)) .^ 2);
    bji = sum(real(Z(i, :)) .* real(Z(j, :))) / sum(real(Z(j, :)) .^ 2);
    o = (corr(abs(Z(i, :))', abs(Z(j, :) - bij * Z(i, :))') + corr(abs(Z(j, :))', abs(Z(i, :) - bji * Z(j, :))')) / 2;
    verifyEqual(tc, M(3, i, j), o, 'AbsTol', 1e-10);
    w = 1:60;                                      % the first window of 0.6 s
    sw = s(w);  qw = mean(exp(1i * angle(sw)));
    verifyEqual(tc, double(pw(1, i, j)), abs(mean(sign(imag(sw)))), 'AbsTol', 1e-6);
    verifyEqual(tc, double(cw(1, i, j)), abs(imag(qw)) / sqrt(1 - real(qw) ^ 2), 'AbsTol', 1e-5);
    w = 1:30000;                                   % the first window of 300 s
    o = (corr(abs(Z(i, w))', abs(Z(j, w) - bij * Z(i, w))') + corr(abs(Z(j, w))', abs(Z(i, w) - bji * Z(j, w))')) / 2;
    verifyEqual(tc, double(ow(1, i, j)), o, 'AbsTol', 1e-6);
end
end

function test_samples_without_amplitude(tc)
% analytic signals with samples of zero amplitude (no phase): ciPLV, PLV,
% wPLI and PLI stay finite in every window (the cumulative sums of ciPLV
% would otherwise carry a NaN into every later window)
Z = tc.TestData.Z(1:4, 1:4000);
Z(:, 1501:2500) = 0;
seg = [1, 4000];
c = connectivity.ciplv_sliding_window(Z, 128, 0.5, 0.25, seg);
[p, ~, w] = connectivity.plv_wpli_sliding_window(Z, 128, 0.5, 0.25, seg);
q = connectivity.pli_sliding_window(Z, 128, 0.5, 0.25, seg);
verifyTrue(tc, all(isfinite(c(:))) && all(isfinite(p(:))) && all(isfinite(w(:))) && all(isfinite(q(:))));
verifyEqual(tc, c(1, 1, 2), single(static_connectivity_value(Z(:, 1:64), 'ciplv', 1, 2)), 'AbsTol', 1e-6);
M = static_connectivity(Z, seg, {'ciplv', 'plv'}, 128, 8, 12);
verifyTrue(tc, all(isfinite(M(:))));
end

function test_debiased_wpli(tc)
% dwPLI (FieldTrip's debiased wPLI, epochs of 6 cycles as trials) in
% windows of epochs, 20 s every 10 s: near 0 without coupling, where the
% wPLI of the 6-cycle windows is well above 0 (the bias it removes); near
% wPLI^2 with strong coupling; windows and epochs inside each piece; windows
% of fewer than 10 epochs: none, with a note (not an error)
rng(5);
fs = 256;
n = fs * 120;
nb = @(m) bandpass(randn(n, m), [8 12], fs).';
x = nb(7);
x(3, :) = 0.3 * x(3, :) + circshift(x(1, :), 6);   % 1-3 coupled with a lag; the other pairs independent
Z = hilbert(x.').';
e = connectivity.dwpli_epoch(fs, 8, 12);
verifyEqual(tc, [e.samples, e.min_epochs], [round(0.6 * fs), 10]);
verifyEqual(tc, e.freqs, (5:7) * fs / e.samples, 'AbsTol', 1e-12);   % 8.3, 10.0 and 11.6 Hz
[d, sd, full_d, win] = connectivity.dwpli_sliding_window(Z, fs, 20, 10, [], e);
verifyEqual(tc, [win.samples, win.step, win.epochs, win.epoch_samples], [20 * fs, 10 * fs, 33, e.samples]);
verifyEqual(tc, size(d), [11, 7, 7]);
verifyEqual(tc, double(squeeze(d(5, :, :))), squeeze(connectivity.dwpli_epochs(Z, fs, ...
    win.first_sample(5) + (0:32) * e.samples, e)), 'AbsTol', 1e-6);
verifyEqual(tc, sd, squeeze(mean(double(d), 1)), 'AbsTol', 1e-6);
verifyEqual(tc, double(squeeze(d(5, 3, 1))), double(squeeze(d(5, 1, 3))));
verifyEqual(tc, double(diag(squeeze(d(5, :, :)))), zeros(7, 1));
[~, ~, w, ~, ~, full_w] = connectivity.plv_wpli_sliding_window(Z, fs, 0.6, 0.3);
null = setdiff(find(tril(true(7), -1)), sub2ind([7 7], 3, 1));
D = reshape(d, size(d, 1), []);  Wp = reshape(w, size(w, 1), []);
verifyGreaterThan(tc, mean(Wp(:, null), 'all'), 0.3);    % wPLI of 6-cycle windows: biased without coupling
verifyLessThan(tc, abs(mean(D(:, null), 'all')), 0.02);  % dwPLI: not
verifyLessThan(tc, max(abs(full_d(null))), 0.05);
verifyEqual(tc, full_d(1, 3), full_w(1, 3) ^ 2, 'AbsTol', 0.05);
verifyGreaterThan(tc, full_d(1, 3), 0.8);
% pieces: the windows and epochs of the first piece are those of the piece alone
seg = [1, fs * 50; fs * 50 + 1, n];
[dp, ~, fp, wp] = connectivity.dwpli_sliding_window(Z, fs, 20, 10, seg, e);
d1 = connectivity.dwpli_sliding_window(Z(:, 1:fs * 50), fs, 20, 10, [], e);
verifyEqual(tc, size(d1, 1), 4);
verifyEqual(tc, dp(1:4, :, :), d1, 'AbsTol', 1e-6);
verifyEqual(tc, wp.first_sample, [0:10:30, 50:10:100] * fs + 1);
L = e.samples;
verifyEqual(tc, fp, squeeze(connectivity.dwpli_epochs(Z, fs, [1:L:fs * 50 - L + 1, fs * 50 + 1:L:n - L + 1], e)), ...
    'AbsTol', 1e-12);
% windows of 5 s hold 8 epochs: no windows, NaN window mean, the note
[d0, s0, ~, w0, note] = connectivity.dwpli_sliding_window(Z, fs, 5, 2.5, [], e);
verifyEqual(tc, size(d0), [0, 7, 7]);
verifyTrue(tc, all(isnan(s0(~eye(7)))) && all(diag(s0) == 0) && isempty(w0.first_sample));
verifySubstring(tc, note, 'at least 10');
% through get_all_connectivity_mats: dwPLI's own windows (60 cycles of the
% band centre every 30: 6 s every 3 s), not the common ones (2 s here); each
% holds 10 epochs although an epoch is 153.6 samples, rounded to 154 (6 s
% would be 1536 samples, 9 epochs; the window is 1540 long)
[c, ~, ~, info] = get_all_connectivity_mats(x, fs, 8, 12, [], struct('conn', {{'wpli', 'dwpli'}}, 'window_length', 2));
d6 = connectivity.dwpli_sliding_window(Z, fs, 6, 3, [], e);
verifyEqual(tc, c{2}, d6, 'AbsTol', 1e-6);
verifyEqual(tc, [info.windows{1}.samples, info.windows{2}.samples, info.windows{2}.step, info.windows{2}.epochs], ...
    [2 * fs, 10 * e.samples, 3 * fs, 10]);
end

function v = static_connectivity_value(Z, m, i, j)
M = static_connectivity(Z, [], {m}, 128, 8, 12);
v = M(1, i, j);
end

%% Known couplings
function test_lagged_and_zero_lag_coupling(tc)
% region 2 follows region 1 with a lag: high wPLI; region 3 is a zero-lag
% copy of region 1 (as leakage): high PLV and AEC, low wPLI, PLI and
% ciPLV; regions 4 and 5 are independent
M = static_connectivity(tc.TestData.Z, [], tc.TestData.all, tc.TestData.fs, 8, 12);
v = @(m, i, j) M(strcmp(tc.TestData.all, m), i, j);
verifyGreaterThan(tc, v('wpli', 1, 2), 0.6);
verifyGreaterThan(tc, v('dwpli', 1, 2), 0.4);
verifyGreaterThan(tc, v('plv', 1, 3), 0.3);
verifyLessThan(tc, v('wpli', 1, 3), 0.2);
verifyLessThan(tc, v('pli', 1, 3), 0.2);
verifyLessThan(tc, v('ciplv', 1, 3), 0.2);
verifyLessThan(tc, v('aec_orth', 1, 3), v('aec', 1, 3));
verifyLessThan(tc, v('wpli', 4, 5), 0.2);
verifyLessThan(tc, abs(v('dwpli', 4, 5)), 0.05);
for k = 1:numel(tc.TestData.all)                       % symmetric, zero diagonal
    A = squeeze(M(k, :, :));
    verifyEqual(tc, A, A.', 'AbsTol', 1e-12);
    verifyEqual(tc, diag(A), zeros(5, 1));
end
end

function test_identical_signals(tc)
% no imaginary part at all: wPLI 0, not NaN
Z = tc.TestData.Z([1 1], :);
M = static_connectivity(Z, [], {'wpli', 'pli', 'plv'}, tc.TestData.fs, 8, 12);
verifyEqual(tc, M(:, 1, 2), [0; 0; 1], 'AbsTol', 1e-12);
verifyError(tc, @() static_connectivity(tc.TestData.Z, [], {'coherence'}, tc.TestData.fs, 8, 12), ...
    'static_connectivity:measure');
end

%% The windowed functions and get_all_connectivity_mats
function test_windowed_functions_give_the_same_whole_recording(tc)
% the whole-recording output of each windowed function of +connectivity
% equals static_connectivity (dwPLI: the epochs of the whole recording,
% not of its windows). Pieces back to back, as bandpass_pieces gives them (the
% windowed AEC correlates every sample of the data, also between pieces)
Z = tc.TestData.Z;  fs = tc.TestData.fs;
seg = [1, 900; 901, 2700; 2701, size(Z, 2)];
win = 6 / 10;  step = win / 2;
[~, ~, ~, ~, fp, fw] = connectivity.plv_wpli_sliding_window(Z, fs, win, step, seg);
[~, ~, fq] = connectivity.pli_sliding_window(Z, fs, win, step, seg);
[~, ~, fc] = connectivity.ciplv_sliding_window(Z, fs, win, step, seg);
[~, ~, fa] = connectivity.aec_sliding_window(Z, fs, 6, 0.5, seg);
[~, ~, fo] = connectivity.aec_sliding_window_corrected_pairwise(Z, fs, 6, 0.5, seg);
[~, ~, fd] = connectivity.dwpli_sliding_window(Z, fs, 20, 10, seg, connectivity.dwpli_epoch(fs, 8, 12));
M = static_connectivity(Z, seg, {'plv', 'wpli', 'pli', 'ciplv', 'aec', 'aec_orth', 'dwpli'}, fs, 8, 12);
F = {fp, fw, fq, fc, fa, fo, fd};
for k = 1:7
    verifyEqual(tc, squeeze(M(k, :, :)), double(F{k}), 'AbsTol', 1e-10, sprintf('measure %d', k));
end
end

function test_whole_recording_with_and_without_windows(tc)
% get_all_connectivity_mats: static_full_cmats and static_power are
% static_connectivity's, with windows (S.dynamic) or without; without,
% no windows at all
Z = tc.TestData.Z;  seg = tc.TestData.seg;  fs = tc.TestData.fs;
m = {'plv', 'wpli', 'ciplv', 'dwpli', 'aec', 'aec_orth'};
opts = struct('conn', {m}, 'window_cycles', 6, 'window_length', [], 'window_overlap', 0.5, ...
    'aec_window', 6, 'aec_step', 0.5);
[c1, s1, f1, i1, p1] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
opts.dynamic = false;
[c0, s0, f0, i0, p0] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
[M, P] = static_connectivity(Z, seg, m, fs, 8, 12);
verifyEqual(tc, f1, M);  verifyEqual(tc, f0, M);
verifyEqual(tc, p1, P);  verifyEqual(tc, p0, P);
verifyEqual(tc, numel(c1), numel(m));  verifyEqual(tc, size(s1), [numel(m), 5, 5]);
verifyTrue(tc, isfield(i1, 'plv_wpli_windows') && isfield(i1, 'window_s'));
verifyEqual(tc, c0, {});  verifyEqual(tc, s0, []);  verifyEqual(tc, i0.windows, {});
verifyFalse(tc, isfield(i0, 'plv_wpli_windows') || isfield(i0, 'window_s'));
end

function test_window_mean_without_states(tc)
% a static-only run that asks for the window mean (S.static_matrices with
% 'window_mean') computes the windows as S.dynamic does: the same windows,
% window mean and whole recording
Z = tc.TestData.Z;  seg = tc.TestData.seg;  fs = tc.TestData.fs;
m = {'plv', 'wpli', 'ciplv', 'aec_orth'};
opts = struct('conn', {m}, 'window_cycles', 6, 'window_length', [], 'window_overlap', 0.5, ...
    'aec_window', 6, 'aec_step', 0.5);
[c1, s1, f1, i1, p1] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
opts.dynamic = false;  opts.static_matrices = {'full', 'window_mean'};
[c0, s0, f0, i0, p0] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
verifyEqual(tc, c0, c1);  verifyEqual(tc, s0, s1);
verifyEqual(tc, f0, f1);  verifyEqual(tc, p0, p1);
verifyEqual(tc, i0.windows, i1.windows);
verifyEqual(tc, i0.no_windows, struct());
for k = 1:numel(m)                               % the window mean is the mean of the windows
    verifyEqual(tc, squeeze(s0(k, :, :)), squeeze(mean(double(c0{k}), 1)), 'AbsTol', 1e-6, m{k});
end
opts.static_matrices = {'full'};                 % the whole recording only: no windows
[c2, s2] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
verifyEqual(tc, c2, {});  verifyEqual(tc, s2, []);
end

function test_window_sizes(tc)
% phase windows: cycles of the band centre, or a fixed length for every
% band; step = (1 - overlap) x window; AEC: its own length and step
Z = tc.TestData.Z;  fs = tc.TestData.fs;
opts = struct('conn', {{'wpli', 'aec'}}, 'dynamic', false, 'static_matrices', {{'window_mean'}});
[~, ~, ~, i] = get_all_connectivity_mats(Z, fs, 8, 12, [], opts);       % defaults
verifyEqual(tc, [i.window_s, i.step_s, i.aec_window_s, i.aec_step_s], [0.6, 0.3, 6, 0.5], 'AbsTol', 1e-12);
verifyEqual(tc, i.windows{1}.samples, ceil(0.6 * fs));
opts.window_cycles = 4;  opts.window_overlap = 0.75;  opts.aec_window = 4;  opts.aec_step = 1;
[~, ~, ~, i] = get_all_connectivity_mats(Z, fs, 8, 12, [], opts);
verifyEqual(tc, [i.window_s, i.step_s, i.aec_window_s, i.aec_step_s], [0.4, 0.1, 4, 1], 'AbsTol', 1e-12);
opts.window_length = 2;  opts.window_overlap = 0.5;                       % fixed length: the cycles unused
[c, ~, ~, i] = get_all_connectivity_mats(Z, fs, 8, 12, [], opts);
verifyEqual(tc, [i.window_s, i.step_s], [2, 1], 'AbsTol', 1e-12);
verifyEqual(tc, i.windows{1}.samples, 2 * fs);
verifyEqual(tc, size(c{1}, 1), floor((size(Z, 2) - 2 * fs) / fs) + 1);
[~, ~, ~, i] = get_all_connectivity_mats(Z, fs, 20, 30, [], opts);       % another band: the same length
verifyEqual(tc, i.window_s, 2);
end

function test_windows_of_one_measure(tc)
% S.measure_windows: a measure with windows of its own is
% computed in them; the other measures keep the common windows (PLV
% unchanged when wPLI has its own); window_cycles alone means cycles
Z = tc.TestData.Z;  seg = tc.TestData.seg;  fs = tc.TestData.fs;
m = {'plv', 'wpli', 'ciplv', 'aec', 'aec_orth'};
opts = struct('conn', {m}, 'dynamic', false, 'static_matrices', {{'window_mean'}});
[c0, s0, f0, i0] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
opts.measure_windows = struct('wpli', struct('window_length', 2), 'aec_orth', struct('aec_window', 4, 'aec_step', 1));
[c1, s1, f1, i1] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
for k = [1 3 4]                                          % plv, ciplv, aec: the common windows
    verifyEqual(tc, c1{k}, c0{k}, m{k});
    verifyEqual(tc, s1(k, :, :), s0(k, :, :), m{k});
    verifyEqual(tc, i1.windows{k}, i0.windows{k}, m{k});
end
[~, ~, w2, sw2, ~, ~, win2] = connectivity.plv_wpli_sliding_window(Z, fs, 2, 1, seg);   % 2 s windows, 50% overlap
verifyEqual(tc, c1{2}, w2);
verifyEqual(tc, squeeze(s1(2, :, :)), sw2);
verifyEqual(tc, i1.windows{2}, win2);
[a4, sa4, ~, wa4] = connectivity.aec_sliding_window_corrected_pairwise(Z, fs, 4, 1, seg);   % 4 s every 1 s
verifyEqual(tc, c1{5}, a4);
verifyEqual(tc, squeeze(s1(5, :, :)), sa4);
verifyEqual(tc, i1.windows{5}, wa4);
verifyEqual(tc, [i1.windows{2}.samples, i1.windows{5}.samples], [2, 4] * fs);
verifyEqual(tc, f1, f0);                                 % the whole recording: no windows
verifyEqual(tc, [i1.window_s, i1.step_s, i1.aec_window_s, i1.aec_step_s], [0.6, 0.3, 6, 0.5], 'AbsTol', 1e-12);   % plv, aec
verifyEqual(tc, i1.measure_windows, opts.measure_windows);
verifyFalse(tc, isfield(i0, 'measure_windows'));         % none given: no measure_windows in the info
opts.window_length = 2;                                  % fixed windows for all, cycles for ciPLV
opts.measure_windows = struct('ciplv', struct('window_cycles', 4));
[~, ~, ~, i2] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
o3 = opts;  o3.window_length = [];  o3.window_cycles = 4;  o3.measure_windows = struct();
[~, ~, ~, i3] = get_all_connectivity_mats(Z, fs, 8, 12, seg, o3);
verifyEqual(tc, i2.windows{3}, i3.windows{3});
verifyEqual(tc, i2.windows{2}.samples, 2 * fs);
end

function test_dwpli_windows_of_epochs(tc)
% dwPLI's windows: 60 cycles of the band centre every 30 by default,
% each cut into consecutive
% epochs of 6 cycles from its start, the trials of its own debiased wPLI
% (Vinck's formula here), 10 although each epoch is rounded to whole
% samples (76.8 to 77 here: the window is 770 samples, not 768); the
% window mean is their mean;
% longer epochs make the windows hold 10 at least; S.measure_windows.dwpli
% sets them; windows of fewer than 10 epochs: none, and why. The other
% measures keep their windows
Z = tc.TestData.Z;  fs = tc.TestData.fs;  N = size(Z, 2);  n = size(Z, 1);
L = round(0.6 * fs);  f = (5:7) * fs / L;
opts = struct('conn', {{'wpli', 'dwpli'}}, 'dynamic', false, 'static_matrices', {{'window_mean'}});
[c, s, ~, i] = get_all_connectivity_mats(Z, fs, 8, 12, [], opts);
w = i.windows{2};
verifyEqual(tc, [w.samples, w.step, w.epochs, w.epoch_samples], [10 * L, 3 * fs, 10, L]);
verifyEqual(tc, w.first_sample, 1:3 * fs:N - 10 * L + 1);
verifyEqual(tc, size(c{2}, 1), numel(w.first_sample));
verifyEqual(tc, i.windows{1}.samples, ceil(0.6 * fs));             % wPLI: its 6-cycle windows
for k = 1:numel(w.first_sample)
    R = vinck(real(Z), w.first_sample(k) + (0:w.epochs - 1) * L, L, f, fs);
    verifyEqual(tc, squeeze(double(c{2}(k, :, :))), R, 'AbsTol', 1e-6);   % windows in single precision
end
verifyEqual(tc, squeeze(s(2, :, :)), squeeze(mean(double(c{2}), 1)), 'AbsTol', 1e-6);
o = opts;  o.dwpli_epoch_length = 3;                               % 10 epochs of 3 s: windows of 30 s
[~, ~, ~, i3] = get_all_connectivity_mats(Z, fs, 8, 12, [], o);
verifyEqual(tc, [i3.windows{2}.samples, i3.windows{2}.step, i3.windows{2}.epochs], [30 * fs, 15 * fs, 10]);
o = opts;  o.measure_windows = struct('dwpli', struct('window_length', 10, 'window_overlap', 0.75));
[c4, ~, ~, i4] = get_all_connectivity_mats(Z, fs, 8, 12, [], o);
verifyEqual(tc, [i4.windows{2}.samples, i4.windows{2}.step, i4.windows{2}.epochs], [10 * fs, 2.5 * fs, 16]);
verifyEqual(tc, c4{1}, c{1});
o.measure_windows.dwpli.window_length = 5;                        % 8 epochs: no windows
[c5, s5, f5, i5] = get_all_connectivity_mats(Z, fs, 8, 12, [], o);
verifyEqual(tc, size(c5{2}, 1), 0);
A = squeeze(s5(2, :, :));
verifyTrue(tc, all(isnan(A(~eye(n)))) && all(diag(A) == 0));
verifyEqual(tc, fieldnames(i5.no_windows), {'dwpli'});
verifySubstring(tc, i5.no_windows.dwpli, 'at least 10');
verifyTrue(tc, all(isfinite(f5), 'all'));
end

function test_no_window_fits(tc)
% pieces shorter than the windows: no window mean for that measure (NaN,
% zero diagonal) and why; the other measures keep theirs
Z = tc.TestData.Z;  fs = tc.TestData.fs;
seg = [1, 4 * fs; 4 * fs + 101, 9 * fs];                % pieces of 4 s and 4.2 s; AEC windows of 6 s
opts = struct('conn', {{'wpli', 'aec'}}, 'dynamic', false, 'static_matrices', {{'window_mean'}});
[c, s, f, i] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
A = squeeze(s(2, :, :));
verifyTrue(tc, all(isnan(A(~eye(5)))) && all(diag(A) == 0));
verifyTrue(tc, all(isfinite(s(1, :, :)), 'all'));
verifyTrue(tc, all(isfinite(f), 'all'));                  % the whole recording is there
verifyEqual(tc, size(c{2}, 1), 0);
verifyEqual(tc, fieldnames(i.no_windows), {'aec'});
verifySubstring(tc, i.no_windows.aec, 'no window of 6 s');
end

%% Split halves
function test_split_halves(tc)
% blocks of 10 s counted from the start of the recording; the halves
% cover the pieces exactly, without overlap
fs = 100;
seg = [1, 2500; 3001, 7000];                    % pieces starting at 2 s and 37 s of the recording
onset_s = [2, 37];
h = split_halves(seg, onset_s, fs, 10);
mark = zeros(1, 7000);
for k = 1:2
    for r = 1:size(h{k}, 1)
        verifyEqual(tc, mark(h{k}(r, 1):h{k}(r, 2)), zeros(1, diff(h{k}(r, :)) + 1));
        mark(h{k}(r, 1):h{k}(r, 2)) = k;
    end
end
in = false(1, 7000);  in(1:2500) = true;  in(3001:7000) = true;
verifyTrue(tc, all(mark(in) > 0) && all(mark(~in) == 0));
t = nan(1, 7000);                               % time in the recording of each sample
t(1:2500) = 2 + (0:2499) / fs;  t(3001:7000) = 37 + (0:3999) / fs;
verifyEqual(tc, mark(in), mod(floor(t(in) / 10), 2) + 1);
verifyEqual(tc, h{1}(1, :), [1, 800]);          % 2-10 s: the first block
verifyEqual(tc, h{2}(1, :), [801, 1800]);       % 10-20 s
end

function test_window_mean_halves_blocks_per_measure(tc)
% the window mean on each half (window_mean_halves): each measure in
% blocks of max(10 s, 3 x its window), so 0.6 s windows keep 10 s blocks,
% 6 s windows (AEC) get 18 s and 20 s windows 60 s; the windows
% used are those wholly inside one block of the half (found here from
% their times in the recording), at least about a third of the windows per
% half, none in both halves
fs = 100;
seg = [1, 30000; 30001, 60000];                 % 0-300 s and 310-610 s of the recording
onset = [0, 310];
len = [0.6, 6, 20];  step = [0.3, 0.5, 10];
cmats = cell(1, 3);  windows = cell(1, 3);
for k = 1:3
    w = round(len(k) * fs);  first = [];
    for p = 1:2
        first = [first, seg(p, 1):round(step(k) * fs):seg(p, 2) - w + 1]; %#ok<AGROW>
    end
    windows{k} = struct('samples', w, 'first_sample', first);
    i = (1:numel(first))';
    cmats{k} = zeros(numel(first), 2, 2);
    cmats{k}(:, 1, 2) = i;  cmats{k}(:, 2, 1) = i .^ 2;
end
[H, block] = window_mean_halves(cmats, windows, seg, onset, fs, 10, [3, 2, 2]);
verifyEqual(tc, block, [10, 18, 60], 'AbsTol', 1e-12);
for k = 1:3
    first = windows{k}.first_sample(:);
    p = 1 + (first > seg(1, 2));
    t0 = onset(p)' + (first - seg(p, 1)) / fs;  % start and last sample of each window (s)
    t1 = t0 + (windows{k}.samples - 1) / fs;
    b0 = floor(t0 / block(k) + 1e-9);
    one_block = b0 == floor(t1 / block(k) + 1e-9);
    i = (1:numel(first))';
    for h = 1:2
        in = one_block & mod(b0, 2) == h - 1;
        verifyEqual(tc, H(h, k, 1, 2), mean(i(in)), 'RelTol', 1e-12);
        verifyEqual(tc, H(h, k, 2, 1), mean(i(in) .^ 2), 'RelTol', 1e-12);
        verifyGreaterThan(tc, nnz(in) / numel(i), 0.3, sprintf('measure %d, half %d', k, h));
    end
end
end

%% Direct formulas
function R = direct(Z)
% PLV, wPLI, PLI, ciPLV, AEC, AEC (pairwise orthogonalised), pair by pair
n = size(Z, 1);
R = repmat({zeros(n)}, 1, 6);
env = abs(Z);
for i = 1:n
    for j = 1:n
        if i == j, continue; end
        s = Z(i, :) .* conj(Z(j, :));
        p = mean(exp(1i * angle(s)));
        R{1}(i, j) = abs(p);
        R{2}(i, j) = abs(sum(imag(s))) / sum(abs(imag(s)));
        R{3}(i, j) = abs(mean(sign(imag(s))));
        R{4}(i, j) = abs(imag(p)) / sqrt(1 - real(p) ^ 2);
        R{5}(i, j) = corr(env(i, :)', env(j, :)');
        x = real(Z(i, :));
        b = sum(real(Z(j, :)) .* x) / sum(x .^ 2);
        R{6}(i, j) = corr(env(i, :)', abs(Z(j, :) - b * Z(i, :))');
    end
end
R{6} = (R{6} + R{6}') / 2;
end

function D = vinck(X, first, L, f, fs)
% debiased wPLI (Vinck et al. 2011) of the epochs X(:, first(e) + (0:L-1)),
% each demeaned and Hanning-tapered, at the frequencies f, then their mean
n = size(X, 1);
D = zeros(n);
for fk = f
    w = exp(-2i * pi * fk * (0:L - 1) / fs) .* hanning(L)';
    F = zeros(numel(first), n);
    for e = 1:numel(first)
        x = X(:, first(e) + (0:L - 1));
        F(e, :) = ((x - mean(x, 2)) * w.').';
    end
    for i = 1:n
        for j = [1:i - 1, i + 1:n]
            a = imag(F(:, i) .* conj(F(:, j)));
            D(i, j) = D(i, j) + (sum(a) ^ 2 - sum(a .^ 2)) / (sum(abs(a)) ^ 2 - sum(a .^ 2)) / numel(f);
        end
    end
end
end
