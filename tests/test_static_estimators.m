function tests = test_static_estimators
%TEST_STATIC_ESTIMATORS The whole-recording measures against their formulas, and the split halves.
%   From the pipeline folder:  results = runtests('tests');  table(results)
%   Synthetic analytic signals only (no data): every measure of
%   static_connectivity.m computed again with a direct, slow formula, on
%   one continuous piece and on several pieces; known couplings; the
%   halves of split_halves.m; the whole-recording values of the windowed
%   functions of +connectivity and of get_all_connectivity_mats with and
%   without windows (S.dynamic).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, '1- EEG source connectivity'));      % static_connectivity, +connectivity
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

function test_dwpli_matches_formula(tc)
% debiased squared wPLI over blocks of B samples, pairs of blocks less than
% three blocks apart left out
Z = tc.TestData.Z;  fs = tc.TestData.fs;  n = size(Z, 1);
B = round(fs / (12 - 8) / 2);
M = static_connectivity(Z, [], {'dwpli'}, fs, 8, 12);
first = 1:B:(size(Z, 2) - B + 1);
R = zeros(n);
for i = 1:n
    for j = 1:n
        if i == j, continue; end
        s = imag(Z(i, :) .* conj(Z(j, :)));
        x = arrayfun(@(f) sum(s(f:f + B - 1)), first);
        num = 0;  den = 0;
        for a = 1:numel(x)
            for c = 1:numel(x)
                if abs(a - c) > 2
                    num = num + x(a) * x(c);
                    den = den + abs(x(a)) * abs(x(c));
                end
            end
        end
        R(i, j) = num / den;
    end
end
verifyEqual(tc, squeeze(M), R, 'AbsTol', 1e-10);
% fewer than 4 blocks (a split half too short): NaN and a note, not an
% error (which stopped the recording); the other measures are computed
[M3, ~, notes] = static_connectivity(Z(:, 1:3 * B), [], {'plv', 'dwpli'}, fs, 8, 12);
verifyTrue(tc, all(isfinite(M3(1, :, :)), 'all'));
D = squeeze(M3(2, :, :));
verifyTrue(tc, all(isnan(D(~eye(n)))) && all(diag(D) == 0));
verifyEmpty(tc, notes{1});
verifySubstring(tc, notes{2}, 'at least 4 blocks');
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
% (it was an error, which lost every measure of the recording)
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
% is not an independent signal. MATLAB's rank() counted it (tolerance at
% machine precision) and the orthogonalisation amplified it; the rank is
% now counted as data_rank counts the channels' (singular values above
% 1e-6 of the largest)
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

function test_long_recordings_in_chunks(tc)
% PLI, ciPLV and pairwise-orthogonalised AEC take their edges (signals) in
% chunks of about 200 MB, as wPLI (128 edges, and all signals at once for
% aec_orth, needed 7 GB for an hour of 68 regions at 1000 Hz): 68 regions
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
% wPLI and PLI stay finite in every window (ciPLV was NaN in every window
% from the first such sample on: the cumulative sums carried the NaN)
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
% (what static_full_cmats was before 3 November 2026) equals
% static_connectivity, except dwPLI, whose blocks they fit into windows.
% Pieces back to back, as bandpass_pieces gives them (the windowed AEC
% correlates every sample of the data, also between pieces)
Z = tc.TestData.Z;  fs = tc.TestData.fs;
seg = [1, 900; 901, 2700; 2701, size(Z, 2)];
win = 6 / 10;  step = win / 2;
[~, ~, ~, ~, fp, fw] = connectivity.plv_wpli_sliding_window(Z, fs, win, step, seg);
[~, ~, fq] = connectivity.pli_sliding_window(Z, fs, win, step, seg);
[~, ~, fc] = connectivity.ciplv_sliding_window(Z, fs, win, step, seg);
[~, ~, fa] = connectivity.aec_sliding_window(Z, fs, 6, 0.5, seg);
[~, ~, fo] = connectivity.aec_sliding_window_corrected_pairwise(Z, fs, 6, 0.5, seg);
M = static_connectivity(Z, seg, {'plv', 'wpli', 'pli', 'ciplv', 'aec', 'aec_orth'}, fs, 8, 12);
F = {fp, fw, fq, fc, fa, fo};
for k = 1:6
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
% S.measure_windows (6 Nov 2026): a measure with windows of its own is
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
verifyFalse(tc, isfield(i0, 'measure_windows'));         % none given: the info as before
opts.window_length = 2;                                  % fixed windows for all, cycles for ciPLV
opts.measure_windows = struct('ciplv', struct('window_cycles', 4));
[~, ~, ~, i2] = get_all_connectivity_mats(Z, fs, 8, 12, seg, opts);
o3 = opts;  o3.window_length = [];  o3.window_cycles = 4;  o3.measure_windows = struct();
[~, ~, ~, i3] = get_all_connectivity_mats(Z, fs, 8, 12, seg, o3);
verifyEqual(tc, i2.windows{3}, i3.windows{3});
verifyEqual(tc, i2.windows{2}.samples, 2 * fs);
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
