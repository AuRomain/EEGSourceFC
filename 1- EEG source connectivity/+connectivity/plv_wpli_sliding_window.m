function [plv, static_plv, wpli, static_wpli, full_plv, full_wpli, win] = plv_wpli_sliding_window(data, srate, window, step, seg)
%PLV_WPLI_SLIDING_WINDOW PLV and wPLI in sliding windows, and over the whole signal.
%   [plv, static_plv, wpli, static_wpli, full_plv, full_wpli, win] =
%   connectivity.plv_wpli_sliding_window(data, srate, window, step, seg)
%   data: nSignals x nSamples, band-passed (real; its analytic signal is
%   taken once per continuous piece) or already analytic (complex).
%   seg: nPieces x 2 continuous pieces (default: one piece); windows of
%   ceil(window*srate) samples every ceil(step*srate) samples, never
%   across a cut.
%   plv, wpli: nWindows x nSignals x nSignals (single), zero diagonal
%   static_plv, static_wpli: mean of the windows ("static")
%   full_plv, full_wpli: over all the samples (whole recording)
%   win: window length, step, first samples and centre samples
%   PLV = |mean(exp(1i*(phi_x - phi_y)))|; wPLI = |sum Im(Sxy)| / sum |Im(Sxy)|
%   (Vinck et al. 2011), as connectivity.compute_plv and calc_wpli, but:
%   the analytic signal of the whole recording is used in every window
%   instead of a Hilbert transform of each window; the sums of
%   cross-products (PLV, wPLI numerator) of each window are one matrix
%   product of its signals (multithreaded), the sums of |Im(Sxy)| (wPLI
%   denominator) cumulative sums (one pass over the data). Checked against
%   compute_plv and calc_wpli in tests/test_stage1.m.
[nb_signals, nb_samples] = size(data);
if nargin < 5 || isempty(seg), seg = [1, nb_samples]; end
Z = connectivity.analytic_signal(data, seg);
win_samples = ceil(window * srate);
nb_shifts = ceil(step * srate);
starts = connectivity.window_starts(seg, win_samples, nb_shifts);
in_pieces = false(1, nb_samples);
for k = 1:size(seg, 1), in_pieces(seg(k, 1):seg(k, 2)) = true; end
[ii, jj] = find(tril(true(nb_signals), -1));
ne = numel(ii);
U = Z ./ abs(Z);
U(abs(Z) == 0) = 0;                             % no phase: adds nothing (else NaN for every later sum)
nw = numel(starts);
lin = sub2ind([nb_signals, nb_signals], ii, jj);
plv_e = zeros(ne, nw, 'single');
num_e = zeros(ne, nw);
for w = 1:nw                                    % sums of cross-products: one product per window
    t = starts(w):starts(w) + win_samples - 1;
    P = U(:, t) * U(:, t)';                     % sum of u_i conj(u_j)
    C = Z(:, t) * Z(:, t)';                     % sum of the cross-spectrum S_ij
    plv_e(:, w) = abs(P(lin)) / win_samples;
    num_e(:, w) = imag(C(lin));
end
P = U(:, in_pieces) * U(:, in_pieces)';
C = Z(:, in_pieces) * Z(:, in_pieces)';
full_plv_e = abs(P(lin)) / nnz(in_pieces);
full_num = imag(C(lin));
% sums of |Im(S_ij)|: not a product, so by cumulative sums, edges in chunks
den_e = zeros(ne, nw);
full_den = zeros(ne, 1);
X = real(Z);  Y = imag(Z);
chunk = max(1, min(256, floor(2e8 / (8 * nb_samples))));    % about 200 MB per array
for c0 = 1:chunk:ne
    e = c0:min(c0 + chunk - 1, ne);
    a = abs(Y(ii(e), :) .* X(jj(e), :) - X(ii(e), :) .* Y(jj(e), :));     % |Im(z_i conj(z_j))|
    [den_e(e, :), full_den(e)] = connectivity.windowed_sums(a, starts, win_samples, in_pieces);
end
wpli_e = single(wpli_ratio(num_e, den_e));
full_wpli_e = wpli_ratio(full_num, full_den);
plv = connectivity.to_matrices(plv_e, ii, jj, nb_signals);
wpli = connectivity.to_matrices(wpli_e, ii, jj, nb_signals);
static_plv = squeeze(mean(double(plv), 1));
static_wpli = squeeze(mean(double(wpli), 1));
full_plv = squeeze(connectivity.to_matrices(full_plv_e, ii, jj, nb_signals));
full_wpli = squeeze(connectivity.to_matrices(full_wpli_e, ii, jj, nb_signals));
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2);
end

function w = wpli_ratio(num, den)
% |sum Im| / sum |Im|, in [0, 1]: 0 when there is no imaginary part at all
% (den = 0, e.g. a window of zero signal: Inf or NaN would stop the
% states), and at most 1 (num and den are summed in different orders)
w = min(abs(num) ./ den, 1);
w(den == 0) = 0;
end
