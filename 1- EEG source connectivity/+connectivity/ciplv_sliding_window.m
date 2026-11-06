function [ciplv, static_ciplv, full_ciplv, win] = ciplv_sliding_window(data, srate, window, step, seg)
%CIPLV_SLIDING_WINDOW Corrected imaginary PLV in sliding windows (Bruna et al. 2018).
%   [ciplv, static_ciplv, full_ciplv, win] = connectivity.ciplv_sliding_window(data, srate, window, step, seg)
%   ciPLV = |Im(mean(e^{i dphi}))| / sqrt(1 - Re(mean(e^{i dphi}))^2), the
%   definition of connectivity.compute_ciPLV, in the windows of
%   plv_wpli_sliding_window (same outputs). Insensitive to zero-lag coupling.
%   A sample without amplitude (|z| = 0) has no phase and adds nothing, as in
%   plv_wpli_sliding_window and static_connectivity (it was NaN, and the
%   cumulative sums made every later window of its edges NaN).
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
U(abs(Z) == 0) = 0;                             % no phase: adds nothing (NaN made every later window NaN)
c_e = zeros(ne, numel(starts), 'single');
full_e = zeros(ne, 1);
ci = @(P) abs(imag(P)) ./ sqrt(max(1 - real(P) .^ 2, eps));
chunk = max(1, min(256, floor(2e8 / (8 * nb_samples))));    % edges per pass, as plv_wpli_sliding_window
for c0 = 1:chunk:ne
    e = c0:min(c0 + chunk - 1, ne);
    [W, T] = connectivity.windowed_sums(U(ii(e), :) .* conj(U(jj(e), :)), starts, win_samples, in_pieces);
    c_e(e, :) = ci(W / win_samples);
    full_e(e) = ci(T / nnz(in_pieces));
end
ciplv = connectivity.to_matrices(c_e, ii, jj, nb_signals);
static_ciplv = squeeze(mean(double(ciplv), 1));
full_ciplv = squeeze(connectivity.to_matrices(full_e, ii, jj, nb_signals));
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2);
end
