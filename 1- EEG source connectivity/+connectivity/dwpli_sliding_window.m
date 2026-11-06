function [dwpli, static_dwpli, full_dwpli, win] = dwpli_sliding_window(data, srate, window, step, seg, tau)
%DWPLI_SLIDING_WINDOW Debiased squared wPLI in sliding windows (Vinck et al. 2011).
%   [dwpli, static_dwpli, full_dwpli, win] =
%   connectivity.dwpli_sliding_window(data, srate, window, step, seg, tau)
%   Same windows and outputs as plv_wpli_sliding_window; tau (s) is the
%   band's decorrelation time, 1 / bandwidth (get_all_connectivity_mats).
%   Each window is cut into blocks of tau / 2: the imaginary part of the
%   cross-spectrum summed over a block, X_b, is one observation, as a trial
%   is in FieldTrip's ft_connectivity_wpli with debiasing, and
%     dwPLI = sum_{b,c far} X_b X_c / sum_{b,c far} |X_b| |X_c|,
%   the sums over the pairs of blocks whose closest samples are at least tau
%   apart (block indices 3 or more apart). Vinck's debiased estimator leaves
%   out only b = c, which is enough for independent trials; neighbouring
%   blocks of one signal are correlated, and their products keep most of
%   the bias (2 s windows of independent 7-13 Hz noise: wPLI^2 0.12, with
%   b = c left out 0.07, with the pairs closer than tau also left out
%   0.005). dwPLI estimates wPLI^2: about 0 without coupling (it can be
%   negative), up to 1. Meant for long windows (S.window_length = 2 or more:
%   24 blocks in 7-13 Hz); a window of fewer than 4 blocks is an error.
%   full_dwpli: over every block of the recording (blocks never cross a cut).
[nb_signals, nb_samples] = size(data);
if nargin < 5 || isempty(seg), seg = [1, nb_samples]; end
Z = connectivity.analytic_signal(data, seg);
win_samples = ceil(window * srate);
nb_shifts = ceil(step * srate);
starts = connectivity.window_starts(seg, win_samples, nb_shifts);
nb = round(win_samples / (tau / 2 * srate));               % blocks per window
if nb < 4
    error('connectivity:dwpli', ['dwPLI needs at least 4 blocks of %.3f s per window (window %.3f s): ' ...
        'use a longer window (S.window_length)'], tau / 2, window);
end
block_samples = floor(win_samples / nb);
first = starts(:) + (0:nb - 1) * block_samples;            % nWindows x nb first samples
all_blocks = connectivity.window_starts(seg, block_samples, block_samples);   % the recording's blocks
[ii, jj] = find(tril(true(nb_signals), -1));
ne = numel(ii);
nw = numel(starts);
d_e = zeros(ne, nw, 'single');
full_e = zeros(ne, 1);
chunk = max(1, min(256, floor(2e8 / (8 * nb_samples))));   % about 200 MB per array
for c0 = 1:chunk:ne
    e = c0:min(c0 + chunk - 1, ne);
    cs = cumsum([zeros(numel(e), 1), imag(Z(ii(e), :) .* conj(Z(jj(e), :)))], 2);
    X = cs(:, first + block_samples) - cs(:, first);       % block sums, edges x (nWindows * nb)
    d_e(e, :) = single(debiased(reshape(X, numel(e), nw, nb)));
    full_e(e) = debiased(reshape(cs(:, all_blocks + block_samples) - cs(:, all_blocks), numel(e), 1, []));
end
dwpli = connectivity.to_matrices(d_e, ii, jj, nb_signals);
static_dwpli = squeeze(mean(double(dwpli), 1));
full_dwpli = squeeze(connectivity.to_matrices(full_e, ii, jj, nb_signals));
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2, 'blocks', nb, 'block_samples', block_samples);
end

function d = debiased(X)
% X: edges x windows x blocks. Products of every pair of blocks (sum^2 minus
% squares) minus those of blocks 1 or 2 apart, over the same for |X|; 0
% when no pair far enough apart has an imaginary part to weigh
A = abs(X);
num = sum(X, 3) .^ 2 - sum(X .^ 2, 3);
den = sum(A, 3) .^ 2 - sum(A .^ 2, 3);
for lag = 1:2
    num = num - 2 * sum(X(:, :, 1:end - lag) .* X(:, :, 1 + lag:end), 3);
    den = den - 2 * sum(A(:, :, 1:end - lag) .* A(:, :, 1 + lag:end), 3);
end
d = num ./ den;
d(den <= 0) = 0;
end
