function [pli, static_pli, full_pli, win] = pli_sliding_window(data, srate, window, step, seg)
%PLI_SLIDING_WINDOW Phase lag index in sliding windows (Stam et al. 2007).
%   [pli, static_pli, full_pli, win] = connectivity.pli_sliding_window(data, srate, window, step, seg)
%   PLI = |mean(sign(Im(Sxy)))|, insensitive to zero-lag coupling. Same
%   windows and outputs as plv_wpli_sliding_window: every window (single,
%   zero diagonal), their mean and the value over the whole signal, from
%   one analytic signal per continuous piece.
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
pli_e = zeros(ne, numel(starts), 'single');
full_e = zeros(ne, 1);
chunk = max(1, min(256, floor(2e8 / (8 * nb_samples))));    % edges per pass: about 200 MB per array, as
for c0 = 1:chunk:ne                                            % plv_wpli_sliding_window (128 edges: 7 GB
    e = c0:min(c0 + chunk - 1, ne);                            % complex for an hour at 1000 Hz)
    [W, T] = connectivity.windowed_sums(sign(imag(Z(ii(e), :) .* conj(Z(jj(e), :)))), starts, win_samples, in_pieces);
    pli_e(e, :) = abs(W) / win_samples;
    full_e(e) = abs(T) / nnz(in_pieces);
end
pli = connectivity.to_matrices(pli_e, ii, jj, nb_signals);
static_pli = squeeze(mean(double(pli), 1));
full_pli = squeeze(connectivity.to_matrices(full_e, ii, jj, nb_signals));
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2);
end
