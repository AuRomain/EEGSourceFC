function [aec, static_aec, full_aec, win] = aec_sliding_window(data, srate, window, step, seg)
%AEC_SLIDING_WINDOW Amplitude envelope correlation in sliding windows.
%   [aec, static_aec, full_aec, win] = connectivity.aec_sliding_window(data, srate, window, step, seg)
%   data: band-passed signals (real; envelope by Hilbert transform of each
%   continuous piece of seg) or their analytic signal (complex).
%   aec: nWindows x nSignals x nSignals (single), Pearson correlation of the
%   envelopes in each window; static_aec: mean of the windows; full_aec:
%   over all the samples. Windows of 6 s every 0.5 s by default
%   (get_all_connectivity_mats). The sign is kept (with abs(),
%   anti-correlated envelopes would look coupled), the diagonal is 0 and
%   windows never cross a cut. No leakage correction: see aec_sliding_window_corrected
%   (symmetric) and aec_sliding_window_corrected_pairwise.
[nb_signals, nb_samples] = size(data);
if nargin < 5 || isempty(seg), seg = [1, nb_samples]; end
env = abs(connectivity.analytic_signal(data, seg));
win_samples = ceil(window * srate);
nb_shifts = ceil(step * srate);
starts = connectivity.window_starts(seg, win_samples, nb_shifts);
nb_windows = numel(starts);
aec = zeros(nb_windows, nb_signals, nb_signals, 'single');
off = ~eye(nb_signals);
for i = 1:nb_windows
    R = corr(env(:, starts(i):starts(i) + win_samples - 1).');
    aec(i, :, :) = R .* off;
end
static_aec = squeeze(mean(double(aec), 1));
full_aec = corr(env.') .* off;
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2);
end
