function [aec, static_aec, full_aec, win] = aec_sliding_window_corrected_pairwise(data, srate, window, step, seg)
%AEC_SLIDING_WINDOW_CORRECTED_PAIRWISE AEC after pairwise leakage correction.
%   [aec, static_aec, full_aec, win] =
%   connectivity.aec_sliding_window_corrected_pairwise(data, srate, window, step, seg)
%   data: band-passed signals (real), or their analytic signal (complex,
%   computed before the filter edges were trimmed), nSignals x nSamples;
%   seg: continuous pieces. For every pair, the zero-lag linear leakage of x into y is
%   regressed out (y - x * (x \ y), Brookes et al. 2012), the envelopes of x
%   and of the corrected y are correlated, and the two directions (x -> y,
%   y -> x) are averaged, as the original. Corrections:
%   - the original recomputed the whole-recording correction inside the
%     window loop and never used the windows, so every "window" held the
%     same value; here the correction is estimated once on the whole
%     recording and the envelopes are correlated in each window (6 s every
%     0.5 s, as aec_sliding_window), plus over the whole recording;
%   - the sign is kept (the original took abs(), E11);
%   - vectorised: one regression per signal and no further Hilbert
%     transform: the analytic signal of y - b x is Z_y - b Z_x (the
%     transform is linear), so the corrected envelopes keep the analytic
%     signal computed before the edges were trimmed (the original ran pinv
%     and hilbert for every pair and every window);
%   - the corrected envelopes of the signals y are made in chunks of about
%     200 MB (all of them at once up to about 180,000 samples, as before):
%     68 complex signals of an hour at 1000 Hz are 4 GB.
[n, nb_samples] = size(data);
if nargin < 5 || isempty(seg), seg = [1, nb_samples]; end
win_samples = ceil(window * srate);
nb_shifts = ceil(step * srate);
starts = connectivity.window_starts(seg, win_samples, nb_shifts);
nw = numel(starts);
Z = connectivity.analytic_signal(data, seg);          % unchanged when already complex
data = real(Z);
env = abs(Z);
A = zeros(nw, n, n);                     % A(:, i, j): x = i, corrected y = j
F = zeros(n);
chunk = max(1, min(n, floor(2e8 / (16 * nb_samples))));    % signals y per pass: about 200 MB per complex
for i = 1:n                                                 % array (all at once up to ~180,000 samples)
    x = data(i, :);
    beta = (data * x') / (x * x');                      % regression of every y on x
    envx = env(i, :);
    for c0 = 1:chunk:n
        c = c0:min(c0 + chunk - 1, n);
        envY = abs(Z(c, :) - beta(c) * Z(i, :));        % envelope of y - beta x
        for w = 1:nw
            idx = starts(w):starts(w) + win_samples - 1;
            A(w, i, c) = corr(envx(idx)', envY(:, idx)');
        end
        F(i, c) = corr(envx', envY');
    end
end
aec = single((A + permute(A, [1 3 2])) / 2);
for w = 1:nw, aec(w, 1:n + 1:end) = 0; end
full_aec = (F + F') / 2;
full_aec(1:n + 1:end) = 0;
static_aec = squeeze(mean(double(aec), 1));
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts, ...
    'centre_sample', starts + (win_samples - 1) / 2);
end
