function [P, f, n_seg] = psd_welch(X, fs, segments, seg_s, overlap, frange)
%PSD_WELCH Power spectral density of each row of X (MATLAB's pwelch), inside the continuous pieces.
%   [P, f, n_seg] = psd_welch(X, fs, segments, seg_s, overlap, frange)
%   X         signals x samples (channels, or regions' time series)
%   fs        sampling rate (Hz)
%   segments  continuous pieces, rows [first last] sample (get_eeg_data):
%             no segment crosses a cut (a rejected stretch, an epoch limit)
%   seg_s     segment length (s); the frequency resolution is 1 / seg_s
%   overlap   fraction of a segment shared with the next (0.5 = 50%)
%   frange    [fmin fmax] Hz kept
%   Each piece, its mean removed, goes to pwelch (Signal Processing
%   Toolbox) with a Hann window of seg_s (hann, symmetric), the overlap,
%   an FFT of the segment's length and fs: the mean of its segments'
%   modified periodograms, a one-sided density in units^2/Hz (microvolts^2/Hz
%   for the channels) whose integral over 0 to fs/2 is the variance of the
%   signal. The pieces' spectra are averaged with the number of segments
%   of each as weight, so P is the mean over all segments of all pieces,
%   as pwelch gives for one piece. Pieces shorter than a segment are not
%   used; n_seg is the number of segments averaged (pwelch's: fix((samples
%   - overlap) / (segment - overlap)) per piece).
L = round(seg_s * fs);
noverlap = L - max(1, round(L * (1 - overlap)));
w = hann(L);
X = double(X);
acc = 0;  n_seg = 0;  f = [];
for p = 1:size(segments, 1)
    x = X(:, segments(p, 1):segments(p, 2))';    % samples x signals, as pwelch takes them
    if size(x, 1) < L, continue; end
    k = fix((size(x, 1) - noverlap) / (L - noverlap));
    [Pp, f] = pwelch(x - mean(x, 1), w, noverlap, L, fs, 'onesided');
    acc = acc + k * Pp;
    n_seg = n_seg + k;
end
if n_seg == 0
    error('psd_welch:pieces', 'No continuous piece of %g s (the segment length)', seg_s);
end
keep = f >= frange(1) - 1e-9 & f <= frange(2) + 1e-9;
P = (acc(keep, :) / n_seg)';
f = f(keep)';
end
