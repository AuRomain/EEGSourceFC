function [C, f, n_seg] = psd_cross_spectra(X, fs, segments, seg_s, overlap, frange)
%PSD_CROSS_SPECTRA Cross-spectral density of every pair of rows of X (Welch's method, as cpsd), inside the continuous pieces.
%   [C, f, n_seg] = psd_cross_spectra(X, fs, segments, seg_s, overlap, frange)
%   The arguments are those of psd_welch: X signals x samples (the
%   channels), fs the sampling rate, segments the continuous pieces (rows
%   [first last] sample), seg_s the segment length (s), overlap the
%   fraction shared by consecutive segments, frange [fmin fmax] Hz kept.
%   C   signals x signals x frequencies: the real part of the one-sided
%       cross-spectral density, the same as MATLAB's cpsd(x, x, hann(L),
%       noverlap, L, fs, 'mimo') (Signal Processing Toolbox) gives for
%       each piece, so diag(C(:, :, k)) is psd_welch's (pwelch's) P(:, k);
%       tests/test_psd.m checks both to 1e-12
%   n_seg  the number of segments averaged, as psd_welch's
%   Each piece, its mean removed, is cut into the segments pwelch and cpsd
%   use (fix((samples - noverlap) / (L - noverlap)) of L samples, L - noverlap
%   apart, from the first sample), each multiplied by the Hann window and
%   Fourier-transformed (FFT of L points); C is the mean over all segments
%   of all pieces of Re(F' * F) / (fs * sum(w.^2)), doubled except at 0 Hz and
%   the Nyquist frequency (one-sided). The products are summed here instead
%   of calling cpsd, which is far slower for many channels.
%   Only the real part is kept: for real weights w (the inverse filters of
%   a source) the power of w' * X is w' * C * w, and the imaginary part of
%   a Hermitian matrix adds nothing to it. This is how source power is
%   computed from the channels' cross-spectra (LORETA-KEY's spectra of
%   current density; FieldTrip's ft_sourceanalysis on a 'freq' input with
%   'csd'; the PSD of every source in Brainstorm), without the sources'
%   time series.
L = round(seg_s * fs);
noverlap = L - max(1, round(L * (1 - overlap)));
step = L - noverlap;
w = hann(L);
X = double(X);
m = size(X, 1);
fall = (0:floor(L / 2))' * fs / L;                 % one-sided frequencies (pwelch's)
keep = find(fall >= frange(1) - 1e-9 & fall <= frange(2) + 1e-9);
scale = 2 * ones(numel(fall), 1);                  % one-sided: twice, but 0 Hz and Nyquist once
scale(1) = 1;
if mod(L, 2) == 0, scale(end) = 1; end
scale = scale(keep) / (fs * sum(w .^ 2));
nk = numel(keep);
acc = zeros(m, m, nk);
n_seg = 0;
for p = 1:size(segments, 1)
    x = X(:, segments(p, 1):segments(p, 2))';    % samples x signals
    if size(x, 1) < L, continue; end
    x = x - mean(x, 1);
    k = fix((size(x, 1) - noverlap) / step);
    starts = (0:k - 1) * step;
    for b = 1:200:k                                % 200 segments at a time
        sb = starts(b:min(b + 199, k));
        seg = zeros(L, numel(sb), m);
        for j = 1:numel(sb)
            seg(:, j, :) = reshape(x(sb(j) + (1:L), :) .* w, L, 1, m);
        end
        F = fft(seg, L, 1);                        % L x segments x signals
        F = F(keep, :, :);
        for q = 1:nk
            A = reshape(F(q, :, :), numel(sb), m); % segments x signals
            acc(:, :, q) = acc(:, :, q) + real(A)' * real(A) + imag(A)' * imag(A);
        end
    end
    n_seg = n_seg + k;
end
if n_seg == 0
    error('psd_cross_spectra:pieces', 'No continuous piece of %g s (the segment length)', seg_s);
end
C = acc .* reshape(scale, 1, 1, nk) / n_seg;
C = (C + permute(C, [2 1 3])) / 2;                  % exactly symmetric
f = fall(keep)';
end
