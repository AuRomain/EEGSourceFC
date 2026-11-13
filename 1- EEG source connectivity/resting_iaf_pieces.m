function [pSum, pChans, f] = resting_iaf_pieces(data, pieces, cmin, fRange, Fs, w, Fw, k)
%RESTING_IAF_PIECES restingIAF (Corcoran et al. 2018) on a recording cut into continuous pieces.
%   [pSum, pChans, f] = resting_iaf_pieces(data, pieces, cmin, fRange, Fs, w, Fw, k)
%   data    channels x samples (rows: channels)
%   pieces  continuous pieces, rows [first last] sample of data, each at
%           least 4 s long
%   cmin, fRange, Fs, w, Fw, k: as restingIAF (minimum channels, frequency
%           range, sampling rate, alpha search window, Savitzky-Golay frame
%           width and order); restingIAF's defaults for the rest (mpow 1,
%           mdiff 0.20, Hamming windows of 4 s, 50% overlap, FFT length the
%           next power of 2, normalised spectra)
%   The steps of restingIAF.m (github.com/corcorana/restingIAF, GPL-3, (c)
%   Andrew W. Corcoran 2016-2018; the version in Toolboxes/restingIAF), with
%   one change: each channel's Welch spectrum is the mean of the Welch
%   spectra of the pieces, weighted by their numbers of segments, i.e. the
%   average over the 4 s segments that lie inside a piece. restingIAF reads
%   one continuous signal: given the joined pieces, its segments would cross
%   the joins, whose steps add broadband power. With a single piece this is
%   restingIAF exactly
%   (tests/test_stage1.m). The smoothing, peak, gravity and averaging steps
%   are restingIAF's own functions (sgfDiff, peakBounds, chanGravs,
%   chanMeans), called as restingIAF calls them.
%   Reference: Corcoran AW, Alday PM, Schlesewsky M, Bornkessel-Schlesewsky
%   I (2018). Toward a reliable, automated method of individual alpha
%   frequency (IAF) quantification. Psychophysiology 55:e13064.
mpow = 1;  mdiff = 0.20;
tlen = round(Fs * 4);                            % restingIAF's defaults
tover = [];                                      % pwelch: 50%
nfft = [];                                       % pwelch: max(256, next power of 2)
nchan = size(data, 1);
pChans = struct('pxx', [], 'minPow', [], 'd0', [], 'd1', [], 'd2', [], ...
    'peaks', [], 'pos1', [], 'pos2', [], 'f1', [], 'f2', [], 'inf1', [], 'inf2', [], ...
    'Q', [], 'Qf', [], 'gravs', [], 'selP', [], 'selG', []);
nseg = floor((diff(pieces, 1, 2) + 1 - tlen / 2) / (tlen / 2));   % segments of each piece (50% overlap)
nseg = max(nseg, 0);
for kx = 1:nchan
    if any(isnan(data(kx, :)))
        warning('Row #%s contains NaNs, skipping channel...', num2str(kx))
        nchan = nchan - 1;
        continue
    end
    use = find(nseg > 0)';
    if isscalar(use)                             % one piece: restingIAF's spectrum, bit for bit
        [pxx, f] = pwelch(data(kx, pieces(use, 1):pieces(use, 2)), hamming(tlen), tover, nfft, Fs);
    else
        pxx = 0;
        for p = use
            [pk, f] = pwelch(data(kx, pieces(p, 1):pieces(p, 2)), hamming(tlen), tover, nfft, Fs);
            pxx = pxx + nseg(p) * pk;
        end
        pxx = pxx / sum(nseg(use));
    end
    frex = dsearchn(f, fRange(1)):dsearchn(f, fRange(2));
    f = f(frex);
    pxx = pxx(frex);
    pChans(kx).pxx = pxx / mean(pxx);
    [pfit, sig] = polyfit(f, log10(pChans(kx).pxx), 1);
    [yval, del] = polyval(pfit, f, sig);
    pChans(kx).minPow = yval + (mpow * del);
    [pChans(kx).d0, pChans(kx).d1, pChans(kx).d2] = sgfDiff(pChans(kx).pxx, Fw, k, Fs, tlen);
    fres = Fs / 2 .^ nextpow2(tlen);
    [pChans(kx).peaks, pChans(kx).pos1, pChans(kx).pos2, pChans(kx).f1, pChans(kx).f2, ...
        pChans(kx).inf1, pChans(kx).inf2, pChans(kx).Q, pChans(kx).Qf] ...
        = peakBounds(pChans(kx).d0, pChans(kx).d1, pChans(kx).d2, f, w, pChans(kx).minPow, mdiff, fres);
end
[gravs, selG, iaw] = chanGravs([pChans(:).d0], f, [pChans(:).f1], [pChans(:).f2]);
[selP, pSum] = chanMeans(gravs, selG, [pChans(:).peaks], [pChans(:).d0], [pChans(:).Qf], cmin);
for cx = 1:nchan
    pChans(cx).gravs = gravs(cx);
    pChans(cx).selP = selP(cx);
    pChans(cx).selG = selG(cx);
end
pSum.pSel = sum(selP);
pSum.gSel = sum(selG);
pSum.iaw = iaw;
end
