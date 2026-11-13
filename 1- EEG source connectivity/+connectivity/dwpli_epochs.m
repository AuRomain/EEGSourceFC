function [d, n_epochs, note] = dwpli_epochs(Z, srate, first, e)
%DWPLI_EPOCHS Debiased wPLI with epochs as trials (Vinck et al. 2011), computed by FieldTrip.
%   [d, n_epochs, note] = connectivity.dwpli_epochs(Z, srate, first, e)
%   Z         nROI x nSamples: band-passed region signals (real), or their
%             analytic signals (complex: the real part, the band-passed
%             signal, is used)
%   first     first sample of every epoch: a vector (one estimate over all
%             its epochs) or a cell of vectors (one estimate per cell, e.g.
%             per window)
%   e         the epochs, connectivity.dwpli_epoch: e.samples per epoch,
%             e.freqs the Fourier frequencies of the band, e.min_epochs
%   d         nEstimates x nROI x nROI, zero diagonal: the mean over e.freqs
%             of the debiased weighted phase lag index; NaN for an estimate
%             with fewer than e.min_epochs epochs, or when the band holds no
%             Fourier frequency of the epochs (note says why)
%   n_epochs  epochs of each estimate
%   Every epoch is one trial, as FieldTrip computes dwPLI:
%   ft_freqanalysis (method 'mtmfft', Hanning taper, no padding, output
%   'fourier', one row per epoch) at the frequencies e.freqs, then
%   ft_connectivity_wpli with 'debias' (Vinck et al. 2011, by its first
%   author; 'isunivariate': the cross-spectra computed epoch by epoch, as
%   ft_connectivityanalysis with method 'wpli_debiased' does for Fourier
%   input) at each frequency, then the mean over the frequencies. dwPLI
%   estimates the squared wPLI: about 0 without lagged coupling (it can be
%   negative), at most 1. Stage 1's band-pass filter is zero-phase and has
%   the same gain in every epoch, so at each frequency the band-passed
%   signals give the dwPLI of the broadband ones (up to the leakage of the
%   Hanning taper between neighbouring frequencies).
if ~iscell(first), first = {first}; end
first = cellfun(@(f) double(f(:)'), first, 'UniformOutput', false);
n = size(Z, 1);
ne = numel(first);
n_epochs = cellfun(@numel, first);
d = nan(ne, n, n);
note = '';
if isempty(e.freqs)
    note = sprintf(['dwPLI: no Fourier frequency of %.3g s epochs (every %.3g Hz) inside the band: use longer ' ...
        'epochs (S.dwpli_epoch_cycles or S.dwpli_epoch_length)'], e.samples / srate, srate / e.samples);
    return
end
use = n_epochs >= e.min_epochs;
if ~any(use)
    note = sprintf(['dwPLI needs at least %d epochs of %.3g s per estimate (trials of the debiased wPLI); ' ...
        'the data hold %d'], e.min_epochs, e.samples / srate, max([0, n_epochs]));
    return
end
all_first = unique([first{use}]);
X = double(real(Z));
L = e.samples;
data = struct('label', {compose('r%d', (1:n)')}, 'fsample', srate, ...
    'trial', {arrayfun(@(f) X(:, f:f + L - 1), all_first, 'UniformOutput', false)}, ...
    'time', {repmat({(0:L - 1) / srate}, 1, numel(all_first))}, ...
    'sampleinfo', [all_first(:), all_first(:) + L - 1]);
cfg = struct('method', 'mtmfft', 'taper', 'hanning', 'foi', e.freqs, 'pad', 'maxperlen', 'output', 'fourier', ...
    'keeptrials', 'yes', 'feedback', 'no', 'showcallinfo', 'no');
txt = evalc('F = ft_freqanalysis(cfg, data);');   % rpttap x chan x freq: one row per epoch (one taper);
for w = regexp(txt, 'Warning: [^\n]*', 'match')     % FieldTrip's lines silenced, its warnings kept
    warning('connectivity:fieldtrip', 'ft_freqanalysis: %s', w{1}(10:end));
end
if numel(F.freq) ~= numel(e.freqs) || max(abs(F.freq(:) - e.freqs(:))) > 1e-6 * max(e.freqs)
    error('connectivity:dwpli', 'FieldTrip computed other frequencies (%s Hz) than asked (%s Hz)', ...
        mat2str(F.freq, 4), mat2str(e.freqs, 4));
end
for k = find(use)
    [~, rows] = ismember(first{k}, all_first);
    W = ft_connectivity_wpli(F.fourierspctrm(rows, :, :), 'debias', true, 'isunivariate', true);
    D = mean(real(W), 3);                           % n x n x freq -> the band's mean
    D(1:n + 1:end) = 0;                             % a region with itself: 0 / 0 in FieldTrip
    d(k, :, :) = D;
end
if ~all(use)
    note = sprintf('dwPLI: %d of %d estimates with fewer than %d epochs of %.3g s are NaN', nnz(~use), ne, ...
        e.min_epochs, e.samples / srate);
end
end
