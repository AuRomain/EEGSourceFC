function e = dwpli_epoch(srate, fmin, fmax, opts)
%DWPLI_EPOCH The epochs (trials) of the debiased wPLI in a band, and the fewest epochs of one estimate.
%   e = connectivity.dwpli_epoch(srate, fmin, fmax, opts)
%   opts.dwpli_epoch_cycles  epoch length in cycles of the band centre
%                            (fmin + fmax) / 2 (default 6: 0.57 s in 8-13 Hz,
%                            1 s in 4-8 Hz); in a band anchored on the IAF,
%                            of its centre at IAF 10 Hz (opts.band_ref,
%                            cycles_centre)
%   opts.dwpli_epoch_length  epoch length in seconds instead, the same in
%                            every band (default []: cycles)
%   e.seconds, e.samples (round(seconds * srate)), e.min_epochs (10: an
%   estimate over fewer epochs is NaN), e.freqs (the Fourier frequencies
%   k * srate / samples inside [fmin, fmax], those dwPLI averages).
%   e = connectivity.dwpli_epoch() gives only e.min_epochs.
%   Debiased wPLI needs trials (Vinck et al. 2011): the epochs are the
%   trials of a continuous recording (connectivity.dwpli_epochs).
e = struct('min_epochs', 10);
if nargin == 0, return; end
if nargin < 4, opts = struct(); end
cycles = 6;  len = [];
if isfield(opts, 'dwpli_epoch_cycles') && ~isempty(opts.dwpli_epoch_cycles), cycles = opts.dwpli_epoch_cycles; end
if isfield(opts, 'dwpli_epoch_length') && ~isempty(opts.dwpli_epoch_length), len = opts.dwpli_epoch_length; end
if isempty(len)
    e.seconds = cycles / cycles_centre(fmin, fmax, opts);   % anchored bands: the centre at IAF 10 Hz
else
    e.seconds = len;
end
e.samples = round(e.seconds * srate);
k = ceil(fmin * e.samples / srate - 1e-9):floor(fmax * e.samples / srate + 1e-9);
e.freqs = k(k > 0) * srate / e.samples;
end
