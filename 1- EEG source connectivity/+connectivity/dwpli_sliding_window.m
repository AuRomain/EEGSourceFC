function [dwpli, static_dwpli, full_dwpli, win, note] = dwpli_sliding_window(data, srate, window, step, seg, e)
%DWPLI_SLIDING_WINDOW Debiased wPLI in sliding windows of epochs (Vinck et al. 2011, computed by FieldTrip).
%   [dwpli, static_dwpli, full_dwpli, win, note] =
%   connectivity.dwpli_sliding_window(data, srate, window, step, seg, e)
%   data    nROI x nSamples band-passed signals (real) or their analytic
%           signals (complex)
%   window, step  in seconds; windows as plv_wpli_sliding_window (from the
%           start of each piece of seg, never across a cut)
%   e       the epochs (connectivity.dwpli_epoch): each window is cut into
%           floor(window / epoch) consecutive epochs from its start, the
%           trials of its debiased wPLI (connectivity.dwpli_epochs:
%           ft_freqanalysis, ft_connectivity_wpli with 'debias'). The
%           epochs are counted from their length in seconds as well as in
%           whole samples, whichever gives more: a window of 10 epochs
%           (60 cycles of 6-cycle epochs) holds 10 although each epoch is
%           rounded to whole samples, the window then a few samples longer
%   dwpli         nWindows x nROI x nROI (single), static_dwpli their mean
%   full_dwpli    over the epochs of the whole recording (each piece cut
%                 into epochs from its start), as static_connectivity
%   win           samples, step, first_sample, centre_sample (as the other
%                 measures), epochs (per window), epoch_samples
%   note          '' or why there are no windows: a window holding fewer
%                 than e.min_epochs epochs has no windows at all (dwpli
%                 0 x nROI x nROI, static_dwpli NaN, first_sample empty)
%   Debiased wPLI needs several trials per estimate: the windows of the
%   other phase measures (6 cycles) hold none, so dwPLI has its own
%   (get_all_connectivity_mats: 60 cycles of the band centre, half
%   overlapping, by default; S.measure_windows.dwpli).
if nargin < 5 || isempty(seg), seg = [1, size(data, 2)]; end
n = size(data, 1);
win_samples = ceil(window * srate);
nb_shifts = ceil(step * srate);
per_window = floor(win_samples / e.samples);
if isfield(e, 'seconds')                           % not one epoch fewer for the epochs' rounding
    per_window = max(per_window, floor(window / e.seconds + 1e-9));
end
win_samples = max(win_samples, per_window * e.samples);
note = '';
starts = connectivity.window_starts(seg, win_samples, nb_shifts);
if per_window < e.min_epochs
    note = sprintf(['dwPLI windows of %.3g s hold %d epochs of %.3g s; the debiased wPLI needs at least %d ' ...
        'per window (S.measure_windows.dwpli): no windows'], window, per_window, e.samples / srate, e.min_epochs);
    starts = zeros(0, 1);
end
first = arrayfun(@(s) s + (0:per_window - 1) * e.samples, starts(:)', 'UniformOutput', false);
if isempty(first)
    dwpli = zeros(0, n, n, 'single');
    static_dwpli = nan(n);  static_dwpli(1:n + 1:end) = 0;
else
    [d, ~, note_w] = connectivity.dwpli_epochs(data, srate, first, e);
    if isempty(note), note = note_w; end
    dwpli = single(d);
    static_dwpli = squeeze(mean(d, 1));
end
full_dwpli = [];
if nargout > 2
    full_dwpli = squeeze(connectivity.dwpli_epochs(data, srate, connectivity.window_starts(seg, e.samples, e.samples), e));
end
win = struct('samples', win_samples, 'step', nb_shifts, 'first_sample', starts(:)', ...
    'centre_sample', starts(:)' + (win_samples - 1) / 2, 'epochs', per_window, 'epoch_samples', e.samples);
end
