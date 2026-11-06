function [aec, static_aec, full_aec, win, note] = aec_sliding_window_corrected(data, srate, window, step, seg)
%AEC_SLIDING_WINDOW_CORRECTED AEC after symmetric (multivariate) leakage correction.
%   [aec, static_aec, full_aec, win, note] = connectivity.aec_sliding_window_corrected(data, srate, window, step, seg)
%   data: band-passed signals (real), or their analytic signal (complex,
%   computed before the filter edges were trimmed); the signals are made mutually
%   orthogonal over the whole recording (closest orthogonal matrix, Colclough
%   et al. 2015; connectivity.symmetric_orthogonalisation), then AEC is
%   computed in windows as aec_sliding_window. The original used
%   ROInets.symmetric_orthogonalise without keeping each signal's
%   magnitude, and abs() on the correlations (E11). Needs as many
%   independent signals as regions (at most the data rank of the channels
%   minus 1 after the average reference): otherwise every edge is NaN (the
%   diagonal 0) and note says why (it stopped the whole recording with an
%   error; use 'aec_orth', the pairwise correction, or fewer regions).
if nargin < 5 || isempty(seg), seg = [1, size(data, 2)]; end
Z = connectivity.analytic_signal(data, seg);          % unchanged when already complex
[Zo, note] = connectivity.symmetric_orthogonalisation(Z, '''aec_orth''');
if isempty(note)
    [aec, static_aec, full_aec, win] = connectivity.aec_sliding_window(Zo, srate, window, step, seg);
    return
end
% not computed: the windows of the uncorrected signals, every edge NaN
[aec, static_aec, full_aec, win] = connectivity.aec_sliding_window(Z, srate, window, step, seg);
[aec, static_aec, full_aec] = connectivity.not_computed(aec, static_aec, full_aec);
end
