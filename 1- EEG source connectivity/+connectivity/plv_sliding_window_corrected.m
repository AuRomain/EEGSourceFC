function [plv, static_plv, full_plv, win, note] = plv_sliding_window_corrected(data, srate, window, step, seg)
%PLV_SLIDING_WINDOW_CORRECTED PLV after symmetric (multivariate) leakage correction.
%   [plv, static_plv, full_plv, win, note] = connectivity.plv_sliding_window_corrected(data, srate, window, step, seg)
%   The signals (band-passed: real, or their analytic signal, complex,
%   computed before the filter edges were trimmed) are made mutually orthogonal over the
%   whole recording (connectivity.symmetric_orthogonalisation, Colclough et
%   al. 2015; each signal keeps its magnitude), then PLV is computed in
%   the windows of plv_wpli_sliding_window. Needs as many
%   independent signals as regions: otherwise every edge is NaN and note
%   says why (see aec_sliding_window_corrected).
if nargin < 5 || isempty(seg), seg = [1, size(data, 2)]; end
Z = connectivity.analytic_signal(data, seg);          % unchanged when already complex
[Zo, note] = connectivity.symmetric_orthogonalisation(Z, '''wpli'' or ''ciplv''');
if isempty(note)
    [plv, static_plv, ~, ~, full_plv, ~, win] = connectivity.plv_wpli_sliding_window(Zo, srate, window, step, seg);
    return
end
% not computed: the windows of the uncorrected signals, every edge NaN
[plv, static_plv, ~, ~, full_plv, ~, win] = connectivity.plv_wpli_sliding_window(Z, srate, window, step, seg);
[plv, static_plv, full_plv] = connectivity.not_computed(plv, static_plv, full_plv);
end
