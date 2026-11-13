function fit = psd_specparam(f, P, fit_range, power_line)
%PSD_SPECPARAM Aperiodic (1/f) and periodic parts of power spectra (specparam / FOOOF).
%   fit = psd_specparam(f, P, fit_range, power_line)
%   f          1 x F frequencies (Hz); P units x F power spectral densities
%   fit_range  [fmin fmax] Hz fitted (default of the PSD step: 2-40 Hz)
%   power_line S.psd.power_line: '50' (default), '60' or 'none'. The fit
%              leaves out 1, 2 and 3 times this frequency +/- 2 Hz where
%              they lie inside fit_range (process_fooof's mask; with 2-40
%              Hz none does; 'none' leaves out nothing). A line elsewhere
%              stays in the fit: psd_line_check flags it (docs/psd.md,
%              "Line noise")
%   Each spectrum is fitted with specparam (Donoghue et al. 2020, Nat
%   Neurosci 23:1655) as FieldTrip's ft_freqanalysis does for its outputs
%   'fooof*' (and spectral_profile_features, of the statistics): Brainstorm's MATLAB
%   implementation (process_fooof, FOOOF_matlab, in FieldTrip's
%   external/brainstorm) with Brainstorm's defaults and FieldTrip's
%   settings: aperiodic mode 'fixed' (no knee: log10 P = offset - exponent
%   log10 f), at most 3 Gaussian peaks of 0.5-12 Hz width (2 SD) and at
%   least 0.3 (log10) above the aperiodic fit, 2 SD peak threshold;
%   constrained peak fitting when the Optimization Toolbox is present.
%   fit (units in rows):
%     offset, exponent   the aperiodic parameters (log10 units^2/Hz; 1/f^exponent).
%                        A spectrum that rises over the fit range has a
%                        negative exponent: process_fooof returns the
%                        absolute value (its line 283) after building its
%                        aperiodic fit with the signed one, and the sign is
%                        taken from that fit (else such a spectrum would get
%                        a positive exponent and an aperiodic curve of the
%                        wrong slope)
%     r2, error          R^2 and mean squared error of the model, log10 power
%     peaks              cell: each unit's peaks [centre (Hz), height above
%                        the aperiodic fit (log10), Gaussian SD (Hz; the
%                        bandwidth specparam reports in Python is 2 SD)]
%     ap                 units x F: the aperiodic fit (log10 power) at the
%                        frequencies f inside fit_range (NaN outside)
%     optimisation       'fmincon' or 'fminsearch'
ft_hastoolbox('brainstorm', 1);                   % FieldTrip's copy of Brainstorm's specparam, as ft_freqanalysis
if nargin < 4 || isempty(power_line), power_line = '50'; end
o = getfield(process_fooof('GetDescription'), 'options'); %#ok<GFLD>
opt = struct('freq_range', fit_range, 'peak_width_limits', o.peakwidth.Value{1}, ...
    'max_peaks', o.maxpeaks.Value{1}, 'min_peak_height', o.minpeakheight.Value{1} / 10, ...
    'aperiodic_mode', 'fixed', 'peak_threshold', 2, 'return_spectrum', 0, 'border_threshold', 1, ...
    'power_line', char(power_line), 'peak_type', o.peaktype.Value, 'proximity_threshold', o.proxthresh.Value{1}, ...
    'guess_weight', o.guessweight.Value, 'thresh_after', true);
hot = exist('fmincon', 'file') ~= 0;              % as ft_freqanalysis
n = size(P, 1);
TF = reshape(double(P), n, 1, []);
[~, fg] = process_fooof('FOOOF_matlab', TF, f, opt, hot);
in = f >= fit_range(1) & f <= fit_range(2);
fit = struct('offset', nan(n, 1), 'exponent', nan(n, 1), 'r2', nan(n, 1), 'error', nan(n, 1), ...
    'peaks', {cell(n, 1)}, 'ap', nan(n, numel(f)), 'optimisation', 'fminsearch');
if hot, fit.optimisation = 'fmincon'; end
for c = 1:n
    a = fg(c).aperiodic_params;
    ex = a(end);
    apf = fg(c).ap_fit;                            % 10^(the aperiodic fit), built with the signed exponent
    if numel(apf) > 1 && apf(end) > apf(1), ex = -abs(ex); end   % rising over the fit range
    fit.offset(c) = a(1);
    fit.exponent(c) = ex;
    fit.r2(c) = fg(c).r_squared;
    fit.error(c) = fg(c).error;
    pk = reshape(fg(c).peak_params, [], 3);
    fit.peaks{c} = pk(pk(:, 1) > 0 & pk(:, 3) > 0, :);   % without a peak, FOOOF gives a row of zeros
    fit.ap(c, in) = a(1) - ex * log10(f(in));
end
end
