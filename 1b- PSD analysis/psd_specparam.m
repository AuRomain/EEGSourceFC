function fit = psd_specparam(f, P, fit_range)
%PSD_SPECPARAM Aperiodic (1/f) and periodic parts of power spectra (specparam / FOOOF).
%   fit = psd_specparam(f, P, fit_range)
%   f          1 x F frequencies (Hz); P units x F power spectral densities
%   fit_range  [fmin fmax] Hz fitted (default of the PSD step: 2-40 Hz)
%   Each spectrum is fitted with specparam (Donoghue et al. 2020, Nat
%   Neurosci 23:1655) as FieldTrip's ft_freqanalysis does for its outputs
%   'fooof*' (and spectral_profile_features, stage 3): Brainstorm's MATLAB
%   implementation (process_fooof, FOOOF_matlab, in FieldTrip's
%   external/brainstorm) with Brainstorm's defaults and FieldTrip's
%   settings: aperiodic mode 'fixed' (no knee: log10 P = offset - exponent
%   log10 f), at most 3 Gaussian peaks of 0.5-12 Hz width and at least
%   0.3 (log10) above the aperiodic fit, 2 SD peak threshold;
%   constrained peak fitting when the Optimization Toolbox is present.
%   fit (units in rows):
%     offset, exponent   the aperiodic parameters (log10 units^2/Hz; 1/f^exponent)
%     r2, error          R^2 and mean squared error of the model, log10 power
%     peaks              cell: each unit's peaks [centre (Hz), height above
%                        the aperiodic fit (log10), bandwidth (Hz)]
%     ap                 units x F: the aperiodic fit (log10 power) at the
%                        frequencies f inside fit_range (NaN outside)
%     optimisation       'fmincon' or 'fminsearch'
ft_hastoolbox('brainstorm', 1);                   % FieldTrip's copy of Brainstorm's specparam, as ft_freqanalysis
o = getfield(process_fooof('GetDescription'), 'options'); %#ok<GFLD>
opt = struct('freq_range', fit_range, 'peak_width_limits', o.peakwidth.Value{1}, ...
    'max_peaks', o.maxpeaks.Value{1}, 'min_peak_height', o.minpeakheight.Value{1} / 10, ...
    'aperiodic_mode', 'fixed', 'peak_threshold', 2, 'return_spectrum', 0, 'border_threshold', 1, ...
    'power_line', '50', 'peak_type', o.peaktype.Value, 'proximity_threshold', o.proxthresh.Value{1}, ...
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
    fit.offset(c) = a(1);
    fit.exponent(c) = a(end);
    fit.r2(c) = fg(c).r_squared;
    fit.error(c) = fg(c).error;
    pk = reshape(fg(c).peak_params, [], 3);
    fit.peaks{c} = pk(pk(:, 1) > 0 & pk(:, 3) > 0, :);   % without a peak, FOOOF gives a row of zeros
    fit.ap(c, in) = a(1) - a(end) * log10(f(in));
end
end
