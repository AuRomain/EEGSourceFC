function M = psd_measures(f, P, fit, o, iaf, iaf_opt)
%PSD_MEASURES Band power and spectral measures of each channel or region of one recording.
%   M = psd_measures(f, P, fit, o, iaf, iaf_opt)
%   f, P     frequencies (1 x F) and PSD (units x F, units^2/Hz: psd_welch)
%   fit      their specparam fit (psd_specparam); [] for spectra without a
%            fit (the sources' PSDs of the vertex maps): absolute and
%            relative only, the measures of the fit NaN
%   o        S.psd: bands, total, fit_range
%   iaf      the recording's individual alpha frequency (NaN if none) for
%            bands anchored on it (band_edges); iaf_opt = S.iaf
%   Per unit (rows) and band (columns of o.bands):
%   absolute   log10 of the band power: the integral of the PSD from the
%              band's lower to its upper edge (trapezoids, the PSD linear
%              between grid points, so adjacent bands add up to the
%              power of both); log10 units^2 (microvolts^2 at the scalp)
%   relative   the band power in % of the power over o.total (1-30 Hz)
%   periodic   the power above the aperiodic fit: the mean of log10 PSD
%              minus the fit over the band's part inside o.fit_range, from
%              edge to edge (its integral by trapezoids, the edges
%              interpolated as for absolute, over its width; log10 units:
%              0.3 = twice the 1/f background; NaN when less than 1 Hz of
%              the band is inside the fit range). Not the plain mean of the
%              grid points inside the band, whose number changes with where
%              an anchored band's edges fall on the grid
%   Per unit:
%   exponent, offset   the aperiodic fit (log10 P = offset - exponent log10 f)
%   fit_r2             R^2 of the specparam model
%   alpha_frequency    the centre of the strongest specparam peak inside
%              iaf_opt.search (7-13 Hz); NaN without one. The centre of
%              gravity (alpha_cog) does not replace a missing peak: on
%              spectra without alpha it lies near the middle of the range,
%              and the measure would mix two kinds of value
%   alpha_peak         true when there is such a peak
%   alpha_cog          the centre of gravity of the power above the
%              aperiodic fit in iaf_opt.search, when that power is at least
%              5% of the range's power (NaN otherwise), with or without a
%              peak: reported (values_<level>_alpha_cog.csv), not tested
%   edges              bands x 2: the edges used (anchored: this IAF's)
nb = numel(o.bands);
n = size(P, 1);
P = double(P);
if isempty(fit)                                    % no fit: its measures NaN
    fit = struct('exponent', nan(n, 1), 'offset', nan(n, 1), 'r2', nan(n, 1), 'ap', nan(size(P)), ...
        'peaks', {repmat({zeros(0, 3)}, n, 1)});
end
edges = nan(nb, 2);
for b = 1:nb
    [edges(b, 1), edges(b, 2)] = band_edges(o.bands(b), iaf, iaf_opt);
end
total = band_power(f, P, o.total(1), o.total(2));
M = struct('absolute', nan(n, nb), 'relative', nan(n, nb), 'periodic', nan(n, nb), ...
    'exponent', fit.exponent, 'offset', fit.offset, 'fit_r2', fit.r2, 'alpha_frequency', nan(n, 1), ...
    'alpha_peak', false(n, 1), 'alpha_cog', nan(n, 1), 'edges', edges);
flat = log10(P) - fit.ap;                          % NaN outside the fit range
infit = find(f >= o.fit_range(1) - 1e-9 & f <= o.fit_range(2) + 1e-9);
fit_lo = f(infit(1));  fit_hi = f(infit(end));     % the fit's first and last frequencies
for b = 1:nb
    lo = edges(b, 1);  hi = edges(b, 2);
    if ~(hi > lo) || lo < f(1) - 1e-9 || hi > f(end) + 1e-9, continue; end   % no IAF, or outside the spectra
    bp = band_power(f, P, lo, hi);
    M.absolute(:, b) = log10(bp);
    M.relative(:, b) = 100 * bp ./ total;
    a = max(lo, fit_lo);  z = min(hi, fit_hi);      % the band's part inside the fit
    if z - a >= 1
        M.periodic(:, b) = band_power(f, flat, a, z) / (z - a);   % its mean from edge to edge
    end
end
s = iaf_opt.search;
in = f >= s(1) & f <= s(2);
above = max(P(:, in) - 10 .^ fit.ap(:, in), 0);   % the power above the 1/f fit
for c = 1:n
    pk = fit.peaks{c};
    if ~isempty(pk), pk = pk(pk(:, 1) >= s(1) & pk(:, 1) <= s(2), :); end
    if ~isempty(pk)
        [~, j] = max(pk(:, 2));
        M.alpha_frequency(c) = pk(j, 1);
        M.alpha_peak(c) = true;
    end
    if sum(above(c, :)) >= 0.05 * sum(P(c, in))
        M.alpha_cog(c) = sum(f(in) .* above(c, :)) / sum(above(c, :));
    end
end
end

function bp = band_power(f, P, lo, hi)
% the integral of P (units x F) from lo to hi, P linear between grid points
% (P: the PSD, or log10 PSD minus the fit for the periodic power)
in = find(f > lo + 1e-9 & f < hi - 1e-9);
x = [lo, f(in), hi];
y = [interp1(f, P', lo)', P(:, in), interp1(f, P', hi)'];
bp = trapz(x, y, 2);
end
