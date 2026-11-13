function [hz, ratio, flagged, n_units] = psd_line_check(f, P)
%PSD_LINE_CHECK The strongest narrow line in a recording's spectra (mains, a stimulator, a device).
%   [hz, ratio, flagged, n_units] = psd_line_check(f, P)
%   f   1 x F frequencies (Hz), evenly spaced (psd_grid: 1 / S.psd.segment Hz)
%   P   units x F PSDs (linear), e.g. the channels of one recording
%   A narrow peak is a local maximum whose excess over the background (the
%   median power 4 to 8 bins away, 1-2 Hz with 4 s segments) falls to at
%   most 10% of itself 2 bins away, where a sinusoid's spectrum in Hann
%   segments has its first zero (under 4% when the line lies between two
%   bins); alpha and other rhythms, which are wider, keep much more of it
%   in real recordings. On the mean of the units' spectra: hz, the
%   frequency of the narrow peak with the highest ratio of its power to
%   the background (NaN and 1 without a narrow peak); flagged = ratio >= 2
%   (the line at least doubles its bin's power). n_units: the units whose
%   own spectrum has such a line (ratio >= 2, at any frequency). A line is
%   not removed: it enters the band powers, the total of relative power,
%   the specparam fit and the periodic power (docs/psd.md, "Line noise").
P = double(P);
[hz, ratio] = strongest(f, mean(P, 1));
flagged = ratio >= 2;
n_units = 0;
for u = 1:size(P, 1)
    [~, r] = strongest(f, P(u, :));
    n_units = n_units + (r >= 2);
end
end

function [hz, best] = strongest(f, p)
hz = NaN;  best = 1;
F = numel(p);
for i = 3:F - 2
    if p(i) < p(i - 1) || p(i) < p(i + 1), continue; end           % a local maximum
    d = abs((1:F) - i);
    bg = median(p(d >= 4 & d <= 8));
    e0 = p(i) - bg;
    if ~(e0 > 0) || max(p(i - 2), p(i + 2)) - bg > 0.1 * e0, continue; end   % narrow
    if p(i) / bg > best
        best = p(i) / bg;  hz = f(i);
    end
end
end
