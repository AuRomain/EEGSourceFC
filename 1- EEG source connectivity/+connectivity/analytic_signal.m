function Z = analytic_signal(data, seg)
%ANALYTIC_SIGNAL Hilbert transform of each continuous piece separately.
%   Z = connectivity.analytic_signal(data, seg) with data nSignals x
%   nSamples (real, band-passed) and seg nPieces x 2. Complex data are
%   returned unchanged (already analytic). One transform per piece and per
%   recording instead of one per window.
if ~isreal(data)
    Z = data;
    return
end
if nargin < 2 || isempty(seg), seg = [1, size(data, 2)]; end
Z = complex(zeros(size(data)));
for k = 1:size(seg, 1)
    idx = seg(k, 1):seg(k, 2);
    Z(:, idx) = hilbert(data(:, idx).').';
end
end
