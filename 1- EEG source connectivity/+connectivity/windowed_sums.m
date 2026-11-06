function [W, total] = windowed_sums(X, starts, win_samples, in_pieces)
%WINDOWED_SUMS Sum of each row of X over every window, by cumulative sums.
%   [W, total] = connectivity.windowed_sums(X, starts, win_samples, in_pieces)
%   X: nRows x nSamples; W(:, w) = sum(X(:, starts(w) + (0:win_samples-1)), 2);
%   total = sum over the samples flagged in in_pieces (all when omitted).
%   One pass over the data whatever the number of windows (the original
%   recomputed every window from scratch).
cs = cumsum([zeros(size(X, 1), 1, 'like', X), X], 2);
W = cs(:, starts + win_samples) - cs(:, starts);
if nargin < 4 || isempty(in_pieces)
    total = cs(:, end);
else
    total = sum(X(:, in_pieces), 2);
end
end
