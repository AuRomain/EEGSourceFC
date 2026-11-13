function G = gt_window_metrics(W)
%GT_WINDOW_METRICS Strength, clustering and path length of every window.
%   G = gt_window_metrics(W) with W nROI x nROI x nWindows
%   (weights in [0, 1]; negative weights of a signed measure such as AEC are
%   set to 0: the BCT definitions hold for non-negative weights only).
%   Debiased wPLI (dwpli) is centred on 0 without coupling, so its positive
%   part is biased upwards (half-normal noise: about 0.4 x its SD): its
%   graph metrics mostly measure estimator noise. The static features do
%   not call this function for dwPLI (static_features: the signed
%   strength only); the states' graph metrics
%   do, with a warning when dwPLI is the state measure (check_settings):
%   use wPLI for graph metrics. Returns nWindows x 1 vectors strength, clustering,
%   path_length. Definitions of BCT (strengths_und,
%   clustering_coef_wu, weight_conversion 'lengths' + distance_wei +
%   charpath), computed for many windows at once: clustering with
%   pagemtimes, shortest paths with a batched Floyd-Warshall, which gives
%   the same distances as distance_wei (to 1e-12: tests/test_stage3.m)
%   many times faster (distance_wei takes about 20 ms per 68-node window).
%   Small-worldness is not computed: on real recordings the weighted index
%   against weight-shuffled null graphs is about 1 in every window, so it
%   carries no information, and its null graphs would be most of the cost.
%   Used once per subject and shared by the raw and smoothed labels; also
%   the graph metrics of the static matrices (static_features: one
%   "window" per recording), hence at the top of the pipeline folder.
[n, ~, nw] = size(W);
W = max(double(W), 0);                                  % weights must be >= 0: a signed measure (AEC)
W(repmat(logical(eye(n)), 1, 1, nw)) = 0;               % enters with its positive part only
low = find(tril(true(n), -1));
X = reshape(W, n * n, nw);
X = X(low, :)';                                         % windows x edges
G.strength = 2 * sum(X, 2) / n;
[G.clustering, G.path_length] = batch(X, n, low);
end

function [C, L] = batch(X, n, low)
nw = size(X, 1);
C = zeros(nw, 1);  L = C;
off = ~eye(n);
d = 1:n + 1:n * n;
for b0 = 1:500:nw
    r = b0:min(b0 + 499, nw);
    m = numel(r);
    A = zeros(n, n, m);
    A(low + (0:m - 1) * n * n) = X(r, :)';
    A = A + permute(A, [2 1 3]);
    A3 = A .^ (1 / 3);
    cyc = reshape(pagemtimes(pagemtimes(A3, A3), A3), n * n, m);
    cyc = cyc(d, :);
    K = reshape(sum(A ~= 0, 2), n, m);
    K(cyc == 0) = inf;
    C(r) = mean(cyc ./ (K .* (K - 1)), 1)';
    D = 1 ./ A;
    D(A == 0) = inf;
    D(repmat(logical(eye(n)), 1, 1, m)) = 0;
    for k = 1:n
        D = min(D, D(:, k, :) + D(k, :, :));
    end
    D = reshape(D, n * n, m);
    D = D(off(:), :);
    Lm = mean(D, 1)';
    Lm(any(isinf(D), 1)) = NaN;
    L(r) = Lm;
end
end
