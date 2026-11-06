function [G, cb, df] = welch_G(Y, X, c, grp)
%WELCH_G Welch-type statistic of one contrast of a linear model with two variance groups.
%   [G, cb, df] = welch_G(Y, X, c, grp)
%   Y    n x m outcomes (one model per column: features, edges or permutations)
%   X    n x p design (full column rank), c  p x 1 contrast
%   grp  n x 1 logical: the two variance groups (the compared groups)
%   G    1 x m statistic, cb  1 x m contrast of the OLS coefficients (c'b)
%   df   1 x m Aspin-Welch degrees of freedom of G (t_to_z converts G to z
%        with them): 1 / sum_g (1 - W_g / sum W)^2 / sum_g R_ii, with W_g
%        the group's weight times its size, as PALM computes them for its v
%        statistic (palm_core.m, fastv, lines 2821-2857; Winkler et al.
%        2014, the Aspin-Welch v of a rank-1 contrast). Without covariates
%        they are the Welch-Satterthwaite degrees of freedom.
%   The G statistic of Winkler et al. (2014, NeuroImage 92:381) for a
%   contrast of rank 1: G = c'b / sqrt(c' inv(X'WX) c), b the OLS
%   coefficients, W diagonal with, for every subject of variance group g,
%   W_g = sum_g R_ii / sum_g r_i^2 (R = I - X pinv(X), r the residuals).
%   Each group's residual variance is estimated on its own, so a group
%   with more variance (often the smaller one) no longer makes the test
%   liberal, as the pooled-variance t of the ANCOVA did (audit V3-3, M3).
%   Without covariates (X = [1, group]) G is Welch's t exactly. A group
%   without residual variance gives an undefined (NaN or Inf) value.
c = c(:);
pX = pinv(X);
b = pX * Y;
r = Y - X * b;
h = 1 - sum(X .* pX', 2);                       % R_ii: 1 - leverage
w0 = sum(h(~grp)) ./ sum(r(~grp, :) .^ 2, 1);   % 1 x m weights of the two groups
w1 = sum(h(grp)) ./ sum(r(grp, :) .^ 2, 1);
% c' inv(w0 A0 + w1 A1) c for every column at once: with A0 + A1 = X'X,
% V'(A0 + A1)V = I and V'A0V = diag(lam), it is sum_k (V'c)_k^2 / (w0 lam_k + w1 (1 - lam_k))
A0 = X(~grp, :)' * X(~grp, :);
A1 = X(grp, :)' * X(grp, :);
[V, D] = eig((A0 + A0') / 2, (A0 + A1 + A0' + A1') / 2, 'chol');
lam = min(max(diag(D), 0), 1);
u = (V' * c) .^ 2;
v = sum(u ./ (lam * w0 + (1 - lam) * w1), 1);
cb = c' * b;
G = cb ./ sqrt(v);
if nargout > 2
    d0 = sum(h(~grp));  d1 = sum(h(grp));
    W0 = w0 * nnz(~grp);  W1 = w1 * nnz(grp);
    df = 1 ./ ((1 - W0 ./ (W0 + W1)) .^ 2 / d0 + (1 - W1 ./ (W0 + W1)) .^ 2 / d1);
end
end
