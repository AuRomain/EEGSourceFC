function [R, D] = perm_ancova_maxT(Y, isB, Z, opts)
%PERM_ANCOVA_MAXT Two-group permutation test per column (Welch t or Welch-type ANCOVA), max-statistic across columns.
%   R = perm_ancova_maxT(Y, isB, Z, opts)
%   [R, D] = perm_ancova_maxT(...) also returns the permutation statistics,
%   for a cluster test over the columns (cluster_perm_test) and the family
%   test: D.t_obs (1 x K) and D.T (n_perm x K, the same permutations in
%   every column) are z values (below), D.df (1 x K) Inf (NaN when not
%   tested), so their users read them as normal deviates
%   Y     n x K outcome (one column per state; NaN = missing for that state)
%   isB   n x 1 logical, the group whose difference is reported: effects are
%         (isB) - (~isB)
%   Z     n x q covariates ([] for none)
%   opts  .n_perm, .seed, .n_boot (bootstrap of Hedges' g; 0 = no CI)
%
%   Without covariates: Welch t; group labels are permuted.
%   With covariates: ANCOVA y = b0 + b1*group + Z*c (both groups in one
%   model, S2), b1 tested with the Welch-type G statistic (welch_G: each
%   group's residual variance estimated on its own; the pooled-variance t
%   was liberal when the groups differ in size and variance, audit V3-3,
%   M3); Freedman-Lane permutation: residuals of the reduced model y ~ Z are
%   permuted and the full model refitted (S14).
%   The statistic of the data and of every permutation is converted to z
%   with its own Aspin-Welch (Welch-Satterthwaite without covariates)
%   degrees of freedom (welch_G, t_to_z) before the p-values are counted,
%   as PALM does for its v statistic (Winkler et al. 2014, NeuroImage
%   92:381; palm_core.m: fastv, then palm_gtoz, lines 1579 and 1603-1605
%   with -twotail): the degrees of freedom change with the variances of
%   each permutation and with the column's subjects, so the raw statistic
%   is not on one scale across permutations or columns, its z is
%   (5 November 2026; before, the raw t or G, whose maximum gave columns
%   with fewer subjects more weight). tests/test_reference_statistics
%   compares G, df, p_unc and p_maxT with PALM itself.
%   Missing values are dropped per column. Every column uses the same
%   random permutation of the n subjects, restricted to its own subjects by
%   following the permutation's cycles (identical for columns without
%   missing values), so the max-statistic distribution is the joint one.
%   An undefined statistic (fewer than 2 subjects per group, no variance, a
%   rank-deficient design, or a model that fits the outcome perfectly:
%   residual sum of squares <= 1e-12 times the total, as when the outcome
%   is one of its own covariates) gives NaN, never a p-value.
%
%   R (table, one row per column): nA, nB, meanA, meanB, diff (B - A, or
%   the adjusted difference b1), t (Welch t, or G with covariates), df (its
%   Welch-Satterthwaite or Aspin-Welch degrees of freedom), p_unc, p_maxT
%   (FWER across the K columns: the largest |z| of each permutation),
%   g (Hedges' g of the raw values on the same subjects as p,
%   S16), g_ci_low, g_ci_high (percentile bootstrap), d_adjusted (b1 over
%   the pooled residual SD of the ANCOVA, an effect size; NaN without
%   covariates).

[n, K] = size(Y);
if nargin < 3, Z = []; end
isB = logical(isB(:));
rs = RandStream('mrg32k3a', 'Seed', opts.seed);
[~, P] = sort(rand(rs, n, opts.n_perm), 1);       % one permutation of all subjects per column of P

t_obs = nan(1, K);  dfw = nan(1, K);  Zp = nan(opts.n_perm, K);
nA = zeros(1, K); nB = nA; mA = nan(1, K); mB = mA; dif = mA; g = mA; glo = mA; ghi = mA; dadj = mA;
for j = 1:K
    ok = ~isnan(Y(:, j));
    if ~isempty(Z), ok = ok & all(~isnan(Z), 2); end
    y = Y(ok, j);  b = isB(ok);
    nA(j) = nnz(~b);  nB(j) = nnz(b);
    if nA(j) < 2 || nB(j) < 2 || numel(y) < size(Z, 2) + 4, continue; end   % >= 2 residual df
    mA(j) = mean(y(~b));  mB(j) = mean(y(b));
    if var(y) == 0, continue; end                 % constant: nothing to test
    ord = restrict_permutations(P, ok);
    if isempty(Z)
        dif(j) = mB(j) - mA(j);
        [t_obs(j), dfw(j)] = welch_t(y, b);
        [Tj, dfj] = welch_t(y(ord), b);
    else
        Zj = Z(ok, :);
        X0 = [ones(numel(y), 1), Zj];
        X1 = [ones(numel(y), 1), double(b), Zj];
        if rank(X1) < size(X1, 2), continue; end  % a covariate constant or collinear in these subjects
        c = [0; 1; zeros(size(Zj, 2), 1)];
        res = y - X1 * (X1 \ y);
        if sum(res .^ 2) <= 1e-12 * sum((y - mean(y)) .^ 2)   % a perfect fit (e.g. y among its own covariates):
            continue;                             % residuals of rounding size, G meaningless
        end
        [t_obs(j), dif(j), dfw(j)] = welch_G(y, X1, c, b);
        dadj(j) = dif(j) / sqrt(sum(res .^ 2) / (numel(y) - size(X1, 2)));
        fit0 = X0 * (X0 \ y);
        r0 = y - fit0;
        Ystar = fit0 + r0(ord);                  % Freedman-Lane
        [Tj, ~, dfj] = welch_G(Ystar, X1, c, b);
    end
    Zp(:, j) = t_to_z(Tj(:), dfj(:));            % each permutation on its own degrees of freedom (PALM)
    [g(j), glo(j), ghi(j)] = hedges_g_boot(y(~b), y(b), opts.n_boot, opts.seed + j);
end
t_obs(~isfinite(t_obs)) = NaN;
dfw(isnan(t_obs)) = NaN;
z_obs = t_to_z(t_obs, dfw);
valid = ~isnan(z_obs);
p_unc = nan(1, K);  p_max = nan(1, K);
maxZ = max(abs(Zp(:, valid)), [], 2);
for j = find(valid)
    thr = abs(z_obs(j)) * (1 - 1e-10);        % a permutation equal to the data counts despite rounding
    p_unc(j) = (1 + sum(abs(Zp(:, j)) >= thr)) / (opts.n_perm + 1);
    p_max(j) = (1 + sum(maxZ >= thr)) / (opts.n_perm + 1);
end
R = table((1:K)', nA', nB', mA', mB', dif', t_obs', dfw', p_unc', p_max', g', glo', ghi', dadj', ...
    'VariableNames', {'state', 'nA', 'nB', 'meanA', 'meanB', 'diff', 't', 'df', 'p_unc', ...
                      'p_maxT', 'g', 'g_ci_low', 'g_ci_high', 'd_adjusted'});
if nargout > 1
    Zp(:, ~valid) = NaN;
    df = inf(1, K);  df(~valid) = NaN;
    D = struct('t_obs', z_obs, 'T', Zp, 'df', df);
end
end

function ord = restrict_permutations(P, ok)
% P: n x nPerm permutations of 1..n; returns the permutations of the
% subjects in OK induced by P: each kept subject goes where P sends it,
% following P further when it lands on a left-out subject (cycle walking)
if all(ok), ord = P; return; end
n = size(P, 1);
keep = find(ok);
local = zeros(n, 1);  local(keep) = 1:numel(keep);
J = P(keep, :);
col = repmat(1:size(P, 2), numel(keep), 1);
bad = ~ok(J);
while any(bad(:))
    J(bad) = P(sub2ind(size(P), J(bad), col(bad)));
    bad = ~ok(J);
end
ord = local(J);
end

function [t, df] = welch_t(Y, b)
% Y: n x P (columns = permutations or one outcome), b: n x 1 logical; t and
% its Welch-Satterthwaite degrees of freedom (P x 1 each)
a = ~b;
ma = mean(Y(a, :), 1);  mb = mean(Y(b, :), 1);
ua = var(Y(a, :), 0, 1) / nnz(a);  ub = var(Y(b, :), 0, 1) / nnz(b);
t = (mb - ma) ./ sqrt(ua + ub);
df = (ua + ub) .^ 2 ./ (ua .^ 2 / (nnz(a) - 1) + ub .^ 2 / (nnz(b) - 1));
t = t(:);  df = df(:);
end

function [g, lo, hi] = hedges_g_boot(a, b, nboot, seed)
% Hedges' g of b - a (hedges_g) and its percentile bootstrap CI
g = hg(a, b);
lo = NaN;  hi = NaN;
if nboot > 0
    rs = RandStream('mrg32k3a', 'Seed', seed);
    gb = zeros(nboot, 1);
    for i = 1:nboot
        gb(i) = hg(a(randi(rs, numel(a), numel(a), 1)), b(randi(rs, numel(b), numel(b), 1)));
    end
    gb = gb(isfinite(gb));
    if ~isempty(gb)
        q = prctile(gb, [2.5, 97.5]);
        lo = q(1);  hi = q(2);
    end
end
end

function g = hg(a, b)
g = hedges_g(b, a);                           % b - a: the effect of the group isB
end
