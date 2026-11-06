function [R, D] = perm_partial_corr_maxT(Y, x, Z, opts)
%PERM_PARTIAL_CORR_MAXT Permutation (partial) correlation per column, max-statistic across columns.
%   R = perm_partial_corr_maxT(Y, x, Z, opts)
%   [R, D] = perm_partial_corr_maxT(...) also returns the permutation
%   statistics as t values, for a cluster test over the columns
%   (cluster_perm_test): D.t_obs (1 x K), D.T (n_perm x K, the same
%   permutations in every column), D.df (1 x K: N - 2 - covariates; NaN
%   when not tested), t = r sqrt(df / (1 - r^2))
%   Y     n x K features (one column per state), NaN = missing
%   x     n x 1 clinical or cognitive score
%   Z     n x q covariates for a partial correlation ([] for none)
%   opts  .n_perm, .seed, .type ('Pearson' or 'Spearman')
%
%   Without covariates x is permuted. With covariates the partial
%   correlation r(y, x | Z) is tested with Freedman-Lane: the residuals of
%   x ~ Z are permuted, added back to the fitted values and the partial
%   correlation recomputed (S14). Spearman ranks y, x and the covariates
%   (among the column's subjects) first, the partial Spearman correlation
%   of MATLAB's partialcorr(..., 'Type', 'Spearman') (before 5 November
%   2026 the covariates were not ranked; the pipeline uses Pearson). The
%   observed r equals partialcorr's (tests/test_reference_statistics). Every
%   column uses the same random permutation of the n subjects, restricted
%   to its own subjects by following the permutation's cycles, so p_maxT
%   controls the FWER across the K columns. The maximum is taken over |r|
%   when every column has the same N, and over |z| when the columns' N
%   differ (missing values in some states): the t = r sqrt(df / (1 - r^2))
%   of each column (the t of the score's coefficient in y ~ x + Z, PALM's
%   default statistic) converted to z on its own df (t_to_z, PALM's
%   palm_gtoz; Winkler et al. 2014; Alberton et al. 2020, NeuroImage
%   216:116760), so that a column with fewer subjects, whose |r| and |t|
%   vary more under the null, does not dominate the maximum (the same p as
%   over |r| when N is equal; before 5 November 2026 the maximum was over
%   |t|, closer to one scale than |r| but not on it).
%   p_unc is per column and the same either way. N is reported per column (S8).
%   An undefined correlation (no variance, too few subjects, a covariate
%   constant in these subjects) gives NaN, never a p-value.

[n, K] = size(Y);
x = x(:);
if nargin < 3, Z = []; end
rs = RandStream('mrg32k3a', 'Seed', opts.seed);
[~, P] = sort(rand(rs, n, opts.n_perm), 1);
r_obs = nan(1, K);  Rp = nan(opts.n_perm, K);  N = zeros(1, K);
for j = 1:K
    ok = ~isnan(Y(:, j)) & ~isnan(x);
    if ~isempty(Z), ok = ok & all(~isnan(Z), 2); end
    N(j) = nnz(ok);
    if N(j) < size(Z, 2) + 4, continue; end         % at least 2 residual degrees of freedom
    y = Y(ok, j);  xx = x(ok);
    if var(y) == 0 || var(xx) == 0, continue; end
    if strcmpi(opts.type, 'Spearman'), y = tiedrank(y); xx = tiedrank(xx); end
    ord = restrict_permutations(P, ok);
    if isempty(Z)
        r_obs(j) = corr(y, xx);
        Rp(:, j) = corr(y, xx(ord)).';
    else
        Zj = Z(ok, :);
        if strcmpi(opts.type, 'Spearman'), Zj = tiedrank(Zj); end   % every variable ranked, as partialcorr
        Z1 = [ones(N(j), 1), Zj];
        if rank(Z1) < size(Z1, 2), continue; end
        H = @(v) v - Z1 * (Z1 \ v);              % residual-forming
        ry = H(y);
        fx = Z1 * (Z1 \ xx);  rx = xx - fx;
        if norm(ry) < 1e-12 * norm(y) || norm(rx) < 1e-12 * norm(xx), continue; end
        r_obs(j) = corr(ry, rx);
        Xs = fx + rx(ord);                        % Freedman-Lane
        Rp(:, j) = corr(ry, H(Xs)).';
    end
end
r_obs(~isfinite(r_obs)) = NaN;
valid = ~isnan(r_obs);
df = N - 2 - size(Z, 2);
df(~valid) = NaN;
Rp(:, ~valid) = NaN;
to_t = @(r) r .* sqrt(df ./ max(1 - r .^ 2, eps));
if numel(unique(N(valid))) > 1                % columns of different N: the max over |z| (the same |r| is
    s_obs = abs(t_to_z(to_t(r_obs), df));     % stronger evidence with more subjects; PALM)
    s_perm = abs(t_to_z(to_t(Rp), repmat(df, size(Rp, 1), 1)));
else
    s_obs = abs(r_obs);  s_perm = abs(Rp);    % the same N everywhere: |r| (the same order as |t|)
end
p_unc = nan(1, K);  p_max = nan(1, K);
mx = max(s_perm(:, valid), [], 2);
for j = find(valid)
    thr = abs(r_obs(j)) * (1 - 1e-10);        % a permutation equal to the data counts despite rounding
    p_unc(j) = (1 + sum(abs(Rp(:, j)) >= thr)) / (opts.n_perm + 1);
    p_max(j) = (1 + sum(mx >= s_obs(j) * (1 - 1e-10))) / (opts.n_perm + 1);
end
R = table((1:K)', N', r_obs', p_unc', p_max', ...
    'VariableNames', {'state', 'n', 'r', 'p_unc', 'p_maxT'});
if nargout > 1
    D = struct('t_obs', to_t(r_obs), 'T', to_t(Rp), 'df', df);
end
end

function ord = restrict_permutations(P, ok)
% permutations of the subjects in OK induced by P (cycle walking), as in perm_ancova_maxT
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
