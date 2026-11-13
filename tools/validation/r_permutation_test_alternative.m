function [R, D] = r_permutation_test_alternative(Y, x, Z, opts)
%R_PERMUTATION_TEST_ALTERNATIVE Permutations of r: a correlation test that the pipeline does not use.
%   Used only to reproduce the column "r, permutations" of the table of
%   the section "Correlations: the studentised slope" of docs/validation.md
%   (correlation_test_table.m) and the same comparison in
%   correlation_test_check.m. The pipeline does not use it:
%   perm_partial_corr_maxT permutes the slope over its HC3 standard error,
%   because the permutations of r reject too often when the spread of the
%   feature changes with the score.
%   [R, D] = r_permutation_test_alternative(Y, x, Z, opts), with the
%   inputs of perm_partial_corr_maxT (opts .n_perm, .seed, .type 'Pearson'
%   or 'Spearman'; no .method):
%   Y     n x K features (one column per state), NaN = missing
%   x     n x 1 score
%   Z     n x q covariates for a partial correlation ([] for none)
%   Without covariates x is permuted. With covariates the partial
%   correlation r(y, x | Z) is tested by permuting the residuals of x ~ Z,
%   adding them back to the fitted values and recomputing the partial
%   correlation. Spearman ranks y, x and the covariates (among the
%   column's subjects) first. Every column uses the same random
%   permutation of the n subjects, restricted to its own subjects by
%   following the permutation's cycles. The maximum over the columns is
%   over |r| when every column has the same N, over |z| (t = r sqrt(df /
%   (1 - r^2)) on df = N - 2 - q, t_to_z) when the N differ.
%   R (table, one row per column): state, n, r, p_unc, p_maxT
%   D: t_obs, T (n_perm x K) as t values, df (1 x K)

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
        Xs = fx + rx(ord);                        % the score's residuals permuted, added back
        Rp(:, j) = corr(ry, H(Xs)).';
    end
end
r_obs(~isfinite(r_obs)) = NaN;
valid = ~isnan(r_obs);
df = N - 2 - size(Z, 2);
df(~valid) = NaN;
Rp(:, ~valid) = NaN;
to_t = @(r) r .* sqrt(df ./ max(1 - r .^ 2, eps));
if numel(unique(N(valid))) > 1                % columns of different N: the max over |z|
    s_obs = abs(t_to_z(to_t(r_obs), df));
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
