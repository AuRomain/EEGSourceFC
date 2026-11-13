function [R, D] = perm_partial_corr_maxT(Y, x, Z, opts)
%PERM_PARTIAL_CORR_MAXT Permutation test of a (partial) correlation per column, max-statistic across columns.
%   R = perm_partial_corr_maxT(Y, x, Z, opts)
%   [R, D] = perm_partial_corr_maxT(...) also returns the permutation
%   statistics, for a cluster test over the columns (cluster_perm_test) and
%   the family test: D.t_obs (1 x K) and D.T (n_perm x K, the same
%   permutations in every column) are z values (below), D.df (1 x K) Inf
%   (NaN when not tested), as those of perm_ancova_maxT
%   Y     n x K features (one column per state), NaN = missing
%   x     n x 1 clinical or cognitive score
%   Z     n x q covariates for a partial correlation ([] for none)
%   opts  .n_perm (permutations, or bootstrap draws), .seed, .type
%         ('Pearson' or 'Spearman', case ignored), .method ('permutation',
%         the default, or 'wild_bootstrap'; S.correlation_test), .se ('HC3',
%         the default; 'HC0' only to compare with references)
%
%   The statistic is the slope b of the feature on the score, in the model
%   y = b0 + b x + Z c, over its heteroscedasticity-consistent standard
%   error HC3 (MacKinnon & White 1985, J Econometrics 29:305; Long & Ervin
%   2000, Am Stat 54:217):
%       t_hc3 = b sum(rx.^2) / sqrt(sum(rx.^2 e.^2 ./ (1 - h).^2))
%   with rx the residuals of the score on the covariates (and intercept), e
%   the residuals of the model and h its leverages. A permutation test of r
%   itself tests that feature and score are independent, not that they are
%   uncorrelated: when the spread of the feature changes with the score it
%   rejects more often than alpha although r is 0 (24% at alpha 0.05 with
%   the SD of y 0.2 + 2|x|, 40 subjects; docs/validation.md). The
%   studentised statistic keeps the test exact under independence and
%   valid, asymptotically, when only the correlation is 0 (DiCiccio &
%   Romano 2017, J Am Stat Assoc 112:1211). t_hc3 equals statsmodels'
%   OLS(...).fit(cov_type=
%   'HC3') (tests/test_reference_statistics).
%   Without covariates the score is permuted. With covariates the residuals
%   of the score on the covariates are permuted, added back to its fitted
%   values and the model refitted (the Smith method of Winkler et al.
%   2014, NeuroImage 92:381: the score residualised on the covariates is
%   permuted, the features are not, so that they keep their joint
%   distribution across the columns; Freedman-Lane would permute the
%   residuals of each feature);
%   every permutation's t is studentised in the same way. Spearman
%   ranks y, x and the covariates (among the column's subjects) first, the
%   partial Spearman correlation of MATLAB's partialcorr(..., 'Type',
%   'Spearman'). The reported r is the (partial) correlation, partialcorr's
%   (tests/test_reference_statistics).
%   With .method 'wild_bootstrap' the null distribution comes from a wild
%   bootstrap instead of permutations: the model without the score fitted
%   (y on the covariates and intercept: the null imposed), its residuals
%   multiplied by Rademacher signs, one per subject, the full model refitted
%   and its t_hc3 computed (the restricted wild bootstrap; Davidson &
%   Flachaire 2008, J Econometrics 146:162; MacKinnon 2013), each subject's
%   sign the same in every column so that p_maxT and the family test keep
%   the joint distribution. It keeps the link between the score and the
%   spread of the feature, which a permutation breaks: with 20 to 60
%   subjects and a strongly heteroscedastic feature it is closer to alpha
%   than the permutations (docs/validation.md). Its p equals that of
%   wildboottest's WildboottestHC (Python port of fwildclusterboot and
%   boottest; Rademacher, null imposed, type 11) within Monte Carlo error
%   (tests/test_reference_statistics, with .se 'HC0', whose p is that of
%   its HC1 statistic).
%   Every column uses the same random permutation of the n subjects,
%   restricted to its own subjects by following the permutation's cycles,
%   so p_maxT controls the FWER across the K columns. The t of the data and
%   of every permutation is converted to z on the column's N - 2 - q
%   degrees of freedom (t_to_z, PALM's palm_gtoz; Winkler et al. 2014;
%   Alberton et al. 2020, NeuroImage 216:116760) and p_unc and p_maxT are
%   counted on |z|, so that a column with fewer subjects (a state missing
%   in some recordings), whose t varies more under the null, does not
%   dominate the maximum; with the same N in every column the order is that
%   of |t_hc3|. N is reported per column.
%   An undefined statistic gives NaN, never a p-value: no variance, too few
%   subjects, a covariate constant or collinear in these subjects, a model
%   that fits the feature exactly, or a subject of leverage 1 (fitted
%   exactly whatever its value, e.g. the only one at 1 on a 0/1 score:
%   HC3 is not defined). A subject fitted exactly by the covariates alone
%   (the only one of a category) adds nothing to b and is left out of the
%   standard error.
%
%   R (table, one row per column): state, n, r, t_hc3, p_unc, p_maxT

[n, K] = size(Y);
x = x(:);
if nargin < 3, Z = []; end
q = size(Z, 2);
spearman = strcmpi(opts.type, 'Spearman');
wild = isfield(opts, 'method') && strcmpi(opts.method, 'wild_bootstrap');
se = 'HC3';
if isfield(opts, 'se') && ~isempty(opts.se), se = upper(opts.se); end
rs = RandStream('mrg32k3a', 'Seed', opts.seed);
if wild
    P = 2 * (rand(rs, n, opts.n_perm) > 0.5) - 1;   % Rademacher signs, one per subject and draw
else
    [~, P] = sort(rand(rs, n, opts.n_perm), 1);     % permutations of the subjects
end
r_obs = nan(1, K);  t_obs = nan(1, K);  Tp = nan(opts.n_perm, K);  N = zeros(1, K);
M = struct('ok', {});                            % what the columns of the same subjects share
for j = 1:K
    ok = ~isnan(Y(:, j)) & ~isnan(x);
    if ~isempty(Z), ok = ok & all(~isnan(Z), 2); end
    N(j) = nnz(ok);
    if N(j) < q + 4, continue; end               % at least 2 residual degrees of freedom
    y = Y(ok, j);  xx = x(ok);
    if var(y) == 0 || var(xx) == 0, continue; end
    if isempty(M) || ~isequal(M.ok, ok)
        M = score_model(xx, Z, ok, P, spearman, wild);
    end
    if ~M.defined, continue; end
    if spearman, y = tiedrank(y); end
    ry = M.H(y);
    if norm(ry) < 1e-12 * norm(y), continue; end
    e = ry - M.rx * ((M.rx' * ry) / sum(M.rx .^ 2));
    if sum(e .^ 2) <= 1e-12 * sum((y - mean(y)) .^ 2), continue; end   % an exact fit: t infinite
    r_obs(j) = corr(ry, M.rx);
    t_obs(j) = hc_t(ry, M.rx, M.hz, se);
    if wild
        Tp(:, j) = hc_t_wild(ry, P(ok, :), M, se);
    else
        Tp(:, j) = hc_t(ry, M.RX, M.hz, se);
    end
end
valid = isfinite(t_obs);
r_obs(~valid) = NaN;  t_obs(~valid) = NaN;
df = N - 2 - q;  df(~valid) = NaN;
z_obs = t_to_z(t_obs, df);
Zp = t_to_z(Tp, repmat(df, opts.n_perm, 1));
Zp(:, ~valid) = NaN;
p_unc = nan(1, K);  p_max = nan(1, K);
mx = max(abs(Zp(:, valid)), [], 2);
for j = find(valid)
    thr = abs(z_obs(j)) * (1 - 1e-10);           % a permutation equal to the data counts despite rounding
    p_unc(j) = (1 + sum(abs(Zp(:, j)) >= thr)) / (opts.n_perm + 1);
    p_max(j) = (1 + sum(mx >= thr)) / (opts.n_perm + 1);
end
R = table((1:K)', N', r_obs', t_obs', p_unc', p_max', ...
    'VariableNames', {'state', 'n', 'r', 't_hc3', 'p_unc', 'p_maxT'});
if nargout > 1
    dz = inf(1, K);  dz(~valid) = NaN;
    D = struct('t_obs', z_obs, 'T', Zp, 'df', dz);
end
end

function M = score_model(xx, Z, ok, P, spearman, wild)
% what the columns with the subjects OK share: the residual-forming map of
% the covariates and intercept (H), their leverages (hz), the score's
% residuals (rx) and, for permutations, its permuted residuals (RX, N x
% n_perm: the residuals permuted, added back to the fitted values, the
% covariates removed again; H(fitted) = 0, so H of the permuted
% residuals); defined = false when the statistic is not
M = struct('ok', ok, 'defined', false, 'H', [], 'hz', [], 'rx', [], 'RX', []);
if spearman, xx = tiedrank(xx); end
Z1 = ones(numel(xx), 1);
if ~isempty(Z)
    Zj = Z(ok, :);
    if spearman, Zj = tiedrank(Zj); end         % every variable ranked, as partialcorr
    Z1 = [Z1, Zj];
    if rank(Z1) < size(Z1, 2), return; end      % a covariate constant or collinear in these subjects
end
[Q, ~] = qr(Z1, 0);
M.H = @(v) v - Q * (Q' * v);
M.hz = sum(Q .^ 2, 2);
M.rx = M.H(xx);
if norm(M.rx) < 1e-12 * norm(xx), return; end  % the score is a function of the covariates
h = M.hz + M.rx .^ 2 / sum(M.rx .^ 2);
if any(h(M.hz < 1 - 1e-10) > 1 - 1e-8), return; end   % leverage 1: HC3 undefined
if ~wild, M.RX = M.H(M.rx(restrict_permutations(P, ok))); end
M.defined = true;
end

function t = hc_t(ry, RX, hz, se)
% t of the slope of RY on each column of RX (residuals on the covariates),
% HC3 standard errors (HC0 with se 'HC0'); HZ the covariates' leverages.
% Subjects fitted exactly by the covariates (hz = 1, rx = 0) add nothing
% and are left out.
use = hz < 1 - 1e-10;
ry = ry(use);  hz = hz(use);
t = zeros(size(RX, 2), 1);
for c = 1:5000:size(RX, 2)                       % 5000 permutations at a time (memory)
    i = c:min(c + 4999, size(RX, 2));
    X = RX(use, i);
    sxx = sum(X .^ 2, 1);
    b = (ry' * X) ./ sxx;
    U = ry - X .* b;                             % e
    if strcmp(se, 'HC3'), U = U ./ (1 - hz - X .^ 2 ./ sxx); end   % e / (1 - h)
    t(i) = (b .* sxx ./ sqrt(sum(X .^ 2 .* U .^ 2, 1)))';
end
end

function t = hc_t_wild(ry, V, M, se)
% the t of every wild bootstrap draw: y* = fitted values without the score
% + ry .* v (V: N x draws of signs), the full model refitted; with the
% covariates' residual-forming map H, H(y*) = H(ry .* v), the slope
% rx' (ry .* v) / sum(rx.^2) and the residuals H(ry .* v) - rx b
use = M.hz < 1 - 1e-10;
rx = M.rx(use);
sxx = sum(rx .^ 2);
a = ones(nnz(use), 1);
if strcmp(se, 'HC3'), a = 1 ./ (1 - M.hz(use) - rx .^ 2 / sxx); end   % 1 / (1 - h)
t = zeros(size(V, 2), 1);
for c = 1:5000:size(V, 2)
    i = c:min(c + 4999, size(V, 2));
    HU = M.H(ry .* V(:, i));
    HU = HU(use, :);
    b = (rx' * HU) / sxx;
    U = (HU - rx .* b) .* a;
    t(i) = (b * sxx ./ sqrt(sum(rx .^ 2 .* U .^ 2, 1)))';
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
