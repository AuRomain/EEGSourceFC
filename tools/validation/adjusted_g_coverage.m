function adjusted_g_coverage(n_sim, n_boot)
%ADJUSTED_G_COVERAGE Bias and CI coverage of the covariate-adjusted g (docs/validation.md, "Effect size with covariates").
%   adjusted_g_coverage()                the documented sizes
%   adjusted_g_coverage(n_sim)           fewer datasets, for a quick check
%   adjusted_g_coverage(n_sim, n_boot)   also fewer bootstrap resamples
%
%   Reproduces the table of the section "Effect size with covariates"
%   of docs/validation.md, with the current
%   perm_ancova_maxT of this repository (its g with covariates: the
%   adjusted difference over the pooled within-group SD of the raw values,
%   times Hedges' J, with a percentile bootstrap CI; g_unadjusted: Hedges'
%   g of the raw values).
%   Data: age as covariate, uniform over 30 years (20-50) plus 5 years in
%   group B, raising the feature by 0.05 per year; a group difference
%   delta; within-group SD 1. Three cases: 20 + 20 with delta 0.5, 10 + 30
%   with delta 0.8, 30 + 30 with delta 0. The target is delta over the
%   within-group SD of the raw values, sqrt(1 + 0.05^2 30^2 / 12) = 1.09.
%   Defaults (those of the documented table): n_sim = 1000 datasets per
%   case, n_boot = 1000 bootstrap resamples (perm_ancova_maxT's n_boot;
%   19 permutations, the p-values are not used).
%   Seeds: dataset i of case c from RandStream('mrg32k3a', 'Seed', 7000 c
%   + i); perm_ancova_maxT's seed is i.
%   Prints one line per case: the target, the mean and SD of the adjusted
%   g, the coverage of the target by its 95% CI and the CI's mean width,
%   and the mean raw g.
%   Uses an open parallel pool (parfor over the datasets; the results do
%   not depend on it); without one it runs serially.
%   Run time with the defaults: about 1 minute on one core.
if nargin < 1 || isempty(n_sim), n_sim = 1000; end
if nargin < 2 || isempty(n_boot), n_boot = 1000; end
[~, n_workers] = validation_setup();
cases = [20 20 0.5; 10 30 0.8; 30 30 0];
fprintf('Adjusted g: age as covariate (group B 5 years older), %d datasets per case, %d bootstrap resamples\n', ...
    n_sim, n_boot);
for c = 1:size(cases, 1)
    nA = cases(c, 1);  nB = cases(c, 2);  delta = cases(c, 3);
    target = delta / sqrt(1 + 0.05 ^ 2 * 30 ^ 2 / 12);   % in units of the within-group SD of the raw values
    covered = false(n_sim, 1);  g = nan(n_sim, 1);  gu = g;  w = g;
    parfor (i = 1:n_sim, n_workers)
        rs = RandStream('mrg32k3a', 'Seed', 7000 * c + i);
        isB = [false(nA, 1); true(nB, 1)];
        age = 20 + 30 * rand(rs, nA + nB, 1) + 5 * isB;   % group B older: age differs between the groups
        y = 0.05 * age + delta * isB + randn(rs, nA + nB, 1);
        R = perm_ancova_maxT(y, isB, age, struct('n_perm', 19, 'seed', i, 'n_boot', n_boot));
        covered(i) = R.g_ci_low <= target && target <= R.g_ci_high;
        g(i) = R.g;  gu(i) = R.g_unadjusted;  w(i) = R.g_ci_high - R.g_ci_low;
    end
    fprintf(['%d + %d, delta %.1f (target %.3f): adjusted g mean %.3f (SD %.3f), 95%% CI covers the target %.3f, ' ...
        'mean width %.2f; raw g mean %.3f\n'], nA, nB, delta, target, mean(g), std(g), mean(covered), mean(w), mean(gu));
end
end
