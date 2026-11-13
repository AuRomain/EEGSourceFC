function correlation_test_check(n_sim, n_perm)
%CORRELATION_TEST_CHECK The correlation test against the permutations of r under the null, 40 subjects, with and without a heteroscedastic feature.
%   correlation_test_check()               the sizes of the check
%   correlation_test_check(n_sim)          fewer datasets, for a quick check
%   correlation_test_check(n_sim, n_perm)  also fewer permutations
%
%   A check behind the section "Correlations: the studentised slope"
%   of docs/validation.md (no table of its own; the
%   table is correlation_test_table.m): with 40 subjects and a correlation
%   of 0, the permutations of r (an alternative the pipeline does not use,
%   r_permutation_test_alternative.m) reject too often when the
%   spread of the feature changes with the score (SD 0.2 + 2|x|), the
%   current perm_partial_corr_maxT (the slope over its HC3 standard error,
%   permutations) much less (with the defaults: 0.225 and 0.085 at p <
%   0.05); with a constant spread both are near 0.05 and their p-values
%   nearly the same. Each test without a covariate and with
%   one (z, normal, mixed into both the feature and the score with weight
%   0.5).
%   Defaults: n_sim = 2000 datasets per case, n_perm = 499 permutations.
%   Seeds: rng(7) before the datasets of each case (heteroscedastic, then
%   constant spread); the tests' seed is the dataset's number.
%   Prints, per case, the share of p_unc < 0.05 of the current test and of
%   the permutations of r, without and with the covariate; for the
%   constant spread also the median and largest |p difference| of the two
%   tests without the covariate.
%   Run time with the defaults: under 1 minute (one core).
%   (The values of t_hc3 against statsmodels are a test of the repository:
%   tests/test_reference_statistics.m, test_hc3_t_against_statsmodels_itself.)
if nargin < 1 || isempty(n_sim), n_sim = 2000; end
if nargin < 2 || isempty(n_perm), n_perm = 499; end
validation_setup();
n = 40;
fprintf('Correlation 0, %d subjects, %d datasets per case, %d permutations: rejection at p < 0.05\n', n, n_sim, n_perm);
for het = [true false]
    p_new = nan(n_sim, 1);  p_old = p_new;  p_newc = p_new;  p_oldc = p_new;
    rng(7);
    for i = 1:n_sim
        xs = randn(n, 1);  zs = randn(n, 1);
        if het, y = (0.2 + 2 * abs(xs)) .* randn(n, 1); else, y = randn(n, 1); end
        y2 = y + 0.5 * zs;  x2 = xs + 0.5 * zs;
        oo = struct('n_perm', n_perm, 'seed', i, 'type', 'Pearson');
        a = perm_partial_corr_maxT(y, xs, [], oo);
        b = r_permutation_test_alternative(y, xs, [], oo);
        c = perm_partial_corr_maxT(y2, x2, zs, oo);
        d = r_permutation_test_alternative(y2, x2, zs, oo);
        p_new(i) = a.p_unc;  p_old(i) = b.p_unc;  p_newc(i) = c.p_unc;  p_oldc(i) = d.p_unc;
    end
    if het, name = 'SD 0.2 + 2|x|'; else, name = 'constant SD'; end
    fprintf('%-14s studentised slope %.3f, r %.3f; with a covariate: studentised slope %.3f, r %.3f\n', name, ...
        mean(p_new < 0.05), mean(p_old < 0.05), mean(p_newc < 0.05), mean(p_oldc < 0.05));
    if ~het
        fprintf('%-14s |p difference| of the two tests (no covariate): median %.4f, largest %.4f\n', name, ...
            median(abs(p_new - p_old)), max(abs(p_new - p_old)));
    end
end
end
