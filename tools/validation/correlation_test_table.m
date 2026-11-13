function correlation_test_table(n_sim, n_perm)
%CORRELATION_TEST_TABLE Null rejection and power of the correlation tests (docs/validation.md, "Correlations: the studentised slope").
%   correlation_test_table()               the documented sizes
%   correlation_test_table(n_sim)          fewer datasets, for a quick check
%   correlation_test_table(n_sim, n_perm)  also fewer permutations or draws
%
%   Reproduces the table of the section "Correlations: the studentised
%   slope" of docs/validation.md (rejection at p < 0.05
%   under the null) and the power figures of the paragraph after it (a
%   slope of 0.4 SD with a constant spread), with the current
%   perm_partial_corr_maxT of this repository, three tests on the same
%   datasets:
%     perm   the studentised slope (HC3), permutations (the default)
%     wild   the same statistic, restricted wild bootstrap
%            (opts.method = 'wild_bootstrap', S.correlation_test)
%     old    the permutations of the partial r, an alternative the pipeline
%            does not use (r_permutation_test_alternative.m, in this folder)
%   Data: n = 20, 40, 60, 120 subjects; age (uniform over 20-50 years) and
%   sex as covariates; the score x normal, entered as x + 0.03 age; the
%   feature 0.02 age + 0.3 sex + beta x + SD(x) e, with SD(x) = 0.2 + 2|x|,
%   1 + 0.5|x| or 1 (constant); beta = 0 (null) for every SD, and 0.4
%   (power) with the constant SD.
%   Defaults (those of the documented table): n_sim = 2000 datasets per
%   row, n_perm = 499 permutations or bootstrap draws.
%   Seeds: dataset i of n subjects and SD number s from
%   RandStream('mrg32k3a', 'Seed', 1000 n + 100 s + i) (the same datasets
%   for beta 0 and 0.4); the tests' seed is i.
%   Prints one line per row: n, the SD, beta, and the share of p_unc < 0.05
%   of the three tests (perm, wild, old).
%   Uses an open parallel pool (parfor over the datasets; the results do
%   not depend on it); without one it runs serially.
%   Run time with the defaults: about 3 minutes on one core.
if nargin < 1 || isempty(n_sim), n_sim = 2000; end
if nargin < 2 || isempty(n_perm), n_perm = 499; end
[~, n_workers] = validation_setup();
names = ["0.2+2|x|", "1+0.5|x|", "constant"];
fprintf(['n   SD          effect  perm   wild   old    (rejection at p < 0.05; age and sex as covariates, %d ' ...
    'permutations or draws, %d datasets)\n'], n_perm, n_sim);
for n = [20 40 60 120]
    for s = 1:3
        for beta = [0 0.4]
            if beta > 0 && s ~= 3, continue; end
            pn = nan(n_sim, 1);  po = pn;  pw = pn;
            parfor (i = 1:n_sim, n_workers)
                rr = RandStream('mrg32k3a', 'Seed', 1000 * n + 100 * s + i);
                age = 20 + 30 * rand(rr, n, 1);  sex = double(rand(rr, n, 1) > 0.5);
                x = randn(rr, n, 1);
                y = 0.02 * age + 0.3 * sex + beta * x + feature_sd(s, x) .* randn(rr, n, 1);
                xs = x + 0.03 * age;
                o = struct('n_perm', n_perm, 'seed', i, 'type', 'Pearson');
                a = perm_partial_corr_maxT(y, xs, [age sex], o);
                b = r_permutation_test_alternative(y, xs, [age sex], o);
                o.method = 'wild_bootstrap';
                w = perm_partial_corr_maxT(y, xs, [age sex], o);
                pn(i) = a.p_unc;  po(i) = b.p_unc;  pw(i) = w.p_unc;
            end
            fprintf('%-4d%-12s%-8.1f%-7.3f%-7.3f%-7.3f\n', n, names(s), beta, mean(pn < 0.05), mean(pw < 0.05), ...
                mean(po < 0.05));
        end
    end
end
end

function sd = feature_sd(s, x)
% the SD of the feature given the score: 0.2 + 2|x|, 1 + 0.5|x| or constant
switch s
    case 1, sd = 0.2 + 2 * abs(x);
    case 2, sd = 1 + 0.5 * abs(x);
    otherwise, sd = ones(size(x));
end
end
