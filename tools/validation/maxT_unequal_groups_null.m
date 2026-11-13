function maxT_unequal_groups_null(n_sim, n_perm)
%MAXT_UNEQUAL_GROUPS_NULL Family-wise error of perm_ancova_maxT's p_maxT under the null with unequal groups and missing values.
%   maxT_unequal_groups_null()               the documented sizes
%   maxT_unequal_groups_null(n_sim)          fewer datasets, for a quick check
%   maxT_unequal_groups_null(n_sim, n_perm)  also fewer permutations
%
%   Reproduces the rows "z ... (the pipeline)" of the table "Family-wise error
%   (p_maxT <= 0.05) | p_unc < 0.05, per column" at the end of the section
%   "Where the methods come from" of docs/validation.md (after the PALM
%   comparison), with the current perm_ancova_maxT of this repository.
%   Data: groups of 10 (SD 3) and 25 (SD 1) recordings, no difference,
%   K = 5 columns, columns 4 and 5 with 8 and 12 missing values (so other
%   subjects and degrees of freedom than columns 1-3), every column raised
%   by 0.02 per year of age (uniform over 20-60 years); without and with
%   age as covariate. Both get the same datasets.
%   Defaults (those of the documented table): n_sim = 2000 datasets,
%   n_perm = 499 permutations each.
%   Seeds: rng(1) before the datasets of each row; perm_ancova_maxT's seed
%   is the dataset's number.
%   Prints, without and with the covariate: the share of datasets with any
%   p_maxT <= 0.05 (family-wise error over the 5 columns) with its 95%
%   Wilson interval, and the share of p_unc <= 0.05 in each column.
%   Run time with the defaults: under 1 minute (one core).
%   The documented table also has rows "raw t" and "raw G", an alternative
%   the pipeline does not use (the raw Welch t or G counted, without the
%   conversion of each permutation's statistic to z on its own degrees of
%   freedom); this script runs the pipeline's code only.
if nargin < 1 || isempty(n_sim), n_sim = 2000; end
if nargin < 2 || isempty(n_perm), n_perm = 499; end
validation_setup();
nA = 10;  nB = 25;  n = nA + nB;  K = 5;
isB = [false(nA, 1); true(nB, 1)];
fprintf('Unequal groups, max-T over %d columns: %d (SD 3) against %d (SD 1), 8 and 12 missing values in columns 4 and 5; %d datasets, %d permutations\n', ...
    K, nA, nB, n_sim, n_perm);
for cov = [false, true]
    rng(1);
    fw = false(n_sim, 1);  pu = false(n_sim, K);
    for s = 1:n_sim
        Y = randn(n, K) .* (1 + 2 * ~isB);            % the smaller group 3 x the SD
        age = 20 + 40 * rand(n, 1);
        Y = Y + 0.02 * age;
        Y(randperm(n, 8), 4) = NaN;  Y(randperm(n, 12), 5) = NaN;
        Z = [];  if cov, Z = age; end
        R = perm_ancova_maxT(Y, isB, Z, struct('n_perm', n_perm, 'seed', s, 'n_boot', 0));
        fw(s) = any(R.p_maxT <= 0.05);
        pu(s, :) = R.p_unc' <= 0.05;
    end
    [lo, hi] = wilson(sum(fw), n_sim);
    if cov, name = 'z, age as covariate'; else, name = 'z, no covariate'; end
    fprintf('%-22s family-wise error %.3f (95%% CI %.3f-%.3f)   p_unc <= 0.05 per column: %s\n', ...
        name, mean(fw), lo, hi, mat2str(mean(pu, 1), 3));
end
end

function [lo, hi] = wilson(k, n)
% Wilson score interval of a binomial proportion k / n (95%)
z = 1.96;  p = k / n;
c = (p + z ^ 2 / (2 * n)) / (1 + z ^ 2 / n);
h = z * sqrt(p * (1 - p) / n + z ^ 2 / (4 * n ^ 2)) / (1 + z ^ 2 / n);
lo = max(0, c - h);  hi = min(1, c + h);
end
