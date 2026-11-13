function unequal_groups_table(n_sim, n_perm)
%UNEQUAL_GROUPS_TABLE Null p-values of the group test with very unequal groups (docs/validation.md, "Unequal groups").
%   unequal_groups_table()               the documented sizes
%   unequal_groups_table(n_sim)          fewer datasets, for a quick check
%   unequal_groups_table(n_sim, n_perm)  also fewer permutations
%
%   Reproduces the rows "z on its df (the pipeline)" of the table of the section
%   "Unequal groups" (under "Statistics of stage 4 and of the PSD
%   analysis") of docs/validation.md, with the current perm_ancova_maxT of
%   this repository. Data: 8 against 30 recordings, the smaller group with
%   3 times the SD, no difference between the groups; one test per dataset:
%     Welch t   no covariate (group labels permuted)
%     ANCOVA    age as covariate (uniform over 20-60 years, raising the
%               outcome by 0.02 per year): Welch-type G, Freedman-Lane
%   Both rows get the same datasets (the same random draws in the same
%   order; age is drawn in both).
%   Defaults (those of the documented table): n_sim = 8000 datasets,
%   n_perm = 1999 permutations each.
%   Seeds: rng(2) before the datasets of each row; perm_ancova_maxT's seed
%   is the dataset's number.
%   Prints, per row, the share of datasets with p_unc below 0.05, 0.01,
%   0.005 and 0.001, and the Monte Carlo SE of a share of 0.05.
%   Run time with the defaults: about 1 minute (one core).
%   The documented table also has rows of alternatives the pipeline does
%   not use (the raw Welch t and G counted, without the conversion of each
%   permutation's statistic to z on its own degrees of freedom, and a
%   sign-flip variant); this script runs the pipeline's code only.
if nargin < 1 || isempty(n_sim), n_sim = 8000; end
if nargin < 2 || isempty(n_perm), n_perm = 1999; end
validation_setup();
lv = [0.05 0.01 0.005 0.001];
nA = 8;  nB = 30;  n = nA + nB;
isB = [false(nA, 1); true(nB, 1)];
fprintf('Unequal groups: %d against %d, the smaller group 3 x the SD, no difference; %d datasets, %d permutations\n', ...
    nA, nB, n_sim, n_perm);
for cov = [false, true]
    rng(2);
    pv = nan(n_sim, 1);
    for s = 1:n_sim
        y = randn(n, 1) .* (1 + 2 * ~isB);
        age = 20 + 40 * rand(n, 1);
        Z = [];
        if cov, Z = age;  y = y + 0.02 * age; end
        R = perm_ancova_maxT(y, isB, Z, struct('n_perm', n_perm, 'seed', s, 'n_boot', 0));
        pv(s) = R.p_unc;
    end
    if cov, name = 'ANCOVA (G, age), z on its df'; else, name = 'Welch t, z on its df'; end
    fprintf('%-30s %s\n', name, strjoin(arrayfun(@(a) sprintf('p < %g: %.4f', a, mean(pv < a)), lv, ...
        'UniformOutput', false), '   '));
end
fprintf('Monte Carlo SE of a share of 0.05 with %d datasets: %.4f\n', n_sim, sqrt(0.05 * 0.95 / n_sim));
end
