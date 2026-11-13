# Validation scripts

The simulations behind the tables of `docs/validation.md`, run against the
code of the repository they sit in (`validation_setup` puts the repository
root and `4- Analysis/fun` at the top of the MATLAB path). Each script is a
MATLAB function whose help says which table it reproduces, its default
sizes (those of the documented numbers), its seeds, its run time and what
it prints. From any folder:

    addpath('<repository>/tools/validation')
    unequal_groups_table        % the documented sizes
    unequal_groups_table(50)    % a quick check: 50 simulations

The first argument of every script is the number of simulations
(`palm_comparison`, which has one dataset: the number of permutations).
The scripts print their results; nothing is written to disk.

Needs MATLAB with the Statistics and Machine Learning Toolbox;
`cluster_test_null` also needs the toolboxes of `setup_toolboxes`
(FieldTrip for the channel neighbours). The external references are
optional arguments: PALM for `palm_comparison`, PERMUTOOLS for check A of
`family_maxT_checks`; without them those parts are skipped with a note.
`correlation_test_table` and `adjusted_g_coverage` use `parfor` on a
parallel pool that is already open (`parpool`); they never open one, and
without one (or without the Parallel Computing Toolbox) they run serially,
with the same results.

| Script | Reproduces (`docs/validation.md`) | Defaults | Run time (one core) |
|---|---|---|---|
| `unequal_groups_table.m` | "Unequal groups": the rows "z on its df (the pipeline)" of the table | 8,000 datasets, 1,999 permutations | about 1 min |
| `correlation_test_table.m` | "Correlations: the studentised slope": the table (null rejection, n 20–120) and the power figures after it | 2,000 datasets per row, 499 permutations or draws | about 3 min |
| `correlation_test_check.m` | the same section: the permutations of r against the studentised slope at n = 40, with and without a heteroscedastic feature and a covariate (no table) | 2,000 datasets per case, 499 permutations | under 1 min |
| `r_permutation_test_alternative.m` | the column "r, permutations" of that table: the permutation test of r, an alternative the pipeline does not use (called by the two scripts above) | | |
| `adjusted_g_coverage.m` | "Effect size with covariates": the table | 1,000 datasets per case, 1,000 bootstrap resamples | about 1 min |
| `family_maxT_checks.m` | "Multiple testing: the max-statistic family test": the checks A (against PERMUTOOLS), B1–B3 (calibration under the null) and C (the Simes–FDR rule against the max-statistic rule) | A 5 datasets, 19,999 permutations; B1, B2 1,000 studies, 499 permutations; B3 400 runs, 199 permutations; C 20 datasets per scenario, 9,999 permutations | about 25 min (C alone about 23) |
| `palm_comparison.m` | "Where the methods come from": G, df and p of `perm_ancova_maxT` against PALM ("On the same data PALM and the pipeline agree ...") | 20,000 permutations | about 5 min |
| `maxT_unequal_groups_null.m` | "Where the methods come from": the rows "z ... (the pipeline)" of the table of family-wise error with unequal groups | 2,000 datasets, 499 permutations | under 1 min |
| `cluster_test_null.m` | "Where the methods come from": the cluster test with FieldTrip's rules under the null | 2,000 datasets, 1,000 permutations | about 3 min |
| `validation_setup.m` | (helper: the paths, an open pool) | | |

The scripts run the pipeline's code. Rows of the documented tables that
describe alternatives the pipeline does not use (the raw t or G, sign
flips) are not reproduced, except the permutation test of r and the
Simes–FDR rule, which the scripts implement for the comparison. With the
defaults the scripts give the documented numbers.
