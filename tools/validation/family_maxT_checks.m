function family_maxT_checks(n_sim, permutools_dir, parts)
%FAMILY_MAXT_CHECKS Checks of the max-statistic family test (docs/validation.md, "Multiple testing: the max-statistic family test").
%   family_maxT_checks()                               the documented sizes, check A skipped
%   family_maxT_checks(n_sim)                          fewer simulations, for a quick check
%   family_maxT_checks(n_sim, permutools_dir)          with check A
%   family_maxT_checks(n_sim, permutools_dir, parts)   only some checks
%   n_sim           null studies of B1 and B2 (default 1000); B3 runs
%                   round(0.4 n_sim) null runs (400) and C max(1,
%                   round(n_sim / 50)) datasets per scenario (20); A has 5
%                   datasets whatever n_sim
%   permutools_dir  the folder of PERMUTOOLS (Crosse et al. 2024,
%                   github.com/mickcrosse/PERMUTOOLS), or a folder above
%                   it: permuttest2.m is looked for inside. Optional:
%                   without it check A is skipped with a note
%   parts           a cell of 'A', 'B1', 'B2', 'B3', 'C' (default all)
%
%   Reproduces the list "Checks (...)" of the section "Multiple testing:
%   the max-statistic family test" of docs/validation.md, with the current
%   perm_ancova_maxT, perm_partial_corr_maxT and family_adjust of this
%   repository; the family's maximum and p_maxT_family are computed as
%   permu_ttest_permu_corr_BNS computes them (its local functions
%   add_to_family_max and family_wise, copied below).
%   A   Against PERMUTOOLS' permuttest2 (Welch t, max correction; complete
%       data, no covariates): 5 datasets of 25 + 35 subjects, 2 features x
%       10 states correlated 0.5, effects of several sizes in 3 columns,
%       19,999 permutations; family p of both, the largest difference of
%       the 20 tests' p and their correlation. Seeds: rng(100 + dataset);
%       ours 7, permuttest2 3. Different random permutations, so they agree
%       within Monte Carlo error.
%   B1  Calibration of group comparisons under the null: 20 against 30
%       subjects, the smaller group twice as variable, age as covariate,
%       10% missing values in one feature, 2 features x 8 states correlated
%       0.5, 499 permutations: share of studies with family p <= 0.05 and
%       <= 0.01. Seeds: rng(1000 + study); permutations seed = study. (The
%       "6.1% and 1.4% with the raw t or G" beside the documented values
%       come from an alternative the pipeline does not use; not run here.)
%   B2  Calibration of partial correlations under the null: N 40, two
%       scores with 10% missing values, age partialled out, 2 features x 8
%       states, 499 permutations. Seeds: rng(5000 + study).
%   B3  The whole adjustment (family_adjust) over 10 null families of 5
%       tests (40 subjects, correlated 0.3), 199 permutations: share of
%       runs with any significant test, with the max-statistic rule and
%       with the Simes-FDR rule (Simes family p, Benjamini-Hochberg over
%       the families, q_bb; local function old_rule).
%       Seeds: rng(20000 + 10 run + family); permutations seed = run.
%   C   The two rules on 15 families of 10 features x 10 states (correlated
%       0.3), 40 + 40 subjects, 9,999 permutations: share of the planted
%       tests found significant and datasets with a false positive, for one
%       strong test (d = 1.5) and for 30 tests of one family (d = 0.6,
%       0.8). Seeds: rng(90000 + 100 scenario + dataset); permutations seed
%       = dataset.
%   Prints the results of each check as Markdown tables and lists.
%   Run time with the defaults: about 25 minutes (one core): A, B1, B2
%   and B3 about 1 minute together, C about 23 minutes.
%   With the defaults it gives the documented values: in A the family p
%   of ours 0.2051, 0.0044, 0.1881, 0.0391, 0.0001 (largest difference
%   0.0084, r >= 0.9999); in B2 0.045 (<= 0.01: 0.009); in B3 0.060 (the
%   Simes-FDR rule: 0).
if nargin < 1 || isempty(n_sim), n_sim = 1000; end
if nargin < 2, permutools_dir = ''; end
if nargin < 3 || isempty(parts), parts = {'A', 'B1', 'B2', 'B3', 'C'}; end
parts = upper(cellstr(parts));
validation_setup();
fprintf('# Max-statistic family test: checks\n\n');

%% A: PERMUTOOLS
if ismember('A', parts)
    fprintf('## A. Cross-check with PERMUTOOLS (permuttest2, Welch t, max correction)\n\n');
    pt = '';
    if ~isempty(permutools_dir), pt = find_file(permutools_dir, 'permuttest2.m'); end
    if isempty(pt)
        fprintf('PERMUTOOLS not given (or permuttest2.m not found in the folder given): check A skipped.\n\n');
    else
        p0 = path;
        restore = onCleanup(@() path(p0));
        addpath(pt);
        if isfolder(fullfile(pt, 'utils')), addpath(fullfile(pt, 'utils')); end
        np = 19999;
        fprintf('| dataset | n_perm | family p (ours) | family p (PERMUTOOLS) | max abs difference over the 20 tests | r of the 20 p-values |\n|---|---|---|---|---|---|\n');
        for ds = 1:5
            rng(100 + ds);
            nA = 25;  nB = 35;
            Y = corr_data(nA + nB, 20, 0.5);
            Y(nA + 1:end, [1 2 11]) = Y(nA + 1:end, [1 2 11]) + [0.5 0.7 0.6] * (ds / 3);   % effects of several sizes
            isB = [false(nA, 1); true(nB, 1)];
            o = struct('n_perm', np, 'seed', 7, 'n_boot', 0);
            fmax = -inf(np, 1);  z = [];
            for f = 1:2
                [~, D] = perm_ancova_maxT(Y(:, (f - 1) * 10 + (1:10)), isB, [], o);
                [fmax, zf] = add_to_family_max(fmax, D);
                z = [z, zf]; %#ok<AGROW>
            end
            ours = family_wise(z, fmax);
            [~, ptp] = permuttest2(Y(isB, :), Y(~isB, :), 'vartype', 'unequal', 'nperm', np, 'correct', 1, 'seed', 3);
            fprintf('| %d | %d | %.4f | %.4f | %.4f | %.4f |\n', ds, np, min(ours), min(ptp), ...
                max(abs(ours(:) - ptp(:))), corr(ours(:), ptp(:)));
        end
        fprintf('\nMonte Carlo SD of a p of 0.05 with %d permutations: %.4f per estimate.\n\n', np, ...
            sqrt(0.05 * 0.95 / (np + 1)));
        clear restore
    end
end

if any(ismember({'B1', 'B2', 'B3'}, parts))
    fprintf('## B. Calibration under the null\n\n');
end
%% B1: null group comparisons
if ismember('B1', parts)
    np = 499;
    fp = nan(n_sim, 1);
    for s = 1:n_sim
        rng(1000 + s);
        nA = 20;  nB = 30;  n = nA + nB;
        isB = [true(nA, 1); false(nB, 1)];               % B = the small group, twice as variable
        age = 30 + 8 * randn(n, 1);
        Y = corr_data(n, 16, 0.5) .* (1 + isB) + 0.05 * (age - 30);
        Y2 = Y(:, 9:16);  Y2(rand(n, 8) < 0.1) = NaN;  Y(:, 9:16) = Y2;   % 10% missing in the second feature
        o = struct('n_perm', np, 'seed', s, 'n_boot', 0);
        fmax = -inf(np, 1);  z = [];
        for f = 1:2
            [~, D] = perm_ancova_maxT(Y(:, (f - 1) * 8 + (1:8)), isB, age, o);
            [fmax, zf] = add_to_family_max(fmax, D);
            z = [z, zf]; %#ok<AGROW>
        end
        fp(s) = min(family_wise(z, fmax));
    end
    fprintf(['- B1 group comparison (20 vs 30, the smaller group twice as variable, age as covariate, 10%% missing ' ...
        'in one feature, 2 features x 8 states correlated 0.5), %d null studies, %d permutations: family p <= 0.05 ' ...
        'in %.3f, <= 0.01 in %.3f (Monte Carlo SD %.3f at 0.05).\n'], ...
        n_sim, np, mean(fp <= 0.05), mean(fp <= 0.01), sqrt(0.05 * 0.95 / n_sim));
end

%% B2: null partial correlations
if ismember('B2', parts)
    np = 499;
    fp = nan(n_sim, 1);
    for s = 1:n_sim
        rng(5000 + s);
        n = 40;
        age = 30 + 8 * randn(n, 1);
        Y = corr_data(n, 16, 0.5) + 0.05 * (age - 30);
        X = [randn(n, 1) + 0.03 * age, randn(n, 1)];
        X(rand(n, 2) < 0.1) = NaN;
        o = struct('n_perm', np, 'seed', s, 'type', 'Pearson');
        fmax = -inf(np, 1);  z = [];
        for f = 1:2
            for k = 1:2
                [~, D] = perm_partial_corr_maxT(Y(:, (f - 1) * 8 + (1:8)), X(:, k), age, o);
                [fmax, zf] = add_to_family_max(fmax, D);
                z = [z, zf]; %#ok<AGROW>
            end
        end
        fp(s) = min(family_wise(z, fmax));
    end
    fprintf(['- B2 partial correlations (N 40, two scores with 10%% missing, age partialled out, 2 features x 8 ' ...
        'states), %d null studies, %d permutations: family p <= 0.05 in %.3f, <= 0.01 in %.3f.\n'], ...
        n_sim, np, mean(fp <= 0.05), mean(fp <= 0.01));
end

%% B3: the whole adjustment over 10 null families
if ismember('B3', parts)
    nrun = max(1, round(0.4 * n_sim));  np = 199;  anysig = 0;  anyold = 0;
    for r = 1:nrun
        tabs = cell(1, 10);  fam = tabs;
        for f = 1:10
            rng(20000 + 10 * r + f);
            n = 40;  isB = [true(20, 1); false(20, 1)];
            Y = corr_data(n, 5, 0.3);
            [R, D] = perm_ancova_maxT(Y, isB, [], struct('n_perm', np, 'seed', r, 'n_boot', 0));
            [fmax, zf] = add_to_family_max(-inf(np, 1), D);
            pw = family_wise(zf, fmax);
            tabs{f} = table(R.p_unc, fdr_bh(R.p_unc), pw(:), repmat(min(pw), 5, 1), 'VariableNames', ...
                {'p_unc', 'q_fdr', 'p_maxT_family', 'family_p'});
            fam{f} = repmat(string(f), 5, 1);
        end
        T = family_adjust(tabs, fam, 0.05);
        anysig = anysig + any(cellfun(@(t) any(t.significant), T));
        anyold = anyold + any(cellfun(@any, old_rule(tabs, 0.05)));
    end
    fprintf(['- B3 the adjustment over families (family_adjust), 10 null families of 5 tests, %d runs, %d ' ...
        'permutations: runs with any significant test %.3f (max-statistic rule), %.3f (Simes rule, ' ...
        'q_bb).\n\n'], nrun, np, anysig / nrun, anyold / nrun);
end

%% C: the two rules, 15 families of 100 tests
if ismember('C', parts)
    nds = max(1, round(n_sim / 50));  np = 9999;
    fprintf('## C. The Simes rule (Simes, FDR within) against the max-statistic rule\n\n');
    fprintf(['15 families of 10 features x 10 states (correlated 0.3), 40 + 40 subjects, %d permutations; %d ' ...
        'datasets per scenario.\n\n'], np, nds);
    fprintf(['| scenario | planted tests significant, before | max-statistic | datasets with a false positive, ' ...
        'before | max-statistic |\n|---|---|---|---|---|\n']);
    scen = {'one strong test (d = 1.5)', 'spread: 30 tests of one family (d = 0.6)', ...
        'spread: 30 tests of one family (d = 0.8)'};
    for sc = 1:numel(scen)
        hit_old = 0;  hit_new = 0;  fp_old = 0;  fp_new = 0;  n_planted = 0;
        for ds = 1:nds
            rng(90000 + 100 * sc + ds);
            n = 80;  isB = [true(40, 1); false(40, 1)];
            tabs = cell(1, 15);  fam = tabs;  planted = tabs;
            for f = 1:15
                fmax = -inf(np, 1);  z = [];  pl = false(100, 1);
                Rs = cell(1, 10);
                for ft = 1:10
                    Y = corr_data(n, 10, 0.3);
                    if f == 1
                        if sc == 1 && ft == 1, Y(isB, 1) = Y(isB, 1) + 1.5;  pl(1) = true; end
                        if sc > 1 && ft <= 3
                            d = 0.6 + 0.2 * (sc == 3);
                            Y(isB, :) = Y(isB, :) + d;  pl((ft - 1) * 10 + (1:10)) = true;
                        end
                    end
                    [R, D] = perm_ancova_maxT(Y, isB, [], struct('n_perm', np, 'seed', ds, 'n_boot', 0));
                    [fmax, zf] = add_to_family_max(fmax, D);
                    z = [z, zf]; %#ok<AGROW>
                    Rs{ft} = R;
                end
                R = vertcat(Rs{:});
                pw = family_wise(z, fmax);
                tabs{f} = table(R.p_unc, fdr_bh(R.p_unc), pw(:), repmat(min(pw), 100, 1), 'VariableNames', ...
                    {'p_unc', 'q_fdr', 'p_maxT_family', 'family_p'});
                fam{f} = repmat(string(f), 100, 1);  planted{f} = pl;
            end
            Tn = family_adjust(tabs, fam, 0.05);
            pl = vertcat(planted{:});
            sn = cell2mat(cellfun(@(t) t.significant, Tn(:), 'UniformOutput', false));
            so = cell2mat(old_rule(tabs, 0.05));
            n_planted = n_planted + nnz(pl);
            hit_new = hit_new + nnz(sn & pl);  hit_old = hit_old + nnz(so & pl);
            fp_new = fp_new + any(sn & ~pl);  fp_old = fp_old + any(so & ~pl);
        end
        fprintf('| %s | %.2f | %.2f | %d / %d | %d / %d |\n', scen{sc}, hit_old / n_planted, hit_new / n_planted, ...
            fp_old, nds, fp_new, nds);
    end
    fprintf('\n');
end
end

function Y = corr_data(n, k, rho)
% n x k normal values, every pair of columns correlated rho
Y = sqrt(rho) * randn(n, 1) + sqrt(1 - rho) * randn(n, k);
end

function [fmax, z_obs] = add_to_family_max(fmax, D)
% as add_to_family_max in permu_ttest_permu_corr_BNS: the tests of one
% feature (D of perm_ancova_maxT or perm_partial_corr_maxT) added to the
% family's largest |z| in each permutation; z_obs their |z| in the data
z_obs = abs(t_to_z(D.t_obs, D.df));
fmax = max(fmax, max(abs(t_to_z(D.T, repmat(D.df, size(D.T, 1), 1))), [], 2));
end

function p = family_wise(z, fmax)
% p_maxT_family of each test, as family_wise in permu_ttest_permu_corr_BNS:
% (1 + the permutations whose largest |z| reaches the test's) / (n_perm + 1)
p = nan(size(z));
for i = find(~isnan(z))
    p(i) = (1 + sum(fmax >= z(i) * (1 - 1e-10))) / (numel(fmax) + 1);
end
end

function sig = old_rule(tabs, alpha)
% the Simes-FDR rule (what family_adjust applies to tables without
% p_maxT_family), one family per table: the family p is the
% Simes combination of its p_unc, the families are selected by
% Benjamini-Hochberg over their family p (F tested, R selected), and a test
% of a selected family is significant when q_bb = min(1, q_fdr F / R) <=
% alpha (q_fdr: Benjamini-Hochberg within the family). sig: a column cell,
% one logical column per table
fp = cellfun(@(t) simes(t.p_unc), tabs(:));
fq = fdr_bh(fp);
F = nnz(~isnan(fp));
sel = fq <= alpha;
R = nnz(sel);
sig = cell(numel(tabs), 1);
for i = 1:numel(tabs)
    sig{i} = false(height(tabs{i}), 1);
    if sel(i), sig{i} = min(1, tabs{i}.q_fdr * F / R) <= alpha; end
end
end

function ps = simes(p)
% Simes combination: min over i of m p_(i) / i (NaN when no p-value)
p = sort(p(~isnan(p)));
m = numel(p);
if m == 0, ps = NaN; return; end
ps = min(1, min(m * p(:) ./ (1:m)'));
end

function d = find_file(d, file)
% the folder of FILE: the folder given, or one below it (not a copy for
% Octave or in a test folder)
if isfile(fullfile(d, file)), return; end
f = dir(fullfile(d, '**', file));
f = f(~contains({f.folder}, {'octave', 'test'}));
if isempty(f), d = ''; else, d = f(1).folder; end
end
