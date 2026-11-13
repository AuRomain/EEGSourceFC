function tests = test_reference_statistics
%TEST_REFERENCE_STATISTICS The statistics of the states and the statistics stage against verified implementations.
%   From the pipeline folder:  results = runtests('tests/test_reference_statistics.m')
%   - perm_ancova_maxT without covariates: the Welch t equals MATLAB's
%     ttest2 (Vartype unequal); p_unc and p_maxT agree with PERMUTOOLS'
%     permuttest2 (Crosse et al. 2024; Welch t, max correction) within
%     Monte Carlo error;
%   - perm_ancova_maxT against PALM (Winkler et al. 2014): without and
%     with a covariate, the Welch t / Welch-type G and their degrees of
%     freedom equal PALM's v statistic and -savedof output, and p_unc and
%     p_maxT (z with each permutation's own df, -twotail) agree with PALM's
%     within Monte Carlo error; also for a column with missing values;
%     t_to_z equals PALM's palm_gtoz;
%   - perm_ancova_maxT with covariates: the adjusted group difference
%     equals fitlm's coefficient, and the adjusted g is that coefficient
%     over the pooled SD of the raw values times Hedges' J on nA + nB - 2
%     (hedges_g's g times the ratio of the differences); g_unadjusted is
%     hedges_g of the raw values;
%   - perm_partial_corr_maxT: the observed (partial) r, Pearson and
%     Spearman, equals MATLAB's corr / partialcorr and, without
%     covariates, PERMUTOOLS' permucorr; t_hc3 equals statsmodels'
%     OLS(...).fit(cov_type='HC3') t (values computed with statsmodels
%     0.14.5, and statsmodels itself when python3 has it); p_unc, p_maxT
%     and the permutation statistics equal those of the permutations redone
%     one at a time with the textbook HC3 sandwich of the full model;
%     with the wild bootstrap (S.correlation_test 'wild_bootstrap'), the
%     draws redone one at a time, and p against wildboottest's
%     WildboottestHC (Rademacher, null imposed, type 11) within Monte Carlo
%     error;
%   - hedges_g equals PERMUTOOLS' booteffectsize (Cohen's d, pooled SD,
%     bias-corrected), and perm_ancova_maxT's g uses it;
%   - the duration checks (duration_rank_summary) with one stratum: the
%     rank test's p agrees with MATLAB's ranksum (exact) within Monte Carlo
%     error, Cliff's delta equals PERMUTOOLS' booteffectsize 'cliff';
%   - the adjusted Rand index of k_stability (the choice of the number of
%     states) equals rand_index(..., 'adjusted') (C. McComb's function, in
%     ClusteringPD; Hubert & Arabie 1985).
%   PERMUTOOLS, PALM (github.com/andersonwinkler/PALM, commit 7c00e18)
%   and ClusteringPD are looked for in <toolboxes>/GitHub resources (or on
%   the path), statsmodels and wildboottest in
%   <toolboxes>/python-references (a Python virtual environment; else
%   python3); a test whose reference is missing is skipped (Incomplete),
%   not failed.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
S = static_test_settings();                   % with or without the states
S.parallel = false;
tc.TestData.S = setup_toolboxes(S);
tc.TestData.root = root;
end

function test_welch_t_equals_ttest2(tc)
rs = RandStream('mrg32k3a', 'Seed', 21);
n = 45;  isB = [false(20, 1); true(25, 1)];
Y = randn(rs, n, 6) .* [ones(20, 6); 2 * ones(25, 6)];      % unequal variances
Y(isB, 2) = Y(isB, 2) + 1;
R = perm_ancova_maxT(Y, isB, [], struct('n_perm', 199, 'seed', 1, 'n_boot', 0));
[~, ~, ~, st] = ttest2(Y(isB, :), Y(~isB, :), 'Vartype', 'unequal');
verifyEqual(tc, R.t, st.tstat', 'RelTol', 1e-10);
verifyEqual(tc, R.diff, (mean(Y(isB, :)) - mean(Y(~isB, :)))', 'AbsTol', 1e-12);
end

function test_permutation_p_against_permuttest2(tc)
pt = reference_folder(tc, 'permuttest2.m', 'PERMUTOOLS-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pt, 'IncludingSubfolders', true));
rs = RandStream('mrg32k3a', 'Seed', 22);
isB = [false(18, 1); true(24, 1)];
Y = randn(rs, 42, 8);
Y(isB, 1:2) = Y(isB, 1:2) + [0.9 0.5];
N = 20000;
R = perm_ancova_maxT(Y, isB, [], struct('n_perm', N, 'seed', 3, 'n_boot', 0));
[t1, p1] = permuttest2(Y(isB, :), Y(~isB, :), 'vartype', 'unequal', 'nperm', N, 'correct', 0, 'seed', 4);
[~, pm] = permuttest2(Y(isB, :), Y(~isB, :), 'vartype', 'unequal', 'nperm', N, 'correct', 1, 'seed', 4);
verifyEqual(tc, R.t, t1', 'RelTol', 1e-10);
% two Monte Carlo estimates of the same p: SD of their difference <= 0.005 at N = 20000
% (permuttest2 takes the maximum over the raw Welch t, perm_ancova_maxT over
% their z, as PALM: the same here, every column having the same subjects)
verifyEqual(tc, R.p_unc, p1', 'AbsTol', 0.02);
verifyEqual(tc, R.p_maxT, pm', 'AbsTol', 0.02);
end

function test_welch_G_against_palm(tc)
% G and its Aspin-Welch df equal PALM's v statistic and -savedof output (to
% the 4 decimals of PALM's CSV files); p_unc and p_maxT agree with PALM's
% uncp and fwep (-twotail: |z| with each permutation's own df, the maximum
% over the columns) within Monte Carlo error; without covariates (Welch t,
% Welch-Satterthwaite df) and with one (Freedman-Lane, PALM's default)
pf = reference_folder(tc, 'palm.m', 'PALM-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pf));
rs = RandStream('mrg32k3a', 'Seed', 24);
nA = 12;  nB = 18;  n = nA + nB;
isB = [false(nA, 1); true(nB, 1)];
age = 20 + 40 * rand(rs, n, 1);
Y = randn(rs, n, 5) .* (1 + 1.5 * isB) + 0.03 * age;      % group B 2.5 x the SD
Y(:, 1) = Y(:, 1) + 1.2 * isB;
Y(:, 2) = 4 * Y(:, 2);                                     % a column on another scale
N = 10000;
for cov = [false, true]
    X = [ones(n, 1), double(isB)];  Z = [];
    if cov, X = [X, age]; Z = age; end %#ok<AGROW>
    R = perm_ancova_maxT(Y, isB, Z, struct('n_perm', N, 'seed', 5, 'n_boot', 0));
    P = run_palm(tc, Y, X, isB, N);
    verifyEqual(tc, R.t, P.G, 'AbsTol', 1e-4);
    verifyEqual(tc, R.df, P.df, 'AbsTol', 1e-4);
    verifyEqual(tc, R.p_unc, P.p_unc, 'AbsTol', mc_tol(P.p_unc, N));
    verifyEqual(tc, R.p_maxT, P.p_fwe, 'AbsTol', mc_tol(P.p_fwe, N));
end
% a column with missing values: its permutations restricted to its subjects
% (cycle walking) against PALM on those subjects
Y([1:4, nA + (1:6)], 3) = NaN;
R = perm_ancova_maxT(Y, isB, age, struct('n_perm', N, 'seed', 5, 'n_boot', 0));
ok = ~isnan(Y(:, 3));
P = run_palm(tc, Y(ok, 3), [ones(nnz(ok), 1), double(isB(ok)), age(ok)], isB(ok), N);
verifyEqual(tc, R.t(3), P.G, 'AbsTol', 1e-4);
verifyEqual(tc, R.df(3), P.df, 'AbsTol', 1e-4);
verifyEqual(tc, R.p_unc(3), P.p_unc, 'AbsTol', mc_tol(P.p_unc, N));
end

function test_t_to_z_equals_palm_gtoz(tc)
pf = reference_folder(tc, 'palm_gtoz.m', 'PALM-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pf));
[t, df] = ndgrid([-12 -5 -2.1 -0.3 0 0.3 1.96 4 9 30], [2.5 7 19.3 60 400]);
verifyEqual(tc, t_to_z(t, df), palm_gtoz(t, 1, df), 'AbsTol', 1e-10);
verifyEqual(tc, t_to_z([-1.5 2 NaN], Inf), [-1.5 2 NaN], 'AbsTol', 1e-12);
end

function test_ancova_difference_equals_fitlm(tc)
rs = RandStream('mrg32k3a', 'Seed', 23);
n = 50;  isB = [false(22, 1); true(28, 1)];
age = 20 + 30 * rand(rs, n, 1);  sex = double(rand(rs, n, 1) > 0.5);
Y = 0.05 * age + 0.3 * sex + randn(rs, n, 3);
Y(isB, 1) = Y(isB, 1) + 0.8;
Y(3, 2) = NaN;                                              % a missing value in one column
R = perm_ancova_maxT(Y, isB, [age sex], struct('n_perm', 199, 'seed', 1, 'n_boot', 0));
for j = 1:3
    ok = ~isnan(Y(:, j));
    m = fitlm([double(isB(ok)), age(ok), sex(ok)], Y(ok, j));
    verifyEqual(tc, R.diff(j), m.Coefficients.Estimate(2), 'AbsTol', 1e-10);
    verifyEqual(tc, R.nA(j) + R.nB(j), nnz(ok));
end
end

function test_partial_r_equals_partialcorr(tc)
rs = RandStream('mrg32k3a', 'Seed', 24);
n = 40;
age = 20 + 30 * rand(rs, n, 1);  sex = double(rand(rs, n, 1) > 0.5);
x = 0.04 * age + randn(rs, n, 1);
Y = 0.03 * age + 0.5 * x .* [1 0 0.3] + randn(rs, n, 3);
Y(5, 3) = NaN;  x(9) = NaN;
for type = ["Pearson", "Spearman"]
    o = struct('n_perm', 99, 'seed', 1, 'type', char(type));
    R0 = perm_partial_corr_maxT(Y, x, [], o);
    R1 = perm_partial_corr_maxT(Y, x, [age sex], o);
    for j = 1:3
        ok = ~isnan(Y(:, j)) & ~isnan(x);
        verifyEqual(tc, R0.r(j), corr(Y(ok, j), x(ok), 'Type', char(type)), 'AbsTol', 1e-12);
        verifyEqual(tc, R1.r(j), partialcorr(Y(ok, j), x(ok), [age(ok) sex(ok)], 'Type', char(type)), ...
            'AbsTol', 1e-12, sprintf('%s partial r, column %d', type, j));
    end
end
end

function test_correlation_r_equals_permucorr(tc)
% r equals permucorr's; the p-values are those of the studentised slope
% (permucorr permutes r): checked below against the permutations redone
% one at a time
pt = reference_folder(tc, 'permucorr.m', 'PERMUTOOLS-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pt, 'IncludingSubfolders', true));
rs = RandStream('mrg32k3a', 'Seed', 25);
n = 36;  x = randn(rs, n, 1);
Y = randn(rs, n, 6) + x .* [0.6 0.3 0 0 0 0];
R = perm_partial_corr_maxT(Y, x, [], struct('n_perm', 99, 'seed', 3, 'type', 'Pearson'));
r1 = permucorr(repmat(x, 1, 6), Y, 'nperm', 99, 'correct', 0, 'seed', 4);
verifyEqual(tc, R.r, r1', 'AbsTol', 1e-12);
end

function test_hc3_t_equals_statsmodels(tc)
% statsmodels 0.14.5: OLS(y, X).fit(cov_type='HC3').tvalues[1] with X =
% [1 x] and [1 x age sex]; Spearman: the same on the ranks (scipy's
% rankdata, ties averaged as tiedrank)
x = [-0.8 -1.32 -0.25 0.42 1.14 0.11 -0.55 -0.78 0.75 1.63 0.27 -1.23 -0.96 1.6]';
age = [56 54 36 40 47 22 42 31 55 23 47 55 29 56]';
sex = [1 0 1 0 1 0 0 0 1 0 0 1 0 1]';
y = [-0.84 1.53 0.37 1.39 1.93 -0.11 0.84 0.34 2.96 -0.94 1.94 -0.52 0.07 0.72]';
o = struct('n_perm', 99, 'seed', 1, 'type', 'Pearson');
R0 = perm_partial_corr_maxT(y, x, [], o);
R1 = perm_partial_corr_maxT(y, x, [age sex], o);
o.type = 'Spearman';
R2 = perm_partial_corr_maxT(y, x, [age sex], o);
verifyEqual(tc, [R0.t_hc3, R1.t_hc3, R2.t_hc3], [0.5995348595349611, 0.9095807434206368, 0.7651971660525754], ...
    'AbsTol', 1e-10);
end

function test_hc3_t_against_statsmodels_itself(tc)
% three columns (one with a missing value), heteroscedastic, with and
% without covariates, against statsmodels run here
py = reference_python(tc, 'statsmodels.api');
tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
rs = RandStream('mrg32k3a', 'Seed', 41);
n = 40;  age = 20 + 30 * rand(rs, n, 1);  sex = double(rand(rs, n, 1) > 0.5);
x = 0.04 * age + randn(rs, n, 1);
Y = 0.03 * age + 0.4 * x .* [1 0 0.3] + (0.3 + abs(x)) .* randn(rs, n, 3);
Y(5, 3) = NaN;
o = struct('n_perm', 99, 'seed', 1, 'type', 'Pearson');
R0 = perm_partial_corr_maxT(Y, x, [], o);
R1 = perm_partial_corr_maxT(Y, x, [age sex], o);
t = [R0.t_hc3, R1.t_hc3];
writematrix([Y x age sex], fullfile(tmp, 'data.csv'));
fid = fopen(fullfile(tmp, 'hc3.py'), 'w');
fprintf(fid, '%s\n', 'import sys, numpy as np, statsmodels.api as sm', ...
    'd = np.loadtxt(sys.argv[1], delimiter=",")', ...
    'Y, x, age, sex = d[:, :3], d[:, 3], d[:, 4], d[:, 5]', ...
    'for k in range(2):', ...
    '    for j in range(3):', ...
    '        ok = ~np.isnan(Y[:, j])', ...
    '        X = np.column_stack([np.ones(ok.sum()), x[ok]] + ([age[ok], sex[ok]] if k else []))', ...
    '        print(repr(float(sm.OLS(Y[ok, j], X).fit(cov_type="HC3").tvalues[1])))');
fclose(fid);
[st, out] = system(sprintf('"%s" -I "%s" "%s"', py, fullfile(tmp, 'hc3.py'), fullfile(tmp, 'data.csv')));
verifyEqual(tc, st, 0, out);
ref = reshape(str2double(splitlines(strtrim(string(out)))), 3, 2);
verifyEqual(tc, t, ref, 'AbsTol', 1e-10);
end

function test_correlation_permutations_brute_force(tc)
% perm_partial_corr_maxT's permutation statistics, p_unc and p_maxT against
% the permutations redone one at a time: the score's fitted values on
% [1 z] plus its residuals permuted (the column's subjects following the
% permutation's cycles), the full model [1 x* z] refitted, its HC3
% sandwich (X'X)^-1 X' diag(e.^2 ./ (1 - h).^2) X (X'X)^-1, t to z on N - 3
rs = RandStream('mrg32k3a', 'Seed', 43);
n = 25;  z = randn(rs, n, 1);  x = 0.5 * z + randn(rs, n, 1);
Y = 0.3 * x + (0.5 + abs(x)) .* randn(rs, n, 3);
Y(4, 2) = NaN;  Y([7 19], 3) = NaN;
o = struct('n_perm', 300, 'seed', 9, 'type', 'Pearson');
[R, D] = perm_partial_corr_maxT(Y, x, z, o);
rp = RandStream('mrg32k3a', 'Seed', o.seed);
[~, P] = sort(rand(rp, n, o.n_perm), 1);
Zb = nan(o.n_perm, 3);  z0 = nan(1, 3);
for j = 1:3
    keep = find(~isnan(Y(:, j)));
    m = numel(keep);
    y = Y(keep, j);  Z1 = [ones(m, 1), z(keep)];
    fx = Z1 * (Z1 \ x(keep));  rx = x(keep) - fx;
    z0(j) = t_to_z(sandwich_t(y, [ones(m, 1), x(keep), z(keep)]), m - 3);
    for p = 1:o.n_perm
        ord = zeros(m, 1);
        for k = 1:m                                 % where the permutation sends keep(k), skipping left-out subjects
            to = P(keep(k), p);
            while ~ismember(to, keep), to = P(to, p); end
            ord(k) = find(keep == to);
        end
        Zb(p, j) = t_to_z(sandwich_t(y, [ones(m, 1), fx + rx(ord), z(keep)]), m - 3);
    end
end
verifyEqual(tc, D.t_obs, z0, 'AbsTol', 1e-9);
verifyEqual(tc, D.T, Zb, 'AbsTol', 1e-8);
mx = max(abs(Zb), [], 2);
for j = 1:3
    verifyEqual(tc, R.p_unc(j), (1 + sum(abs(Zb(:, j)) >= abs(z0(j)) * (1 - 1e-10))) / (o.n_perm + 1), 'AbsTol', 1e-12);
    verifyEqual(tc, R.p_maxT(j), (1 + sum(mx >= abs(z0(j)) * (1 - 1e-10))) / (o.n_perm + 1), 'AbsTol', 1e-12);
end
end

function test_wild_bootstrap_brute_force(tc)
% the wild bootstrap's draws, p_unc and p_maxT against the draws redone one
% at a time: y* = the fit of y on [1 z] + its residuals times the draw's
% Rademacher signs (the same random numbers: rand > 0.5), the full model
% [1 x z] refitted, the textbook HC3 sandwich, t to z on N - 3
rs = RandStream('mrg32k3a', 'Seed', 44);
n = 22;  z = randn(rs, n, 1);  x = 0.5 * z + randn(rs, n, 1);
Y = 0.3 * x + (0.5 + abs(x)) .* randn(rs, n, 3);
Y(4, 2) = NaN;  Y([7 19], 3) = NaN;
o = struct('n_perm', 300, 'seed', 9, 'type', 'Pearson', 'method', 'wild_bootstrap');
[R, D] = perm_partial_corr_maxT(Y, x, z, o);
rv = RandStream('mrg32k3a', 'Seed', o.seed);
V = 2 * (rand(rv, n, o.n_perm) > 0.5) - 1;
Zb = nan(o.n_perm, 3);  z0 = nan(1, 3);
for j = 1:3
    keep = ~isnan(Y(:, j));
    m = nnz(keep);
    y = Y(keep, j);  Z1 = [ones(m, 1), z(keep)];  X = [ones(m, 1), x(keep), z(keep)];
    f0 = Z1 * (Z1 \ y);  u0 = y - f0;
    z0(j) = t_to_z(sandwich_t(y, X), m - 3);
    for b = 1:o.n_perm
        Zb(b, j) = t_to_z(sandwich_t(f0 + u0 .* V(keep, b), X), m - 3);
    end
end
verifyEqual(tc, D.t_obs, z0, 'AbsTol', 1e-9);
verifyEqual(tc, D.T, Zb, 'AbsTol', 1e-8);
mx = max(abs(Zb), [], 2);
for j = 1:3
    verifyEqual(tc, R.p_unc(j), (1 + sum(abs(Zb(:, j)) >= abs(z0(j)) * (1 - 1e-10))) / (o.n_perm + 1), 'AbsTol', 1e-12);
    verifyEqual(tc, R.p_maxT(j), (1 + sum(mx >= abs(z0(j)) * (1 - 1e-10))) / (o.n_perm + 1), 'AbsTol', 1e-12);
end
end

function test_wild_bootstrap_against_wildboottest(tc)
% p against wildboottest 0.3.2's WildboottestHC (Rademacher, null imposed,
% type 11: HC1 statistic, residuals unchanged) on three heteroscedastic
% columns with a covariate, 20,000 draws each: within Monte Carlo error;
% with .se 'HC0' (the p of an HC1 statistic is that of HC0, a constant
% factor apart), whose t times sqrt((N - k) / (N - 1)) is wildboottest's t
py = reference_python(tc, 'wildboottest.wildboottest');
tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
rs = RandStream('mrg32k3a', 'Seed', 45);
n = 30;  age = 20 + 30 * rand(rs, n, 1);  x = 0.03 * age + randn(rs, n, 1);
Y = 0.02 * age + [0.5 0.25 0] .* x + (0.3 + abs(x)) .* randn(rs, n, 3);
N = 20000;
R = perm_partial_corr_maxT(Y, x, age, struct('n_perm', N, 'seed', 5, 'type', 'Pearson', ...
    'method', 'wild_bootstrap', 'se', 'HC0'));
writematrix([Y x age], fullfile(tmp, 'data.csv'));
fid = fopen(fullfile(tmp, 'wb.py'), 'w');
fprintf(fid, '%s\n', 'import sys, numpy as np', 'from wildboottest.wildboottest import WildboottestHC', ...
    'd = np.loadtxt(sys.argv[1], delimiter=",")', ...
    'X = np.column_stack([np.ones(d.shape[0]), d[:, 3], d[:, 4]])', ...
    'for j in range(3):', ...
    '    wb = WildboottestHC(X=X, Y=d[:, j], R=np.array([0, 1, 0]), r=0, B=int(sys.argv[2]), seed=11 + j)', ...
    '    wb.get_adjustments(bootstrap_type="11")', ...
    '    wb.get_uhat(impose_null=True)', ...
    '    wb.get_tboot(weights_type="rademacher")', ...
    '    wb.get_tstat()', ...
    '    wb.get_pvalue(pval_type="two-tailed")', ...
    '    print(float(np.ravel(wb.pvalue)[0]), float(np.ravel(wb.t_stat)[0]))');
fclose(fid);
[st, out] = system(sprintf('"%s" -I "%s" "%s" %d', py, fullfile(tmp, 'wb.py'), fullfile(tmp, 'data.csv'), N));
verifyEqual(tc, st, 0, out);
ref = reshape(str2double(split(strtrim(string(out)))), 2, 3)';
verifyEqual(tc, R.t_hc3 * sqrt((n - 3) / (n - 1)), ref(:, 2), 'AbsTol', 1e-10);
verifyEqual(tc, R.p_unc, ref(:, 1), 'AbsTol', mc_tol(ref(:, 1), N));
end

function t = sandwich_t(y, X)
% t of the coefficient of column 2, HC3 (MacKinnon & White 1985)
B = inv(X' * X);
e = y - X * (B * X' * y);
h = sum((X * B) .* X, 2);
V = B * X' * diag(e .^ 2 ./ (1 - h) .^ 2) * X * B;
b = B * X' * y;
t = b(2) / sqrt(V(2, 2));
end

function test_adjusted_g_equals_fitlm(tc)
% with covariates g = J b1 / s_pooled (fitlm's group coefficient over the
% pooled within-group SD of the raw values, the J of hedges_g): hedges_g
% times fitlm's b1 over the raw difference of the means; g_unadjusted =
% hedges_g of the raw values; without covariates g = g_unadjusted; the CIs
% bracket g
rs = RandStream('mrg32k3a', 'Seed', 27);
n = 44;  isB = [false(20, 1); true(24, 1)];
age = 20 + 30 * rand(rs, n, 1);  sex = double(rand(rs, n, 1) > 0.5);
Y = 0.05 * age + 0.3 * sex + randn(rs, n, 2);
Y(isB, 1) = Y(isB, 1) + 0.9;
Y(6, 2) = NaN;
R = perm_ancova_maxT(Y, isB, [age sex], struct('n_perm', 99, 'seed', 1, 'n_boot', 500));
for j = 1:2
    ok = ~isnan(Y(:, j));
    m = fitlm([double(isB(ok)), age(ok), sex(ok)], Y(ok, j));
    gu = hedges_g(Y(ok & isB, j), Y(ok & ~isB, j));
    raw = mean(Y(ok & isB, j)) - mean(Y(ok & ~isB, j));
    verifyEqual(tc, R.g(j), gu * m.Coefficients.Estimate(2) / raw, 'RelTol', 1e-10);
    verifyEqual(tc, R.g_unadjusted(j), gu, 'RelTol', 1e-12);
    verifyTrue(tc, R.g_ci_low(j) < R.g(j) && R.g(j) < R.g_ci_high(j));
    verifyTrue(tc, R.g_unadjusted_ci_low(j) < R.g_unadjusted(j) && R.g_unadjusted(j) < R.g_unadjusted_ci_high(j));
end
R0 = perm_ancova_maxT(Y, isB, [], struct('n_perm', 99, 'seed', 1, 'n_boot', 200));
verifyEqual(tc, [R0.g, R0.g_ci_low, R0.g_ci_high], [R0.g_unadjusted, R0.g_unadjusted_ci_low, R0.g_unadjusted_ci_high]);
end

function test_hedges_g_equals_booteffectsize(tc)
pt = reference_folder(tc, 'booteffectsize.m', 'PERMUTOOLS-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pt, 'IncludingSubfolders', true));
rs = RandStream('mrg32k3a', 'Seed', 26);
for n = [3 5 12 40]
    a = randn(rs, n, 1);  b = 0.7 + 1.5 * randn(rs, n + 2, 1);
    d = booteffectsize(b, a, 'paired', 0, 'effect', 'cohen', 'vartype', 'equal', 'correct', 1);
    verifyEqual(tc, hedges_g(b, a), d, 'RelTol', 1e-12, sprintf('n = %d', n));
end
isB = [false(15, 1); true(17, 1)];
Y = randn(rs, 32, 2);  Y(isB, 1) = Y(isB, 1) + 1;
R = perm_ancova_maxT(Y, isB, [], struct('n_perm', 99, 'seed', 1, 'n_boot', 0));
d = booteffectsize(Y(isB, :), Y(~isB, :), 'paired', 0, 'effect', 'cohen', 'vartype', 'equal', 'correct', 1);
verifyEqual(tc, R.g, d', 'RelTol', 1e-12, 'perm_ancova_maxT''s g: B - A');
end

function test_duration_checks_against_ranksum(tc)
assumeNotEmpty(tc, which('duration_rank_summary'), 'no states in this pipeline');
pt = reference_folder(tc, 'booteffectsize.m', 'PERMUTOOLS-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pt, 'IncludingSubfolders', true));
rs = RandStream('mrg32k3a', 'Seed', 28);
inA = [true(8, 1); false(9, 1)];
for shift = [0 0.6 1.2]
    dur = exp(randn(rs, 17, 1) * 0.3) + shift * inA;
    cov = rand(rs, 17, 1);
    R = duration_rank_summary(dur, cov, inA, [], struct('n_perm', 20000, 'seed', 5, 'outlier_iqr', 3));
    p = ranksum(dur(inA), dur(~inA), 'method', 'exact');
    verifyEqual(tc, R.rank_p, p, 'AbsTol', 0.02, sprintf('shift %g: permutation p against the exact rank sum p', shift));
    d = booteffectsize(dur(inA), dur(~inA), 'paired', 0, 'effect', 'cliff');
    verifyEqual(tc, R.cliff_delta, d, 'AbsTol', 1e-12);
    verifyEqual(tc, R.hl_shift, median(dur(inA) - dur(~inA)', 'all'), 'AbsTol', 1e-12);
end
end

function test_adjusted_rand_equals_rand_index(tc)
ri = reference_folder(tc, 'rand_index.m', 'ClusteringPD-main');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(ri));
% k_stability's adjusted_rand is a local function: copied, unchanged, to a temporary file
ks = fullfile(stage_folder(tc.TestData.root, 'Dynamic connectivity'), 'fun', 'Clustering', 'k_stability.m');
assumeTrue(tc, isfile(ks), 'no states in this pipeline');
src = fileread(ks);
i = strfind(src, 'function a = adjusted_rand(x, y)');
verifyNumElements(tc, i, 1);
rest = src(i:end);
e = strfind(rest, [newline 'function ']);                   % up to the next local function, if any
if ~isempty(e), rest = rest(1:e(1)); end
tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
fid = fopen(fullfile(tmp, 'adjusted_rand.m'), 'w');  fwrite(fid, rest);  fclose(fid);
tc.applyFixture(matlab.unittest.fixtures.PathFixture(tmp));
rs = RandStream('mrg32k3a', 'Seed', 27);
for k = 1:20
    n = 50 + 10 * k;
    x = randi(rs, 2 + mod(k, 5), n, 1);
    y = x;  flip = rand(rs, n, 1) < k / 25;
    y(flip) = randi(rs, 3 + mod(k, 4), nnz(flip), 1);       % from identical to unrelated
    verifyEqual(tc, adjusted_rand(x, y), rand_index(x, y, 'adjusted'), 'AbsTol', 1e-12, sprintf('case %d', k));
end
verifyEqual(tc, adjusted_rand([1 1 2 2 3 3], [5 5 7 7 9 9]), 1, 'AbsTol', 1e-12);   % relabelled: identical
end

function py = reference_python(tc, module)
% the Python of <toolboxes>/python-references (a virtual environment made
% for these tests), else python3, when it imports MODULE
py = fullfile(tc.TestData.S.toolboxes, 'python-references', 'bin', 'python');
if ~isfile(py), py = 'python3'; end
[st, ~] = system(sprintf('"%s" -I -c "import %s"', py, module));
assumeEqual(tc, st, 0, sprintf('Python with %s not found: reference check skipped', module));
end

function folder = reference_folder(tc, file, repo)
% the folder of a reference function: on the path, or under <toolboxes>/GitHub resources/<repo>
w = which(file);
if ~isempty(w), folder = fileparts(w); return; end
folder = '';
base = fullfile(tc.TestData.S.toolboxes, 'GitHub resources', repo);
if isfolder(base)
    d = dir(fullfile(base, '**', file));
    d = d(~contains({d.folder}, {'octave', 'test'}));
    if ~isempty(d), folder = d(1).folder; end
end
assumeNotEmpty(tc, folder, sprintf('%s (%s) not found: reference check skipped', file, repo));
end

function P = run_palm(tc, Y, X, isB, N)
% PALM (the downloaded copy) on CSV files: the contrast of the group column,
% the two groups as variance groups (so its v statistic), -twotail,
% -savedof; PALM's defaults otherwise (Freedman-Lane, Ridgway partitioning)
d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
writematrix(Y, fullfile(d, 'y.csv'));
writematrix(X, fullfile(d, 'x.csv'));
writematrix([0, 1, zeros(1, size(X, 2) - 2)], fullfile(d, 'c.csv'));
writematrix(1 + double(isB), fullfile(d, 'vg.csv'));
o = fullfile(d, 'r'); %#ok<NASGU>
p0 = path;  restore = onCleanup(@() path(p0)); %#ok<NASGU>   % palm adds its own folders
evalc(['palm(''-i'', fullfile(d, ''y.csv''), ''-d'', fullfile(d, ''x.csv''), ''-t'', fullfile(d, ''c.csv''), ', ...
       '''-vg'', fullfile(d, ''vg.csv''), ''-n'', N, ''-twotail'', ''-savedof'', ''-o'', o, ''-seed'', 11, ''-quiet'')']);
rd = @(suffix) readmatrix(fullfile(d, ['r_dat_vstat' suffix]));
P.G = reshape(rd('.csv'), [], 1);
P.p_unc = reshape(rd('_uncp.csv'), [], 1);
P.p_fwe = reshape(rd('_fwep.csv'), [], 1);
if isfile(fullfile(d, 'r_dat_vstat_dof.csv'))           % one df per column
    P.df = reshape(rd('_dof.csv'), [], 1);
else                                                    % one column: df1 and df2 in a text file
    v = readmatrix(fullfile(d, 'r_dat_vstat_dof.txt'), 'FileType', 'text');
    P.df = v(end);
end
end

function tol = mc_tol(p, N)
% 4.5 SD of the difference of two Monte Carlo estimates of p (N permutations
% each), plus PALM's rounding to 4 decimals; p kept 3 permutations away from
% 0 and 1, where the SD of a count of a few permutations is not near 0
p = min(max(p(:), 3 / N), 1 - 3 / N);
tol = 4.5 * sqrt(2 * p .* (1 - p) / N) + 1e-4;
end
