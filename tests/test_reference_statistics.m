function tests = test_reference_statistics
%TEST_REFERENCE_STATISTICS The statistics of stages 2-3 against verified implementations.
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
%     equals fitlm's coefficient;
%   - perm_partial_corr_maxT: the observed (partial) r, Pearson and
%     Spearman, equals MATLAB's corr / partialcorr; without covariates
%     p_unc and p_maxT agree with PERMUTOOLS' permucorr;
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
%   the path); a test whose reference is missing is skipped (Incomplete),
%   not failed.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
S = static_test_settings();                   % in both pipelines (v_static shares this test)
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

function test_correlation_p_against_permucorr(tc)
pt = reference_folder(tc, 'permucorr.m', 'PERMUTOOLS-master');
tc.applyFixture(matlab.unittest.fixtures.PathFixture(pt, 'IncludingSubfolders', true));
rs = RandStream('mrg32k3a', 'Seed', 25);
n = 36;  x = randn(rs, n, 1);
Y = randn(rs, n, 6) + x .* [0.6 0.3 0 0 0 0];
N = 20000;
R = perm_partial_corr_maxT(Y, x, [], struct('n_perm', N, 'seed', 3, 'type', 'Pearson'));
[r1, p1] = permucorr(repmat(x, 1, 6), Y, 'nperm', N, 'correct', 0, 'seed', 4);
[~, pm] = permucorr(repmat(x, 1, 6), Y, 'nperm', N, 'correct', 1, 'seed', 4);
verifyEqual(tc, R.r, r1', 'AbsTol', 1e-12);
verifyEqual(tc, R.p_unc, p1', 'AbsTol', 0.02);
verifyEqual(tc, R.p_maxT, pm', 'AbsTol', 0.02);
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
assumeNotEmpty(tc, which('duration_rank_summary'), 'no states in this pipeline (v_static)');
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
ks = fullfile(tc.TestData.root, '2- Dynamic connectivity', 'fun', 'Clustering', 'k_stability.m');
assumeTrue(tc, isfile(ks), 'no states in this pipeline (v_static)');
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
