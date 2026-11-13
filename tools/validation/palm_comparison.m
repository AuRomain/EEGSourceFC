function palm_comparison(n_perm, palm_dir)
%PALM_COMPARISON The group test of perm_ancova_maxT against PALM on one simulated dataset (docs/validation.md, "Where the methods come from").
%   palm_comparison()                   our side only (PALM not given)
%   palm_comparison(n_perm, palm_dir)   with PALM (Winkler et al. 2014)
%   n_perm    permutations of each test (default 20000, as documented;
%             fewer for a quick check). The script has one dataset, so
%             this takes the place of the number of simulations.
%   palm_dir  the folder of PALM (github.com/andersonwinkler/PALM; the
%             documented check used commit 7c00e18), or a folder above it:
%             palm.m is looked for inside. Optional: without it, or when
%             palm.m is not found there, only the pipeline's values are
%             printed, with a note.
%
%   Reproduces the sentence after the table of the section "Where the
%   methods come from" of docs/validation.md: "On the same data PALM and
%   the pipeline agree: G and df equal, p_maxT of the strongest column
%   0.040 and 0.042 (20,000 permutations each ...)", with the current
%   perm_ancova_maxT of this repository.
%   Data: 12 against 18 recordings, group B 2.5 times the SD, a difference
%   of 1.2 in column 1, column 2 on 4 times the scale, all 5 columns raised
%   by 0.03 per year of age (uniform over 20-60 years), age as covariate:
%   ANCOVA, Welch-type G, Freedman-Lane permutations. Scenario 1: no
%   missing values; scenario 2: 5 and 8 missing values in columns 3 and 4.
%   PALM: its v statistic with the two groups as variance groups (-vg),
%   -twotail, -savedof, its default Freedman-Lane, -seed 11; each column
%   alone for the uncorrected p (PALM's CSV input takes no NaN), and in
%   scenario 1 all columns together for the family-wise p.
%   Seeds: rng(7) for the data; perm_ancova_maxT seed 3; PALM seed 11.
%   Prints, per scenario and column: G (ours) and v (PALM), their degrees
%   of freedom, p_unc of both; in scenario 1 also p_maxT (ours) and PALM's
%   family-wise p over the 5 columns, and PALM's df when all columns are
%   tested together. PALM writes 4 decimals; the p-values come from
%   different random permutations, so they agree within Monte Carlo error.
%   Run time with the defaults: about 5 minutes, almost all in PALM.
if nargin < 1 || isempty(n_perm), n_perm = 20000; end
if nargin < 2, palm_dir = ''; end
validation_setup();
palm_ok = false;
if ~isempty(palm_dir)
    palm_dir = find_palm(palm_dir);
    palm_ok = ~isempty(palm_dir);
end
if ~palm_ok
    fprintf('PALM not given (or palm.m not found in the folder given): the PALM columns are skipped.\n');
end
rng(7);
nA = 12;  nB = 18;  n = nA + nB;  K = 5;
isB = [false(nA, 1); true(nB, 1)];
age = 20 + 40 * rand(n, 1);
Y = randn(n, K) .* (1 + 1.5 * isB) + 0.03 * age;     % group B 2.5 x the SD
Y(:, 1) = Y(:, 1) + 1.2 * isB;
Y(:, 2) = Y(:, 2) * 4;                                % a column on another scale
fprintf('perm_ancova_maxT against PALM: %d against %d, age as covariate, %d permutations each\n', nA, nB, n_perm);
for scenario = 1:2
    Yc = Y;
    if scenario == 2                                  % missing values: other subjects, other df
        Yc([1 2 3 13 14], 3) = NaN;  Yc([4 5 15 16 17 18 19 20], 4) = NaN;
        fprintf('\n== scenario 2 (missing values in columns 3 and 4)\n');
    else
        fprintf('\n== scenario 1 (no missing values)\n');
    end
    R = perm_ancova_maxT(Yc, isB, age, struct('n_perm', n_perm, 'seed', 3, 'n_boot', 0));
    pal = nan(K, 4);                                  % G, p_unc, p_fwe, df of PALM, column by column
    if palm_ok
        for j = 1:K
            ok = ~isnan(Yc(:, j));
            P = run_palm(palm_dir, Yc(ok, j), isB(ok), age(ok), n_perm);
            pal(j, :) = [P.G, P.p_unc, P.p_fwe, P.df];
        end
    end
    fprintf('col   G_ours     v_PALM     df_ours   df_PALM   p_unc ours  p_unc PALM\n');
    for j = 1:K
        fprintf('%d  %9.5f  %9.5f  %8.4f  %8.4f   %.4f      %.4f\n', j, R.t(j), pal(j, 1), R.df(j), pal(j, 4), ...
            R.p_unc(j), pal(j, 2));
    end
    if scenario == 1
        fprintf('p_maxT ours (max |z|):            %s\n', mat2str(R.p_maxT', 4));
        if palm_ok
            P = run_palm(palm_dir, Yc, isB, age, n_perm);
            fprintf('family-wise p PALM (max |z|):     %s\n', mat2str(P.p_fwe', 4));
            fprintf('df PALM, all columns together:    %s\n', mat2str(P.df', 6));
        end
    end
end
end

function d = find_palm(d)
% the folder of palm.m: the folder given, or one below it
if isfile(fullfile(d, 'palm.m')), return; end
f = dir(fullfile(d, '**', 'palm.m'));
if isempty(f), d = ''; else, d = f(1).folder; end
end

function P = run_palm(palm_dir, Y, isB, age, n_perm) %#ok<INUSD>   (n_perm: used inside evalc)
% PALM on CSV files in a temporary folder: the contrast of the group
% column, the two groups as variance groups (its v statistic), -twotail,
% -savedof; PALM's defaults otherwise. The MATLAB path is restored after
% (PALM adds its own folders).
p0 = path;
restore = onCleanup(@() path(p0));
addpath(palm_dir);
d = tempname;  mkdir(d);
remove = onCleanup(@() rmdir(d, 's'));
writematrix(Y, fullfile(d, 'y.csv'));
writematrix([ones(numel(isB), 1), double(isB), age], fullfile(d, 'x.csv'));
writematrix([0 1 0], fullfile(d, 'c.csv'));
writematrix(1 + double(isB), fullfile(d, 'vg.csv'));
o = fullfile(d, 'r'); %#ok<NASGU>   % used inside evalc
evalc(['palm(''-i'', fullfile(d, ''y.csv''), ''-d'', fullfile(d, ''x.csv''), ''-t'', fullfile(d, ''c.csv''), ', ...
       '''-vg'', fullfile(d, ''vg.csv''), ''-n'', n_perm, ''-twotail'', ''-savedof'', ''-o'', o, ''-seed'', 11, ''-quiet'')']);
rd = @(suffix) readmatrix(fullfile(d, ['r_dat_vstat' suffix]));
P.G = reshape(rd('.csv'), [], 1);
P.p_unc = reshape(rd('_uncp.csv'), [], 1);
P.p_fwe = reshape(rd('_fwep.csv'), [], 1);
if isfile(fullfile(d, 'r_dat_vstat_dof.csv'))         % one df per column
    P.df = reshape(rd('_dof.csv'), [], 1);
else                                                  % one column: df1 and df2 in a text file
    v = readmatrix(fullfile(d, 'r_dat_vstat_dof.txt'), 'FileType', 'text');
    P.df = v(end);
end
clear remove restore
end
