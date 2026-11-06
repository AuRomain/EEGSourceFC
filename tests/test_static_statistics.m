function tests = test_static_statistics
%TEST_STATIC_STATISTICS The static analyses on synthetic matrices: features, reliability, enrichment, NBS.
%   From the pipeline folder:  results = runtests('tests/test_static_statistics.m')
%   Matrices of the 68 Desikan-Killiany regions made here (no stage 1):
%   20 recordings, the second group with stronger connectivity inside the
%   default mode network (DMN); the features are checked against their
%   definitions and the edge tests must find the planted network.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
S = setup_toolboxes(static_test_settings());
S.nbs = struct('thresh', 3.1, 'perms', 500, 'size', 'extent');
S.enrichment = struct('n_perm', 2000, 'seed', 1);   % p >= 1 / 2001: q of 28 blocks can be below 0.05
A = load(fullfile(root, '1- EEG source connectivity', 'inputs', 'scout_Desikan-Killiany_68.mat'));
labels = {A.Scouts.Label};
nodes = load_network_table(S, labels);
tc.TestData.S = S;
tc.TestData.nodes = nodes;
tc.TestData.dmn = network_assignment(nodes) == "DMN";
tc.TestData.D = synthetic(labels, tc.TestData.dmn, 0.08, 7);
tc.TestData.D0 = synthetic(labels, tc.TestData.dmn, 0, 8);
tc.TestData.Gatlas = atlas_regions(S, labels, nodes);
tc.TestData.out = tempname;
end

function teardownOnce(tc)
if isfolder(tc.TestData.out), rmdir(tc.TestData.out, 's'); end
end

%% Features
function test_feature_definitions(tc)
D = tc.TestData.D;  nodes = tc.TestData.nodes;
F = static_features(D, 1, nodes, {'global', 'networks', 'regions', 'power'});
n = 68;  r = 3;
W = squeeze(D.conn(1, :, :, r));
Wp = max(W, 0);                                             % graph metrics: the positive part
low = tril(true(n), -1);
v = @(f) F.value(F.id == D.ids(r) & F.feature == f);
verifyEqual(tc, v("Strength"), sum(Wp(:)) / n, 'AbsTol', 1e-12);          % mean node strength
G = gt_window_metrics(W);
verifyEqual(tc, v("ClustCoeff"), G.clustering, 'AbsTol', 1e-12);
verifyEqual(tc, v("CPL"), G.path_length, 'AbsTol', 1e-12);
dmn = tc.TestData.dmn;
B = W(dmn, dmn);
verifyEqual(tc, v("within_DMN"), mean(B(tril(true(nnz(dmn)), -1))), 'AbsTol', 1e-12);
lab = network_assignment(nodes);
B = W(lab == "VIS", lab == "DMN");
verifyEqual(tc, v("between_VIS_DMN"), mean(B(:)), 'AbsTol', 1e-12);
verifyEqual(tc, v("strength_" + regexprep(D.roi_labels{5}, '[^A-Za-z0-9]+', '_')), sum(Wp(5, :)), 'AbsTol', 1e-12);
verifyEqual(tc, v("logpower_all"), mean(log10(D.power(:, r))), 'AbsTol', 1e-12);
verifyEqual(tc, v("logpower_DMN"), mean(log10(D.power(dmn, r))), 'AbsTol', 1e-12);
fam = F(F.id == D.ids(1), :);
verifyEqual(tc, nnz(fam.family == "global"), 3);
verifyEqual(tc, nnz(fam.family == "networks"), 7 + 21);      % Yeo-7: 7 within, 21 between
verifyEqual(tc, nnz(fam.family == "regions"), 68);
verifyEqual(tc, nnz(fam.family == "power"), 8);
verifyTrue(tc, all(F.state == 0));
F1 = static_features(D, 1, nodes, {'global'}, 1);           % a half: the odd blocks
W1 = max(squeeze(D.halves(1, 1, :, :, r)), 0);
verifyEqual(tc, F1.value(F1.id == D.ids(r) & F1.feature == "Strength"), sum(W1(:)) / n, 'AbsTol', 1e-12);
verifyEqual(tc, v("Strength"), mean(Wp(low)) * (n - 1), 'AbsTol', 1e-12);   % = (n - 1) x mean connectivity
F0 = static_features(D, 1, table(), {'global', 'networks', 'regions', 'power'});   % no network table
verifyFalse(tc, any(F0.family == "networks"));
verifyEqual(tc, nnz(F0.family == "power" & F0.id == D.ids(1)), 1);
end

function test_graph_metrics_equal_bct(tc)
% the definitions of the Brain Connectivity Toolbox (not needed by the
% static analyses; the test runs where BCT is installed, in S.toolboxes/BCT)
bct = fullfile(tc.TestData.S.toolboxes, 'BCT', '2019_03_03_BCT');
assumeTrue(tc, isfolder(bct), 'BCT not installed: graph metrics not compared with it');
added = isempty(which('clustering_coef_wu'));
if added
    addpath(bct);
    c = onCleanup(@() rmpath(bct));
end
rs = RandStream('mrg32k3a', 'Seed', 9);
n = 12;  W = zeros(n, n, 3);
for k = 1:3
    A = rand(rs, n);  A = (A + A') / 2;  A(1:n + 1:end) = 0;  W(:, :, k) = A;
end
W(1, 2, 3) = -0.2;  W(2, 1, 3) = -0.2;                      % a negative weight: its positive part (0)
G = gt_window_metrics(W);
for k = 1:3
    A = max(W(:, :, k), 0);
    verifyEqual(tc, G.strength(k), mean(strengths_und(A)), 'AbsTol', 1e-12);
    verifyEqual(tc, G.clustering(k), mean(clustering_coef_wu(A)), 'AbsTol', 1e-12);
    verifyEqual(tc, G.path_length(k), charpath(distance_wei(weight_conversion(A, 'lengths')), 0, 0), 'AbsTol', 1e-12);
end
end

%% Reliability
function test_reliability(tc)
% halves = the matrix plus independent noise: the split-half r of a
% feature is var(signal) / (var(signal) + var(noise)), and identical halves
% give 1
D = tc.TestData.D;  nodes = tc.TestData.nodes;
R = static_reliability(D, 1, nodes, {'global', 'networks'});
verifyEqual(tc, R.Properties.VariableNames, {'family', 'feature', 'n', 'r', 'r_sb', 'note'});
w = R(R.feature == "within_DMN", :);
verifyGreaterThan(tc, w.r, 0.8);
verifyEqual(tc, w.r_sb, 2 * w.r / (1 + w.r), 'AbsTol', 1e-12);
verifyEqual(tc, R.r(R.feature == "identification"), 1);      % every recording has its own pattern
D.halves(2, :, :, :, :) = D.halves(1, :, :, :, :);
R = static_reliability(D, 1, nodes, {'global'});
verifyEqual(tc, R.r(R.feature == "Strength"), 1, 'AbsTol', 1e-12);
verifyEqual(tc, R.r_sb(R.feature == "edge r_sb, median"), 1, 'AbsTol', 1e-12);
D.halves(2, :, :, :, :) = 1 - D.halves(1, :, :, :, :);      % halves that disagree: r < 0, r_sb = r
R = static_reliability(D, 1, nodes, {'global'});
r = R.r(R.feature == "Strength");
verifyLessThan(tc, r, -0.99);
verifyEqual(tc, R.r_sb(R.feature == "Strength"), r);         % 2r / (1 + r) gave NaN or values below -1
verifyLessThan(tc, R.r_sb(R.feature == "edge r_sb, median"), 0);
end

%% Edges: enrichment and NBS
function test_enrichment_finds_the_planted_network(tc)
S = tc.TestData.S;  D = tc.TestData.D;
n = 68;  W = reshape(D.conn(1, :, :, :), n, n, []);
grp = D.group == 2;
X = [ones(numel(grp), 1), double(grp)];
[P, stat] = network_enrichment(W, X, grp, tc.TestData.nodes, tc.TestData.Gatlas, S.enrichment);
verifyEqual(tc, stat, stat.', 'AbsTol', 1e-12);
d = P(P.network_a == "DMN" & P.network_b == "DMN", :);
verifyGreaterThan(tc, d.z_random, 3);
verifyLessThan(tc, d.q_random, 0.05);
verifyLessThan(tc, d.p_spin, 0.01);                          % rotations of 20 of 68 regions overlap them:
verifyEqual(tc, max(P.z_spin), d.z_spin);                    % a wider null than random labels
verifyEqual(tc, d.n_edges, nnz(tc.TestData.dmn) * (nnz(tc.TestData.dmn) - 1) / 2);
verifyEqual(tc, height(P), 28);
other = P(~(P.network_a == "DMN" & P.network_b == "DMN"), :);
verifyFalse(tc, any(other.q_random < 0.05 & other.z_random > 0 & other.network_a ~= "DMN" & other.network_b ~= "DMN"));
D0 = tc.TestData.D0;                                         % no difference: nothing stands out
W0 = reshape(D0.conn(1, :, :, :), n, n, []);
P0 = network_enrichment(W0, X, grp, tc.TestData.nodes, tc.TestData.Gatlas, S.enrichment);
verifyLessThanOrEqual(tc, nnz(P0.q_spin < 0.05), 1);
end

function test_nbs_finds_the_planted_network(tc)
S = tc.TestData.S;  D = tc.TestData.D;
T = table(D.ids(:), string(tc.TestData.S.groups(D.group))', 'VariableNames', {'ID', 'group'});
out = fullfile(tc.TestData.out, 'planted');
[R, E] = compute_NBS_analysis(S, D, 1, T, tc.TestData.nodes, tc.TestData.Gatlas, out, 'full');
verifyTrue(tc, isfile(fullfile(out, 'NBS_full.csv')) && isfile(fullfile(out, 'NBS_full.mat')) && ...
    isfile(fullfile(out, 'enrichment_full.csv')) && isfile(fullfile(out, 'edges_patient_control_full.csv')));
up = R(R.contrast == "patient > control", :);
verifyNotEmpty(tc, up);
[~, i] = max(up.n_edges);
verifyLessThan(tc, up.p_both(i), 0.05);
verifyEqual(tc, up.p_both(i), min(1, 2 * up.p(i)), 'AbsTol', 1e-12);
regs = strtrim(split(up.regions(i), ';'));
verifyGreaterThan(tc, mean(ismember(regs, string(D.roi_labels(tc.TestData.dmn)))), 0.9);
down = R(R.contrast == "control > patient", :);
verifyTrue(tc, isempty(down) || all(down.p_both > 0.05));
verifyTrue(tc, any(E.network_a == "DMN" & E.network_b == "DMN" & E.q_random < 0.05));
Te = readtable(fullfile(out, 'edges_patient_control_full.csv'), 'TextType', 'string');
verifyEqual(tc, height(Te), 68 * 67 / 2);
verifyEqual(tc, Te.difference, Te.mean_A - Te.mean_B, 'AbsTol', 1e-12);
D0 = tc.TestData.D0;                                         % no difference
R0 = compute_NBS_analysis(S, D0, 1, T, tc.TestData.nodes, tc.TestData.Gatlas, fullfile(tc.TestData.out, 'null'), 'full');
verifyTrue(tc, isempty(R0) || all(R0.p_both > 0.05));
end

function test_nbs_figures(tc)
% with S.figures: the edge statistics and every component with p_both < alpha
S = tc.TestData.S;  D = tc.TestData.D;
S.figures = true;
S.nbs.perms = 200;
T = table(D.ids(:), string(S.groups(D.group))', 'VariableNames', {'ID', 'group'});
out = fullfile(tc.TestData.out, 'figures');
R = compute_NBS_analysis(S, D, 1, T, tc.TestData.nodes, tc.TestData.Gatlas, out, 'full');
verifyTrue(tc, isfile(fullfile(out, 'patient_control_full_edges.png')));
verifyNotEmpty(tc, dir(fullfile(out, 'patient_control_full_NBS_1_component_*.png')));
verifyTrue(tc, any(R.p_both < S.alpha));
end

%% Synthetic data
%% Small or unusual datasets (the NIRP and Geneva test datasets, 4 November 2026)
function test_temporal_channels_with_old_names(tc)
% T3 and T4 (kept by the 10-05 cap under these names) are temporal channels:
% the muscle index of a 27-channel cap labelled T3/T4 was NaN
opt = tc.TestData.S.iaf;
f = 1:0.25:45;  ff = 2:0.25:40;
L = {'Fp1', 'Fpz', 'Fp2', 'F3', 'Fz', 'F4', 'T3', 'T4', 'O1', 'Oz', 'O2', 'Pz'};
nc = numel(L);
q = struct('offset', ones(nc, 1), 'exponent', ones(nc, 1), 'alpha_cf', nan(nc, 1), 'alpha_pw', nan(nc, 1), ...
    'r2', ones(nc, 1));
iaf = struct('value', NaN, 'paf', NaN, 'cog', NaN, 'irasa', NaN, 'specparam', NaN, 'source', "", 'flag', "no peak");
R = struct('labels', {L}, 'iaf', iaf, 'freq', f, 'pow', single(repmat(10 ./ f, nc, 1)), 'freq_fit', ff, ...
    'aperiodic', single(repmat(10 ./ ff, nc, 1)), 'fooof', q);
[M, sets] = spectral_profile_summary(R, opt);
verifyEqual(tc, sets.temp, {'T3', 'T4'});
s = 10 ./ f;
verifyEqual(tc, M.hf_temporal, log10(sum(s(f >= 30 & f <= 45)) / sum(s(f >= 1 & f <= 45))), 'AbsTol', 1e-6);
end

function test_groups_from_the_start_of_the_ids(tc)
% S.group_column = '': the table gets a column group from the IDs, so that
% the comparisons can name the groups (the Geneva table has no group
% column: every family stopped on "Column group is not in the subject table")
d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
f = fullfile(d, 'subjects.csv');
writetable(table(["C047"; "P120"; "X001"], [35; 42; 30], 'VariableNames', {'Code', 'age'}), f);
S = tc.TestData.S;
S.subs_info = f;  S.id_column = 'Code';  S.group_column = '';  S.groups = {'C', 'P'};
S.spectral_profile = false;  S.iaf_used = false;
T = read_subject_table(S);
verifyEqual(tc, T.group(1:2), ["C"; "P"]);
verifyTrue(tc, ismissing(T.group(3)));
verifyEqual(tc, select_subjects(T, {'group', '==', 'P'}), [false; true; false]);
writetable(table(["C047"; "P120"], ["a"; "b"], 'VariableNames', {'Code', 'group'}), f);   % a column group of its own
verifyWarning(tc, @() read_subject_table(S), 'read_subject_table:group');
w = warning('off', 'read_subject_table:group');  cw = onCleanup(@() warning(w));
T = read_subject_table(S);
verifyEqual(tc, T.group, ["a"; "b"]);
end

function test_sample_sizes_are_reported(tc)
% 2 + 2 recordings (the NIRP test dataset): with a covariate the
% comparison and the correlations cannot be tested, and nothing said so
% (rows of NaN, no correlation rows); a warning now names them. Without the
% covariate both can be tested. Recordings whose stage 1 run_all reused
% (status "skipped (up to date)") count: the first version counted only
% "ok", so statistics run after a reused stage 1 saw no recordings.
d = tempname;  mkdir(fullfile(d, 'Connectivity'));  c = onCleanup(@() rmdir(d, 's'));
writetable(table(["s1_rest.set"; "s2_rest.set"; "s3_rest.set"; "s4_rest.set"; "s5_rest.set"], ...
    ["ok"; "ok"; "skipped (up to date)"; "skipped (up to date)"; "failed"], 'VariableNames', {'file', 'status'}), ...
    fullfile(d, 'Connectivity', 'stage1_summary.csv'));                      % s3, s4: stage 1 reused
f = fullfile(d, 'subjects.csv');
writetable(table(["s1"; "s2"; "s3"; "s4"; "s5"], ["HV"; "HV"; "SZ"; "SZ"; "SZ"], [25; 27; 30; 37; 40], ...
    [15; 16; 13; 13; 12], 'VariableNames', {'ID', 'grp', 'Age', 'EDU'}), f);
S = tc.TestData.S;
S.results = d;  S.subs_info = f;  S.id_column = 'ID';  S.id_regexp = '^([^_]+)';
S.group_column = 'grp';  S.groups = {'HV', 'SZ'};  S.spectral_profile = false;  S.iaf_used = false;
S.comparisons = struct('name', 'SZ - HV', 'A', {{'grp', '==', 'SZ'}}, 'B', {{'grp', '==', 'HV'}});
S.covariates = {'Age'};
S.correlations = struct('name', {'Education', 'Symptoms'}, 'columns', {'^EDU$', '^PANSS'}, 'subjects', {{}, {}}, ...
    'covariates', {{'Age'}, {}}, 'zero_as_missing', {false, false});
S.covariate_sets = struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {});
verifyWarning(tc, @() check_sample_sizes(S), 'run_statistics:sample_size');
w = warning('off', 'run_statistics:sample_size');  cw = onCleanup(@() warning(w));
[~, r] = evalc('check_sample_sizes(S)');
verifyEqual(tc, r.name', ["SZ - HV", "Education", "Symptoms"]);
verifyEqual(tc, r.n', ["2 + 2", "4", "0"]);                % s5 failed in stage 1: not counted
verifyEqual(tc, r.testable', [false, false, false]);
S.covariates = {};
S.correlations(1).covariates = {};
[~, r] = evalc('check_sample_sizes(S)');
verifyEqual(tc, r.testable', [true, true, false]);
verifyEqual(tc, r.needed(1:2)', [4, 4]);
end

function test_channel_sets_of_a_cap_with_other_names(tc)
% A cap whose channels have no 10-05 names (the Geneva ANT cap: Z1 ...
% RR14; EGI: E1 ... E256): every posterior, frontal, frontopolar and
% temporal measure of the spectral profile was NaN, and the IAF could not be
% found. The sets are now matched by position (cap_channels); a cap with
% 10-20 names keeps its channels and S.iaf.
T = get_ftElec('ASA_10-05_343');
keep = {'Fp1', 'Fpz', 'Fp2', 'F3', 'Fz', 'F4', 'FT7', 'T7', 'T8', 'FT8', 'C3', 'Cz', 'C4', 'P3', 'Pz', 'P4', ...
    'PO7', 'PO8', 'O1', 'Oz', 'O2'};
[~, i] = ismember(keep, T.label);
names = cellstr("E" + (1:numel(keep)));
elec = struct('label', {names'}, 'unit', 'm', 'chanpos', T.chanpos(i, :) + [0.003 0 0]);   % 3 mm off
e = @(lab) names(cellfun(@(x) find(strcmp(keep, x)), lab));                % the E name of 10-05 names
% by name: as before (any case, the spelling of the labels), the cap not used
verifyEqual(tc, cap_channels({'O1', 'Oz', 'O2', 'PO3'}, {'o1', 'Oz', 'X1'}, elec), {'o1', 'Oz'});
verifyEqual(tc, cap_channels({'O1', 'Oz'}, {'X1', 'X2'}, []), cell(1, 0));
% by position: each name's nearest channel within 20 mm, each channel once
[ch, how] = cap_channels({'O1', 'Oz', 'O2', 'T7', 'T3'}, names, elec);
verifyEqual(tc, ch, e({'O1', 'Oz', 'O2', 'T7'}));                         % T3: at T7, already taken
o1 = e({'O1'});
verifyTrue(tc, contains(how, [o1{1} ' for O1 (3 mm)']));
verifyEmpty(tc, cap_channels({'O1', 'Oz'}, names, elec, 2));               % farther than 2 mm
% the spectral profile's sets on that cap
opt = tc.TestData.S.iaf;
f = 1:0.25:45;  ff = 2:0.25:40;  nc = numel(names);
q = struct('offset', ones(nc, 1), 'exponent', ones(nc, 1), 'alpha_cf', nan(nc, 1), 'alpha_pw', nan(nc, 1), ...
    'r2', ones(nc, 1));
iaf = struct('value', NaN, 'paf', NaN, 'cog', NaN, 'irasa', NaN, 'specparam', NaN, 'source', "", 'flag', "no peak");
R = struct('labels', {names}, 'iaf', iaf, 'freq', f, 'pow', single(repmat(10 ./ f, nc, 1)), 'freq_fit', ff, ...
    'aperiodic', single(repmat(10 ./ ff, nc, 1)), 'fooof', q);
[~, sets0] = spectral_profile_summary(R, opt);                             % without the cap: by name, none
verifyEmpty(tc, [sets0.post, sets0.front, sets0.fp, sets0.temp]);
[M, sets] = spectral_profile_summary(R, opt, elec);
verifyEqual(tc, sort(sets.post), sort(e({'O1', 'Oz', 'O2', 'PO7', 'PO8', 'P3', 'Pz', 'P4'})));
verifyEqual(tc, sort(sets.front), sort(e({'F3', 'Fz', 'F4'})));
verifyEqual(tc, sort(sets.fp), sort(e({'Fp1', 'Fpz', 'Fp2'})));
verifyEqual(tc, sort(sets.temp), sort(e({'FT7', 'FT8', 'T7', 'T8'})));
verifyTrue(tc, all(isfinite([M.log_alpha_post, M.hf_temporal, M.eye_frontopolar, M.alpha_front_post])));
% S.iaf.channels on a cap file with other names: matched by position, once
d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
f1 = fullfile(d, 'cap_E21.csv');
writetable(table(string(names'), elec.chanpos(:, 1), elec.chanpos(:, 2), elec.chanpos(:, 3), ...
    'VariableNames', {'label', 'x', 'y', 'z'}), f1);
S = tc.TestData.S;  S.montages = {f1};
[~, S2] = evalc('setup_toolboxes(S)');
verifyEqual(tc, S2.iaf.channels(1:numel(S.iaf.channels)), cellstr(S.iaf.channels(:)'));   % the names kept
verifyEqual(tc, sort(S2.iaf.channels(numel(S.iaf.channels) + 1:end)), ...
    sort(e({'O1', 'Oz', 'O2', 'PO7', 'PO8', 'P3', 'Pz', 'P4'})));
[~, S3] = evalc('setup_toolboxes(S2)');
verifyEqual(tc, S3.iaf.channels, S2.iaf.channels);                          % once
% a cap file with 10-20 names: S.iaf unchanged (the keys of stage 1 too)
f2 = fullfile(d, 'cap_1020.csv');
writetable(table(string(keep'), elec.chanpos(:, 1), elec.chanpos(:, 2), elec.chanpos(:, 3), ...
    'VariableNames', {'label', 'x', 'y', 'z'}), f2);
S.montages = {f2};
[~, S4] = evalc('setup_toolboxes(S)');
verifyEqual(tc, S4.iaf, S.iaf);
end

function D = synthetic(labels, dmn, effect, seed)
% 10 + 10 recordings; connectivity = a common pattern + a recording's own
% pattern (so recordings differ and are identifiable) + edge noise; the
% second group's DMN edges raised by EFFECT; halves = the matrix plus
% independent noise each
rs = RandStream('mrg32k3a', 'Seed', seed);
n = numel(labels);  nr = 20;
sym = @(A) triu(A, 1) + triu(A, 1).';
base = sym(0.3 + 0.05 * randn(rs, n));
conn = zeros(1, n, n, nr);  halves = zeros(2, 1, n, n, nr);
for r = 1:nr
    own = sym(0.04 * randn(rs, n)) + 0.01 * randn(rs);       % pattern and overall level
    W = base + own;
    if r > nr / 2, W(dmn, dmn) = W(dmn, dmn) + effect; end
    W(1:n + 1:end) = 0;
    conn(1, :, :, r) = W;
    for h = 1:2
        H = W + sym(0.01 * randn(rs, n));
        H(1:n + 1:end) = 0;
        halves(h, 1, :, :, r) = H;
    end
end
power = 10 .^ (1 + 0.2 * randn(rs, n, nr));
D = struct('ids', "sub-" + compose("%02d", (1:nr)'), 'group', [ones(nr / 2, 1); 2 * ones(nr / 2, 1)], ...
    'measures', {{'wpli'}}, 'conn', conn, 'halves', halves, 'power', power, ...
    'power_halves', permute(cat(3, power, power .* 1.01), [1 3 2]), 'roi_labels', {labels}, ...
    'band_hz', repmat([7 13], nr, 1));
end
