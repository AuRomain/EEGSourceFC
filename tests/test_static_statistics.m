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
mkdir(tc.TestData.out);                                     % for the tests that write into it directly
end

function teardownOnce(tc)
if ~isfield(tc.TestData, 'out'), return; end                % setupOnce stopped before its end
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

function test_dwpli_features_use_signed_weights(tc)
% dwPLI is centred on 0 without coupling (negative half the time): its
% strengths are the sums of all its weights (BCT strengths_und_sign:
% positive minus negative strength), not of the positive part, which the
% noise biases upwards; no clustering coefficient or
% path length (they need non-negative weights)
D = tc.TestData.D0;  nodes = tc.TestData.nodes;  n = 68;  r = 4;
D.measures = {'dwpli'};
D.conn = D.conn - 0.3;                                      % centred on 0: about half the weights negative
F = static_features(D, 1, nodes, {'global', 'networks', 'regions'});
W = squeeze(D.conn(1, :, :, r));
W(1:n + 1:end) = 0;
verifyGreaterThan(tc, nnz(W < 0), 0.3 * n * (n - 1));
v = @(f) F.value(F.id == D.ids(r) & F.feature == f);
verifyEqual(tc, sort(unique(F.feature(F.family == "global"))), "Strength");
verifyEqual(tc, v("Strength"), sum(W(:)) / n, 'AbsTol', 1e-12);
verifyEqual(tc, v("strength_" + regexprep(D.roi_labels{5}, '[^A-Za-z0-9]+', '_')), sum(W(5, :)), 'AbsTol', 1e-12);
Wp = max(W, 0);
verifyGreaterThan(tc, sum(Wp(:)) / n, v("Strength") + 1);  % the positive part: far above
D.measures = {'wpli'};                                      % other measures: the positive part
F = static_features(D, 1, nodes, {'global'});
verifyEqual(tc, F.value(F.id == D.ids(r) & F.feature == "Strength"), sum(Wp(:)) / n, 'AbsTol', 1e-12);
end

function test_dwpli_strength_equals_bct_signed(tc)
% the dwPLI strength against BCT's strengths_und_sign (positive minus
% negative strength), where BCT is installed in S.toolboxes/BCT (a test of
% its own, reported as skipped without BCT)
bct = fullfile(tc.TestData.S.toolboxes, 'BCT', '2019_03_03_BCT');
assumeTrue(tc, isfolder(bct), 'BCT not installed: the signed strength not compared with it');
if isempty(which('strengths_und_sign')), addpath(bct);  c = onCleanup(@() rmpath(bct)); end
D = tc.TestData.D0;  nodes = tc.TestData.nodes;  n = 68;  r = 4;
D.measures = {'dwpli'};
D.conn = D.conn - 0.3;
F = static_features(D, 1, nodes, {'global'});
W = squeeze(D.conn(1, :, :, r));
W(1:n + 1:end) = 0;
[Spos, Sneg] = strengths_und_sign(W);
verifyEqual(tc, F.value(F.id == D.ids(r) & F.feature == "Strength"), mean(Spos - Sneg), 'AbsTol', 1e-12);
end

function test_recordings_without_the_measure(tc)
% a recording whose matrix is NaN (stage 1 could not compute the measure,
% e.g. the symmetric orthogonalisation with fewer independent signals than
% regions): its connectivity features are NaN, missing in the tests (not
% 0), its power kept; the NBS leaves it out with a warning (its edges do
% not enter the statistic as 0); the reliability uses the
% other recordings
S = tc.TestData.S;  S.nbs.perms = 50;  S.enrichment.n_perm = 50;  S.figures = false;
D = tc.TestData.D;  nodes = tc.TestData.nodes;  n = 68;
A = nan(n);  A(1:n + 1:end) = 0;
D.conn(1, :, :, 3) = A;  D.halves(:, 1, :, :, 3) = repmat(reshape(A, [1 1 n n]), 2, 1);
F = static_features(D, 1, nodes, {'global', 'networks', 'regions', 'power'});
x = F(F.id == D.ids(3), :);
verifyTrue(tc, all(isnan(x.value(x.family ~= "power"))));
verifyTrue(tc, all(isfinite(x.value(x.family == "power"))));
verifyTrue(tc, all(isfinite(F.value(F.id ~= D.ids(3)))));
R = static_reliability(D, 1, nodes, {'global'});
verifyEqual(tc, R.n(R.feature == "Strength"), 19);
T = table(D.ids, repmat("control", 20, 1), 'VariableNames', {S.id_column, S.group_column});
T.(S.group_column)(D.group == 2) = "patient";
out = fullfile(tc.TestData.out, 'nan_recording');
Rn = verifyWarning(tc, @() compute_NBS_analysis(S, D, 1, T, nodes, tc.TestData.Gatlas, out, 'full'), ...
    'compute_NBS_analysis:not_computed');
verifyTrue(tc, istable(Rn));
E = readtable(fullfile(out, 'edges_patient_control_full.csv'), 'TextType', 'string');
verifyTrue(tc, all(isfinite(E.mean_A)) && all(isfinite(E.mean_B)));
end

function test_reliability_needs_half_of_the_recordings(tc)
% a split-half reliability only from at least max(3, half of the
% recordings) recordings with both halves: 9 of 20 give NaN with the
% reason and a warning, 10 of 20 a value
D = tc.TestData.D;  nodes = tc.TestData.nodes;  n = 68;
D.halves(2, 1, :, :, 1:11) = NaN;
R = verifyWarning(tc, @() static_reliability(D, 1, nodes, {'global'}), 'static_reliability:few');
verifyTrue(tc, all(isnan(R.r(R.family ~= "edges"))) && all(isnan(R.r_sb)));
verifyEqual(tc, unique(R.n), 9);
verifyTrue(tc, all(contains(R.note, "9 of 20 recordings with both halves, fewer than 10")));
D.halves(2, 1, :, :, 11) = D.halves(1, 1, :, :, 11);
R = verifyWarningFree(tc, @() static_reliability(D, 1, nodes, {'global'}));
verifyEqual(tc, R.n(R.feature == "Strength"), 10);
verifyTrue(tc, all(isfinite(R.r_sb(R.family == "edges"))) && all(isfinite(R.r(R.family == "global"))));
verifyEqual(tc, size(D.halves, 3), n);
end

function test_no_feature_families(tc)
% S.static_families empty (the edges only): an empty feature table with
% its columns, and the reliability of the edges and patterns
D = tc.TestData.D;  nodes = tc.TestData.nodes;
F = static_features(D, 1, nodes, {});
verifyEqual(tc, height(F), 0);
verifyEqual(tc, F.Properties.VariableNames, {'id', 'family', 'feature', 'state', 'value'});
verifyEqual(tc, height(static_features(D, 1, table(), {'networks'})), 0);   % no network table: none either
R = static_reliability(D, 1, nodes, {});
verifyEqual(tc, unique(R.family), ["edges"; "pattern"]);
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
verifyEqual(tc, R.Properties.VariableNames, {'family', 'feature', 'n', 'r', 'r_sb', 'r_within_groups', 'note'});
w = R(R.feature == "within_DMN", :);
verifyGreaterThan(tc, w.r, 0.8);
verifyEqual(tc, w.r_sb, 2 * w.r / (1 + w.r), 'AbsTol', 1e-12);
verifyTrue(tc, all(isnan(R.r_within_groups)));             % without the groups
% with the groups: the planted group difference in the DMN
% raises the pooled r; within the groups it is r of the recordings' own
% variation, the halves of each group centred
g = ["A"; "B"];  g = g(D.group);
Rg = static_reliability(D, 1, nodes, {'global', 'networks'}, g);
wg = Rg(Rg.feature == "within_DMN", :);
verifyEqual(tc, wg.r, w.r, 'AbsTol', 1e-12);
verifyLessThan(tc, wg.r_within_groups, wg.r);
F1 = static_features(D, 1, nodes, {'networks'}, 1);  F2 = static_features(D, 1, nodes, {'networks'}, 2);
a = F1.value(F1.feature == "within_DMN");  b = F2.value(F2.feature == "within_DMN");
for k = 1:2
    m = D.group == k;
    a(m) = a(m) - mean(a(m));  b(m) = b(m) - mean(b(m));
end
verifyEqual(tc, wg.r_within_groups, corr(a, b), 'AbsTol', 1e-12);
verifyEqual(tc, R.r(R.feature == "identification"), 1);      % every recording has its own pattern
D.halves(2, :, :, :, :) = D.halves(1, :, :, :, :);
R = static_reliability(D, 1, nodes, {'global'});
verifyEqual(tc, R.r(R.feature == "Strength"), 1, 'AbsTol', 1e-12);
verifyEqual(tc, R.r_sb(R.feature == "edge r_sb, median"), 1, 'AbsTol', 1e-12);
D.halves(2, :, :, :, :) = 1 - D.halves(1, :, :, :, :);      % halves that disagree: r < 0, r_sb = r
R = static_reliability(D, 1, nodes, {'global'});
r = R.r(R.feature == "Strength");
verifyLessThan(tc, r, -0.99);
verifyEqual(tc, R.r_sb(R.feature == "Strength"), r);         % 2r / (1 + r) would give NaN or values below -1
verifyLessThan(tc, R.r_sb(R.feature == "edge r_sb, median"), 0);
end

function test_leverage_and_skipped_bands(tc)
% a score with a recording of high leverage warns (the HC3
% permutation test can reject too often with an extreme score and a
% spread that changes with it), not with the wild bootstrap or Spearman;
% recordings left out of an anchored band are counted per group
d = tempname;  mkdir(fullfile(d, 'Connectivity'));  c = onCleanup(@() rmdir(d, 's'));
ids = "s" + (1:30)';
skipped = strings(30, 1);  skipped([1 2 4 6]) = "alpha_iaf";  skipped(28) = "alpha_iaf theta_iaf";
summary = fullfile(d, 'Connectivity', 'stage1_summary.csv');
writetable(table(ids + "_rest.set", repmat("ok", 30, 1), skipped, 'VariableNames', {'file', 'status', 'bands_skipped'}), ...
    summary);
grp = repmat(["HV"; "SZ"], 15, 1);                   % s1 HV, s2 SZ, ...
score = (1:30)' / 30;  score(5) = 40;                % one extreme score
f = fullfile(d, 'subjects.csv');
writetable(table(ids, grp, score, 'VariableNames', {'ID', 'grp', 'PANSS'}), f);
S = tc.TestData.S;
S.results = d;  S.subs_info = f;  S.id_column = 'ID';  S.id_regexp = '^([^_]+)';
S.group_column = 'grp';  S.groups = {'HV', 'SZ'};  S.spectral_profile = false;  S.iaf_used = false;
S.comparisons = struct('name', 'SZ - HV', 'A', {{'grp', '==', 'SZ'}}, 'B', {{'grp', '==', 'HV'}});
S.covariates = {};
S.correlations = struct('name', 'Symptoms', 'columns', '^PANSS$', 'subjects', {{}}, 'covariates', {{}}, ...
    'zero_as_missing', false);
S.covariate_sets = struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {});
S.correlation_type = 'pearson';  S.correlation_test = 'permutation';
verifyWarning(tc, @() check_sample_sizes(S), 'run_statistics:leverage');
lastwarn('');  evalc('check_sample_sizes(S)');  [msg, id] = lastwarn();   % the correlations come last
verifyEqual(tc, id, 'run_statistics:leverage');
verifySubstring(tc, msg, 'PANSS (h = ');
for o = {{'correlation_test', 'wild_bootstrap'}, {'correlation_type', 'spearman'}}
    So = S;  So.(o{1}{1}) = o{1}{2};
    lastwarn('');  evalc('check_sample_sizes(So)');  [~, id] = lastwarn();
    verifyNotEqual(tc, id, 'run_statistics:leverage');
end
S.correlations = S.correlations([]);
verifyWarning(tc, @() check_sample_sizes(S), 'run_statistics:bands_skipped');
lastwarn('');  evalc('check_sample_sizes(S)');  msg = lastwarn();       % the last band named: theta_iaf
verifySubstring(tc, msg, 'Band theta_iaf: left out of 1 recording(s)');
verifySubstring(tc, msg, 'SZ - HV: 1 of 15 and 0 of 15');               % s28
skipped(28) = "alpha_iaf";
writetable(table(ids + "_rest.set", repmat("ok", 30, 1), skipped, 'VariableNames', {'file', 'status', 'bands_skipped'}), ...
    summary);
lastwarn('');  evalc('check_sample_sizes(S)');  msg = lastwarn();
verifySubstring(tc, msg, 'Band alpha_iaf: left out of 5 recording(s)');
verifySubstring(tc, msg, 'SZ - HV: 4 of 15 and 1 of 15');               % s2, s4, s6, s28; s1
end

function test_subject_table_text_and_missing_codes(tc)
% a covariate of numbers with one text value ("33 years") is an error
% naming it (read as text it would enter as one indicator per value);
% a usual missing-value code far from its column's values warns, once
d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
f = fullfile(d, 'subjects.csv');
age = string(20 + (1:12)');  age(3) = "33 years";
edu = 10 + mod((1:12)', 5);  edu(7) = -999;
writetable(table("s" + (1:12)', repmat(["HV"; "SZ"], 6, 1), age, edu, repmat(["F"; "M"; "M"], 4, 1), ...
    'VariableNames', {'ID', 'grp', 'Age', 'EDU', 'sex'}), f);
S = tc.TestData.S;
S.subs_info = f;  S.id_column = 'ID';  S.group_column = 'grp';
S.correlations = struct('name', {}, 'columns', {}, 'subjects', {}, 'covariates', {}, 'zero_as_missing', {});
S.covariate_sets = struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {});
S.covariates = {'Age', 'sex'};
verifyError(tc, @() read_subject_table(S), 'read_subject_table:text_covariate');
try
    read_subject_table(S);
catch err
    verifySubstring(tc, err.message, '"33 years"');
end
S.categorical_covariates = {'Age'};                       % named as categories: no error
S.covariates = {'sex'};                                   % (Age not a covariate: no error either)
verifyWarning(tc, @() read_subject_table(S), 'read_subject_table:missing_codes');
lastwarn('');
T = read_subject_table(S);                               % once per file and session
[~, id] = lastwarn();
verifyNotEqual(tc, id, 'read_subject_table:missing_codes');
verifyEqual(tc, T.EDU(T.ID == "s7"), -999);              % used as a number: the warning says so
age(3) = "-";                                             % no digits: a missing value, no error
writetable(table("s" + (1:12)', repmat(["HV"; "SZ"], 6, 1), age, ...
    'VariableNames', {'ID', 'grp', 'Age'}), f);
S.categorical_covariates = {};  S.covariates = {'Age'};
T = read_subject_table(S);
verifyTrue(tc, isnumeric(T.Age) && isnan(T.Age(3)) && T.Age(4) == 24);
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

function test_group_matrices_without_the_measure(tc)
% recordings whose matrix is NaN, one in each group, are left out of the
% group-matrix figure (they would make both group means NaN and the colour
% scale invalid: "Error setting property 'CLim'"); with no recording left: no
% figure and a warning, no error
S = tc.TestData.S;  S.figures = true;
D = tc.TestData.D;  n = 68;
A = nan(n);  A(1:n + 1:end) = 0;
D.conn(1, :, :, 3) = A;  D.conn(1, :, :, 13) = A;
verifyEqual(tc, D.group([3 13])', [1 2]);
T = table(D.ids(:), string(S.groups(D.group))', 'VariableNames', {'ID', 'group'});
file = fullfile(tc.TestData.out, 'group_matrices_nan.png');
verifyWarningFree(tc, @() group_matrices(S, D, 1, T, tc.TestData.nodes, file));
verifyTrue(tc, isfile(file));
D.conn(:) = NaN;
file0 = fullfile(tc.TestData.out, 'group_matrices_all_nan.png');
verifyWarning(tc, @() group_matrices(S, D, 1, T, tc.TestData.nodes, file0), 'group_matrices:no_values');
verifyFalse(tc, isfile(file0));
end

function test_failure_summary_figures_apart(tc)
% the end of a run: figure steps that failed are listed after the others;
% when only figures failed, the run is not called "Done with errors"
F = ["eloreta_alpha_wpli, window_mean: group matrices figure", "CLim"; ...
     "PSD: PSD figures", "no axes"];
[head, lines, n_fig] = failure_summary(F, 12.34);
verifyEqual(tc, n_fig, 2);
verifySubstring(tc, head, 'Done in 12.3 min: every result written, but 2 figure step(s) failed');
verifyFalse(tc, contains(head, 'with errors'));
verifyEqual(tc, lines, ["eloreta_alpha_wpli, window_mean: group matrices figure (figure only): CLim"; ...
    "PSD: PSD figures (figure only): no axes"]);
F = [F; "eloreta_alpha_wpli: microfeatures_raw", "covariate missing"];
[head, lines, n_fig] = failure_summary(F, 5);
verifyEqual(tc, n_fig, 2);
verifyEqual(tc, head, 'Done with errors in 5.0 min: 3 step(s) failed (2 of them figures only)');
verifyEqual(tc, lines(1), "eloreta_alpha_wpli: microfeatures_raw: covariate missing");
[head, ~, n_fig] = failure_summary(F(3, :), 5);
verifyEqual(tc, n_fig, 0);
verifyEqual(tc, head, 'Done with errors in 5.0 min: 1 step(s) failed');
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
%% Small or unusual datasets
function test_temporal_channels_with_old_names(tc)
% T3 and T4 (kept by the 10-05 cap under these names) are temporal channels:
% the muscle index of a cap labelled T3/T4 is not NaN
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
% the comparisons can name the groups (without it, every family would stop
% on "Column group is not in the subject table")
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
% 2 + 2 recordings: with a covariate the comparison and the correlations
% cannot be tested (rows of NaN, no correlation rows), and a warning names
% them. Without the covariate both can be tested. Recordings whose stage 1
% run_all reused (status "skipped (up to date)") count too, so statistics
% run after a reused stage 1 see them.
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
% A cap whose channels have no 10-05 names (an ANT cap: Z1 ...
% RR14; EGI: E1 ... E256): matched by name, every posterior, frontal,
% frontopolar and temporal measure of the spectral profile would be NaN,
% and the IAF could not be found, so the sets are matched by position
% (cap_channels); a cap with 10-20 names is matched by name.
T = get_ftElec('ASA_10-05_343');
keep = {'Fp1', 'Fpz', 'Fp2', 'F3', 'Fz', 'F4', 'FT7', 'T7', 'T8', 'FT8', 'C3', 'Cz', 'C4', 'P3', 'Pz', 'P4', ...
    'PO7', 'PO8', 'O1', 'Oz', 'O2'};
[~, i] = ismember(keep, T.label);
names = cellstr("E" + (1:numel(keep)));
elec = struct('label', {names'}, 'unit', 'm', 'chanpos', T.chanpos(i, :) + [0.003 0 0]);   % 3 mm off
e = @(lab) names(cellfun(@(x) find(strcmp(keep, x)), lab));                % the E name of 10-05 names
% by name (any case, the spelling of the labels), the cap not used
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
