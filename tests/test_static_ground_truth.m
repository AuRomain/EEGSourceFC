function tests = test_static_ground_truth
%TEST_STATIC_GROUND_TRUTH The static-only pipeline recovers a planted coupling difference, and finds nothing without one.
%   From the pipeline folder:  results = runtests('tests/test_static_ground_truth.m')
%   Simulated datasets (run_static_simulation, in a temporary folder kept
%   between runs): 16 + 16 recordings of 60 s whose default mode network
%   (DMN) regions share a lagged 8-12 Hz driver, more strongly in the
%   second group, with a score correlated with the coupling in that group;
%   and a null dataset, the same with equal coupling and no correlation.
%   Stage 1 in the alpha band, wPLI; statistics with 2000 permutations.
%   About 10 minutes with a parallel pool.
%   What is checked is what EEG source connectivity can claim: the planted
%   difference is found, in its direction, and is largest where the planted
%   network is involved. It is not confined to that network: source leakage
%   carries the DMN regions' lagged driver into neighbouring regions, so
%   most network blocks differ too (3 Nov 2026: within_DMN g 1.8,
%   between_SAL_DMN g 2.3, 17 other blocks significant) and an NBS
%   component spans several networks.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'));
S = static_test_settings();
S.conn = {'wpli'};
base = fullfile(tempdir, 'static_ground_truth');
tc.TestData.base = base;
tc.TestData.run = S.run_name;
evalc(['tc.TestData.planted = run_static_simulation(fullfile(base, ''planted''), ' ...
    'struct(''n_per_group'', [16 16], ''duration_s'', 60), S);']);
evalc(['tc.TestData.null = run_static_simulation(fullfile(base, ''null''), ' ...
    'struct(''n_per_group'', [16 16], ''duration_s'', 60, ''coupling'', [0.375 0.375], ''planted_r'', 0), S);']);
end

function test_planted_coupling_is_measured(tc)
% the measured within-DMN wPLI follows each recording's planted coupling
verifyGreaterThan(tc, tc.TestData.planted.coupling_r, 0.6);
end

function test_group_difference_is_found(tc)
% the planted block differs, in the planted direction, and is among the
% largest differences; blocks involving the DMN differ more than the others
rep = tc.TestData.planted;
g = rep.group_effect;
verifyGreaterThan(tc, g.g, 0.8);
verifyLessThan(tc, g.p_maxT_family, 0.05);
G = readtable(fullfile(tc.TestData.base, 'planted', 'results', 'Statistics', tc.TestData.run, ...
    'Statistical analysis', 'eloreta_alpha_wpli', 'static_networks_full_group_comparison.csv'), 'TextType', 'string');
[~, order] = sort(G.g, 'descend');
verifyLessThanOrEqual(tc, find(G.Parameter(order) == "within_DMN"), 3);
dmn = contains(G.Parameter, "DMN");
verifyGreaterThan(tc, mean(G.g(dmn)) - mean(G.g(~dmn)), 0.3);
end

function test_edges_differ_in_the_planted_direction(tc)
% an NBS component of stronger connectivity in the group with the stronger
% coupling, none in the other direction; the DMN edges' mean statistic is
% positive and above that of random networks of the same size
rep = tc.TestData.planted;
up = rep.nbs(rep.nbs.contrast == "patient > control", :);
verifyNotEmpty(tc, up);
verifyEmpty(tc, rep.nbs(rep.nbs.contrast == "control > patient", :));
e = rep.enrichment;
verifyNotEmpty(tc, e);
verifyGreaterThan(tc, e.mean_stat, 0);
verifyGreaterThan(tc, e.z_random, 0);
end

function test_correlation_is_found(tc)
c = tc.TestData.planted.correlation;
verifyNotEmpty(tc, c);
verifyGreaterThan(tc, c.r(1), 0.3);
end

function test_null_dataset_gives_nothing(tc)
rep = tc.TestData.null;
verifyGreaterThan(tc, rep.group_effect.p_maxT_family, 0.05);
verifyEmpty(tc, rep.other_blocks);
verifyTrue(tc, isempty(rep.nbs) || height(rep.nbs) == 0);
Fm = readtable(fullfile(tc.TestData.base, 'null', 'results', 'Statistics', tc.TestData.run, ...
    'Statistical analysis', 'families.csv'), 'TextType', 'string');
verifyLessThanOrEqual(tc, sum(Fm.n_significant), 1);    % at most one false positive in the run
end
