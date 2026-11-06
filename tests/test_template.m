function tests = test_template
%TEST_TEMPLATE The settings template of this folder and the examples.
%   From v_static/:  results = runtests('tests/test_template.m')
%   pipeline_settings.m is this folder's own (tests/test_same_code_as_v_dynamic.m):
%   static only, six bands with the individual-band option, wPLI and
%   orthogonalised AEC, every feature family, the PSD step with its eight
%   bands; every example is valid.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'examples'));
tc.TestData.root = root;
end

function test_template_is_valid(tc)
S = check_settings(pipeline_settings());
verifyFalse(tc, S.dynamic);
verifyEqual(tc, {S.bands.name}, {'delta', 'theta', 'alpha', 'beta1', 'beta2', 'gamma'});
verifyEqual(tc, [S.bands.fmin; S.bands.fmax], [1 4 7 13 20 30; 4 7 13 20 30 45]);
verifyEqual(tc, S.conn, {'wpli', 'aec_orth'});
verifyEqual(tc, S.static_measures, {'wpli', 'aec_orth'});
verifyEqual(tc, S.static_matrices, {'window_mean'});   % the default (6 Nov 2026)
verifyEqual(tc, S.min_segment, 6);                     % the AEC windows fit in every piece kept
verifyEqual(tc, S.measure_windows, struct());
verifyEqual(tc, S.static_families, {'global', 'networks', 'regions', 'power'});
verifyTrue(tc, isfile(S.network_table));
verifyFalse(tc, isfolder(fullfile(S.root, '2- Dynamic connectivity')));
verifyTrue(tc, S.psd.run);                                   % the PSD step, with its own bands
verifyEqual(tc, S.psd.levels, {'sensor', 'source'});
verifyEqual(tc, {S.psd.bands.name}, {'delta', 'theta', 'alpha', 'alpha1', 'alpha2', 'beta', 'beta1', 'beta2'});
verifyEqual(tc, [S.psd.bands.fmin; S.psd.bands.fmax], [1 4 7 7 10 13 13 20; 4 7 13 10 13 30 20 30]);
verifyTrue(tc, isfolder(fullfile(S.root, '1b- PSD analysis')));
end

function test_individual_bands_switch(tc)
% the six bands anchored on the individual alpha frequency: with IAF = 10 Hz
% the fixed ones
S = pipeline_settings();
S.individual_bands = true;
S = check_settings(S);
verifyTrue(tc, S.iaf_used);
edges = zeros(numel(S.bands), 2);
for b = 1:numel(S.bands)
    [edges(b, 1), edges(b, 2)] = band_edges(S.bands(b), 10, S.iaf);
end
verifyEqual(tc, edges, [1 4; 4 8; 8 12; 12 20; 20 30; 30 45]);
end

function test_examples_are_valid(tc)
d = dir(fullfile(tc.TestData.root, 'examples', 'settings_*.m'));
verifyNotEmpty(tc, d);
for i = 1:numel(d)
    [~, name] = fileparts(d(i).name);
    S = check_settings(feval(name));
    verifyFalse(tc, S.dynamic, name);
    verifyWarningFree(tc, @() check_settings(feval(name)), name);   % no setting of the windows
end
verifyWarningFree(tc, @() check_settings(pipeline_settings()));
end
