function S = static_test_settings()
%STATIC_TEST_SETTINGS Settings of the static-only tests (S.dynamic = false).
%   S = static_test_settings(): pipeline_settings() with every field that
%   the static tests depend on set here, so that the tests (test_static_*)
%   give the same results with a template with windows and states and with
%   a static-only one: eLORETA on the 68 Desikan-Killiany regions of
%   the 10-05 template cap, the alpha band, wPLI and orthogonalised AEC, the
%   four feature families, split halves of 10-s blocks, Brainstorm's filter
%   widths, the whole-recording matrices (S.static_matrices = {'full'}; the
%   default is the window mean) with the template's windows and pieces of
%   at least 4 s. The tests change the data, results and statistics fields.
S = pipeline_settings();
S.dynamic = false;
S.montages = {'ASA_10-05_343'};
S.channel_subset = {};
S.atlases = {'Desikan_Killiany_68'};
S.anatomy = '';
S.anatomy_individual = false;
S.anatomy_fallback = false;
S.inv = {'eloreta'};
S.orientation = 'normal';
S.scout_function = 'mean';
S.filter_transition = 'brainstorm';
S.edge_trim = 1;
S.bands = struct('name', 'alpha', 'fmin', 7, 'fmax', 13);
S.individual_bands = false;
S.conn = {'wpli', 'aec_orth'};
S.state_conn = 'wpli';
S.save_scouts = false;
S.reliability = true;
S.reliability_block = 10;
S.static_measures = {};
S.static_matrices = {'full'};                      % the whole recording (the default: the window mean)
S.min_segment = 4;
S.window_cycles = 6;
S.window_length = [];
S.window_overlap = 0.5;
S.aec_window = 6;
S.aec_step = 0.5;
S.static_families = {'global', 'networks', 'regions', 'power'};
S.network_table = fullfile(stage_folder(S.root, 'Analysis'), 'DK68_networks.csv');
S.language_file = '';
S.groups = {'control', 'patient'};
S.group_column = 'group';
S.comparisons = struct('name', 'patient - control', 'A', {{'group', '==', 'patient'}}, ...
    'B', {{'group', '==', 'control'}});
S.covariates = {};
S.correlations = struct('name', {}, 'columns', {}, 'subjects', {}, 'covariates', {}, 'zero_as_missing', {});
S.covariate_sets = struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {});
S.n_perm = 2000;
S.alpha = 0.05;
S.seed = 20261015;
S.nbs = struct('thresh', 3.1, 'perms', 1000, 'size', 'extent');
S.enrichment = struct('n_perm', 2000);
S.spectral_profile = false;
S.figures = false;
S.run_name = 'test';
end
