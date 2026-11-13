function rep = run_static_simulation(name, over, S)
%RUN_STATIC_SIMULATION Simulate a dataset with a planted coupling difference, run the static-only pipeline, score the recovery.
%   rep = run_static_simulation(name, over, S)
%   name  folder under simulated/ (simulated/<name>/: the recordings,
%         participants.tsv, truth.mat, results/), or a full path
%   over  struct of fields of sim.static_params to change (e.g.
%         struct('n_per_group', [8 8], 'coupling', [0.3 0.3]) for a null)
%   S     settings (default pipeline_settings() with S.static_matrices =
%         {'full'}); run static only (S.dynamic = false) with the family
%         'networks' among S.static_families; the data, subject table,
%         groups, comparison and correlation are set here for the dataset,
%         and only the band of S.bands that holds the planted band is run.
%         The report reads the whole-recording results: S.static_matrices
%         must hold 'full' (an error before anything runs otherwise)
%   The recordings are simulated once (a second call reuses them and the
%   stage 1 outputs). rep: what was planted and what the statistics found,
%   for the first measure of S.conn (wPLI when it is there):
%     planted        network, coupling per group
%     coupling_r     correlation across recordings of the planted coupling
%                    with the measured within-network connectivity
%     group_effect   the networks family's row of the planted network
%                    (g, p_maxT_family, p_bb, significant)
%     best_block     the block with the largest |g| (should be the planted one)
%     other_blocks   significant blocks other than the planted network's
%     enrichment     the enrichment row of the planted network (z_spin, q_spin)
%     nbs            NBS components with p_both < alpha, and the share of
%                    their regions inside the planted network
%     correlation    the correlation of the planted score with the
%                    network's connectivity in the last group (r, q)
%   and writes results/recovery_report.md.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'));
if nargin < 2 || isempty(over), over = struct(); end
if nargin < 3 || isempty(S), S = pipeline_settings();  S.static_matrices = {'full'}; end
% the report reads the results of the whole recording ('full'): stop here,
% before the simulation and stage 1, when they will not be computed
m = {};
if isfield(S, 'static_matrices'), m = cellstr(S.static_matrices); end
if ~ismember('full', m)
    error('run_static_simulation:full', ['The recovery report reads the results of the whole recording: ' ...
        'S.static_matrices must hold ''full'' ({''full''}, or {''full'', ''window_mean''} for both); ' ...
        'it asks for the window mean only']);
end
S.dynamic = false;
if ~ismember('networks', S.static_families), S.static_families{end + 1} = 'networks'; end
S.spectral_profile = false;                     % scalp measures of simulated noise: not informative (and no restingIAF needed)
S = setup_toolboxes(S);
p = sim.static_params(S);
p.out_dir = fullfile(S.root, 'simulated', name);
if contains(name, filesep), p.out_dir = name; end      % a folder
for f = fieldnames(over)', p.(f{1}) = over.(f{1}); end
if ~isfile(fullfile(p.out_dir, 'truth.mat'))
    fprintf('Simulating %s: %s recordings\n', name, strjoin(string(p.n_per_group), ' + '));
    sim.simulate_static_dataset(S, p);
end
L = load(fullfile(p.out_dir, 'truth.mat'));
truth = L.truth;
if ~isequal(truth.params.n_per_group, p.n_per_group) || ~isequal(truth.params.coupling, p.coupling) || ...
        truth.params.duration_s ~= p.duration_s
    error('run_static_simulation:stale', '%s holds another simulation: delete the folder to simulate again', p.out_dir);
end
% the dataset's settings: BIDS layout, the simulated cap, the planted band
S.data = p.out_dir;
S.file_pattern = {'sub-*/eeg/*_eeg.set', 'sub-*/eeg/*_eeg.edf', 'sub-*/eeg/*_eeg.bdf', 'sub-*/eeg/*_eeg.vhdr'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.results = fullfile(p.out_dir, 'results');
S.subs_info = fullfile(p.out_dir, 'participants.tsv');
S.id_column = 'participant_id';
S.group_column = 'group';
S.groups = p.groups;
S.comparisons = struct('name', sprintf('%s - %s', p.groups{end}, p.groups{1}), ...
    'A', {{'group', '==', p.groups{end}}}, 'B', {{'group', '==', p.groups{1}}});
S.correlations = struct('name', 'Score', 'columns', '^score$', 'subjects', {{'group', '==', p.groups{end}}}, ...
    'covariates', {{'age'}}, 'zero_as_missing', false);
S = check_settings(S);
bands = S.bands;
if S.individual_bands, bands = S.bands_individual; end
b = find(arrayfun(@(B) isnumeric(B.fmin) && isnumeric(B.fmax) && B.fmin <= p.band(1) && B.fmax >= p.band(2), bands), 1);
if isempty(b)
    error('run_static_simulation:band', 'No fixed band of S.bands holds the planted %g-%g Hz', p.band(1), p.band(2));
end
S.bands = bands(b);                             % the band of the planted coupling only
S.individual_bands = false;
S.psd.run = false;                              % connectivity is evaluated (the PSD step: tests/test_psd.m)
S = check_settings(S);
t = tic;
run_full_pipeline(S);                           % stage 1 and the statistics (no stage 2: S.dynamic = false)
minutes = toc(t) / 60;

% what was found
out = analysis_dir(S);
measure = S.static_measures{1};
if ismember('wpli', S.static_measures), measure = 'wpli'; end
tag = sprintf('%s_%s_%s', S.inv{1}, S.bands(1).name, measure);
sdir = fullfile(out, 'Statistical analysis', tag);
ndir = fullfile(out, 'NBS', tag);
G = readtable(fullfile(sdir, 'static_networks_full_group_comparison.csv'), 'TextType', 'string');
feat = "within_" + string(p.network);
rep.planted = struct('network', p.network, 'coupling', p.coupling, 'n_per_group', p.n_per_group, ...
    'band_hz', p.band, 'planted_r', p.planted_r);
rep.measure = measure;
rep.minutes = minutes;
rep.group_effect = G(G.Parameter == feat, {'Parameter', 'g', 't', 'p_unc', 'p_maxT_family', 'p_bb', 'significant'});
[~, ib] = max(abs(G.g));
rep.best_block = G(ib, {'Parameter', 'g', 'q_fdr'});
rep.other_blocks = G(G.significant & G.Parameter ~= feat, {'Parameter', 'g', 'p_bb'});
F = readtable(fullfile(sdir, 'features_full.csv'), 'TextType', 'string');
v = F(F.feature == feat, :);
[ok, row] = ismember(truth.ids, v.id);
rep.coupling_r = corr(truth.coupling(ok), v.value(row(ok)));
E = table();
ef = fullfile(ndir, 'enrichment_full.csv');
if isfile(ef)
    E = readtable(ef, 'TextType', 'string');
    E = E(E.network_a == string(p.network) & E.network_b == string(p.network), :);
end
rep.enrichment = E;
N = readtable(fullfile(ndir, 'NBS_full.csv'), 'TextType', 'string');
N = N(N.p_both < S.alpha, :);
inside = nan(height(N), 1);
for i = 1:height(N)
    r = strtrim(split(N.regions(i), ';'));
    inside(i) = mean(ismember(r, string(truth.network_regions)));
end
N.share_regions_in_network = inside;
rep.nbs = N;
C = table();
cf = fullfile(sdir, 'static_networks_full_correlations.csv');
if isfile(cf)
    C = readtable(cf, 'TextType', 'string');
    C = C(C.Parameter == feat, {'Scale', 'r', 'p_unc', 'p_maxT_family', 'p_bb', 'significant'});
end
rep.correlation = C;
report(rep, fullfile(S.results, 'recovery_report.md'), tag);
end

function report(rep, file, tag)
fid = fopen(file, 'w');
c = onCleanup(@() fclose(fid));
P = rep.planted;
fprintf(fid, '# Recovery of a planted static coupling difference (%s)\n\n', tag);
fprintf(fid, 'Planted: coupling of the %s network in %g-%g Hz, %s (groups %s recordings); score correlated r = %.2f with it in the last group.\n\n', ...
    P.network, P.band_hz(1), P.band_hz(2), mat2str(P.coupling), mat2str(P.n_per_group), P.planted_r);
fprintf(fid, '- Planted coupling vs measured within-network connectivity, across recordings: r = %.2f\n', rep.coupling_r);
if height(rep.group_effect)
    g = rep.group_effect;
    fprintf(fid, '- Group difference of within_%s: g = %.2f, family-wise p = %.4f, over families = %.4f, significant = %d\n', ...
        P.network, g.g, g.p_maxT_family, g.p_bb, g.significant);
end
fprintf(fid, '- Largest group difference among the network blocks: %s (g = %.2f)\n', rep.best_block.Parameter, rep.best_block.g);
fprintf(fid, '- Other significant blocks: %d\n', height(rep.other_blocks));
if height(rep.enrichment)
    fprintf(fid, '- Enrichment of the group difference inside %s: z_spin = %.2f, q_spin = %.4f\n', P.network, ...
        rep.enrichment.z_spin, rep.enrichment.q_spin);
end
fprintf(fid, '- NBS components with p (both directions) < alpha: %d', height(rep.nbs));
if height(rep.nbs)
    fprintf(fid, ' (largest %d edges, %.0f%% of its regions in %s)', max(rep.nbs.n_edges), ...
        100 * rep.nbs.share_regions_in_network(find(rep.nbs.n_edges == max(rep.nbs.n_edges), 1)), P.network);
end
fprintf(fid, '\n');
if height(rep.correlation)
    fprintf(fid, '- Correlation of the score with within_%s (last group): r = %.2f, family-wise p = %.4f\n', P.network, ...
        rep.correlation.r(1), rep.correlation.p_maxT_family(1));
end
fprintf(fid, '\nStage 1 and the statistics took %.1f min.\n', rep.minutes);
end
