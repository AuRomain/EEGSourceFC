function rep = run_static_demo()
%RUN_STATIC_DEMO The static-only pipeline on simulated recordings, and what it recovered.
%   rep = run_static_demo() simulates 16 + 16 recordings of 120 s (64
%   channels of a 10-10 cap, 256 Hz, EEGLAB files in a BIDS layout,
%   participants.tsv) in which the regions of the default mode network
%   share a lagged 8-12 Hz driver, more strongly in the patients, and in
%   the patients a score correlated with that coupling (sim.static_params);
%   runs stage 1 and the statistics with S.dynamic = false (no states) on
%   the files as on real data, in the alpha band (7-13 Hz), with wPLI and
%   orthogonalised AEC over the whole recording (S.static_matrices =
%   {'full'}, below) and the four feature families; and compares what the
%   pipeline found with what was planted (run_static_simulation). No real
%   data are needed, only the toolboxes of the README, sections 1 and 2
%   (FieldTrip and Brainstorm only).
%   Outputs in simulated/static_demo/: the recordings, participants.tsv,
%   truth.mat and results/ (Connectivity/, Statistics/main/) with
%   results/recovery_report.md. A second call reuses the simulated
%   recordings and the stage 1 outputs. About 10 minutes with a parallel
%   pool of 8 workers. The settings are the defaults of pipeline_settings.m
%   with the static-only values of the static template, but the whole
%   recording instead of the window mean (the default): the planted
%   coupling is weak and stable, which the whole recording measures and the
%   mean of 0.6 s alpha windows does not (on these recordings: within-DMN
%   g 0.42 instead of 2.17, nothing significant; docs/static_mode.md,
%   "Whole recording or windows"). The demo's report reads the
%   whole-recording results.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'simulation'));
S = pipeline_settings();
S.dynamic = false;
S.conn = {'wpli', 'aec_orth'};
S.static_matrices = {'full'};          % the whole recording (header; the report reads it)
S.static_families = {'global', 'networks', 'regions', 'power'};
S.filter_transition = 'brainstorm';
S.save_scouts = false;
S.bands = struct('name', 'alpha', 'fmin', 7, 'fmax', 13);
S.individual_bands = false;
S.run_name = 'main';
t = tic;
rep = run_static_simulation('static_demo', struct(), S);
fprintf('\nDemo finished in %.1f min.\n', toc(t) / 60);
fprintf('Planted: stronger coupling of the DMN regions in the patients (%s).\n', mat2str(rep.planted.coupling));
fprintf('Measured within-DMN %s against the planted coupling, across recordings: r = %.2f\n', rep.measure, rep.coupling_r);
if height(rep.group_effect)
    g = rep.group_effect;
    fprintf('within_DMN, patient - control: g = %.2f, family-wise p = %.4f (over all families: %.4f)\n', g.g, ...
        g.p_maxT_family, g.p_bb);
end
fprintf(['Largest block difference: %s (g = %.2f); %d other blocks significant (source leakage spreads ' ...
    'the planted coupling to neighbouring regions: README)\n'], rep.best_block.Parameter, ...
    rep.best_block.g, height(rep.other_blocks));
fprintf('NBS components with p (both directions) < 0.05: %d\n', height(rep.nbs));
if height(rep.correlation)
    fprintf('Score against within_DMN in the patients: r = %.2f, family-wise p = %.4f\n', rep.correlation.r(1), ...
        rep.correlation.p_maxT_family(1));
end
fprintf('Outputs and the recovery report: %s\n', fullfile(root, 'simulated', 'static_demo', 'results'));
end
