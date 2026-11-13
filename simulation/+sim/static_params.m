function p = static_params(S)
%STATIC_PARAMS Settings of a simulated dataset with a planted static coupling difference.
%   p = sim.static_params(S), S from pipeline_settings (the forward model
%   is the one of stage 1: S.montages{1}, S.atlases{1}, S.network_table).
%   For the static-only pipeline (S.dynamic = false: run_static_simulation,
%   run_static_demo; the states are simulated by sim.default_params). The
%   model, per recording (sim.simulate_static_dataset):
%   - every region carries an independent oscillation in p.band and 1/f
%     background activity, on a 1-cm cortical patch around its seed;
%   - the regions of the network p.network (Yeo-7) share a common driver
%     in p.band with fixed, region-specific phase lags (non-zero: wPLI sees
%     it), with weight c, the recording's coupling: drawn around its
%     group's p.coupling with SD p.coupling_sd; nothing else differs
%     between the groups (no amplitude change: p.state_gain = 1);
%   - diffuse background from p.brain_noise_sources random cortical
%     vertices, then white sensor noise.
%   The subject table (participants.tsv) has age and sex, unrelated to the
%   EEG, and in the last group a score correlated with the coupling
%   (p.planted_r). p.coupling equal in the groups and p.planted_r = 0
%   give a null dataset. Nothing here comes from real participants.
p.name = 'static_sim';
p.out_dir = fullfile(S.root, 'simulated', p.name);
p.groups = {'control', 'patient'};
p.n_per_group = [16, 16];
p.srate = 256;
p.channels = {'Fp1', 'Fpz', 'Fp2', 'AF7', 'AF3', 'AFz', 'AF4', 'AF8', 'F7', 'F5', 'F3', 'F1', 'Fz', ...
    'F2', 'F4', 'F6', 'F8', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'T7', ...
    'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', ...
    'CP4', 'CP6', 'TP8', 'TP10', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'PO7', 'PO3', ...
    'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2', 'Iz'};   % a 64-channel 10-10 cap (labels of S.montages{1}); {} = every channel
p.duration_s = 120;
p.format = 'set';                  % 'set', 'edf', 'bdf' or 'vhdr', or a cell (one per recording, cycled)
p.band = [8, 12];
p.network = 'DMN';                 % Yeo-7 network of the planted coupling (DK-68: 20 regions)
p.coupling = [0.30, 0.45];         % group means of the coupling weight (0-1), one per group
p.coupling_sd = 0.08;              % recording-to-recording SD
p.planted_r = 0.6;                 % correlation of the last group's score with its coupling
p.lag_jitter_rad = 0.2;            % recording-to-recording variation of the regions' phase lags
p.state_gain = 1;                  % amplitude of the network unchanged
p.crossfade_s = 0.05;
p.patch_radius_m = 0.01;
p.osc_rms = 1;
p.pink_rms = 1;
p.alpha_rms = 0;
p.brain_noise_sources = 300;
p.brain_noise_rms = 1;
p.sensor_noise = 0.05;             % fraction of the sensor signal RMS
p.scale_uv = 10;                   % sensor RMS in microvolts
p.seed = 1;
end
