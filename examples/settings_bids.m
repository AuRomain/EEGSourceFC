function S = settings_bids()
%SETTINGS_BIDS Template: a BIDS dataset recorded with a standard cap, two groups.
%   S = settings_bids(): copy this file, rename it (settings_mystudy.m) and
%   edit the paths, columns, groups and scores. Every field not set here
%   keeps the default of pipeline_settings.m.
%       addpath('examples');  S = setup_toolboxes(settings_mystudy());
%       run_all(S);  run_statistics(S);
S = pipeline_settings();
S.results = fullfile(S.root, 'results_mystudy');

%% Recordings: cleaned files in sub-*/eeg/ folders
S.data          = '/data/mystudy/derivatives/cleaned';
S.file_pattern  = {'sub-*/eeg/*_eeg.set', 'sub-*/eeg/*_eeg.vhdr', 'sub-*/eeg/*_eeg.edf'};
S.id_regexp     = '^(sub-[^_]+)';   % sub-01_task-rest_eeg -> sub-01
S.resample      = 250;              % recordings of different rates: all resampled
S.core_duration = 240;              % the same duration for everyone (the shortest recording)
S.montages      = {'ASA_10-05_343'};   % channels matched by label (10-20, 10-10, 10-05)

%% Bands and measures (defaults: six bands, delta to gamma; wPLI and orthogonalised AEC)
% every band is tested and adjusted for: keep the bands of your question
% S.bands = struct('name', {'theta', 'alpha'}, 'fmin', {4, 7}, 'fmax', {7, 13});
% S.individual_bands = true;      % bands anchored on each recording's alpha frequency (S.bands_individual)
% S.conn  = {'wpli', 'aec_orth', 'plv'};

%% Subject table (BIDS participants.tsv) and the tests
S.subs_info     = '/data/mystudy/participants.tsv';
S.id_column     = 'participant_id';
S.group_column  = 'diagnosis';
S.groups        = {'control', 'patient'};
S.comparisons   = struct('name', 'patient - control', 'A', {{'diagnosis', '==', 'patient'}}, ...
    'B', {{'diagnosis', '==', 'control'}});
S.covariates    = {'age', 'sex'};   % text columns (sex) become indicators
S.correlations  = struct('name', 'Symptoms', 'columns', '^PANSS_', ...
    'subjects', {{'diagnosis', '==', 'patient'}}, 'covariates', {{'age'}}, 'zero_as_missing', false);
end
