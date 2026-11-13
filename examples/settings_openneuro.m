function S = settings_openneuro(data, subject_table)
%SETTINGS_OPENNEURO Static connectivity of the OpenNeuro first-episode psychosis example.
%   S = settings_openneuro(data, subject_table): the defaults of
%   pipeline_settings.m and what differs for these data.
%   Data: two public datasets of eyes-open resting EEG (Salisbury, Seebold,
%   Coffman; CC0), 60-channel 10-10 caps recorded with an Elekta Neuromag
%   system at 1000 Hz, about 5 minutes each:
%     D1  ds003944, "EEG: First Episode Psychosis vs. Control Resting Task 1"
%         (doi:10.18112/openneuro.ds003944.v1.0.0)
%     D2  ds003947, "EEG: First Episode Psychosis vs. Control Resting Task 2"
%         (doi:10.18112/openneuro.ds003947.v1.0.0)
%   data           folder of the recordings, cleaned by the user beforehand
%                  (OpenNeuro distributes the raw BrainVision files, and
%                  the pipeline does not clean data), saved as EEGLAB .set
%                  named sub-<ID>_task-Rest_eeg_ready.set, directly in it or
%                  one subfolder down (default: <project>/Datasets_OpenNeuro).
%   subject_table  the user's own subject table for these recordings (the
%                  OpenNeuro participant and phenotype files are not read),
%                  one row per recording (default: subjects.csv in data),
%                  with the columns participant_id, type (Control or
%                  Psychosis), age, gender, dataset (D1 or D2), CPZ_at_scan
%                  (chlorpromazine-equivalent dose; 0 = unmedicated), the
%                  SAPS and SANS global items (SAPS_..., SANS_...), the positive, negative and
%                  disorganisation dimensions (POS_SYMPT, NEG_SYMPT,
%                  DISORGANIZATION) and the MATRICS domain T-scores
%                  (SPEEDTSCR ... OVERALLTSCR).
%   The analysis: every channel at its template 10-05 position (the two
%   caps share 49 channels); the first 300 s of every recording; the six
%   bands, wPLI and orthogonalised AEC (the defaults);
%   Psychosis - Control and, in the patients, Medicated - Unmedicated,
%   adjusted for age, sex and dataset; correlations in the patients with
%   symptoms and cognition, age partialled out, then also CPZ. Outputs in
%   results_openneuro/.
%       addpath('examples');  S = setup_toolboxes(settings_openneuro('/data/openneuro_fep'));
%       run_all(S);  run_statistics(S);
S = pipeline_settings();
project = fileparts(S.root);
if nargin < 1 || isempty(data), data = fullfile(project, 'Datasets_OpenNeuro'); end
if nargin < 2 || isempty(subject_table), subject_table = fullfile(data, 'subjects.csv'); end
S.results = fullfile(S.root, 'results_openneuro');

%% Recordings
S.data          = data;
S.file_pattern  = {'*_eeg_ready.set', '*/*_eeg_ready.set'};   % in data or one subfolder down
S.id_regexp     = '^(sub-[0-9]+A?)_';
S.srate         = 1000;
S.core_duration = 300;          % the same duration for every recording
S.montages      = {'ASA_10-05_343'};

%% Subject table, comparisons, correlations
S.subs_info     = subject_table;
S.id_column     = 'participant_id';
S.group_column  = 'type';
S.groups        = {'Control', 'Psychosis'};
S.comparisons   = struct('name', {'Psychosis - Control', 'Medicated - Unmedicated'}, ...
    'A', {{'type', '==', 'Psychosis'}, {'type', '==', 'Psychosis'; 'CPZ_at_scan', '>', 0}}, ...
    'B', {{'type', '==', 'Control'}, {'type', '==', 'Psychosis'; 'CPZ_at_scan', '==', 0}});
S.covariates    = {'age', 'gender', 'dataset'};
S.correlations  = struct( ...
    'name',       {'Symptoms', 'Cognition'}, ...
    'columns',    {'^(SAPS_|SANS_)|^(POS_SYMPT|NEG_SYMPT|DISORGANIZATION)$', ...
                   '^(SPEED|ATT_VIG|WM|VERB|VIS|RPS|SOCCOG|OVERALL)TSCR$'}, ...
    'subjects',   {{'type', '==', 'Psychosis'}, {'type', '==', 'Psychosis'}}, ...
    'covariates', {{'age'}, {'age'}}, ...
    'zero_as_missing', {false, false});
S.covariate_sets = struct('name', {'Symptoms, CPZ adjusted', 'Cognition, CPZ adjusted'}, ...
    'correlations', {'Symptoms', 'Cognition'}, ...
    'covariates', {{'age', 'CPZ_at_scan'}, {'age', 'CPZ_at_scan'}}, ...
    'columns', {'', ''});
S.n_perm        = 50000;        % many families (6 bands x 2 measures x feature families x 6 comparisons
                                % and correlation sets, and more in the PSD statistics): with 50,000
                                % permutations a family is selected alone when a few permutations reach its
                                % largest statistic (at 10,000: only when none does; docs/validation.md,
                                % "Permutation resolution")
end
