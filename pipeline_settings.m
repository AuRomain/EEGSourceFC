function S = pipeline_settings()
%PIPELINE_SETTINGS Paths and parameters of the static connectivity pipeline: the template for a study.
%   S = pipeline_settings() returns the default settings. The fields marked
%   EDIT describe a study (recordings, cap, subject table, groups,
%   comparisons, correlations) and must be set for yours; every other value
%   is the validated default. Everything is relative to this folder, so the
%   folder can be moved or shared.
%
%   A study: write a function that starts from this one and changes only
%   what differs (examples/settings_openneuro.m: the public OpenNeuro
%   datasets ds003944 and ds003947; examples/settings_bids.m: a BIDS dataset), then run
%       S = setup_toolboxes(settings_mystudy());
%       run_all(S);          % stage 1: connectivity of every recording, every band
%       run_statistics(S);   % statistics: NBS, graph metrics, networks, regions, power
%   or both, one after the other, with the PSD analysis between them
%   (section "PSD analysis" below): run_full_pipeline(settings_mystudy()).
%   Each stage called without settings calls this function; a variable S of
%   the workspace is never used. check_settings.m validates the values and
%   fills defaults; a field it does not know (a misspelled name) stops the run.
%   The connectivity of each recording is static (S.dynamic = false): the
%   mean of each measure's windows (the default), or its connectivity over
%   the whole recording (S.static_matrices).

S.root      = fileparts(mfilename('fullpath'));
project     = fileparts(S.root);                       % folder holding Datasets/ and Toolboxes/
S.toolboxes = fullfile(project, 'Toolboxes');
S.results   = fullfile(S.root, 'results');             % outputs of both stages
S.dynamic   = false;                                   % static connectivity only (this folder has no states)

%% Recordings (stage 1)
% The recordings must be cleaned before (artefacts removed or corrected, bad
% channels removed or interpolated): the pipeline only takes the average
% reference, leaves flat channels out and cuts the data at boundary events.
S.data          = fullfile(project, 'Datasets');   % EDIT: folder of the recordings
S.file_pattern  = {'*.set', '*.edf', '*.bdf', '*.vhdr'};   % EDIT: one or more patterns; any format FieldTrip
                                        % reads (.set, .edf, .bdf, .vhdr, .fif, .cnt, .mff, ...); a pattern may
                                        % name subfolders, e.g. 'sub-*/eeg/*_eeg.set' (BIDS) or '**/*.vhdr'
S.id_regexp     = '^([^_]+)';           % EDIT: subject ID = first token of the file name without extension
                                        % (sub-01_task-rest_eeg -> sub-01); whole name when it does not match.
                                        % Each recording needs its own ID (stage 1 stops otherwise): BIDS with
                                        % sessions, '^(sub-[^_]+_ses-[^_]+)' (sub-01_ses-1, one table row per
                                        % session), or one session per run (file_pattern 'sub-*/ses-1/eeg/*.set')
S.srate         = [];           % expected sampling rate (checked), e.g. 512; [] accepts any rate
S.resample      = [];           % Hz: resample every recording (e.g. 256); [] keeps the rate
S.core_start    = 0;            % s: start of the analysed segment
S.core_duration = [];           % s analysed (e.g. 180: the first 3 minutes); [] = whole recording.
                                % Connectivity estimates depend on the amount of data: give every
                                % recording the same duration when the recordings differ in length
S.edge_trim     = 1;            % s discarded at both ends of each continuous piece after filtering
S.filter_transition = 'brainstorm';   % band-pass transition width: 'brainstorm' (Brainstorm's default:
                                % 0.5 Hz at a lower edge at or below 5 Hz, else 1 Hz: sharp bands, long filters, whose
                                % transients are trimmed), 'eeglab' (25% of the lower edge), 'window'
                                % (the narrowest transition whose filter fits two windows of 6 cycles of the
                                % band) or a width in Hz
S.boundary_events = {'boundary', 'New Segment'};  % events marking a discontinuity (EEGLAB rejections, BrainVision pauses)
S.bad_events    = {'bad*'};     % events marking bad data, removed over their duration (case ignored, * = any
                                % text): MNE's BAD_ annotations in EEGLAB/EDF+/BrainVision files, BrainVision
                                % 'Bad Interval'; {} = none. FIF annotations are not read (get_eeg_data)
S.min_segment   = 6;            % s: continuous pieces shorter than this (after trimming) are dropped; 6: the
                                % AEC windows of the window mean (S.aec_window) fit in every piece kept (4 is
                                % enough with S.static_matrices = {'full'})

%% Montage, head model, atlas, inverse (stage 1)
% Montage = electrode positions on the ICBM152 head with a leadfield in inputs/:
%   'ASA_10-05_343'  any cap with 10-20 / 10-10 / 10-05 labels (Fp1, AFz, T7 or T3, ...):
%                    the channels of each recording are matched by label
%   'EEGLAB_61'      a 61-channel ANT Waveguard cap (EEGLAB channel file; positions
%                    corrected: docs/anatomy.md)
%   a file           Brainstorm channel_*.mat or CSV (label,x,y,z in m, ICBM152): the
%                    leadfield is computed once with OpenMEEG (see README)
% Validated with about 60 channels; fewer than 32 give a warning (68 regions from
% so few channels have not been validated), fewer than S.min_channels stop.
S.montages      = {'ASA_10-05_343'};   % EDIT if the cap has no standard labels
S.channel_subset = {};          % labels to keep, e.g. a 32- or 19-channel subset; {} = all
S.min_channels  = 16;           % fewer matched channels stops the run
S.drop_flat_channels = true;    % channels with (near) zero variance are left out (stage1_summary.csv)
S.atlases       = {'Desikan_Killiany_68'};   % or a Brainstorm scout file on tess_cortex_pial_low, or an atlas
                                % of the cortex by name: 'Destrieux', 'DKT', 'Schaefer_100_17net', ...
S.anatomy       = '';           % head model: '' = the shipped ICBM152 (inputs/); a folder made by
                                % tools/make_anatomy.m (Brainstorm) = another head for every recording
                                % (a template: 'Colin27_2016', ...), or one folder per subject (below)
S.anatomy_individual = false;   % true: the folder <S.anatomy>/<subject ID> of each recording (individual
                                % heads from T1 MRIs, CAT12); the atlases are then given by name
S.anatomy_fallback = false;     % true: subjects without a folder use the shipped ICBM152 (stage1_summary.csv,
                                % column anatomy); false: they stop the run
S.inv           = {'eloreta'};  % 'eloreta' (FieldTrip), 'lcmv' (FieldTrip), 'wmne' (+wmne)
S.orientation   = 'normal';     % sources along the cortical normal (one signed filter each); 'free': three
                                % orientations per source and one principal component per region (S.scout_function = 'pca')
S.scout_function = 'mean';      % each region's signal from its sources (Brainstorm's scout functions):
                                % 'mean' (sign-flipped mean), 'pca' (first principal component);
                                % also 'max', 'power', 'rms' (not linear; power and RMS have no phase)
S.eloreta_lambda = 0.05;        % FieldTrip default; relative (5% of the mean eigenvalue it inverts)
S.lcmv_lambda   = '5%';
S.lcmv_kappa    = [];           % LCMV: [] = FieldTrip's rule (the regularised covariance inverted whole);
                                % 'rank' = only its components within the data rank (as MNE-Python's
                                % beamformers); or a number of components
S.wmne          = struct('weightExp', 0.5, 'weightLimit', 10, 'SNR', 3);   % lower SNR = stronger regularisation

%% Connectivity (stage 1)
S.bands         = struct('name', {'delta', 'theta', 'alpha', 'beta1', 'beta2', 'gamma'}, ...
                         'fmin', {1, 4, 7, 13, 20, 30}, 'fmax', {4, 7, 13, 20, 30, 45});
% (every band is analysed and adjusted for in the statistics: keep only those of your question;
% alpha 7-13 Hz: a slower alpha, as in drowsiness, stays in the band)
%% Individual bands (optional): ON/OFF HERE
% false: the fixed bands of S.bands above. true: S.bands_individual is used
% instead: bands whose edges are anchored on each recording's individual
% alpha frequency (IAF, estimated in stage 1 with the settings S.iaf; needs
% restingIAF). An edge is a number (Hz) or 'iaf', 'iaf+2', 'iaf-6',
% 'iaf+bins', 'iaf-bins' (band_edges.m); with IAF = 10 Hz the set below
% gives delta 1-4, theta 4-8, alpha 8-12, beta1 12-20, beta2 20-30 and gamma
% 30-45 Hz (docs/individual_alpha_frequency.md).
S.individual_bands = false;
S.bands_individual = struct('name', {'delta_iaf', 'theta_iaf', 'alpha_iaf', 'beta1_iaf', 'beta2', 'gamma'}, ...
    'fmin', {1, 'iaf-6', 'iaf-bins', 'iaf+bins', 20, 30}, ...
    'fmax', {'iaf-6', 'iaf-bins', 'iaf+bins', 20, 30, 45});
S.iaf = struct( ...                % individual alpha frequency, used only by anchored bands
    'bins', 2, 'bin_hz', 1, ...    % alpha band = IAF +/- bins x bin_hz (IAF 9 Hz -> 7-11 Hz)
    'measure', 'paf', ...          % IAF used: restingIAF's peak ('paf') or centre of gravity ('cog')
    'no_peak', 'exclude', ...      % no IAF: 'exclude' (anchored bands not computed) or 'cog'
    'channels', {{'O1', 'Oz', 'O2', 'PO3', 'POz', 'PO4', 'PO7', 'PO8', 'P3', 'Pz', 'P4'}}, ...
    'search', [7 13], 'cmin', 3, 'sgf_frame', 11, 'sgf_order', 5, ...   % restingIAF's defaults
    'checks', {{'irasa', 'specparam'}}, ...   % FieldTrip checks, reported in stage1_summary.csv
    'max_disagreement', 0.5, ...   % Hz: an IAF further from a check's peak is flagged (kept)
    'edge_check', true);           % a 'peak' on the flank of beta or theta at the search range's
                                   % edge (no IRASA maximum inside the range) is not taken as the IAF
S.conn          = {'wpli', 'aec_orth'};   % measures, each estimated as S.static_matrices below says:
%   'wpli'  weighted phase lag index: phase coupling with a lag (zero-lag coupling, which
%           source leakage produces, is ignored); the measure to read first
%   'aec_orth' amplitude envelope correlation after pairwise orthogonalisation (leakage removed)
%   also: 'plv', 'aec' (both dominated by leakage between neighbouring regions: compare with
%   the others, do not read alone), 'pli', 'ciplv' (corrected imaginary PLV), 'dwpli' (debiased
%   wPLI of Vinck et al. 2011, by FieldTrip, epochs as trials: below), 'aec_symorth', 'plv_symorth'
%   (symmetric orthogonalisation; need fewer regions than channels - 1)
S.static_matrices = {'window_mean'};   % how each measure is estimated (each one analysed on its own):
%   'window_mean'  in the windows below, then averaged
%   'full'         over all the samples of the recording at once (static_connectivity.m), as
%                  Brainstorm, FieldTrip and MNE estimate the connectivity of continuous data;
%   both: {'full', 'window_mean'}.
%   Short windows bias phase measures upwards (wPLI of 0.6 s windows of independent noise: about
%   0.5), the more so where the rhythm is stronger: docs/methods.md, "Whole recording or windows".
%   For phase measures (wPLI, ciPLV, PLV, PLI) 'full' is recommended: their window mean follows
%   band power more than coupling. run_statistics warns once per run; with the window mean, add
%   band power (stage 2) as a covariate
S.window_cycles  = 6;           % windows of the phase measures (window mean): 6 cycles of the band
                                % centre (alpha 7-13 Hz: 0.6 s) ...
S.window_length  = [];          % ... or a fixed length in s for every band (e.g. 2) instead of cycles
S.window_overlap = 0.5;         % overlap of the phase windows (0.5 = 50%)
S.aec_window     = 6;           % s: windows of the AEC measures (window mean; O'Neill et al. 2017) ...
S.aec_step       = 0.5;         % s: ... starting every 0.5 s
S.measure_windows = struct();   % windows of one measure instead of the common ones above, e.g.
                                % struct('wpli', struct('window_length', 2), 'aec_orth', struct('aec_window', 10, 'aec_step', 1));
                                % phase measures: window_cycles, window_length, window_overlap; AEC: aec_window, aec_step.
                                % Several window settings for one measure: separate runs, each with its own S.results
S.dwpli_epoch_cycles = 6;       % dwPLI's epochs (its trials): 6 cycles of the band centre ...
S.dwpli_epoch_length = [];      % ... or a fixed length in s for every band
                                % dwPLI's windows are its own (at least 10 epochs each): 60 cycles of the
                                % band centre every 30 cycles (10 epochs of 6 cycles: 6 s in alpha; 10
                                % epochs when 60 cycles hold fewer), or S.measure_windows.dwpli
S.reliability_block = 10;       % s: connectivity is also computed on the odd and on the even blocks of
                                % this length, for the split-half reliability of every feature and edge
S.save_scouts   = false;        % keep the region time series (results/Connectivity/Scouts time-series/)
S.parallel      = true;         % parfor over recordings when the Parallel Computing Toolbox is present
S.max_workers   = 8;

%% Subjects and groups (statistics)
S.subs_info     = fullfile(S.data, 'participants.tsv');   % EDIT: subject table, .csv, .tsv or .xlsx
S.id_column     = 'participant_id';  % EDIT: column of the subject IDs (as S.id_regexp gives them)
S.group_column  = 'group';      % EDIT: column of the groups; '' = group from the start of the subject ID
S.groups        = {'control', 'patient'};   % EDIT: values of that column, reference group first

%% Statistics (run_statistics: static_analysis.m for every band and measure)
S.run_name      = 'main';       % outputs in results/Statistics/<run_name>/
S.static_measures = {};         % measures analysed; {} = every measure of S.conn
S.static_families = {'global', 'networks', 'regions', 'power'};   % feature families tested:
%   'global'   Strength (mean node strength), ClustCoeff (mean weighted clustering coefficient),
%              CPL (characteristic path length, lengths 1 / weight); the Brain Connectivity
%              Toolbox's definitions, on the positive part of a signed measure such as AEC;
%              dwPLI (centred on 0): Strength only, from its signed weights
%   'networks' mean connectivity inside each network and between each pair (S.network_table)
%   'regions'  node strength of every region (sum of its positive weights; dwPLI: signed)
%   'power'    log10 band power, mean over all regions and over each network
%   (edges: the network-based statistic, S.nbs; spectral profile: S.spectral_profile)
S.network_table = fullfile(S.root, '3- Analysis', 'DK68_networks.csv');   % regions -> networks (Yeo et
                                % al. 2011, 7 networks) and MNI centroids, for DK-68; '' = none
S.n_perm        = 10000;        % smallest p 1/(n_perm + 1): a family can be selected only when it is at
                                % most S.alpha / families (the statistics say how many they have, at the end)
S.alpha         = 0.05;
S.seed          = 20261015;     % every random step
% Group comparisons: subjects of A vs subjects of B, each defined by
% conditions on the subject table {column, operator, value; ...} (all must hold).
% EDIT the comparisons, covariates and correlations.
S.comparisons   = struct('name', {'patient - control'}, ...
    'A', {{'group', '==', 'patient'}}, 'B', {{'group', '==', 'control'}});
%   a subgroup is a further condition, e.g. medicated patients against controls:
%   'medicated - control': A = {'group','==','patient'; 'dose','>',0}, B = {'group','==','control'}
S.covariates    = {'age'};      % group comparisons: ANCOVA with these covariates (columns of the table;
                                % text columns become indicators); {} = Welch t
S.categorical_covariates = {};  % numeric covariates that code categories (e.g. a site coded 1, 2, 3):
                                % indicators, one per value but the smallest, instead of one straight line
% Correlations of every feature with clinical scores, in the subjects
% selected, with covariates partialled out of feature and score. None by
% default; for example, every column whose name starts with PANSS_, in the
% patients, with age partialled out:
%   S.correlations = struct('name', 'Symptoms', 'columns', '^PANSS_', ...
%       'subjects', {{'group', '==', 'patient'}}, 'covariates', {{'age'}}, 'zero_as_missing', false);
S.correlations  = struct('name', {}, 'columns', {}, 'subjects', {}, 'covariates', {}, 'zero_as_missing', {});
S.correlation_type = 'pearson'; % 'pearson', or 'spearman' (feature, score and covariates ranked: for skewed
                                % scores, where a few subjects would drive r); the slope over its HC3
                                % standard error is tested either way (perm_partial_corr_maxT)
S.correlation_test = 'permutation';   % its null: 'permutation', or 'wild_bootstrap' (Rademacher, null
                                % imposed; S.n_perm draws): closer to alpha with fewer than about 100
                                % subjects when the spread of a feature changes with the score
% Extra covariate sets: each re-runs a set of correlations with more covariates.
S.covariate_sets = struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {});
S.nbs           = struct('thresh', 3.1, 'perms', 5000, 'size', 'extent');   % network-based statistic
                                % of the edges, per comparison (Zalesky et al. 2010)
S.enrichment    = struct('n_perm', 10000);   % which networks concentrate a group difference of the
                                % edges, against random and spatially rotated (spin) network labels
S.reliability   = true;         % split-half reliability of every feature and of the edges
S.spectral_profile = true;      % family 'spectral', before the others: each recording's alpha frequency
                                % and peak, 1/f background, muscle and eye indices at the scalp
                                % (needs restingIAF), joined to the subject table so that they can be
                                % covariates (e.g. S.covariates = {'age', 'alpha_pw_post'}; alpha_pw_post
                                % itself is then tested adjusted for age only)
S.figures       = true;         % box and scatter plots of the significant results (drawn for
                                % p_maxT_family <= alpha, kept when significant over the families), group
                                % mean matrices, edge statistics and NBS components

%% PSD analysis (stage 2: run_psd_analysis.m, in "2- PSD analysis")
% Band power of every channel and region, compared between the groups and
% correlated with the scores of the statistics above (S.comparisons,
% S.covariates, S.correlations, S.covariate_sets, S.n_perm, S.alpha), with
% scalp and cortex maps. Its own bands (not S.bands); an edge can be
% anchored on the IAF as in S.bands_individual ('iaf-bins', ...; S.iaf).
% Measures: 'absolute' (log10 band power), 'relative' (% of S.psd.total),
% 'periodic' (log10 power above the 1/f fit), 'aperiodic' (1/f exponent and
% offset of the specparam fit), 'alpha_frequency' (specparam's alpha peak;
% NaN without one; the centre of gravity above the 1/f fit is reported,
% not tested).
S.psd = struct( ...
    'run', true, ...                % run_full_pipeline runs it after stage 1 (run_psd_analysis(S) runs it alone)
    'levels', {{'sensor', 'source'}}, ...   % channels, and regions of S.atlases{1} for every inverse of S.inv
    'bands', struct('name', {'delta', 'theta', 'alpha', 'alpha1', 'alpha2', 'beta', 'beta1', 'beta2'}, ...
                    'fmin', {1, 4, 7, 7, 10, 13, 13, 20}, 'fmax', {4, 7, 13, 10, 13, 30, 20, 30}), ...
    'total', [1 30], ...            % Hz: relative power = % of the power in this range
    'measures', {{'absolute', 'relative', 'periodic', 'aperiodic', 'alpha_frequency'}}, ...
    'segment', 4, 'overlap', 0.5, ...   % Welch: 4 s Hann segments (0.25 Hz resolution), 50% overlap
    'freq_range', [1 45], ...       % Hz: the spectra
    'fit_range', [2 40], ...        % Hz: the specparam fit (aperiodic 'fixed', as the spectral profile above)
    'power_line', '50', ...         % the mains: '50', '60' or 'none'; the fit leaves out 1-3 times it +/- 2 Hz
    ...                             % inside fit_range. Narrow lines anywhere are flagged, not removed
    ...                             % (recordings.csv: line_hz, line_ratio; docs/psd.md, "Line noise")
    'region_spectra', 'signal', ... % a region's spectrum: 'signal', the PSD of its signal (stage 1's: the
    ...                             % sign-flipped mean of its sources); 'vertex_mean', the mean of its
    ...                             % sources' PSDs (levels source_<inverse>_vertexmean). Either way
    ...                             % region_power_check.csv compares the two
    'vertex_maps', true, ...        % maps of every source's band power (descriptive, not tested; the
    ...                             % shipped or one common head only)
    'cluster', struct('run', true, 'alpha', 0.05), ...   % cluster test over neighbouring channels and regions
    'n_perm', []);                  % permutations of its statistics ([] = S.n_perm)
                                    % Cluster test: cluster-forming threshold two-sided p < 0.05 of each one's t.
                                    % n_perm: the PSD has many families (levels x measures x bands x comparisons
                                    % and correlation sets), and a family can be selected only when
                                    % 1/(n_perm + 1) <= S.alpha / families: about 10 x families / S.alpha
                                    % leaves room (the first line of its statistics says how many)
end
