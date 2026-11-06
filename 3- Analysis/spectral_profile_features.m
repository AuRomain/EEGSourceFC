function [P, tested] = spectral_profile_features(S)
%SPECTRAL_PROFILE_FEATURES Alpha frequency, alpha peak, 1/f background and artefact indices of every recording.
%   [P, tested] = spectral_profile_features(S) (stage 3, S.spectral_profile)
%   reads every recording of S.data as stage 1 does (get_eeg_data: the
%   analysed stretch and its continuous pieces; the channels of the first
%   cap of S.montages, flat channels left out, average reference) and
%   measures at the scalp what the connectivity of a band depends on besides
%   coupling: how strong and how fast the alpha rhythm is, the 1/f
%   background, and muscle and eye artefacts. Each recording:
%   - its IAF (individual_alpha_frequency with S.iaf, as stage 1: restingIAF
%     peak and centre of gravity on the posterior channels, IRASA and
%     specparam checks);
%   - each channel's power spectrum, 1-45 Hz in 0.25 Hz steps (Hann taper,
%     4 s pieces inside the continuous pieces, FieldTrip mtmfft);
%   - each channel's specparam fit, 2-40 Hz (FieldTrip's 'fooof_aperiodic',
%     Brainstorm's FOOOF with its defaults, aperiodic 'fixed'): offset,
%     exponent, the strongest peak centred in the alpha search range
%     S.iaf.search, R^2, and the aperiodic fit.
%   These are kept in results/Connectivity/spectral_profile/<recording>.mat
%   with a key (the recording's files by content, the reading settings,
%   S.iaf and the code): a recording is computed again only when one changes.
%   Measures (spectral_profile_summary; alpha = S.iaf.search, 7-13 Hz by
%   default; channel sets: those of each set present in every recording):
%   IAF, IAF_PAF, IAF_CoG, IAF_IRASA, IAF_specparam, IAF_source, IAF_flag
%                      as iaf.csv of stage 1; IAF_cov the IAF, or the CoG
%                      when there is no peak
%   alpha_pw_post      height of the alpha peak above the 1/f fit (log10
%                      power), median of the posterior channels (S.iaf.channels;
%                      0 where a channel has none): how strong the rhythm is
%   alpha_peak_share   share of posterior channels with an alpha peak
%   alpha_cf_post      median centre of those peaks
%   alpha_front_post   the peak height at frontal (F3 F1 Fz F2 F4) minus
%                      posterior channels (alpha moving forward: drowsiness)
%   exponent_all, offset_all, fit_r2   1/f exponent, offset and the fit's R^2,
%                      median of all channels; exponent_post of the posterior ones
%   log_alpha_post     log10 alpha power, posterior; rel_alpha_post its share
%                      of 1-45 Hz; theta_alpha_post log10 of 4 Hz to the alpha
%                      range over the alpha range
%   purity_fixed       share of the alpha power (posterior) above the 1/f fit:
%                      how much of what a fixed alpha band holds is rhythm;
%                      purity_iaf the same for IAF +/- S.iaf.bins x S.iaf.bin_hz
%                      (the alpha band anchored on the IAF)
%   capture_fixed, capture_iaf   share of the rhythmic power of 5-15 Hz in
%                      the alpha range, or in the IAF band
%   hf_temporal        log10 share of 30-45 Hz in 1-45 Hz at temporal channels
%                      (FT7 FT8 T7 T8 TP7 TP8 FT9 FT10 T9 T10 TP9 TP10): muscle
%   eye_frontopolar    log10 of 1-3 Hz power at Fp1 Fpz Fp2 over posterior: eyes
%   A cap without these 10-05 names (an ANT equidistant or EGI cap) uses the
%   channels nearest to their positions (cap_channels; S.iaf.channels is
%   matched so by setup_toolboxes). A set with no channel in every recording gives
%   NaN; so does a recording that fails (e.g. no continuous piece of 4 s;
%   a warning names it, and IAF_flag is 'failed').
%   Writes results/Connectivity/spectral_profile.csv (recording, ID, the
%   measures), which read_subject_table joins to the subject table by ID,
%   so that any measure can be a covariate (S.covariates, S.covariate_sets).
%   tested: the measures stage 3 tests as the family 'spectral' (the IAF
%   variants, specparam's alpha centre and R^2 and the 1/f exponent of the
%   posterior channels are left out as near-duplicates or fit checks; IAF is
%   left out when the family IAF tests it, S.iaf_used). On the 138 OpenNeuro
%   recordings, about 6 min on 14 workers (decision of 28 October 2026).
path_to_save = fullfile(S.results, 'Connectivity');
cache_dir = fullfile(path_to_save, 'spectral_profile');
if ~isfolder(cache_dir), mkdir(cache_dir); end
subs = get_subs(S.data, S.file_pattern);
n = numel(subs);
ftElec = get_ftElec(S.montages{1});
% what the spectra depend on: the reading settings (get_eeg_data,
% get_eeg_montage), S.iaf and the code
V = struct();
for f = {'srate', 'resample', 'core_start', 'core_duration', 'boundary_events', 'drop_flat_channels', 'channel_subset'}
    if isfield(S, f{1}), V.(f{1}) = S.(f{1}); end
end
if isfield(S, 'bad_events') && ~isequal(S.bad_events, {'bad*'}), V.bad_events = S.bad_events; end   % as settings_key
V.montage = S.montages{1};
if isfile(V.montage), V.montage_file = file_fingerprint({V.montage}); end
V.iaf = S.iaf;
V.code = file_fingerprint(cellfun(@which, {'spectral_profile_features', 'individual_alpha_frequency', ...
    'get_eeg_data', 'get_eeg_montage', 'get_ftElec', 'restingIAF', 'ft_freqanalysis'}, 'UniformOutput', false));
keys = cell(n, 1);
todo = false(n, 1);
for s = 1:n
    [folder, stem] = fileparts(fullfile(S.data, subs{s}));
    d = dir(fullfile(folder, [stem '.*']));               % the file and its companions (.fdt, ...), as settings_key
    d = d(~[d.isdir]);
    keys{s} = hash_text(jsonencode(struct('settings', V, 'input', file_fingerprint(fullfile(folder, sort({d.name}))))));
    c = fullfile(cache_dir, [stem '.mat']);
    todo(s) = ~isfile(c) || ~strcmp(getfield(load(c, 'key'), 'key'), keys{s});
end
fprintf('Spectral profile: %d recordings, %d to compute\n', n, nnz(todo));
t0 = tic;
files = fullfile(S.data, subs(todo));
keys = keys(todo);
if S.parallel && numel(files) > 12 && ~isempty(ver('parallel'))
    if isempty(gcp('nocreate')), parpool('Processes', min(S.max_workers, numel(files))); end
    parfor i = 1:numel(files)
        one_recording(files{i}, S, ftElec, cache_dir, keys{i});
    end
else
    for i = 1:numel(files)
        one_recording(files{i}, S, ftElec, cache_dir, keys{i});
    end
end
if any(todo), fprintf('Spectral profile: %d recordings computed (%.1f min)\n', nnz(todo), toc(t0) / 60); end
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
R = cell(n, 1);
for s = 1:n
    w = load(fullfile(cache_dir, [stems{s} '.mat']), 'r');
    R{s} = w.r;
end
bad = cellfun(@(r) isfield(r, 'error'), R);
if all(bad), error('spectral_profile:failed', 'No recording has a spectral profile: %s', R{1}.error); end
if any(bad)
    msg = string(stems(bad)) + ": " + string(cellfun(@(r) char(r.error), R(bad), 'UniformOutput', false));
    warning('spectral_profile:failed', '%d recording(s) without a spectral profile (missing values): %s', ...
        nnz(bad), strjoin(msg, '; '));
end
[M, sets] = spectral_profile_summary([R{~bad}], S.iaf, ftElec);   % a cap with other names: sets by position
rows = cumsum(~bad);
M = M(max(rows, 1), :);                                     % a row for every recording ...
M{bad, varfun(@isnumeric, M, 'OutputFormat', 'uniform')} = NaN;   % ... missing for those that failed
M.IAF_source(bad) = "";  M.IAF_flag(bad) = "failed";
for set = string(fieldnames(sets))'
    fprintf('  %-5s %s\n', set, strjoin(string(sets.(set)), ' '));
end
P = [table(string(stems), subject_ids(S, stems), 'VariableNames', {'recording', 'ID'}), M];
writetable(P, fullfile(path_to_save, 'spectral_profile.csv'));
tested = {'IAF', 'IAF_CoG', 'alpha_pw_post', 'alpha_peak_share', 'alpha_front_post', 'exponent_all', ...
    'offset_all', 'log_alpha_post', 'rel_alpha_post', 'theta_alpha_post', 'purity_fixed', 'purity_iaf', ...
    'capture_fixed', 'hf_temporal', 'eye_frontopolar'};
if S.iaf_used, tested = setdiff(tested, {'IAF'}, 'stable'); end
end

function one_recording(file, S, ftElec, cache_dir, key)
% the spectra, specparam fits and IAF of one recording, saved with its key;
% a recording that fails is saved as its error (the same input and code fail
% again), and gets missing values
[~, stem] = fileparts(file);
try
    r = measure(file, S, ftElec);
catch err
    fprintf(2, '%s: no spectral profile: %s\n', stem, err.message);
    r = struct('error', string(err.message));
end
save(fullfile(cache_dir, [stem '.mat']), 'r', 'key');
end

function r = measure(file, S, ftElec)
eeg = get_eeg_data(file, S, ftElec.label);
[X, labels] = get_eeg_montage(eeg.data, eeg.label, ftElec.label, S.channel_subset);
sd = std(double(X), 0, 2);
flat = sd <= 1e-6 * median(sd);
if any(flat) && S.drop_flat_channels
    X = X(~flat, :);  labels = labels(~flat);
end
X = X - mean(X, 1);                                   % average reference, as stage 1
fs = eeg.srate;
iaf = individual_alpha_frequency(X, labels, fs, eeg.segments, S.iaf);
% 4 s pieces inside the continuous pieces (as the IAF checks)
n4 = round(4 * fs);
keep = eeg.segments(eeg.segments(:, 2) - eeg.segments(:, 1) + 1 >= n4, :);
data = struct('label', {cellstr(labels)}, 'fsample', fs, 'trial', {{}}, 'time', {{}}, 'sampleinfo', zeros(0, 2));
x = double(X);
for k = 1:size(keep, 1)
    for t0 = keep(k, 1):n4:(keep(k, 2) - n4 + 1)
        data.trial{end + 1} = x(:, t0:t0 + n4 - 1);
        data.time{end + 1} = (0:n4 - 1) / fs;
        data.sampleinfo(end + 1, :) = [t0, t0 + n4 - 1];
    end
end
if isempty(data.trial)
    error('spectral_profile:pieces', '%s: no continuous piece of 4 s', file);
end
cfg = struct('method', 'mtmfft', 'taper', 'hanning', 'foilim', [1 45], 'pad', 4, ...
    'output', 'pow', 'keeptrials', 'no', 'feedback', 'no');
P = ft_freqanalysis(cfg, data);
cfg.foilim = [2 40];
cfg.output = 'fooof_aperiodic';
A = ft_freqanalysis(cfg, data);
nc = numel(labels);
fp = struct('offset', nan(nc, 1), 'exponent', nan(nc, 1), 'alpha_cf', nan(nc, 1), 'alpha_pw', nan(nc, 1), ...
    'alpha_bw', nan(nc, 1), 'r2', nan(nc, 1));
for c = 1:nc
    q = A.fooofparams(c);
    fp.offset(c) = q.aperiodic_params(1);  fp.exponent(c) = q.aperiodic_params(end);
    fp.r2(c) = q.r_squared;
    pk = q.peak_params;
    if ~isempty(pk), pk = pk(pk(:, 1) >= S.iaf.search(1) & pk(:, 1) <= S.iaf.search(2), :); end
    if ~isempty(pk)
        [~, j] = max(pk(:, 2));
        fp.alpha_cf(c) = pk(j, 1);  fp.alpha_pw(c) = pk(j, 2);  fp.alpha_bw(c) = pk(j, 3);
    end
end
r = struct('labels', {cellstr(labels)}, 'n_pieces', numel(data.trial), 'iaf', iaf, 'freq', P.freq, ...
    'pow', single(P.powspctrm), 'freq_fit', A.freq, 'aperiodic', single(A.powspctrm), 'fooof', fp);
end
