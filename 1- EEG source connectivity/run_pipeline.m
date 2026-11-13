function qc = run_pipeline(sub, data_path, path_to_save, S)
%RUN_PIPELINE Stage 1 for one recording: EEG -> ROI time series -> connectivity matrices.
%   qc = run_pipeline(sub, data_path, path_to_save, S) with sub the file name
%   (e.g. sub-01_eeg.set), data_path its folder, path_to_save the
%   Connectivity results folder and S from pipeline_settings. Loops over
%   montages -> atlases -> inverse solutions -> bands; outputs:
%     cmats/<sub>/cmats_<inv>_<montage>_<atlas>_<band>.mat         cmats{k}: S.conn{k} (PLV, wPLI, AEC)
%     cmats/<sub>/static_cmats_<inv>_<montage>_<atlas>_<band>.mat  window mean and whole recording
%     Scouts time-series/<sub>/Scouts_<inv>_<montage>_<atlas>.mat
%   static_cmats_*.mat holds static_cmats (the window mean), static_full_cmats
%   (whole recording), static_power (band power of each region),
%   static_halves and static_power_halves (the same on the odd and on the
%   even blocks of S.reliability_block s: split_halves, for the split-half
%   reliability of the statistics) and info. Static only (S.dynamic = false): no
%   windows, so no cmats_*.mat, and static_cmats is []; unless the window
%   mean is asked for (S.static_matrices with 'window_mean'): the windows
%   are computed, their mean saved as static_cmats, the windows themselves
%   not saved (no cmats_*.mat). With windows, static_halves_window is the
%   window mean on each half, for the split-half reliability of the window
%   mean: the windows that lie wholly inside the half's blocks, in blocks
%   of max(S.reliability_block, 3 x the measure's window) per measure
%   (static_halves_block, s), so that at least about 2/3 of the windows
%   are used and the halves share no data (blocks of S.reliability_block s
%   alone would hold few or none of a long window, such as AEC's or
%   dwPLI's).
%   qc: one row of the stage 1 summary (run_all writes them to stage1_summary.csv).
%   Its data_rank is the number of independent signals in the channels used
%   (data_rank, below): ICA components removed, or interpolated channels,
%   lower it; below S.min_channels a warning. Reused outputs without it
%   (made by equivalent code, equivalent_code.txt) get it once: the
%   recording is read again, nothing is recomputed. Its anatomy is the head model used
%   (get_anatomy: the shipped ICBM152, the folder S.anatomy, or the
%   subject's folder with S.anatomy_individual; info.anatomy in the
%   outputs): another head is part of the outputs' input key, so outputs
%   made with another head are recomputed.
%   Its analysed_s, n_pieces and n_windows are those of the band with the
%   least data (the longest filter transient trims most at each cut; each
%   output's info.analysed_s and info.windows give its own band's);
%   nan_removed_s, flat_removed_s
%   and bad_removed_s the seconds of the analysed segment removed as
%   non-finite, flat or marked bad by events (get_eeg_data), unit the
%   unit of the EEG channels in the file (converted to µV), measures_nan
%   the measures that could not be computed for this recording, NaN in its
%   outputs (the symmetric orthogonalisation with fewer independent region
%   signals than regions; a warning names the recording and the reason;
%   the other measures are computed).
%   A band anchored on the IAF that this recording cannot have (no IAF, or
%   a band too narrow for its IAF) is left out (bands_skipped) and counts
%   as done when the IAF of the same file and S.iaf is stored, so the
%   recording is not recomputed at every run.
%   An error is reported in qc (status "failed", message) and printed, so
%   one bad recording does not stop the others; S.stop_on_error = true
%   rethrows it instead.
%   Further: any file format and sampling rate, channels
%   matched by label, discontinuities kept apart (get_eeg_data,
%   get_eeg_montage); corrected electrode positions; fixed normal orientation
%   (one signed filter per source); ROI mean folded
%   into the inverse operator (S.scout_function: also PCA, or max,
%   power and RMS of the source time series); average reference on data and
%   leadfield, over the channels used; filter edges discarded; one
%   Hilbert transform per continuous piece; zero diagonal; signed AEC;
%   bands, windows and measures from the settings.

[~, sub_name] = fileparts(sub);
qc = struct('file', string(sub), 'status', "ok", 'message', "", 'srate', NaN, 'analysed_s', NaN, ...
    'n_pieces', NaN, 'n_boundaries', NaN, 'nan_removed_s', NaN, 'flat_removed_s', NaN, 'bad_removed_s', NaN, ...
    'unit', "", 'n_channels', NaN, 'data_rank', NaN, 'anatomy', "", 'channels_left_out', "", ...
    'subset_missing', "", 'flat_channels', "", 'sd_ratio', NaN, 'n_windows', NaN, ...
    'iaf', NaN, 'iaf_paf', NaN, 'iaf_cog', NaN, 'iaf_irasa', NaN, 'iaf_specparam', NaN, ...
    'iaf_channels', NaN, 'iaf_source', "", 'iaf_flag', "", 'bands_skipped', "", 'measures_nan', "");

try
    qc = process(sub, sub_name, data_path, path_to_save, S, qc);
catch err
    fprintf(2, '%s failed:\n%s\n', sub, err.getReport('extended', 'hyperlinks', 'off'));
    qc.status = "failed";
    qc.message = string(err.message);
    if isfield(S, 'stop_on_error') && S.stop_on_error, rethrow(err); end
end
end


function qc = process(sub, sub_name, data_path, path_to_save, S, qc)
% make directories for connectivity matrices and scout time series
currConditionFolder = fullfile(path_to_save, 'cmats', sub_name);
if exist(currConditionFolder, 'dir') ~= 7, mkdir(currConditionFolder); end
currConditionFolderBis = fullfile(path_to_save, 'Scouts time-series', sub_name);
if S.save_scouts && exist(currConditionFolderBis, 'dir') ~= 7, mkdir(currConditionFolderBis); end   % (none left empty)

for m = 1:numel(S.montages)
    montage = S.montages{m};
    for a = 1:numel(S.atlases)
        atlas = S.atlases{a};
        tag_ma = [file_tag(montage) '_' file_tag(atlas)];
        out = @(inv, band) fullfile(currConditionFolder, sprintf('cmats_%s_%s_%s.mat', inv, tag_ma, band));
        if ~S.dynamic                                 % static only: no windows, no cmats_*.mat
            out = @(inv, band) fullfile(currConditionFolder, sprintf('static_cmats_%s_%s_%s.mat', inv, tag_ma, band));
        end
        % skip only outputs made from the same file with the same settings
        % (the key of each output: its inverse method and band)
        S_ma = setfield(setfield(S, 'montages', {montage}), 'atlases', {atlas}); %#ok<SFLD>
        [~, input_key] = settings_key(S_ma, fullfile(data_path, sub));   % with the head's files (S.anatomy)
        A = get_anatomy(S_ma, sub);                   % the head: shipped ICBM152, S.anatomy or the subject's
        qc.anatomy = string(A.name);
        key = @(inv, band) settings_key(S_ma, '', inv, band);
        accepted = @(inv, band) accepted_keys(S_ma, inv, band);   % key, or that of equivalent code
        % individual alpha frequency of this recording (only for bands
        % anchored on it): a recording without one has no output for those
        % bands, which is up to date when its IAF file is
        iaf_file = fullfile(currConditionFolder, ['iaf_' tag_ma '.mat']);
        anchored = arrayfun(@(b) is_anchored(S.bands(b), S.iaf), 1:numel(S.bands));
        % bands this recording cannot have (no IAF, or too narrow for its
        % IAF): done when its IAF of the same file and S.iaf is stored (else a
        % band too narrow would be recomputed at every run)
        left_out = false(1, numel(S.bands));
        if any(anchored) && isfile(iaf_file)
            w = load(iaf_file);
            if strcmp(w.input_key, input_key) && isequal(w.opts, S.iaf) && isfield(w, 'method') && ...
                    strcmp(w.method, 'pieces')            % IAF of each piece's segments
                left_out = anchored & arrayfun(@(b) band_left_out(S.bands(b), w.iaf.value, S.iaf), 1:numel(S.bands));
            end
        end
        [I, B] = ndgrid(1:numel(S.inv), 1:numel(S.bands));
        done = arrayfun(@(i, b) same_key(out(S.inv{i}, S.bands(b).name), accepted(S.inv{i}, S.bands(b).name), input_key) ...
            || left_out(b), I(:), B(:));
        if S.save_scouts && ~all(left_out)            % the region time series too (else a run without
            scouts = fullfile(currConditionFolderBis, strcat('Scouts_', S.inv, ['_' tag_ma '.mat']));
            done(end + 1) = all(cellfun(@isfile, scouts));   % them would leave none for later runs)
        end
        if all(done)
            fprintf('%s: outputs exist with the same file and settings, skipped\n', sub_name);
            made = find(arrayfun(@(i, b) isfile(out(S.inv{i}, S.bands(b).name)), I(:), B(:)), 1);
            if isempty(made)
                w = load(iaf_file, 'qc');  stored = w.qc;     % every band anchored and left out
            else
                w = load(out(S.inv{I(made)}, S.bands(B(made)).name), 'info');
                stored = struct();
                if isfield(w.info, 'qc'), stored = w.info.qc; end
            end
            for f = intersect(fieldnames(qc), fieldnames(stored))', qc.(f{1}) = stored.(f{1}); end   % the summary row of the run that made them
            qc.status = "skipped (up to date)";
            qc.anatomy = string(A.name);              % this run's head (a fallback subject reuses the shipped head's outputs)
            if ~isfield(stored, 'data_rank')
                % outputs without the data rank (equivalent code): the
                % recording is read once more for it, and it is kept with them
                try
                    X = channels_used(fullfile(data_path, sub), get_ftElec(montage), S);
                    qc.data_rank = data_rank(X);
                    warn_rank(sub_name, qc.data_rank, size(X, 1), S.min_channels);
                    files = arrayfun(@(i, b) out(S.inv{i}, S.bands(b).name), I(:), B(:), 'UniformOutput', false);
                    for f = [files(isfile(files)); {iaf_file}]'
                        add_rank(f{1}, qc.data_rank);
                    end
                    fprintf('%s: data rank %d of %d channels added to the outputs\n', sub_name, qc.data_rank, size(X, 1));
                catch err
                    warning('run_pipeline:data_rank', '%s: data rank not computed: %s', sub_name, err.message);
                end
            end
            continue
        end

        %% electrodes, headmodel and leadfield (FieldTrip format)
        ftElec = get_ftElec(montage);
        [ftHeadmodel, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(montage, atlas, S, A);
        C = load(A.cortex, 'VertNormals');
        VertNormals = C.VertNormals;

        %% load the EEG data, keep the montage's channels, without the flat ones
        [X, labels, idx, rep, eeg, flat_labels, sd] = channels_used(fullfile(data_path, sub), ftElec, S);
        % sd_ratio, the largest channel SD over the median, flags a channel
        % much noisier than the rest
        qc.flat_channels = strjoin(string(flat_labels), ' ');
        qc.sd_ratio = max(sd) / median(sd);
        if ~isempty(flat_labels) && S.drop_flat_channels
            fprintf('%s: flat channel(s) left out: %s\n', sub, qc.flat_channels);
        end
        if numel(labels) < S.min_channels
            error('run_pipeline:channels', ['%s: %d channels match montage %s (minimum %d). Channels of the file: %s.\n' ...
                'Labels that are not 10-20/10-10/10-05 names (EGI E1..E256, BioSemi A1..D32, ...) need the cap''s ' ...
                'own positions: set S.montages to its channel file, a Brainstorm channel_*.mat on the ICBM152 head ' ...
                '(Brainstorm has many in %s: channel_GSN_HydroCel_128_E1.mat, channel_BioSemi_128_A1.mat, ...) ' ...
                'or a CSV label,x,y,z in metres on the ICBM152 head; its leadfield is computed once (docs/settings.md, Caps)'], ...
                sub, numel(labels), montage, S.min_channels, strjoin(eeg.label, ', '), ...
                fullfile(S.toolboxes, 'brainstorm3', 'defaults', 'eeg', 'ICBM152'));
        end
        % outside what the pipeline was validated on (60-61 cleaned channels,
        % no channel far noisier than the others): a warning
        if numel(labels) < 32
            warning('run_pipeline:few_channels', ['%s: %d channels. The pipeline was validated with 60-61 ' ...
                'channels; regions estimated from fewer than 32 channels have not been validated'], sub_name, numel(labels));
        end
        if qc.sd_ratio > 10
            warning('run_pipeline:noisy_channel', ['%s: the noisiest channel has %.0f times the median channel SD: ' ...
                'is the recording cleaned (bad channels removed or interpolated)?'], sub_name, qc.sd_ratio);
        end
        qc.data_rank = data_rank(X);
        warn_rank(sub_name, qc.data_rank, numel(labels), S.min_channels);
        [ftElec, ftLeadfield] = keep_channels(ftElec, ftLeadfield, idx);
        X = X - mean(X, 1);                                   % average reference
        srate = eeg.srate;

        %% individual alpha frequency, when a band is anchored on it (optional)
        iaf = struct('value', NaN);
        if any(anchored)
            iaf = individual_alpha_frequency(X, labels, srate, eeg.segments, S.iaf);
            qc.iaf = iaf.value;  qc.iaf_paf = iaf.paf;  qc.iaf_cog = iaf.cog;
            qc.iaf_irasa = iaf.irasa;  qc.iaf_specparam = iaf.specparam;
            qc.iaf_channels = iaf.n_paf;  qc.iaf_source = iaf.source;  qc.iaf_flag = iaf.flag;
            fprintf('%s: IAF %.2f Hz (%s; PAF %.2f, CoG %.2f, IRASA %.2f, specparam %.2f)%s\n', sub_name, ...
                iaf.value, iaf.source, iaf.paf, iaf.cog, iaf.irasa, iaf.specparam, ...
                char(" " + iaf.flag));
            opts = S.iaf;  method = 'pieces'; %#ok<NASGU>
            save(iaf_file, 'iaf', 'input_key', 'opts', 'method', 'qc');
        end
        skipped = strings(0, 1);
        nan_measures = strings(0, 1);
        qc.analysed_s = NaN;  qc.n_pieces = NaN;  qc.n_windows = NaN;
        noiseCov = eye(numel(labels));                        % identity: no noise recording

        %% solving the inverse problem: one filter row per source and method
        filters = get_inverse_solution(X, srate, noiseCov, ftHeadmodel, ftLeadfield, ftElec, VertNormals, S.inv, S);

        Xc = double(X) - mean(double(X), 2);
        Cx = Xc * Xc' / size(X, 2);                           % channel covariance (PCA, free orientation)
        clear Xc
        for iv = 1:numel(S.inv)
            % ROI operator: Brainstorm scout mean with sign flip (or PCA)
            % applied to the filters, then to the data (same result as on
            % the 15,002 source time series, 200 times less memory);
            % max, power and RMS are not linear: computed from the
            % source time series, region by region (get_scouts_timeseries)
            if ismember(S.scout_function, {'mean', 'pca'})
                [scout_filters, roi_labels] = get_scouts_timeseries(squeeze(filters(iv, :, :)), atlas, Cx, S.scout_function, [], A);
                src_data = scout_filters * X;
            else
                [src_data, roi_labels] = get_scouts_timeseries(squeeze(filters(iv, :, :)), atlas, Cx, S.scout_function, X, A);
                scout_filters = [];
            end
            if S.save_scouts
                segments = eeg.segments; onset_s = eeg.onset_s; channels = labels; %#ok<NASGU>
                save(fullfile(currConditionFolderBis, ['Scouts_' S.inv{iv} '_' tag_ma]), ...
                    'src_data', 'scout_filters', 'srate', 'roi_labels', 'channels', 'segments', 'onset_s');
            end

            %% connectivity in every band
            for b = 1:numel(S.bands)
                % fixed edges, or edges anchored on this recording's IAF (band_edges)
                [fmin, fmax] = band_edges(S.bands(b), iaf.value, S.iaf);
                if anchored(b) && band_left_out(S.bands(b), iaf.value, S.iaf)
                    % no IAF, or a band too narrow for this IAF (e.g. delta below 8 Hz):
                    % left out, and outputs of an earlier run for this band removed
                    % (either edge may be the anchored one: delta_iaf is 1 Hz to
                    % IAF - 6; band_left_out checks both)
                    skipped(end + 1) = string(S.bands(b).name); %#ok<AGROW>
                    tag = [S.inv{iv} '_' tag_ma '_' S.bands(b).name];
                    for old = {['cmats_' tag '.mat'], ['static_cmats_' tag '.mat']}
                        if isfile(fullfile(currConditionFolder, old{1})), delete(fullfile(currConditionFolder, old{1})); end
                    end
                    continue
                end
                if anchored(b) && fmax - fmin < 1.5
                    warning('run_pipeline:narrow_band', '%s: band %s is only %.2f Hz wide (%.2f-%.2f Hz, IAF %.2f Hz)', ...
                        sub, S.bands(b).name, fmax - fmin, fmin, fmax, iaf.value);
                end
                if fmax >= srate / 2
                    error('run_pipeline:band', '%s: band %s (%g-%g Hz) reaches the Nyquist frequency (%g Hz)', ...
                        sub, S.bands(b).name, fmin, fmax, srate / 2);
                end
                % an anchored band: windows, epochs and filter of the band at IAF
                % 10 Hz, the same for every recording (cycles_centre, bandpass_pieces)
                Sb = S;  Sb.band_ref = [];
                if anchored(b)
                    [r1, r2] = band_edges(S.bands(b), 10, S.iaf);
                    Sb.band_ref = [r1, r2];
                end
                [src_data_filtered, seg, onset, tran] = bandpass_pieces(src_data, srate, fmin, fmax, eeg.segments, eeg.onset_s, Sb);
                [cmats, static_cmats, static_full_cmats, info, static_power] = ...
                    get_all_connectivity_mats(src_data_filtered, srate, fmin, fmax, seg, Sb); %#ok<ASGLU>
                for mn = fieldnames(info.not_computed)'     % this measure only: NaN for this recording
                    warning('run_pipeline:not_computed', '%s, band %s, %s: %s not computed (NaN): %s', ...
                        sub_name, S.bands(b).name, S.inv{iv}, mn{1}, info.not_computed.(mn{1}));
                    nan_measures(end + 1) = string(mn{1}); %#ok<AGROW>
                end
                if isfield(info, 'no_windows') && isfield(S, 'static_matrices') && ...
                        any(strcmp(cellstr(S.static_matrices), 'window_mean'))
                    for mn = fieldnames(info.no_windows)'  % its window mean is NaN (the whole recording is computed)
                        warning('run_pipeline:no_windows', '%s, band %s, %s, %s: %s', ...
                            sub_name, S.bands(b).name, S.inv{iv}, mn{1}, info.no_windows.(mn{1}));
                    end
                end
                % the same on the odd and on the even blocks of the recording
                % (split-half reliability of the statistics); a half without samples
                % (a recording shorter than one block) is NaN, a measure the
                % half is too short for (dwPLI) NaN in that half
                halves = split_halves(seg, onset, srate, S.reliability_block);
                static_halves = nan([2, size(static_full_cmats)]);
                static_power_halves = nan(numel(static_power), 2);
                for h = find(~cellfun(@isempty, halves))
                    [static_halves(h, :, :, :), static_power_halves(:, h)] = ...
                        static_connectivity(src_data_filtered, halves{h}, S.conn, srate, fmin, fmax, Sb);
                end
                saved = {'static_cmats', 'static_full_cmats', 'static_power', 'static_halves', 'static_power_halves', 'info'};
                if ~isempty(cmats)                        % with windows: the window mean on each half
                    [static_halves_window, static_halves_block] = window_mean_halves(cmats, info.windows, ...
                        seg, onset, srate, S.reliability_block, size(static_cmats)); %#ok<ASGLU>
                    saved = [saved, {'static_halves_window', 'static_halves_block'}]; %#ok<AGROW>
                end
                % the band with the least data gives the summary row's duration
                if isnan(qc.analysed_s) || size(src_data_filtered, 2) / srate < qc.analysed_s
                    qc.analysed_s = size(src_data_filtered, 2) / srate;
                    qc.n_pieces = size(seg, 1);
                    if ~isempty(cmats), qc.n_windows = size(cmats{min(S.co, numel(cmats))}, 1); end
                end
                info.halves = halves;
                info.halves_s = cellfun(@(p) sum(diff(p, 1, 2) + 1), halves) / srate;
                info.analysed_s = size(src_data_filtered, 2) / srate;
                info.band = S.bands(b).name;
                info.band_hz = [fmin, fmax];              % fixed, or this recording's anchored edges
                info.band_ref = Sb.band_ref;              % anchored: the edges at IAF 10 Hz, whose centre and
                                                          % transition set the windows and filter ([] fixed)
                info.iaf = iaf;                           % its individual alpha frequency (NaN if not used)
                info.inverse = S.inv{iv};
                info.scout_function = S.scout_function;   % how each region's signal was formed
                info.montage = montage;
                info.atlas = atlas;
                info.anatomy = struct('name', A.name, 'folder', A.dir, 'key', A.key);   % the head model
                info.roi_labels = roi_labels;
                info.channels = labels;
                info.file = eeg.info;
                info.piece_onset_s = onset;               % time of each piece's first sample in the recording
                info.settings_key = key(S.inv{iv}, S.bands(b).name);   % stage 1 settings of this output
                info.input_key = input_key;               % the recording's files
                info.filter_transition_hz = tran;
                for k = 1:numel(info.windows)             % window centres in seconds of the recording
                    info.windows{k}.centre_s = sample_to_time(info.windows{k}.centre_sample, seg, onset, srate);
                end
                if isfield(info, 'plv_wpli_windows')
                    info.plv_wpli_windows.centre_s = sample_to_time(info.plv_wpli_windows.centre_sample, seg, onset, srate);
                end
                tag = [S.inv{iv} '_' tag_ma '_' S.bands(b).name];
                if S.dynamic
                    save(fullfile(currConditionFolder, ['cmats_' tag]), 'cmats', 'info', '-v7.3');
                end
                save(fullfile(currConditionFolder, ['static_cmats_' tag]), saved{:});
            end
        end
        qc.srate = srate;
        qc.n_boundaries = eeg.info.n_boundaries;
        qc.nan_removed_s = eeg.info.nan_s;
        qc.flat_removed_s = eeg.info.flat_s;
        qc.bad_removed_s = eeg.info.bad_s;
        qc.unit = string(eeg.info.unit);
        qc.n_channels = numel(labels);
        qc.channels_left_out = strjoin(string([rep.left_out(:); eeg.info.nan_channels(:)]), ' ');
        qc.subset_missing = strjoin(string(rep.subset_missing), ' ');
        qc.bands_skipped = strjoin(unique(skipped, 'stable'), ' ');
        qc.measures_nan = strjoin(unique(nan_measures, 'stable'), ' ');
        for iv = 1:numel(S.inv)                       % keep the summary row with the outputs
            for b = 1:numel(S.bands)
                f = out(S.inv{iv}, S.bands(b).name);
                if ~isfile(f), continue; end          % band left out (no IAF)
                w = load(f, 'info');  info = w.info;  info.qc = qc; %#ok<NASGU>
                save(f, 'info', '-append');
            end
        end
        if any(anchored), save(iaf_file, 'qc', '-append'); end
    end
end
fprintf('%s done\n', sub_name);
end


function tf = same_key(file, keys, input_key)
% true when FILE exists and was made from the same recording with the same
% settings (one of keys: the current key or that of equivalent code)
tf = false;
if ~isfile(file), return; end
try
    w = load(file, 'info');
    tf = isfield(w.info, 'settings_key') && any(strcmp(w.info.settings_key, keys)) && ...
        isfield(w.info, 'input_key') && strcmp(w.info.input_key, input_key);
catch
end
end


%%% ========== CHANNELS ========== %%%
function [X, labels, idx, rep, eeg, flat_labels, sd] = channels_used(file, ftElec, S)
% the recording's channels of the montage (get_eeg_data, get_eeg_montage);
% flat channels (zero variance: a disconnected electrode, or the
% reference) are left out when S.drop_flat_channels, so that they distort
% neither the average reference nor the inverse.
% A reference kept as a channel of zeros is left out too: the average
% reference over the other channels is still consistent for data and
% leadfield, but the reference site is lost as a sensor (docs/settings.md,
% drop_flat_channels). flat_labels: the flat channels; sd: the SD of every
% matched channel
eeg = get_eeg_data(file, S, ftElec.label);
[X, labels, idx, rep] = get_eeg_montage(eeg.data, eeg.label, ftElec.label, S.channel_subset);
sd = std(double(X), 0, 2);
flat = sd <= 1e-6 * median(sd);
flat_labels = labels(flat);
if any(flat) && S.drop_flat_channels
    X = X(~flat, :);  labels = labels(~flat);  idx = idx(~flat);
end
end

function r = data_rank(X)
% the number of independent signals in the channels used, before the
% average reference (which removes one more, unless the recording is
% already average-referenced), by connectivity.signal_rank (MNE-Python's
% rank with tol = 'float32': the rounding of data stored as 32-bit floats
% does not count). Each ICA component removed, each interpolated channel
% (and each pair of bridged channels) lowers it by one. Data stored as
% 16-bit integers (EDF) are rounded above this limit: a file cleaned and
% then saved as EDF has the rank of its number of channels
r = connectivity.signal_rank(X);
end

function warn_rank(sub_name, r, n, min_channels)
% a warning when the recording has fewer independent signals than the
% fewest channels the pipeline accepts (S.min_channels). Recordings
% cleaned by ICA can have a rank far below their number of channels
if r < min_channels
    warning('run_pipeline:low_rank', ['%s: data rank %d of %d channels, fewer independent signals than ' ...
        'S.min_channels (%d): were most ICA components removed? Its regions are mixtures of only %d signals'], ...
        sub_name, r, n, min_channels, r);
end
end

function add_rank(file, r)
% the data rank added to the summary row stored with an output (info.qc)
% or with the individual alpha frequency (qc)
if ~isfile(file), return; end
w = whos('-file', file);
if ismember('info', {w.name})
    s = load(file, 'info');
    if ~isfield(s.info, 'qc'), return; end
    s.info.qc.data_rank = r;
elseif ismember('qc', {w.name})
    s = load(file, 'qc');
    s.qc.data_rank = r;
else
    return
end
save(file, '-struct', 's', '-append');
end

function [ftElec, ftLeadfield] = keep_channels(ftElec, ftLeadfield, idx)
% electrodes and leadfield rows of the channels used, leadfield average-
% referenced over these channels
ftElec.label = ftElec.label(idx);
ftElec.chanpos = ftElec.chanpos(idx, :);
ftElec.elecpos = ftElec.elecpos(idx, :);
ftElec.tra = eye(numel(idx));
ftLeadfield.label = ftLeadfield.label(idx);
for i = 1:numel(ftLeadfield.leadfield)
    g = ftLeadfield.leadfield{i}(idx, :);
    ftLeadfield.leadfield{i} = g - mean(g, 1);
end
end


%%% ========== BAND-PASS FILTER ========== %%%
function t = sample_to_time(s, seg, onset, srate)
% sample index in the concatenated pieces -> seconds in the recording
p = arrayfun(@(v) find(seg(:, 1) <= v, 1, 'last'), s(:)');
t = reshape(onset(p), 1, []) + (s(:)' - reshape(seg(p, 1), 1, [])) / srate;   % row, like the window starts
end

function tf = is_anchored(band, opt)
% true when an edge of the band refers to the individual alpha frequency
[~, ~, tf] = band_edges(band, 10, opt);
end

function tf = band_left_out(band, iaf, opt)
% a band anchored on the IAF that this recording cannot have: no IAF, or
% edges less than 1 Hz apart or not above 0 Hz for its IAF (e.g. delta,
% 1 Hz to IAF - 6, for an IAF below 8 Hz)
[fmin, fmax] = band_edges(band, iaf, opt);
tf = isnan(fmin) || isnan(fmax) || fmax - fmin < 1 || fmin <= 0;
end

function keys = accepted_keys(S, inv, band)
[~, ~, keys] = settings_key(S, '', inv, band);
end
