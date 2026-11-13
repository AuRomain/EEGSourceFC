function [R, stems, files] = psd_spectra(S)
%PSD_SPECTRA The power spectra of every recording, at the scalp and in the regions (kept, with keys).
%   [R, stems, files] = psd_spectra(S) (run_psd_analysis) for the
%   recordings of S.data (S.file_pattern), each:
%   - read as stage 1 reads it (get_eeg_data: the analysed stretch, its
%     continuous pieces; the channels of the first cap S.montages{1},
%     flat channels left out, S.channel_subset; average reference);
%   - 'sensor' (S.psd.levels): each channel's PSD (psd_welch: Hann
%     segments of S.psd.segment s with S.psd.overlap, inside the continuous
%     pieces; S.psd.freq_range on a grid of 1 / S.psd.segment Hz) and its
%     specparam fit over S.psd.fit_range (psd_specparam);
%   - 'source': each region's time series for every inverse method of
%     S.inv on the atlas S.atlases{1}, as stage 1 makes them
%     (psd_source_signals), and their PSD and fit, the same way; and the
%     PSD of every source of the cortex, w' C(f) w with w its inverse
%     filter and C the channels' cross-spectra (psd_cross_spectra, the
%     same segments; with free orientation the sum of its three), whose
%     mean over each region's sources (pow_vertex_mean, with its fit
%     fit_vertex_mean) has none of the cancellation of opposite or
%     out-of-phase sources in the region's signal; the sources' PSDs
%     themselves (pow_vertex, single) with S.psd.vertex_maps, when every
%     recording has the same cortex (not S.anatomy_individual);
%   - its individual alpha frequency (individual_alpha_frequency with
%     S.iaf, as stage 1), only when a band of S.psd.bands is anchored on it.
%   Kept in results/PSD/spectra/<sensor|source>/<recording>.mat with a key:
%   the recording's files (by content), the settings each level depends on
%   (reading, cap, S.psd's spectral settings, S.iaf when used; for the
%   regions also the atlas, head, inverse and scout settings, the shipped
%   inputs and S.psd.vertex_maps) and the code (stage 1's, whose functions read and
%   project the data, and this step's). A level is computed again only
%   when its key changes; changing a band or a test recomputes nothing.
%   A recording that fails (no piece as long as a segment, too few
%   channels, memory) is left out with a warning; its error is saved
%   without a key, so the next run tries it again.
%   R: one struct per recording (the failed ones: .error), stems its names;
%   R(i).file is its file of the regions' level, which keeps pow_vertex
%   (left out of R, read one recording at a time by psd_vertex_maps).
o = S.psd;
cache = fullfile(S.results, 'PSD', 'spectra');
lv = {'sensor', 'source'};
for k = 1:2
    if ~isfolder(fullfile(cache, lv{k})), mkdir(fullfile(cache, lv{k})); end
end
subs = get_subs(S.data, S.file_pattern);
files = fullfile(S.data, subs);
n = numel(subs);
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
if numel(unique(stems)) < n
    error('psd_spectra:names', 'Two recordings with the same file name in %s: their spectra would share a file', S.data);
end
use = [true, ismember('source', o.levels)];      % the channels always (they also give the IAF)
anchored = any(arrayfun(@(b) is_anchored(b, S.iaf), o.bands));
% what each level depends on
V = struct();
for f = {'srate', 'resample', 'core_start', 'core_duration', 'boundary_events', 'drop_flat_channels', ...
        'channel_subset', 'min_channels'}
    V.(f{1}) = S.(f{1});
end
if isfield(S, 'bad_events') && ~isequal(S.bad_events, {'bad*'}), V.bad_events = S.bad_events; end   % as settings_key
V.montage = S.montages{1};
if isfile(V.montage), V.montage_file = file_fingerprint({V.montage}); end
V.psd = rmfield(o, setdiff(fieldnames(o), {'segment', 'overlap', 'freq_range', 'fit_range', 'power_line'}));
if anchored                                          % the IAF: its settings and restingIAF's files
    V.iaf = S.iaf;
    if isfield(S, 'versions') && isfield(S.versions, 'fingerprints') && isfield(S.versions.fingerprints, 'restingiaf')
        V.restingiaf = S.versions.fingerprints.restingiaf;
    end
end
here = fileparts(mfilename('fullpath'));
fooof = which('process_fooof');                      % the specparam that psd_specparam calls, on the
if isempty(fooof)                                    % path before the key (or the key changes once it is)
    clear ft_hastoolbox                              % it repeats its last answer for the same question
    ft_hastoolbox('brainstorm', 1);
    fooof = which('process_fooof');
end
if isempty(fooof)
    error('psd_spectra:fooof', 'process_fooof (Brainstorm''s specparam, in FieldTrip''s external/brainstorm) not found');
end
V.code = hash_text([code_fingerprint(S, 'stage1'), file_fingerprint(fullfile(here, ...
    {'psd_spectra.m', 'psd_welch.m', 'psd_specparam.m', 'psd_source_signals.m', 'psd_cross_spectra.m', ...
    'psd_vertex_spectra.m', 'psd_grid.m'})), ...
    file_fingerprint({fooof})]);
V.optimiser = specparam_optimiser();                 % fmincon (Optimization Toolbox) or fminsearch: the fits differ
W = V;
for f = {'atlases', 'inv', 'eloreta_lambda', 'lcmv_lambda', 'wmne', 'orientation', 'scout_function', ...
        'anatomy', 'anatomy_individual', 'anatomy_fallback'}
    W.(f{1}) = S.(f{1});
end
W.atlases = W.atlases(1);
if isfile(W.atlases{1}), W.atlas_file = file_fingerprint(W.atlases(1)); end
if ismember('lcmv', cellstr(S.inv)) && isfield(S, 'lcmv_kappa') && ~isempty(S.lcmv_kappa)
    W.lcmv_kappa = S.lcmv_kappa;                     % LCMV's truncation, as settings_key
end
W.inputs = builtin_inputs(S);
W.vertex_maps = o.vertex_maps;                       % the sources' PSDs kept or not
keys = cell(n, 2);
todo = false(n, 2);
S_ma = S;  S_ma.montages = S.montages(1);  S_ma.atlases = S.atlases(1);
S_ch = S_ma;  S_ch.anatomy = '';
for s = 1:n
    [~, input_key] = settings_key(S_ch, files{s});        % the recording's files, by content
    keys{s, 1} = hash_text(jsonencode(struct('settings', V, 'input', input_key)));
    try
        [~, input_key] = settings_key(S_ma, files{s});    % with its head, when not the shipped one
    catch err                                            % no head (S.anatomy_individual): fails below
        input_key = ['no head: ' err.message];
    end
    keys{s, 2} = hash_text(jsonencode(struct('settings', W, 'input', input_key)));
    for k = find(use)
        c = fullfile(cache, lv{k}, [stems{s} '.mat']);
        todo(s, k) = ~isfile(c) || ~strcmp(getfield(load(c, 'key'), 'key'), keys{s, k}); %#ok<GFLD>
    end
end
fprintf('PSD: %d recordings; to compute: %d (sensor), %d (source); specparam fits with %s\n', n, ...
    nnz(todo(:, 1)), nnz(todo(:, 2)), V.optimiser);
t0 = tic;
idx = find(any(todo, 2));
if S.parallel && numel(idx) > 4 && ~isempty(ver('parallel'))
    if isempty(gcp('nocreate')), parpool('Processes', min(S.max_workers, numel(idx))); end
    parfor i = 1:numel(idx)
        s = idx(i); %#ok<PFBNS>
        one_recording(files{s}, subs{s}, stems{s}, S, todo(s, :), keys(s, :), cache);
    end
else
    for s = idx'
        one_recording(files{s}, subs{s}, stems{s}, S, todo(s, :), keys(s, :), cache);
    end
end
if ~isempty(idx), fprintf('PSD: spectra of %d recordings computed (%.1f min)\n', numel(idx), toc(t0) / 60); end
R = cell(n, 1);
for s = 1:n
    r = struct();
    for k = find(use)
        c = fullfile(cache, lv{k}, [stems{s} '.mat']);
        w = load(c, 'r');
        if k == 2 && isfield(w.r, 'inverses') && isfield(w.r.inverses, 'pow_vertex')
            w.r.inverses = rmfield(w.r.inverses, 'pow_vertex');   % the sources' PSDs stay in the file
            w.r.file = c;                                       % (psd_vertex_maps reads them one by one)
        end
        r.(lv{k}) = w.r;                             % .error: this level failed
    end
    R{s} = r;
end
R = [R{:}]';
end

function one_recording(file, sub, stem, S, todo, keys, cache)
% the levels to compute of one recording, each saved with its key; an
% error is saved as the recording's result without a key (tried again in
% the next run)
lv = {'sensor', 'source'};
try
    r = measure(file, sub, S, todo);
catch err
    fprintf(2, '%s: no spectra: %s\n', stem, err.message);
    r = struct('sensor', struct('error', string(err.message)), 'source', struct('error', string(err.message)));
end
for k = find(todo)
    key = keys{k};
    if isfield(r.(lv{k}), 'error'), key = ''; end   % not kept: the next run computes it again
    save_level(fullfile(cache, lv{k}, [stem '.mat']), r.(lv{k}), key);
end
end

function save_level(file, r, key) %#ok<INUSD>
save(file, 'r', 'key');
end

function r = measure(file, sub, S, todo)
o = S.psd;
ftElec = get_ftElec(S.montages{1});
eeg = get_eeg_data(file, S, ftElec.label);
[X, labels, idx] = get_eeg_montage(eeg.data, eeg.label, ftElec.label, S.channel_subset);
sd = std(double(X), 0, 2);
flat = sd <= 1e-6 * median(sd);
if any(flat) && S.drop_flat_channels                 % as stage 1 (run_pipeline, channels_used)
    X = X(~flat, :);  labels = labels(~flat);  idx = idx(~flat);
end
if numel(labels) < S.min_channels
    error('psd_spectra:channels', ['%d channels match the cap %s (S.min_channels = %d). Labels that are ' ...
        'not 10-20 names (EGI E1.., BioSemi A1..) need the cap''s channel file in S.montages, e.g. from %s ' ...
        '(docs/settings.md, Caps)'], numel(labels), S.montages{1}, S.min_channels, ...
        fullfile(S.toolboxes, 'brainstorm3', 'defaults', 'eeg', 'ICBM152'));
end
X = X - mean(X, 1);                                  % average reference, as stage 1
fs = eeg.srate;
if o.freq_range(2) >= fs / 2
    error('psd_spectra:nyquist', 'S.psd.freq_range reaches %g Hz, the Nyquist frequency of this recording is %g Hz', ...
        o.freq_range(2), fs / 2);
end
grid = psd_grid(o);
info = struct('srate', fs, 'analysed_s', sum(diff(eeg.segments, 1, 2) + 1) / fs, 'n_pieces', size(eeg.segments, 1));
r = struct();
if todo(1)
    [P, f, nseg] = psd_welch(X, fs, eeg.segments, o.segment, o.overlap, o.freq_range);
    P = on_grid(P, f, grid);
    iaf = struct('value', NaN);
    if any(arrayfun(@(b) is_anchored(b, S.iaf), o.bands))
        iaf = individual_alpha_frequency(X, labels, fs, eeg.segments, S.iaf);
    end
    r.sensor = struct('labels', {cellstr(labels(:)')}, 'freq', grid, 'pow', single(P), ...
        'fit', psd_specparam(grid, P, o.fit_range, o.power_line), 'n_segments', nseg, 'iaf', iaf, 'info', info);
end
if todo(2)                                           % the regions: an error leaves the channels' spectra
    try
        [src, roi_labels, A, filt, vert] = psd_source_signals(X, idx, fs, sub, S);
        [Cx, fc] = psd_cross_spectra(X, fs, eeg.segments, o.segment, o.overlap, o.freq_range);
        keep_vertex = o.vertex_maps && ~S.anatomy_individual;
        lev = struct('inverse', {}, 'labels', {}, 'freq', {}, 'pow', {}, 'fit', {}, 'n_segments', {}, ...
            'pow_vertex_mean', {}, 'fit_vertex_mean', {}, 'pow_vertex', {});
        for iv = 1:numel(S.inv)
            [P, f, nseg] = psd_welch(src{iv}, fs, eeg.segments, o.segment, o.overlap, o.freq_range);
            P = on_grid(P, f, grid);
            [Pm, Pv] = psd_vertex_spectra(filt{iv}, Cx, vert, strcmp(S.orientation, 'free'));
            filt{iv} = [];
            Pm = on_grid(Pm, fc, grid);
            if keep_vertex, Pv = single(on_grid(Pv, fc, grid)); else, Pv = single([]); end
            lev(iv) = struct('inverse', S.inv{iv}, 'labels', {cellstr(roi_labels(:)')}, 'freq', grid, ...
                'pow', single(P), 'fit', psd_specparam(grid, P, o.fit_range, o.power_line), 'n_segments', nseg, ...
                'pow_vertex_mean', single(Pm), 'fit_vertex_mean', psd_specparam(grid, Pm, o.fit_range, o.power_line), ...
                'pow_vertex', Pv);
        end
        r.source = struct('inverses', lev, 'anatomy', A.name, 'channels', {cellstr(labels(:)')}, 'info', info);
    catch err
        [~, stem] = fileparts(file);
        fprintf(2, '%s: no region spectra: %s\n', stem, err.message);
        r.source = struct('error', string(err.message));
    end
end
end

function P = on_grid(P, f, grid)
% the spectra on the common frequency grid (the same for every recording
% when segment x sampling rate is a whole number of samples, as usual)
if numel(f) == numel(grid) && max(abs(f(:) - grid(:))) < 1e-6, return; end
P = interp1(f(:), P', grid(:), 'linear')';
if any(isnan(P(:)))
    error('psd_spectra:grid', 'The spectrum (%g-%g Hz) does not cover S.psd.freq_range', f(1), f(end));
end
end

function tf = is_anchored(band, opt)
[~, ~, tf] = band_edges(band, 10, opt);
end

function o = specparam_optimiser()
% the optimiser of Brainstorm's specparam (process_fooof): fmincon when the
% Optimization Toolbox is installed, fminsearch otherwise (the fits differ)
o = 'fminsearch';
if exist('fmincon') == 2, o = 'fmincon'; end %#ok<EXIST>
end

function h = builtin_inputs(S)
% the files of stage 1's inputs/ that the regions' signals depend on, by
% content: the cortex, the scouts and the files of the cap and the atlas
in = fullfile(S.root, '1- EEG source connectivity', 'inputs');
d = dir(fullfile(in, '*.mat'));
names = {d.name};
use = ismember(names, {'tess_cortex_pial_low.mat', 'scout_Desikan-Killiany_68.mat'});
for c = [S.montages(1), S.atlases(1)]
    use = use | contains(names, file_tag(c{1}), 'IgnoreCase', true);
    if isfile(c{1})                                % a cap given as a file: its leadfield's name (as settings_key)
        [~, stem] = fileparts(c{1});
        use = use | contains(names, ['ftLeadfield_' matlab.lang.makeValidName(stem) '_'], 'IgnoreCase', true);
    end
end
h = file_fingerprint(fullfile(in, sort(names(use))));
end
