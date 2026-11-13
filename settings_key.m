function [key, input_key, accepted] = settings_key(S, file, inv, band)
%SETTINGS_KEY Fingerprints of the stage 1 settings and of one input recording.
%   [key, input_key] = settings_key(S, file)
%   [key, input_key] = settings_key(S, file, inv, band)
%   [key, input_key, accepted] = settings_key(...)
%   key        changes when anything that shapes stage 1 changes: its
%              settings, its code (code_fingerprint 'stage1', with the
%              FieldTrip and Brainstorm versions) and the built-in inputs it
%              reads (cap, leadfield, cortex and atlas files of inputs/, by
%              content). With inv and band, the key of that output only:
%              other bands and inverse methods of S do not enter it (adding
%              or removing a band does not make the other bands' outputs stale)
%   accepted   key and the keys the same settings have with other stage 1
%              code listed in "1- EEG source connectivity/equivalent_code.txt"
%              as computing the same outputs as the current code (a change
%              that only adds a measure, for instance): outputs with any of
%              these keys are up to date; for a static-only run without
%              the window mean, also the keys of the same run with it
%              (S.static_matrices with 'window_mean': those outputs also
%              hold the whole-recording matrices)
%   Runs with dwPLI (computed by FieldTrip with epochs as trials, in
%   windows of 60 cycles by default: V.dwpli) or a symmetric
%   orthogonalisation (its rank rule: V.rank_rule) have keys of their own,
%   as have bands anchored on the IAF (windows, epochs and filter of the
%   band at IAF 10 Hz, the IAF from each piece's segments: V.anchored).
%   input_key  changes when the recording changes: its file and the files
%              next to it with the same name (EEGLAB .fdt, BrainVision .eeg
%              and .vmrk, ...), by content, so a copy on another computer
%              keeps its key ('' when file is ''); a recording that is a
%              folder (EGI .mff) by the files in it; with another head
%              (S.anatomy, get_anatomy) also that head's files: outputs
%              made with another head are not reused; the shipped
%              ICBM152 (also a fallback subject's) adds nothing to the key
%   run_pipeline stores both with the outputs and recomputes a recording
%   whose stored keys differ, instead of reusing results made from another
%   file, with other settings or by other code; the states stop if the
%   recordings it combines were made differently.
if nargin > 2 && ~isempty(inv), S.inv = {inv}; end
if nargin > 3 && ~isempty(band), S.bands = S.bands(strcmp({S.bands.name}, band)); end
f = {'srate', 'resample', 'core_start', 'core_duration', 'edge_trim', 'filter_transition', 'filter_span', 'drop_flat_channels', ...
     'boundary_events', 'min_segment', 'montages', 'channel_subset', 'atlases', 'inv', ...
     'eloreta_lambda', 'lcmv_lambda', 'wmne', 'bands', 'window_cycles', 'window_length', ...
     'window_overlap', 'aec_window', 'aec_step', 'conn', 'orientation', 'dynamic', 'reliability_block'};
V = struct();
for i = 1:numel(f)
    if isfield(S, f{i}), V.(f{i}) = S.(f{i}); end
end
% the scout function, unless it is the mean (the default adds nothing to
% the key)
if isfield(S, 'scout_function') && ~strcmp(S.scout_function, 'mean'), V.scout_function = S.scout_function; end
% LCMV's truncation (S.lcmv_kappa), only when set and LCMV is computed
% (otherwise it adds nothing to the key)
if isfield(S, 'lcmv_kappa') && ~isempty(S.lcmv_kappa) && isfield(S, 'inv') && ismember('lcmv', cellstr(S.inv))
    V.lcmv_kappa = S.lcmv_kappa;
end
% the events marking bad data, unless the default (as the scout function)
if isfield(S, 'bad_events') && ~isequal(S.bad_events, {'bad*'}), V.bad_events = S.bad_events; end
for fld = {'montages', 'atlases'}                % caps and atlases given as files: their content too
    if isfield(S, fld{1})
        for c = cellstr(S.(fld{1}))
            if isfile(c{1})
                [~, name] = fileparts(c{1});
                V.([fld{1} '_files']).(matlab.lang.makeValidName(name)) = file_fingerprint(c(1));
            end
        end
    end
end
% dwPLI computed by FieldTrip with epochs as trials: its estimator, epochs
% and windows enter the key of every run with dwPLI; runs without dwPLI
% have no such part
if isfield(S, 'conn') && ismember('dwpli', cellstr(S.conn))
    V.dwpli = struct('estimator', 'fieldtrip_wpli_debiased_epochs', 'epoch_cycles', field_or(S, 'dwpli_epoch_cycles', 6), ...
        'epoch_length', field_or(S, 'dwpli_epoch_length', []), ...
        'windows', 'cycles60_or_10_epochs');   % default windows, epochs not lost to rounding
    if isfield(S, 'versions') && isfield(S.versions, 'fingerprints') && isfield(S.versions.fingerprints, 'fieldtrip_dwpli')
        V.dwpli.fieldtrip_files = S.versions.fingerprints.fieldtrip_dwpli;   % ft_connectivity_wpli, ft_specest_mtmfft
    end
end
% the rank that decides the symmetric orthogonalisation (connectivity.
% signal_rank, MNE-Python's rule): part of the key of runs with it
if isfield(S, 'conn') && any(ismember({'aec_symorth', 'plv_symorth'}, cellstr(S.conn)))
    V.rank_rule = 'mne_float32';
end
if isfield(S, 'iaf')                               % bands anchored on the IAF: its settings too
    for b = 1:numel(S.bands)
        [~, ~, anchored] = band_edges(S.bands(b), 10, S.iaf);
        if anchored
            V.iaf = S.iaf;
            % windows, epochs and filter of the band at IAF 10 Hz
            % (cycles_centre), the IAF from each piece's segments
            % (resting_iaf_pieces), and restingIAF's files: part of the key
            % of anchored bands only (not of fixed bands)
            V.anchored = 'iaf10_pieces';
            if isfield(S, 'versions') && isfield(S.versions, 'fingerprints') && isfield(S.versions.fingerprints, 'restingiaf')
                V.restingiaf = S.versions.fingerprints.restingiaf;
            elseif isfield(S, 'versions') && isfield(S.versions, 'restingiaf')
                V.restingiaf = S.versions.restingiaf;
            end
            break
        end
    end
end
% a static-only run that asks for the window mean also computes and saves
% the windows' mean (S.dynamic always does): another key. Without it, the
% key has no window-mean part, and outputs made with the window mean are
% accepted too (they hold everything such a run reads)
% windows of single measures (S.measure_windows): only when some are given
% (otherwise they add nothing to the key) and windows are computed (the
% states, or the window mean)
mw = [];
if isfield(S, 'measure_windows') && isstruct(S.measure_windows) && ~isempty(fieldnames(S.measure_windows))
    mw = orderfields(S.measure_windows);
    for m = fieldnames(mw)'
        if isstruct(mw.(m{1})), mw.(m{1}) = orderfields(mw.(m{1})); end
    end
end
dynamic = ~isfield(S, 'dynamic') || S.dynamic;
if dynamic && ~isempty(mw), V.measure_windows = mw; end
variants = {V};
if ~dynamic
    Vw = V;  Vw.window_mean = true;
    if ~isempty(mw), Vw.measure_windows = mw; end
    if ~isfield(S, 'static_matrices') || isempty(S.static_matrices) || ...   % the default, or asked for
            any(strcmp(cellstr(S.static_matrices), 'window_mean'))
        variants = {Vw};
    else
        variants{2} = Vw;
    end
end
[code_key, code, toolboxes] = code_fingerprint(S, 'stage1');
codes = [{code_key}, cellfun(@(e) hash_text([e toolboxes]), equivalent_code(S, code), 'UniformOutput', false)];
inputs = builtin_inputs(S);
accepted = cell(1, 0);
for iv = 1:numel(variants)
    for ic = 1:numel(codes)
        Vk = variants{iv};
        Vk.code = codes{ic};
        Vk.inputs = inputs;
        accepted{end + 1} = hash_text(jsonencode(Vk)); %#ok<AGROW>
    end
end
key = accepted{1};
input_key = '';
if nargin > 1 && ~isempty(file) && (isfile(file) || isfolder(file))
    [folder, stem] = fileparts(file);
    d = dir(fullfile(folder, [stem '.*']));
    d = d(~[d.isdir]);
    names = sort({d.name});
    files = cell(1, 0);
    if ~isempty(names), files = fullfile(folder, names); end
    if isfolder(file)                              % a recording that is a folder (EGI .mff): its files
        c = dir(fullfile(file, '**', '*'));
        c = c(~[c.isdir]);
        files = [files, sort(fullfile({c.folder}, {c.name}))];
    end
    input_key = file_fingerprint(files);
    if isfield(S, 'anatomy') && ~isempty(S.anatomy)
        A = get_anatomy(S, file);
        if ~isempty(A.key), input_key = hash_text([input_key '|anatomy:' A.key]); end
    end
end
end

function h = builtin_inputs(S)
% the files of inputs/ that stage 1 reads, by content: the cortex, the
% scouts and the files of each cap and atlas named by the settings
in = fullfile(S.root, '1- EEG source connectivity', 'inputs');
d = dir(fullfile(in, '*.mat'));
names = {d.name};
use = ismember(names, {'tess_cortex_pial_low.mat', 'scout_Desikan-Killiany_68.mat'});
for c = [cellstr(S.montages), cellstr(S.atlases)]
    t = file_tag(c{1});
    use = use | contains(names, t, 'IgnoreCase', true);
    if isfile(c{1})                                % a cap given as a file: its leadfield on the shipped
        [~, stem] = fileparts(c{1});               % head is named with makeValidName (so a name
        use = use | contains(names, ['ftLeadfield_' matlab.lang.makeValidName(stem) '_'], 'IgnoreCase', true);
    end                                            % with a space or a leading digit is found too)
end
h = file_fingerprint(fullfile(in, sort(names(use))));
end

function v = field_or(S, f, default)
v = default;
if isfield(S, f) && ~isempty(S.(f)), v = S.(f); end
end

function earlier = equivalent_code(S, code)
% the earlier stage 1 codes that equivalent_code.txt lists for the current
% one: rows '<current> <earlier> <note>', lines starting with % are comments
earlier = {};
f = fullfile(S.root, '1- EEG source connectivity', 'equivalent_code.txt');
if ~isfile(f), return; end
rows = strtrim(splitlines(fileread(f)));
rows = rows(~startsWith(rows, '%') & strlength(rows) > 0);
for i = 1:numel(rows)
    w = split(rows{i});
    if numel(w) >= 2 && strcmp(w{1}, code), earlier{end + 1} = w{2}; end %#ok<AGROW>
end
end
