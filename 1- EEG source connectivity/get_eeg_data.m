function eeg = get_eeg_data(file, S, montage_labels)
%GET_EEG_DATA Read one recording, in any format FieldTrip reads, and cut the analysed segment.
%   eeg = get_eeg_data(file, S) with S from pipeline_settings. This is the
%   "load the EEG data" step of the original run_pipeline, which read an
%   EEGLAB .set with load() and kept the first 180 s, without checking the
%   sampling rate, the channel labels or the length (E5).
%   eeg = get_eeg_data(file, S, montage_labels): the channels of the file
%   that match these labels (the cap's, get_ftElec; with S.channel_subset,
%   as get_eeg_montage matches them) are the EEG channels: only they decide
%   the cuts at non-finite samples and flat stretches, the channels left
%   out as entirely non-finite and the units reported. Other channels of
%   the file (EOG, ECG, a status or trigger channel), which stage 1 leaves
%   out, no longer cut the EEG (a NaN stretch in an ECG channel cut every
%   channel). Without montage_labels, or when none matches, every channel.
%   Here:
%   - any format FieldTrip reads (ft_read_header / ft_read_data /
%     ft_read_event): EEGLAB .set (data inside or in a .fdt), EDF/BDF,
%     BrainVision, FIF, EGI .mff (a folder), ...; continuous or epoched;
%   - the data are in microvolts: channels whose unit is V, mV or nV are
%     converted (FIF files are in V). The unit is FieldTrip's hdr.chanunit;
%     for BrainVision the unit of each channel in the .vhdr (FieldTrip
%     always says uV, though the values are in the .vhdr's unit), for BDF
%     the file's PhysDim when FieldTrip says 'unknown'. A unit that is
%     unknown or missing (EEGLAB files often) is taken as microvolts. The
%     units read and the factors applied are in eeg.info (unit,
%     unit_scale); a median channel SD below 0.01 or above 10,000 µV after
%     this gives a warning (get_eeg_data:amplitude: are the units right?).
%     Connectivity does not depend on the unit; band power does;
%   - the sampling rate is checked against S.srate (when set), and the data
%     resampled to S.resample (when set);
%   - discontinuities are kept apart: EEGLAB 'boundary' events (rejected
%     stretches; also the numeric boundary type -99 of EEGLAB files saved
%     with numeric event types) and the other events of S.boundary_events
%     (BrainVision 'New Segment'), and, for epoched files, the epoch
%     limits split the data into continuous pieces, which are filtered and
%     windowed separately by run_pipeline (no window or filter crosses a
%     cut);
%   - stretches marked bad by events are removed: events whose type or
%     value matches S.bad_events (default {'bad*'}, case ignored, * any
%     text: MNE's BAD_ annotations, BrainVision Analyzer's 'Bad Interval'),
%     as FieldTrip reads them from EEGLAB events, EDF+ annotations and
%     BrainVision markers, over their duration (a zero duration: a cut at
%     that sample, like MNE's 'BAD boundary'). FIF files: FieldTrip reads
%     no annotations, so their BAD_ spans are not removed (export the
%     recording from MNE to EEGLAB, EDF or BrainVision, where they become
%     events, or remove the spans before);
%   - samples that are NaN or Inf (stretches marked bad, e.g. by EEGLAB or
%     FieldTrip) are cuts like boundary events: the bad stretch is removed and
%     the data on each side are separate pieces; a channel with no finite
%     sample is left out (V3-7; earlier one NaN sample made a whole piece NaN);
%   - stretches of at least 0.1 s in which every channel stays exactly
%     constant (blanked or held data, e.g. a rejected stretch filled with
%     its last value) are cut the same way: there is no EEG in them, and
%     their zero band-passed signal gave infinite wPLI values that stopped
%     stage 2 (found on the OpenNeuro data, 22 Oct 2026: 10.7 and 16.7 s in
%     two recordings);
%   - the analysed segment is [S.core_start, S.core_start + S.core_duration]
%     s (the whole recording when S.core_duration is empty).
%   eeg.data      nChannels x nSamples (double, µV), pieces one after the other
%   eeg.label     channel labels of the file
%   eeg.srate     sampling rate of eeg.data
%   eeg.segments  nPieces x 2: first and last sample of each continuous piece in eeg.data
%   eeg.onset_s   nPieces x 1: time of each piece's first sample in the recording (s)
%   eeg.info      file, format, original rate, duration, number of cuts,
%                 seconds of the analysed segment removed as non-finite
%                 (nan_s), as flat (flat_s) and as marked bad by events
%                 (bad_s; n_bad_events), channels left out as entirely
%                 non-finite (nan_channels), the units of the EEG channels
%                 (unit) and the factors that made them µV (unit_scale)

if nargin < 3, montage_labels = {}; end
if ~isfield(S, 'bad_events'), S.bad_events = {'bad*'}; end
fmt = ft_filetype(file);
hdr = ft_read_header(file);
fs0 = hdr.Fs;
if ~isempty(S.srate) && abs(fs0 - S.srate) > 1e-6
    error('get_eeg_data:srate', '%s: sampling rate %g Hz, settings say %g Hz (set S.srate = [] to accept any rate)', ...
        file, fs0, S.srate);
end
nt = max(1, hdr.nTrials);                                 % epochs (1 = continuous)
ns = hdr.nSamples;                                        % samples per epoch
dat = ft_read_data(file, 'header', hdr, 'checkboundary', false);
dat = reshape(double(dat), size(dat, 1), []);             % chan x samples, epochs one after the other
label = hdr.label(:);
% the EEG channels: those of the cap (as get_eeg_montage matches them)
use = true(numel(label), 1);
if ~isempty(montage_labels)
    w = warning('off', 'get_eeg_montage:mixed_reference');   % run_pipeline gives it, once
    rows = get_eeg_montage((1:numel(label))', label, montage_labels, S.channel_subset);
    warning(w);
    if ~isempty(rows), use(:) = false;  use(rows) = true; end
end
% units: microvolts
[unit, scale] = channel_units(file, fmt, hdr);
if any(scale ~= 1), dat = dat .* scale; end
% channels without any finite sample are left out (V3-7)
dead = use & all(~isfinite(dat), 2);
nan_channels = label(dead);
if any(dead)
    warning('get_eeg_data:nan_channel', '%s: channel(s) %s have no finite sample and are left out', ...
        file, strjoin(nan_channels, ', '));
    dat = dat(~dead, :);  label = label(~dead);  use = use(~dead);
    unit = unit(~dead);  scale = scale(~dead);
end
n = size(dat, 2);
% the analysed segment, in samples of the file (removed seconds are counted there)
t0 = S.core_start;
t1 = inf;
if ~isempty(S.core_duration), t1 = t0 + S.core_duration; end
tf = (0:n - 1) / fs0;
in_core = tf >= t0 - 1e-9 & tf < t1 - 1e-9;
% continuous pieces: epochs, then cuts at boundary events and at bad
% events without duration; spans of bad events with a duration removed
cuts = (1:nt - 1)' * ns;                                  % last sample of each epoch
try
    ev = ft_read_event(file, 'header', hdr);
catch err
    ev = [];
    if ~isempty(S.boundary_events) || ~isempty(S.bad_events)
        warning('get_eeg_data:events', '%s: events could not be read (%s); no cuts from events', file, err.message);
    end
end
nb_boundaries = 0;
marked = false(1, n);                                     % spans of bad events
n_bad_events = 0;
if ~isempty(ev)
    is_b = arrayfun(@(e) is_boundary(e, S.boundary_events, strcmp(fmt, 'eeglab_set')), ev(:)');
    b = round([ev(is_b).sample]);
    b = b(b > 1 & b < n);
    nb_boundaries = numel(b);
    is_bad = ~is_b & arrayfun(@(e) matches_pattern(e.type, S.bad_events) || ...
        matches_pattern(e.value, S.bad_events), ev(:)');
    for e = find(is_bad)
        s = round(ev(e).sample);
        if ~(isnumeric(s) && isscalar(s) && isfinite(s)), continue; end
        d = 0;
        if isfield(ev, 'duration') && isnumeric(ev(e).duration) && isscalar(ev(e).duration) && ...
                isfinite(ev(e).duration)
            d = round(ev(e).duration);
        end
        if d >= 1
            marked(max(1, s):min(n, s + d - 1)) = true;
        elseif s > 1 && s < n
            b(end + 1) = s; %#ok<AGROW>                   % a marker without duration: a cut
        end
    end
    n_bad_events = nnz(is_bad);
    cuts = unique([cuts; b(:) - 1]);                      % piece ends just before the boundary sample
end
seg = [[1; cuts + 1], [cuts; n]];
seg = seg(seg(:, 2) >= seg(:, 1), :);
% non-finite samples (in an EEG channel) are cut out: each piece is split
% into its runs of finite samples (V3-7)
nonfinite = any(~isfinite(dat(use, :)), 1);
nan_s = nnz(nonfinite & in_core) / fs0;
% flat stretches: every EEG channel unchanged for at least 0.1 s (no EEG there)
same = [false, all(diff(dat(use, :), 1, 2) == 0, 1)];   % equal to the previous sample in every channel
d = diff([false, same, false]);
st = find(d == 1);  en = find(d == -1) - 1;              % runs of repeated samples
flat = false(1, n);
for k = find(en - st + 2 >= round(0.1 * fs0))            % the run plus the sample it repeats
    flat(st(k) - 1:en(k)) = true;
end
flat_s = nnz(flat & ~nonfinite & in_core) / fs0;
if flat_s > 0
    warning('get_eeg_data:flat', '%s: %.2f s in which every channel is constant cut out', file, flat_s);
end
bad_s = nnz(marked & ~nonfinite & ~flat & in_core) / fs0;
if n_bad_events > 0
    warning('get_eeg_data:bad_events', '%s: %d event(s) marking bad data (S.bad_events): %.2f s of the analysed segment cut out', ...
        file, n_bad_events, bad_s);
end
bad = nonfinite | flat | marked;
if any(bad)
    good = {};
    for k = 1:size(seg, 1)
        r = seg(k, 1):seg(k, 2);
        g = ~bad(r);
        st = r(find(diff([false, g]) == 1));
        en = r(find(diff([g, false]) == -1));
        good{end + 1} = [st(:), en(:)]; %#ok<AGROW>
    end
    seg = vertcat(good{:});
    if nan_s > 0
        warning('get_eeg_data:nan', '%s: %.2f s of non-finite samples cut out (%d pieces left)', file, nan_s, size(seg, 1));
    end
end
onset = (seg(:, 1) - 1) / fs0;                            % s from the start of the file

% analysed segment (in seconds of the recording)
pieces = {};
keep_onset = [];
for k = 1:size(seg, 1)
    tk = onset(k) + (0:seg(k, 2) - seg(k, 1)) / fs0;      % time of each sample of the piece
    in = tk >= t0 - 1e-9 & tk < t1 - 1e-9;
    if nnz(in) > 1
        idx = seg(k, 1) - 1 + find(in);
        pieces{end + 1} = dat(:, idx); %#ok<AGROW>
        keep_onset(end + 1, 1) = tk(find(in, 1)); %#ok<AGROW>
    end
end
if isempty(pieces)
    error('get_eeg_data:empty', '%s: no data between %g and %g s', file, t0, t1);
end
available = sum(cellfun(@(p) size(p, 2), pieces)) / fs0;
if ~isempty(S.core_duration) && available < S.core_duration - 1 / fs0
    if nb_boundaries == 0 && nt == 1 && nan_s == 0 && flat_s == 0 && n_bad_events == 0
        error('get_eeg_data:short', '%s is %.1f s long from %g s, shorter than the %g s analysed', ...
            file, available, t0, S.core_duration);
    end
    warning('get_eeg_data:short', '%s: %.1f s of usable data between %g and %g s (cuts, epochs or bad samples)', ...
        file, available, t0, t1);
end
% units: a cleaned EEG has channel SDs of a few to a few hundred µV
sd = median(cellfun(@(p) median(std(p(use, :), 0, 2, 'omitnan')), pieces));
if sd < 0.01 || sd > 1e4
    warning('get_eeg_data:amplitude', ['%s: median channel SD %.3g µV after reading (unit %s): are the units ' ...
        'right? (power depends on them; connectivity does not)'], file, sd, strjoin(unique(unit(use)), ', '));
end

% resampling (each piece separately: the anti-aliasing filter must not cross a cut)
fs = fs0;
if ~isempty(S.resample) && abs(S.resample - fs0) > 1e-6
    [p, q] = rat(S.resample / fs0, 1e-9);
    for k = 1:numel(pieces)
        pieces{k} = resample(pieces{k}', p, q)';
    end
    fs = fs0 * p / q;
end
len = cellfun(@(x) size(x, 2), pieces);
eeg.data = [pieces{:}];
eeg.label = label;
eeg.srate = fs;
eeg.segments = [cumsum([1, len(1:end - 1)])', cumsum(len)'];
eeg.onset_s = keep_onset;
[~, name, ext] = fileparts(file);
eeg.info = struct('file', [name ext], 'format', fmt, 'srate_file', fs0, 'srate', fs, ...
    'duration_file_s', n / fs0, 'analysed_s', size(eeg.data, 2) / fs, 'n_pieces', numel(pieces), ...
    'n_boundaries', nb_boundaries, 'n_epochs', nt, 'nan_s', nan_s, 'flat_s', flat_s, 'bad_s', bad_s, ...
    'n_bad_events', n_bad_events, 'nan_channels', {nan_channels}, ...
    'unit', strjoin(unique(unit(use)), ' '), 'unit_scale', unique(scale(use))');
end

function tf = is_boundary(e, list, eeglab)
% an event of S.boundary_events (type or value, case ignored); in EEGLAB
% files also the numeric boundary type -99 (EEGLAB's boundary when the
% event types are numbers), when 'boundary' is listed
tf = matches_any(e.type, list) || matches_any(e.value, list);
if ~tf && eeglab && any(strcmpi(list, 'boundary'))
    tf = (isnumeric(e.value) && isscalar(e.value) && e.value == -99) || ...
         (ischar(e.value) && strcmp(strtrim(e.value), '-99'));
end
end

function tf = matches_any(v, list)
tf = (ischar(v) || isstring(v)) && any(strcmpi(strtrim(char(v)), list));
end

function tf = matches_pattern(v, list)
% v matches one of the patterns (case ignored, * any text)
tf = false;
if isempty(list) || ~(ischar(v) || isstring(v)) || strlength(string(v)) == 0, return; end
p = strcat('^', regexptranslate('wildcard', cellstr(list)), '$');
tf = any(~cellfun(@isempty, regexpi(strtrim(char(v)), p, 'once')));
end

function [unit, scale] = channel_units(file, fmt, hdr)
% the unit of each channel and the factor that makes it µV (1 when unknown)
nc = numel(hdr.label);
unit = repmat({'unknown'}, nc, 1);
if isfield(hdr, 'chanunit') && numel(hdr.chanunit) == nc
    unit = cellstr(string(hdr.chanunit(:)));
end
if startsWith(fmt, 'brainvision')
    % FieldTrip ignores the unit field of the .vhdr ("Ch1=Fp1,,0.1,µV") and
    % always says uV, but the resolution, which it applies, is in that unit
    [folder, stem] = fileparts(file);
    v = vhdr_units(fullfile(folder, [stem '.vhdr']), nc);
    if ~isempty(v), unit = v; end
elseif strcmp(fmt, 'biosemi_bdf') && isfield(hdr, 'orig') && isfield(hdr.orig, 'PhysDim')
    p = strtrim(cellstr(hdr.orig.PhysDim));
    if numel(p) >= nc
        unknown = strcmpi(unit, 'unknown');
        unit(unknown) = p(unknown);
    end
end
unit = strtrim(unit);
unit(cellfun(@isempty, unit)) = {'unknown'};
scale = cellfun(@microvolt_factor, unit);
scale(isnan(scale)) = 1;                                  % unknown: taken as µV
end

function f = microvolt_factor(u)
% µV per unit; NaN when the unit is not a voltage
u = lower(strtrim(char(u)));
u = regexprep(u, '^[^\x00-\x7f]+v', 'uv');                % µV, μV, in any encoding
u = regexprep(u, '^micro', 'u');
switch u
    case {'uv', 'uvolt', 'uvolts'}, f = 1;
    case {'v', 'volt', 'volts'}, f = 1e6;
    case {'mv', 'millivolt', 'millivolts'}, f = 1e3;
    case {'nv', 'nanovolt', 'nanovolts'}, f = 1e-3;
    otherwise, f = NaN;
end
end

function unit = vhdr_units(vhdr, nc)
% units of the channels of a BrainVision header ([Channel Infos], 4th field;
% microvolts when it is empty, as the format defines)
unit = {};
if ~isfile(vhdr), return; end
txt = strtrim(splitlines(string(fileread(vhdr))));
unit = repmat({'uV'}, nc, 1);
section = "";
for i = 1:numel(txt)
    if startsWith(txt(i), "["), section = lower(txt(i)); continue; end
    if section ~= "[channel infos]", continue; end     % not [Coordinates], whose lines are also Ch<n>=
    t = regexp(char(txt(i)), '^Ch(\d+)=(.*)$', 'tokens', 'once');
    if isempty(t), continue; end
    k = str2double(t{1});
    f = strsplit(t{2}, ',', 'CollapseDelimiters', false);   % Fp1,,0.1,µV: the empty reference is a field
    if k >= 1 && k <= nc && numel(f) >= 4 && ~isempty(strtrim(f{4}))
        unit{k} = strtrim(f{4});
    end
end
end
