function [y, seg, onset, tran] = bandpass_pieces(x, srate, fmin, fmax, seg_in, onset_in, S)
%BANDPASS_PIECES Band-pass and analytic signal of each continuous piece, edges discarded.
%   [y, seg, onset, tran] = bandpass_pieces(x, srate, fmin, fmax, seg_in, onset_in, S)
%   Each continuous piece (rows of seg_in, samples of x; onset_in: its time
%   in the recording) is filtered separately and turned into its analytic
%   signal (Hilbert transform), THEN its edges are discarded (at least
%   S.edge_trim s, and at least the filter's transient; E5), so that the
%   edge distortion of both the filter and the Hilbert transform falls in
%   the discarded part; pieces shorter than S.min_segment s after trimming
%   are dropped. The transient is the half-length of the impulse response
%   that holds 99.9% of its energy (V3-11): Brainstorm's 99% leaves 1% of
%   the energy, about 10% of the amplitude, outside it, so long low-band
%   filters left errors of the analytic signal in short pieces (95th
%   percentile, 20 s pieces: delta 17% -> 6%, theta 7% -> 2%, alpha 1.7% ->
%   0.9%; beta1 and above are trimmed by S.edge_trim, unchanged). y is complex (real(y) is the band-passed signal); seg: the
%   pieces kept, in y; onset: their first sample's time in the recording;
%   tran: the transition width used (Hz; [] = Brainstorm's default).
%   Used by run_pipeline (stage 1) and by the surrogate test of stage 2.
tran = filter_transition(S, fmin, fmax, srate);
pieces = {};
onset = [];
for k = 1:size(seg_in, 1)
    xk = x(:, seg_in(k, 1):seg_in(k, 2));
    if size(xk, 2) / srate < S.min_segment + 2 * S.edge_trim, continue; end
    [yk, transient] = bandpass_filter(xk, srate, fmin, fmax, tran);
    trim = round(max(S.edge_trim, transient) * srate);
    if size(yk, 2) - 2 * trim < S.min_segment * srate, continue; end
    yk = hilbert(yk.').';
    pieces{end + 1} = yk(:, trim + 1:end - trim); %#ok<AGROW>
    onset(end + 1, 1) = onset_in(k) + trim / srate; %#ok<AGROW>
end
if isempty(pieces)
    error('run_pipeline:short', ['No continuous piece of at least %g s in %g-%g Hz after removing %g s, or the ' ...
        'filter''s transient when longer, at each edge; the longest piece is %.1f s (an epoched file is cut at ' ...
        'every epoch: give the continuous recording, or longer epochs)'], S.min_segment, fmin, fmax, S.edge_trim, ...
        max(diff(seg_in, 1, 2) + 1) / srate);
end
len = cellfun(@(p) size(p, 2), pieces);
y = [pieces{:}];
seg = [cumsum([1, len(1:end - 1)])', cumsum(len)'];
end

function [x, transient] = bandpass_filter(x, srate, fmin, fmax, tran)
% Brainstorm band-pass, current design ('bst-hfilter-2019', the default of
% Brainstorm's process_bandpass since 2019): even-order linear-phase FIR,
% 0 dB in [fmin, fmax], 60 dB stop band TRAN Hz beyond each edge, delay
% compensated, mirrored edges. The original called the legacy
% bst_bandpass_filtfilt, whose transition band is fixed at 5% of Nyquist
% (12.8 Hz at 512 Hz): its 13-20 Hz filter (order 68) lets alpha through
% (-3.5 dB at 10 Hz, -7.7 dB at 8 and 25 Hz; measured in development).
% The transition width sets the filter's length, i.e. its time resolution
% (beta1 at 512 Hz; a 6-cycle window is 0.37 s):
%   1 Hz (Brainstorm's default): order 1858 (3.6 s), 99% of its energy
%     within 0.86 s (2.3 windows); -72 dB at 10 Hz
%   1.25 Hz (S.filter_transition 'window', S.filter_span 2, the default):
%     99% within 0.69 s (1.9 windows); 9% of the output power outside 13-20 Hz
%   3.25 Hz (EEGLAB's rule, 25% of the lower edge, 'eeglab'): order 572
%     (1.1 s), 99% within 0.34 s (one window); -38 dB at 10 Hz, -11 dB at
%     11 Hz; 31% of the output power outside 13-20 Hz
%   (every band: docs/settings.md, the table of the filters)
[x, FiltSpec] = bst_bandpass_hfilter(x, srate, fmin, fmax, 1, 0, [], tran, 'bst-hfilter-2019');
b = FiltSpec.b(:);                                % half-length with 99.9% of the energy (V3-11)
h = b(floor(numel(b) / 2) + 1:end) .^ 2;
E = cumsum(h) / sum(h);
transient = (find(E >= 0.999, 1) - 1) / srate;    % seconds
end

function tran = filter_transition(S, fmin, fmax, srate)
% width of the filter's transition bands (Hz) from S.filter_transition:
% 'window' = the narrowest transition whose filter keeps 99% of its impulse
%   response energy within S.filter_span analysis windows of the band (6
%   cycles of the band centre, or S.window_length): as band-specific as the
%   time resolution of the windows allows (docs/settings.md);
% 'eeglab' = 25% of the lower edge, at least 2 Hz, at most fmin and the room
%   below Nyquist (EEGLAB's pop_eegfiltnew rule);
% 'brainstorm' = Brainstorm's default ([]: 0.5 Hz below 5 Hz, else 1 Hz);
% a number = that width
tran = S.filter_transition;
if ischar(tran) || isstring(tran)
    switch lower(char(tran))
        case 'window'
            if isfield(S, 'window_length') && ~isempty(S.window_length)
                win = S.window_length;
            else
                win = S.window_cycles / (0.5 * (fmin + fmax));
            end
            tran = window_matched_transition(srate, fmin, fmax, S.filter_span * win);
        case 'eeglab', tran = min([max(0.25 * fmin, 2), fmin, srate / 2 - fmax]);
        case 'brainstorm', tran = [];
        otherwise, error('run_pipeline:filter', ...
                'S.filter_transition: ''window'', ''eeglab'', ''brainstorm'' or a width in Hz');
    end
end
end

function tran = window_matched_transition(srate, fmin, fmax, span)
% narrowest transition (0.01 Hz steps, at least 0.1 Hz, at most fmin and the
% room below Nyquist) whose filter holds 99% of its energy within SPAN
% seconds; bisection, as the filter shortens when the transition widens.
% Computed once per band and rate (cached).
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%.6g_%.6g_%.6g_%.6g', srate, fmin, fmax, span);
if isKey(cache, key), tran = cache(key); return; end
lo = 0.1;
hi = min(fmin, srate / 2 - fmax);
if energy_span(srate, fmin, fmax, hi) > span
    tran = hi;                                    % even the widest does not fit
elseif energy_span(srate, fmin, fmax, lo) <= span
    tran = lo;
else
    while hi - lo > 0.01
        mid = (lo + hi) / 2;
        if energy_span(srate, fmin, fmax, mid) <= span, hi = mid; else, lo = mid; end
    end
    tran = hi;
end
cache(key) = tran;
end

function s = energy_span(srate, fmin, fmax, tran)
% span (s) of the filter's impulse response holding 99% of its energy
[~, F] = bst_bandpass_hfilter(zeros(1, 8), srate, fmin, fmax, 0, 0, [], tran, 'bst-hfilter-2019');
e = cumsum(F.b(:) .^ 2) / sum(F.b .^ 2);
s = (find(e > 0.995, 1) - find(e > 0.005, 1)) / srate;
end
