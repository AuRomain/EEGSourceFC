function [cmats, static_cmats, static_full_cmats, info, static_power] = get_all_connectivity_mats(src_data, srate, fmin, fmax, seg, opts)
%GET_ALL_CONNECTIVITY_MATS Windowed and whole-recording connectivity of band-passed ROI signals.
%   [cmats, static_cmats, static_full_cmats, info, static_power] =
%   get_all_connectivity_mats(src_data, srate, fmin, fmax, seg, opts)
%   src_data  nROI x nSamples, band-passed: real, or its analytic signal
%             (complex, computed by run_pipeline before the edges were trimmed)
%   seg       nPieces x 2 continuous pieces (default: one piece)
%   opts      conn (default {'plv', 'wpli', 'aec'}, the original's three),
%             window_cycles (6), window_length ([] = cycles), window_overlap
%             (0.5), aec_window (6 s), aec_step (0.5 s), measure_windows
%             (struct(): the same settings for one measure, see below),
%             dynamic (true; false: static only), static_matrices ({}; with
%             'window_mean', the windows of a static-only run): from
%             pipeline_settings
%   Same outputs as the original: cmats{k} = nWindows x nROI x nROI for the
%   k-th measure of opts.conn (cmats{1..3} = PLV, wPLI, AEC by default);
%   static_cmats(k,:,:) = the mean of the windows (the original "static"
%   definition). Windows: 6 cycles of the band centre with 50% overlap for
%   the phase measures, 6 s every 0.5 s for AEC (O'Neill et al. 2017), as the
%   original; or a fixed length for every band (opts.window_length).
%   Measures: 'plv', 'wpli', 'aec' (original), 'pli', 'ciplv', 'aec_orth'
%   (pairwise leakage correction), 'aec_symorth', 'plv_symorth' (symmetric
%   correction), the variants of the original +connectivity folder; 'dwpli'
%   (debiased squared wPLI; the band's decorrelation time 1 / bandwidth sets
%   its blocks: for long windows, S.window_length).
%   Corrections and additions: one Hilbert transform per continuous piece
%   instead of one per window (O4); windows never cross a cut; zero
%   diagonals (E10, D3); signed AEC (E11); static_full_cmats = the same
%   measures over the whole recording (the window mean of wPLI is dominated
%   by the short-window bias, E8), by static_connectivity, and
%   static_power the band power of each region (nROI x 1); single
%   precision (O5); info holds the window sizes and positions.
%   opts.dynamic = false (static only, S.dynamic): no windows; cmats is {},
%   static_cmats [] and info.windows {}; static_full_cmats and static_power
%   are the same as with windows. A static-only run that asks for the
%   window mean (opts.static_matrices with 'window_mean', 6 November 2026)
%   computes the windows as S.dynamic does: cmats, static_cmats and the
%   window sizes as above (run_pipeline does not save cmats then).
%   Windows of one measure (opts.measure_windows, 6 November 2026): a field
%   named after the measure holds the window settings it uses instead of
%   the common ones (window_cycles, window_length, window_overlap for the
%   phase measures; aec_window, aec_step for AEC), e.g.
%   struct('wpli', struct('window_length', 2)); window_cycles given alone
%   means cycles (window_length []). Each measure's windows are in
%   info.windows{k}; info.window_s, step_s, aec_window_s and aec_step_s are
%   those of the first phase and the first AEC measure (as
%   plv_wpli_windows and aec_windows), and info.measure_windows the
%   settings given (only when some are).
%   A measure none of whose windows fits in a piece (every piece shorter
%   than its window, e.g. the 6 s AEC windows in pieces of 4-6 s) has no
%   window mean: static_cmats is NaN (zero diagonal) for it, and
%   info.no_windows.(measure) says why.
%   A measure that cannot be computed for this recording (the symmetric
%   orthogonalisation with fewer independent region signals than regions;
%   dwPLI with fewer than 4 blocks) is NaN (zero diagonal) in its windows
%   and matrices, and info.not_computed.(measure) says why; the other
%   measures are computed (run_pipeline warns, naming the recording). It
%   was an error that lost every measure and band of the recording.
if nargin < 5 || isempty(seg), seg = [1, size(src_data, 2)]; end
if nargin < 6, opts = struct(); end
def = struct('conn', {{'plv', 'wpli', 'aec'}}, 'window_cycles', 6, 'window_length', [], ...
    'window_overlap', 0.5, 'aec_window', 6, 'aec_step', 0.5, 'measure_windows', struct(), 'dynamic', true, ...
    'static_matrices', {{}});
for f = fieldnames(def)', if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end, end

aec = {'aec', 'aec_orth', 'aec_symorth'};
cfg = cell(1, numel(opts.conn));                       % the windows of each measure
for k = 1:numel(opts.conn)
    method = 'default';
    if ismember(opts.conn{k}, aec), method = 'aec'; end
    cfg{k} = cfg_connectivity(srate, fmin, fmax, method, measure_opts(opts, opts.conn{k}));
end
Z = connectivity.analytic_signal(src_data, seg);       % unchanged when already complex; every
n = size(Z, 1);                                        % measure starts from it (Hilbert once)
nc = numel(opts.conn);
cmats = cell(1, nc);
static_cmats = zeros(nc, n, n);
windows = cell(1, nc);
plv_wpli = [];
not_computed = struct();
no_windows = struct();
% windows: for the states (S.dynamic), or for the window mean of a
% static-only run; none otherwise
use_windows = opts.dynamic || any(strcmp(cellstr(opts.static_matrices), 'window_mean'));
for k = 1:nc * use_windows
    c = cfg{k};
    switch opts.conn{k}
        case {'plv', 'wpli'}
            if isempty(plv_wpli) || ~isequal(plv_wpli.size, [c.window, c.step])   % both from one pass
                [plv_wpli.plv, plv_wpli.s_plv, plv_wpli.wpli, plv_wpli.s_wpli, plv_wpli.f_plv, plv_wpli.f_wpli, plv_wpli.win] = ...
                    connectivity.plv_wpli_sliding_window(Z, c.srate, c.window, c.step, seg);
                plv_wpli.size = [c.window, c.step];   % (one pass for each if their windows differ)
            end
            m = opts.conn{k};
            [cmats{k}, s, w] = deal(plv_wpli.(m), plv_wpli.(['s_' m]), plv_wpli.win);
        case 'pli'
            [cmats{k}, s, ~, w] = connectivity.pli_sliding_window(Z, c.srate, c.window, c.step, seg);
        case 'dwpli'
            [cmats{k}, s, ~, w] = connectivity.dwpli_sliding_window(Z, c.srate, c.window, c.step, seg, 1 / (fmax - fmin));
        case 'ciplv'
            [cmats{k}, s, ~, w] = connectivity.ciplv_sliding_window(Z, c.srate, c.window, c.step, seg);
        case 'plv_symorth'
            [cmats{k}, s, ~, w, note] = connectivity.plv_sliding_window_corrected(Z, c.srate, c.window, c.step, seg);
            if ~isempty(note), not_computed.(opts.conn{k}) = note; end
        case 'aec'
            [cmats{k}, s, ~, w] = connectivity.aec_sliding_window(Z, c.srate, c.window, c.step, seg);
        case 'aec_orth'
            [cmats{k}, s, ~, w] = connectivity.aec_sliding_window_corrected_pairwise(Z, c.srate, c.window, c.step, seg);
        case 'aec_symorth'
            [cmats{k}, s, ~, w, note] = connectivity.aec_sliding_window_corrected(Z, c.srate, c.window, c.step, seg);
            if ~isempty(note), not_computed.(opts.conn{k}) = note; end
        otherwise
            error('get_all_connectivity_mats:conn', 'Unknown connectivity measure %s', opts.conn{k});
    end
    static_cmats(k, :, :) = s;
    windows{k} = w;
    if isempty(w.first_sample)                         % no window fits in any piece: no window mean
        A = nan(n);  A(1:n + 1:end) = 0;
        static_cmats(k, :, :) = A;
        no_windows.(opts.conn{k}) = sprintf(['no window of %.3g s fits in the pieces of the recording ' ...
            '(the longest is %.3g s): no window mean'], w.samples / srate, max(diff(seg, 1, 2) + 1) / srate);
    end
end
% the whole recording: every measure over all the samples of the pieces,
% and the band power of each region (with or without windows)
[static_full_cmats, static_power, notes] = static_connectivity(Z, seg, opts.conn, srate, fmin, fmax);
for k = find(~cellfun(@isempty, notes))
    if ~isfield(not_computed, opts.conn{k}), not_computed.(opts.conn{k}) = notes{k}; end
end
if ~use_windows
    cmats = {};  static_cmats = [];  windows = {};
end
% window sizes: those of the first phase and the first AEC measure (the
% common settings when the measures have none of their own; a list without
% phase or AEC measures has the common ones, as before)
iw = find(ismember(opts.conn, {'plv', 'wpli', 'dwpli', 'pli', 'ciplv', 'plv_symorth'}), 1);
ia = find(ismember(opts.conn, aec), 1);
cw = cfg_connectivity(srate, fmin, fmax, 'default', opts);
ca = cfg_connectivity(srate, fmin, fmax, 'aec', opts);
if ~isempty(iw), cw = cfg{iw}; end
if ~isempty(ia), ca = cfg{ia}; end
info = struct('conn', {opts.conn}, 'fmin', fmin, 'fmax', fmax, 'srate', srate, ...
    'window_s', cw.window, 'step_s', cw.step, 'aec_window_s', ca.window, 'aec_step_s', ca.step, ...
    'windows', {windows}, 'segments', seg, 'not_computed', not_computed);
if ~use_windows                                        % no windows: no window sizes
    info = rmfield(info, {'window_s', 'step_s', 'aec_window_s', 'aec_step_s'});
else
    info.no_windows = no_windows;
    if isstruct(opts.measure_windows) && ~isempty(fieldnames(opts.measure_windows))
        info.measure_windows = opts.measure_windows;
    end
end
% original field names, read by stage 2
if ~isempty(iw) && use_windows, info.plv_wpli_windows = windows{iw}; end
if ~isempty(ia) && use_windows, info.aec_windows = windows{ia}; end
end


%%% ========== WINDOW SETTINGS OF ONE MEASURE ========== %%%
function o = measure_opts(opts, measure)
% the common window settings, with those given for this measure instead
% (opts.measure_windows.(measure)); window_cycles alone: cycles
o = opts;
if ~isstruct(opts.measure_windows) || ~isfield(opts.measure_windows, measure), return; end
mw = opts.measure_windows.(measure);
if isfield(mw, 'window_cycles') && ~isfield(mw, 'window_length'), o.window_length = []; end
for f = fieldnames(mw)'
    o.(f{1}) = mw.(f{1});
end
end


%%% ========== GET CONFIGURATION STRUCT ========== %%%
function cfg = cfg_connectivity(srate, fmin, fmax, method, opts)
cfg = [];
cfg.srate = srate;
cfg.fmin = fmin;
cfg.fmax = fmax;
cfg.conn_meth = method;
switch method
    case "aec"
        cfg.window = opts.aec_window;                % s (O'Neill et al. 2017)
        cfg.step = opts.aec_step;
    otherwise
        if isempty(opts.window_length)
            cfg.window = opts.window_cycles / (0.5 * (fmin + fmax));   % 6 cycles of the band centre
        else
            cfg.window = opts.window_length;         % fixed length (s) for every band
        end
        cfg.step = (1 - opts.window_overlap) * cfg.window;   % 50% overlap by default
end
end
