function D = load_static(S, inv, band, which_static)
%LOAD_STATIC The static matrices of one inverse method and band, for every recording in the groups.
%   D = load_static(S, inv, band) returns
%     subs, ids, group   recording (file name without extension), subject ID,
%                        index into S.groups
%     files              static_cmats_<inv>_<montage>_<atlas>_<band>.mat of each
%     measures           S.conn, the order of the first dimension of conn
%     conn               nMeasures x nROI x nROI x nRecordings: every measure
%                        over the whole recording (static_full_cmats)
%     halves             2 x nMeasures x nROI x nROI x nRecordings: the same on
%                        the odd and on the even blocks (static_halves)
%     power              nROI x nRecordings band power (static_power);
%                        power_halves nROI x 2 x nRecordings
%     analysed_s         seconds analysed per recording; halves_s nRecordings x 2
%     band_hz            nRecordings x 2 (individual bands differ per recording)
%     roi_labels         region labels
%     halves_why         nMeasures x nRecordings text: why a recording has no
%                        halves of a measure ("" when it has both): outputs
%                        made before they existed, by earlier code, a half
%                        too short, the measure not computed
%     halves_block       nMeasures x 1 block length of the halves (s)
%   D = load_static(S, inv, band, 'window_mean'): conn is the mean of the
%   windows instead (static_cmats; in a static-only run when
%   S.static_matrices asks for it), halves the window mean on each half
%   (static_halves_window, in blocks of max(S.reliability_block, 3 x the
%   measure's window), static_halves_block), power that of the whole
%   recording. Stage 1 outputs without the window mean's halves have none;
%   in stage 1 outputs without static_halves_block, every measure's
%   windows are in blocks of S.reliability_block: the halves of a measure whose
%   window is longer than a third of a block (AEC's 6 s, dwPLI's 20 s:
%   under a quarter of the windows per half, or none) are left out (NaN),
%   halves_why saying to run stage 1 again. D.matrix: which_static.
%   Recordings and checks as in the states and their statistics (stage1_outputs): those of
%   S.data in the groups with stage 1 outputs, made with the current stage 1
%   settings, code and inputs; but every such recording, also one with few
%   windows of S.state_conn or without it (stage1_outputs(..., 'static')):
%   each measure keeps its own recordings.
if nargin < 4 || isempty(which_static), which_static = 'full'; end
if ~ismember(which_static, {'full', 'window_mean'})
    error('load_static:which', 'which_static is full or window_mean');
end
if strcmp(which_static, 'window_mean') && ~S.dynamic && ...
        ~(isfield(S, 'static_matrices') && any(strcmp(cellstr(S.static_matrices), 'window_mean')))
    error('load_static:which', ['A static-only run (S.dynamic = false) has windows only when ' ...
        'S.static_matrices includes ''window_mean'' (stage 1 computes them then)']);
end
[subs, ids, group, files] = stage1_outputs(S, inv, band, 'static');   % checked: same settings, code, inputs
for i = 1:numel(files)
    [folder, name] = fileparts(files{i});
    if ~startsWith(name, 'static_'), name = ['static_' name]; end   % the static file of cmats_*.mat
    files{i} = fullfile(folder, [name '.mat']);
end
n = numel(files);
for i = 1:n
    w = load(files{i});
    if strcmp(which_static, 'full')
        M = w.static_full_cmats;
        H = w.static_halves;
        why = strings(numel(w.info.conn), 1);
        block = repmat(S.reliability_block, numel(w.info.conn), 1);
    else
        M = w.static_cmats;
        if isempty(M)
            error('load_static:which', '%s has no window mean (made without windows): run stage 1 again', files{i});
        end
        [H, why, block] = window_halves(w, S);
    end
    why = missing_halves(H, why, block, w.info);
    if i == 1
        D = struct('subs', {subs}, 'ids', ids, 'group', group, 'files', {files}, 'matrix', which_static, ...
            'measures', {w.info.conn}, 'conn', zeros([size(M), n]), 'halves', zeros([size(H), n]), ...
            'halves_why', strings(numel(w.info.conn), n), 'halves_block', nan(numel(w.info.conn), 1), ...
            'power', zeros(numel(w.static_power), n), 'power_halves', zeros(numel(w.static_power), 2, n), ...
            'analysed_s', zeros(n, 1), 'halves_s', zeros(n, 2), 'band_hz', zeros(n, 2), ...
            'roi_labels', {w.info.roi_labels});
    elseif ~isequal(w.info.roi_labels, D.roi_labels) || ~isequal(w.info.conn, D.measures)
        error('load_static:regions', '%s has other regions or measures than %s', files{i}, files{1});
    end
    D.conn(:, :, :, i) = M;
    D.halves(:, :, :, :, i) = H;
    D.halves_why(:, i) = why;
    D.halves_block = max(D.halves_block, block(:), 'omitnan');   % the same for every recording, but of individual bands
    D.power(:, i) = w.static_power;
    D.power_halves(:, :, i) = w.static_power_halves;
    D.analysed_s(i) = w.info.analysed_s;
    D.halves_s(i, :) = w.info.halves_s;
    D.band_hz(i, :) = w.info.band_hz;
end
m = numel(D.roi_labels);                           % zero diagonal, as load_static_matrices
C = reshape(D.conn, size(D.conn, 1), m * m, n);
C(:, 1:m + 1:end, :) = 0;
D.conn = reshape(C, size(D.conn));
end

function [H, why, block] = window_halves(w, S)
% the window mean on each half (static_halves_window), the block of each
% measure (static_halves_block), and why a measure has none
nm = numel(w.info.conn);
why = strings(nm, 1);
block = nan(nm, 1);
H = nan(size(w.static_halves));
if ~isfield(w, 'static_halves_window')
    why(:) = "stage 1 outputs made before the window mean's halves existed: run stage 1 again";
    return
end
H = w.static_halves_window;
if isfield(w, 'static_halves_block')
    block = w.static_halves_block(:);
else                                              % no static_halves_block: every measure in blocks of S.reliability_block
    ns = sum(cellfun(@(p) sum(diff(p, 1, 2) + 1), w.info.halves));   % the halves hold every sample
    win = cellfun(@(x) x.samples, w.info.windows(:)) * w.info.analysed_s / ns;   % window length (s)
    block = repmat(S.reliability_block, nm, 1);
    old = 3 * win > S.reliability_block;
    H(:, old, :, :) = NaN;
    why(old) = sprintf("stage 1 outputs of earlier code (windows in blocks of %g s, few or none in a half): " + ...
        "run stage 1 again", S.reliability_block);
end
end

function why = missing_halves(H, why, block, info)
% why a measure without a reason yet has a half without values: not
% computed for this recording, or a half without data for it
for k = find(why(:)' == "")
    if all(arrayfun(@(h) any(isfinite(H(h, k, :, :)), 'all'), 1:2)), continue; end
    if isfield(info, 'not_computed') && isfield(info.not_computed, info.conn{k})
        why(k) = "not computed for this recording";
    else
        why(k) = sprintf("a half without data for it (too short for blocks of %g s)", block(k));
    end
end
end
