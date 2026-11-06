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
%   D = load_static(S, inv, band, 'window_mean'): conn is the mean of the
%   windows instead (static_cmats, the original "static" matrix; in a
%   static-only run when S.static_matrices asks for it, 6 November 2026),
%   halves the window mean on each half (static_halves_window; NaN for
%   outputs made before it existed, which have no reliability of the window
%   mean), power that of the whole recording. D.matrix: which_static.
%   Recordings and checks as in stages 2 and 3 (stage1_outputs): those of
%   S.data in the groups with stage 1 outputs, made with the current stage 1
%   settings, code and inputs.
if nargin < 4 || isempty(which_static), which_static = 'full'; end
if ~ismember(which_static, {'full', 'window_mean'})
    error('load_static:which', 'which_static is full or window_mean');
end
if strcmp(which_static, 'window_mean') && ~S.dynamic && ...
        ~(isfield(S, 'static_matrices') && any(strcmp(cellstr(S.static_matrices), 'window_mean')))
    error('load_static:which', ['A static-only run (S.dynamic = false) has windows only when ' ...
        'S.static_matrices includes ''window_mean'' (stage 1 computes them then)']);
end
[subs, ids, group, files] = stage1_outputs(S, inv, band);   % checked: same settings, code, inputs
for i = 1:numel(files)
    [folder, name] = fileparts(files{i});
    if ~startsWith(name, 'static_'), name = ['static_' name]; end   % the static file of cmats_*.mat
    files{i} = fullfile(folder, [name '.mat']);
end
n = numel(files);
for i = 1:n
    w = load(files{i});
    H = w.static_halves;
    if strcmp(which_static, 'full')
        M = w.static_full_cmats;
    else
        M = w.static_cmats;
        if isempty(M)
            error('load_static:which', '%s has no window mean (made without windows): run stage 1 again', files{i});
        end
        H = nan(size(w.static_halves));           % outputs made before the window mean's halves existed
        if isfield(w, 'static_halves_window'), H = w.static_halves_window; end
    end
    if i == 1
        D = struct('subs', {subs}, 'ids', ids, 'group', group, 'files', {files}, 'matrix', which_static, ...
            'measures', {w.info.conn}, 'conn', zeros([size(M), n]), 'halves', zeros([size(H), n]), ...
            'power', zeros(numel(w.static_power), n), 'power_halves', zeros(numel(w.static_power), 2, n), ...
            'analysed_s', zeros(n, 1), 'halves_s', zeros(n, 2), 'band_hz', zeros(n, 2), ...
            'roi_labels', {w.info.roi_labels});
    elseif ~isequal(w.info.roi_labels, D.roi_labels) || ~isequal(w.info.conn, D.measures)
        error('load_static:regions', '%s has other regions or measures than %s', files{i}, files{1});
    end
    D.conn(:, :, :, i) = M;
    D.halves(:, :, :, :, i) = H;
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
