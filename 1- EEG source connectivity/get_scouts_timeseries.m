function [scouts_data, labels] = get_scouts_timeseries(src_data, atlas, Cx, scout_function, X, A)
%GET_SCOUTS_TIMESERIES Brainstorm scout values: one value per ROI from its sources (mean, PCA, max, power, RMS).
%   [scouts_data, labels] = get_scouts_timeseries(src_data, atlas) with
%   src_data nSources x nTimes, or the inverse filters (nSources x
%   nChannels): the scout mean with orientation-based sign flip is linear
%   and does not depend on the data, so applying it to the filters gives an
%   nROI x nChannels operator (scouts = operator * EEG) and avoids the
%   15,002 x 92,160 source matrix of the original (E4, O1). Same Brainstorm
%   call as the original (bst_scout_value, 'mean', normals, sign flip).
%   The sign flip is now meaningful: sources have fixed normal
%   orientations (E2, E12).
%   atlas: 'Desikan_Killiany_68', a Brainstorm scout file (Scouts with
%   Vertices and Label) defined on the cortex, or the name of one of the
%   cortex's atlases ('Destrieux', 'Schaefer_100_17net', ...): get_atlas.
%   The cortex is the shipped ICBM152's (tess_cortex_pial_low, 15,002
%   vertices) or that of the head A (get_anatomy; S.anatomy), given last:
%   get_scouts_timeseries(src_data, atlas, Cx, scout_function, X, A).
%   labels: the ROI labels (Scouts.Label), in the order of the rows.
%   Free orientation (S.orientation = 'free'): src_data holds three filters
%   per source (3 x 15,002 rows, x/y/z of each source) and Cx is the
%   channel covariance of the recording (broadband, average reference). Each
%   region is then the first principal component of its sources' three
%   orientations (the weights w maximise w' F C F' w, F the region's
%   filters), with the sign that makes the weights point along the cortical
%   normals on average, so that the sign convention matches the normal
%   orientation's; the operator is still one row per region.
%
%   scout_function (S.scout_function; default 'mean'), Brainstorm's scout
%   functions (bst_scout_value, sign flip on), sources along the normal:
%     'mean'   the mean of the region's sources, each source's sign flipped
%              to the region's dominant orientation (above)
%     'pca'    the first principal component of the region's (sign-flipped)
%              sources: Brainstorm's PCA ('pca2023', from the sources'
%              covariance, F C F' with Cx the channel covariance), sign so
%              that its weights sum to a positive number, scaled by
%              1/sqrt(number of sources) like the mean. One set of weights
%              per recording (from the whole analysed recording)
%     'max'    at each sample, the (sign-flipped) source of largest
%              absolute value, with its sign
%     'power'  at each sample, the mean of the squared sources
%     'rms'    the square root of 'power'
%   'mean' and 'pca' are linear: applied to the filters they give the
%   operator (scouts = operator * EEG), as above. 'max', 'power' and 'rms'
%   are not: they need the source time series, which are computed region by
%   region from the filters and X (channels x samples), in pieces of time
%   (never all 15,002 sources at once):
%   scouts_data = get_scouts_timeseries(filters, atlas, Cx, scout_function, X)
%   returns the regions' time series (nROI x nSamples), for any function.
%   As in Brainstorm's connectivity processes (scout function applied
%   "before"), run_pipeline applies the function to the broadband source
%   signals and band-passes the result: for 'power' and 'rms' (no sign) the
%   band-passed signal is the band's part of the fluctuations of the
%   region's power, not the region's oscillation (check_settings warns
%   when phase measures are asked for). With the free orientation only
%   'pca' is defined (the principal component above).
if nargin < 4 || isempty(scout_function), scout_function = 'mean'; end
scout_function = lower(char(scout_function));
if ~ismember(scout_function, {'mean', 'pca', 'max', 'power', 'rms'})
    error('get_scouts_timeseries:function', 'Unknown scout function %s (mean, pca, max, power, rms)', scout_function);
end
linear = ismember(scout_function, {'mean', 'pca'});
with_data = nargin >= 5 && ~isempty(X);
if nargin < 6, A = []; end
if ~linear && ~with_data
    error('get_scouts_timeseries:data', ['Scout function ''%s'' is not linear: give the filters and the data ' ...
        '(get_scouts_timeseries(filters, atlas, Cx, ''%s'', X))'], scout_function, scout_function);
end
[Scouts, C] = get_atlas(atlas, A);
VertNormals = C.VertNormals;
labels = {Scouts.Label};
if size(src_data, 1) == 3 * size(VertNormals, 1)
    if ~strcmp(scout_function, 'pca') && nargin >= 4
        error('get_scouts_timeseries:free', 'Free orientation: each region is its first principal component (scout function ''pca'')');
    end
    if nargin < 3 || isempty(Cx), error('get_scouts_timeseries:free', 'Free orientation: give the channel covariance Cx'); end
    scouts_data = zeros(numel(Scouts), size(src_data, 2));
    for i = 1:numel(Scouts)
        v = Scouts(i).Vertices(:)';
        rows = reshape(3 * v + (-2:0)', 1, []);      % x, y, z of each vertex
        F = src_data(rows, :);
        M = F * Cx * F';
        [w, ~] = eigs((M + M') / 2, 1);
        nrm = VertNormals(v, :)';
        if sum(w .* nrm(:)) < 0, w = -w; end      % along the normals on average
        scouts_data(i, :) = w' * F;
    end
    if with_data, scouts_data = scouts_data * X; end
    return
end
if strcmp(scout_function, 'pca') && (nargin < 3 || isempty(Cx))
    error('get_scouts_timeseries:pca', 'Scout function ''pca'': give the channel covariance Cx');
end
if linear                                          % one operator row per region
    scouts_data = zeros(numel(Scouts), size(src_data, 2));
    for i = 1:numel(Scouts)
        F = src_data(Scouts(i).Vertices, :);
        Orient = VertNormals(Scouts(i).Vertices, :);
        if strcmp(scout_function, 'mean')
            [~, v] = evalc('bst_scout_value(F, ''mean'', Orient, 1, ''none'', 1)');  % quiet "Flipped the sign" messages
            scouts_data(i, :) = v;
        else                                       % Brainstorm's PCA from the sources' covariance
            Covar = F * Cx * F';
            [~, ~, w] = evalc('bst_scout_value([], ''pca2023'', Orient, 1, ''none'', 1, '''', Covar)');
            scouts_data(i, :) = w' * F;            % w includes the sign flip and the 1/sqrt(n) scaling
        end
    end
    if with_data, scouts_data = scouts_data * X; end
    return
end
% max, power, rms: the region's source time series, in pieces of time
nT = size(X, 2);
scouts_data = zeros(numel(Scouts), nT);
for i = 1:numel(Scouts)
    F = src_data(Scouts(i).Vertices, :);
    Orient = VertNormals(Scouts(i).Vertices, :);
    step = max(1, floor(2e7 / size(F, 1)));        % at most 2e7 values (160 MB) at a time
    for t0 = 1:step:nT
        t = t0:min(nT, t0 + step - 1);
        Y = F * double(X(:, t));
        [~, v] = evalc('bst_scout_value(Y, scout_function, Orient, 1, ''none'', 1)');
        scouts_data(i, t) = v;
    end
end
end
