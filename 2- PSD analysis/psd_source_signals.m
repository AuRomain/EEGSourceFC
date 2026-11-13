function [src, roi_labels, A, filt, vert] = psd_source_signals(X, idx, srate, recording, S)
%PSD_SOURCE_SIGNALS The regions' time series of one recording, as stage 1 computes them.
%   [src, roi_labels, A, filt, vert] = psd_source_signals(X, idx, srate, recording, S)
%   X          the recording's channels used, average-referenced: rows idx
%              of the cap S.montages{1} (get_eeg_montage, flat channels
%              left out), as psd_recording reads them
%   recording  its file name (the subject's head with S.anatomy_individual)
%   src        1 x numel(S.inv) cell: regions x samples for each inverse
%              method of S.inv, on the atlas S.atlases{1}
%   roi_labels the regions; A the head (get_anatomy)
%   The steps of stage 1 (run_pipeline), with its functions and settings:
%   the head (get_anatomy), head model and leadfield of the cap
%   (get_ftHeadmodel_and_ftLeadfield, computed once and kept), the
%   electrodes and leadfield rows of the channels used with the leadfield
%   average-referenced over them, the inverse filters (get_inverse_solution:
%   eLORETA, LCMV or wMNE with S.eloreta_lambda, S.lcmv_lambda, S.wmne and
%   S.orientation; identity noise covariance), and each region's signal
%   (get_scouts_timeseries with S.scout_function: the sign-flipped mean by
%   default). The same as stage 1's Scouts_<inverse>_<cap>_<atlas>.mat
%   (src_data; tests/test_psd.m checks it), without its band-pass filters.
%   filt       1 x numel(S.inv) cell: the inverse filters of the sources
%              (sources x channels; with S.orientation = 'free' three rows
%              per source, its x, y and z), so source j's signal is
%              filt{iv}(j, :) * X
%   vert       1 x number of regions cell: the sources (vertices of the
%              cortex) of each region, in the order of roi_labels (get_atlas)
montage = S.montages{1};
atlas = S.atlases{1};
S_ma = S;
S_ma.montages = {montage};
S_ma.atlases = {atlas};
A = get_anatomy(S_ma, recording);
ftElec = get_ftElec(montage);
[ftHeadmodel, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(montage, atlas, S, A);
C = load(A.cortex, 'VertNormals');
[ftElec, ftLeadfield] = keep_channels(ftElec, ftLeadfield, idx);
filters = get_inverse_solution(X, srate, eye(numel(idx)), ftHeadmodel, ftLeadfield, ftElec, C.VertNormals, S.inv, S);
Xc = double(X) - mean(double(X), 2);
Cx = Xc * Xc' / size(X, 2);                           % channel covariance (PCA, free orientation)
clear Xc
src = cell(1, numel(S.inv));
for iv = 1:numel(S.inv)
    if ismember(S.scout_function, {'mean', 'pca'})
        [scout_filters, roi_labels] = get_scouts_timeseries(squeeze(filters(iv, :, :)), atlas, Cx, S.scout_function, [], A);
        src{iv} = scout_filters * X;
    else
        [src{iv}, roi_labels] = get_scouts_timeseries(squeeze(filters(iv, :, :)), atlas, Cx, S.scout_function, X, A);
    end
end
if nargout > 3
    filt = arrayfun(@(iv) reshape(filters(iv, :, :), size(filters, 2), size(filters, 3)), 1:numel(S.inv), ...
        'UniformOutput', false);
    Sc = get_atlas(atlas, A);                     % the regions as get_scouts_timeseries reads them
    if ~isequal({Sc.Label}, cellstr(roi_labels(:)'))
        error('psd_source_signals:regions', 'The atlas regions are not in the order of the region signals');
    end
    vert = arrayfun(@(s) s.Vertices(:)', Sc(:)', 'UniformOutput', false);
end
end

function [ftElec, ftLeadfield] = keep_channels(ftElec, ftLeadfield, idx)
% electrodes and leadfield rows of the channels used, leadfield average-
% referenced over these channels: the same as in run_pipeline
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
