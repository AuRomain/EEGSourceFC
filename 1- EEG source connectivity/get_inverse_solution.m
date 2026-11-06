function filters = get_inverse_solution(eeg, srate, noiseCov, ftHeadmodel, ftLeadfield, ftElec, VertNormals, inv, opts)
%GET_INVERSE_SOLUTION Inverse filters (sources x channels): eLORETA, LCMV, wMNE.
%   filters = get_inverse_solution(eeg, srate, noiseCov, ftHeadmodel,
%   ftLeadfield, ftElec, VertNormals, inv, opts) returns nInv x nSources x
%   nChannels, one slice per method in INV ({'eloreta', 'lcmv', 'wmne'}),
%   as the original. eeg and ftLeadfield must have the same channels, in
%   the same order (run_pipeline), both average-referenced.
%   opts: eloreta_lambda (0.05), lcmv_lambda ('5%'), wmne (weightExp 0.5,
%   weightLimit 10, SNR 3): the original values; orientation 'normal'
%   (default: one filter per source, along the cortical normal) or 'free'
%   (three filters per source, x/y/z, rows 3i-2..3i; needs a leadfield with
%   three orientations, e.g. EEGLAB_61; get_scouts_timeseries then takes
%   one principal component per region).
%   Same calls as the original (ft_timelockanalysis, ft_sourceanalysis with
%   'eloreta' and 'lcmv'; wmne.ComputeWMNE), corrected:
%   - E2/E3: the original computed free-orientation filters (3 per source),
%     reshaped the 3*nSources x nChan matrix as if its rows were ordered
%     source by source (they are x/y/z by source) and took the RMS, which
%     scrambled sources and removed every sign. Here each source's leadfield
%     is projected on its cortical normal (VertNormals), so every method
%     returns one signed filter per source; the scout sign flip of
%     get_scouts_timeseries is then meaningful. Leadfields stored with the
%     normal orientation only (ftLeadfield.ori = 'normal') are used as they are;
%   - E6: data and leadfield are average-referenced (by the caller);
%   - LCMV and wMNE were computed for every recording but never used; they
%     are now computed only when listed in INV (S.inv). LCMV uses the
%     covariance of this recording (broadband, as the original); eLORETA and
%     wMNE do not depend on the data (noiseCov = identity, as the original).
%   This code was originally developed by Sahar Allouch.
if nargin < 8, inv = {'eloreta'}; end
if nargin < 9, opts = struct(); end
if ~isfield(opts, 'eloreta_lambda'), opts.eloreta_lambda = 0.05; end
if ~isfield(opts, 'lcmv_lambda'), opts.lcmv_lambda = '5%'; end
if ~isfield(opts, 'wmne'), opts.wmne = struct('weightExp', 0.5, 'weightLimit', 10, 'SNR', 3); end
if ~isfield(opts, 'orientation'), opts.orientation = 'normal'; end
free = strcmp(opts.orientation, 'free');
nb_sources = size(VertNormals, 1);
nb_chan = size(eeg, 1);
if numel(ftLeadfield.leadfield) ~= nb_sources || size(ftLeadfield.leadfield{1}, 1) ~= nb_chan
    error('get_inverse_solution:size', 'Leadfield (%d sources x %d channels) does not match the data (%d channels) or the cortex (%d)', ...
        numel(ftLeadfield.leadfield), size(ftLeadfield.leadfield{1}, 1), nb_chan, nb_sources);
end
if free && size(ftLeadfield.leadfield{1}, 2) ~= 3
    error('get_inverse_solution:orientation', ['S.orientation = ''free'' needs a leadfield with three ' ...
        'orientations per source (this one has the normal orientation only); for a cap of your own, ' ...
        'build it with tools/build_leadfield_openmeeg.py --free']);
end
n_ori = 1 + 2 * free;                              % filters per source
filters = zeros(numel(inv), n_ori * nb_sources, nb_chan);

% leadfield projected on the normals (fixed orientation), or kept (free)
lf = ftLeadfield;
if size(ftLeadfield.leadfield{1}, 2) == 3 && ~free
    for i = 1:nb_sources
        lf.leadfield{i} = ftLeadfield.leadfield{i} * VertNormals(i, :)';
    end
end
for f = intersect(fieldnames(lf)', {'mom', 'ori', 'cfg'}), lf = rmfield(lf, f{1}); end
lf.label = ftElec.label(:);

ftData.trial = {eeg};
ftData.time = {(0:size(eeg, 2) - 1) / srate};
ftData.elec = ftElec;
ftData.label = ftElec.label(:);
ftData.fsample = srate;
ftData.sampleinfo = [1, size(eeg, 2)];

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = 'all';
cfg.keeptrials = 'no';
[~, timelock] = evalc('ft_timelockanalysis(cfg, ftData)');

for iv = 1:numel(inv)
    switch inv{iv}
        case 'eloreta'
            cfg = [];
            cfg.method = 'eloreta';
            cfg.sourcemodel = lf;
            if ~isempty(ftHeadmodel), cfg.headmodel = ftHeadmodel; end
            cfg.eloreta.keepfilter = 'yes';
            cfg.eloreta.keepmom = 'no';
            cfg.eloreta.lambda = opts.eloreta_lambda;
            [~, src] = evalc('ft_sourceanalysis(cfg, timelock)');   % FieldTrip prints every eLORETA iteration
            filters(iv, :, :) = reshape(cell2mat(src.avg.filter(:)), [1, n_ori * nb_sources, nb_chan]);
        case 'lcmv'
            cfg = [];
            cfg.method = 'lcmv';
            cfg.sourcemodel = lf;
            if ~isempty(ftHeadmodel), cfg.headmodel = ftHeadmodel; end
            cfg.lcmv.keepfilter = 'yes';
            cfg.lcmv.keepmom = 'no';
            cfg.lcmv.fixedori = 'no';            % orientation fixed by the leadfield (normals), or free
            cfg.lcmv.lambda = opts.lcmv_lambda;
            cfg.lcmv.projectnoise = 'no';
            cfg.keepleadfield = 'no';
            [~, src] = evalc('ft_sourceanalysis(cfg, timelock)');
            filters(iv, :, :) = reshape(cell2mat(src.avg.filter(:)), [1, n_ori * nb_sources, nb_chan]);
        case 'wmne'
            Gain = cell2mat(lf.leadfield);                       % nChan x nSources (fixed) or 3 nSources (free)
            orient = ones(nb_sources, 1);
            if free, orient = []; end                            % every orientation kept
            filters(iv, :, :) = reshape(wmne.ComputeWMNE(noiseCov, Gain, lf.pos, orient, ...
                opts.wmne.weightExp, opts.wmne.weightLimit, opts.wmne.SNR), [1, n_ori * nb_sources, nb_chan]);
        otherwise
            error('get_inverse_solution:method', 'Unknown inverse method %s', inv{iv});
    end
end
end
