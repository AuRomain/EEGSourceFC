function [ftHeadmodel, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(montage, atlas, S, A)
%GET_FTHEADMODEL_AND_FTLEADFIELD Head model and leadfield (FieldTrip format), cached in inputs/.
%   [ftHeadmodel, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(montage, atlas, S)
%   Same logic as the original: load the cached leadfield of the montage, or
%   compute it (BEM of the ICBM152 head: scalp, outer and inner skull,
%   conductivities 1 / 0.0125 / 1; sources = the 15,002 vertices of
%   tess_cortex_pial_low, on which every atlas of inputs/ is defined).
%   Montages (get_ftElec):
%     'EEGLAB_61'      ftLeadfield_EEGLab_61_corrected.mat: the corrected
%                      electrode positions (E1), average-referenced (E6)
%     'ASA_10-05_343'  ftLeadfield_ASA_10-05_343_normal.mat: 343 standard
%                      positions, normal orientation only (single, 20 MB)
%     a channel file   ftLeadfield_<file name>_normal.mat, computed the
%                      first time (see below), and again when the file's
%                      labels or positions change
%   ftLeadfield.leadfield{i} is nChannels x 3 (x, y, z) or, for the normal
%   files, nChannels x 1 (dipole along the cortical normal, ftLeadfield.ori
%   = 'normal'); get_inverse_solution handles both. With S.orientation =
%   'free', a cap given as a file gets its three-orientation leadfield
%   (G_free, computed with it; an older file without one is recomputed);
%   the ASA_10-05_343 file has the normal orientation only. The average reference
%   is applied in run_pipeline to the channels actually used.
%   Computing a leadfield needs OpenMEEG: FieldTrip's 'openmeeg' method
%   when the OpenMEEG programs are installed (the original route), or the
%   Python OpenMEEG package through tools/build_leadfield_openmeeg.py
%   (S.openmeeg_python = the Python executable). Files are found from this
%   folder, whatever the current directory (N15). The head model (177 MB)
%   is only needed by FieldTrip's route.
%   [ftHeadmodel, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(montage, atlas, S, A)
%   on another head A (get_anatomy; S.anatomy): the same BEM (conductivities,
%   OpenMEEG) on that head's surfaces and cortex, with its cap
%   channel_<montage>.mat, cached in its folder as ftLeadfield_<montage>_normal.mat
%   (computed the first time, and again when the surfaces or the cap change:
%   their fingerprint is kept in the file, anatomy_key).
if nargin < 3, S = struct(); end
if nargin < 4 || isempty(A), A = get_anatomy(struct()); end
free = isfield(S, 'orientation') && strcmp(S.orientation, 'free');
in = fullfile(fileparts(mfilename('fullpath')), 'inputs');
if ~strcmp(atlas, 'Desikan_Killiany_68') && ~isfile(atlas)
    get_atlas(atlas, A);                            % an atlas of the cortex, by name (else an error)
end
if ~A.shipped
    [ftHeadmodel, ftLeadfield] = other_head(montage, A, S, free);
    return
end
cortex_file = fullfile(in, 'tess_cortex_pial_low.mat');
switch montage
    case 'EEGLAB_61'
        file_channel = fullfile(in, 'channel_EEGLab_61_corrected_bs.mat');
        file_leadfield = fullfile(in, 'ftLeadfield_EEGLab_61_corrected.mat');
    case 'ASA_10-05_343'
        file_channel = fullfile(in, 'channel_ASA_10-05_343_bs.mat');
        file_leadfield = fullfile(in, 'ftLeadfield_ASA_10-05_343_normal.mat');
    otherwise
        if ~isfile(montage)
            error('get_ftHeadmodel_and_ftLeadfield:montage', 'No inputs defined for montage %s', montage);
        end
        file_channel = montage;
        [~, name] = fileparts(montage);
        file_leadfield = fullfile(in, ['ftLeadfield_' matlab.lang.makeValidName(name) '_normal.mat']);
end
file_headmodel = fullfile(in, 'ftHeadmodel_EEGLab_61.mat');     % same BEM for every montage
ftHeadmodel = [];                     % not needed: ft_sourceanalysis uses the precomputed leadfield

if isfile(file_leadfield) && ~any(strcmp(montage, {'EEGLAB_61', 'ASA_10-05_343'}))
    % a cap given as a file: recompute when its labels or positions changed
    % since the leadfield was made (the file name alone would not tell)
    w = whos('-file', file_leadfield);
    elec = BS_to_ft_channels(file_channel, []);
    stale = true;
    if ismember('elec_pos', {w.name})
        L = load(file_leadfield, 'labels', 'elec_pos');
        stale = ~isequal(lower(cellstr(L.labels(:))), lower(elec.label(:))) || ...
            ~isequal(size(L.elec_pos), size(elec.chanpos)) || max(abs(L.elec_pos(:) - elec.chanpos(:))) > 1e-6;
    end
    if free && ~stale && ~ismember('G_free', {w.name})
        fprintf('%s: no three-orientation leadfield (S.orientation = ''free''): recomputed\n', file_leadfield);
        stale = true;
    end
    if stale
        fprintf('%s: electrode labels or positions changed, leadfield recomputed\n', file_channel);
        delete(file_leadfield);
    end
end
if ~isfile(file_leadfield)
    compute_leadfield(file_channel, file_leadfield, file_headmodel, A, S);
end
ftLeadfield = load_leadfield(file_leadfield, free);
elec = get_ftElec(montage);
if ~isequal(lower(elec.label(:)), lower(ftLeadfield.label(:)))
    error('get_ftHeadmodel_and_ftLeadfield:labels', '%s: the leadfield channels differ from the montage''s', file_leadfield);
end
end


function [ftHeadmodel, ftLeadfield] = other_head(montage, A, S, free)
% another head (S.anatomy): its cap and leadfield, in its folder
ftHeadmodel = [];
file_channel = fullfile(A.dir, ['channel_' file_tag(montage) '.mat']);
file_leadfield = fullfile(A.dir, ['ftLeadfield_' file_tag(montage) '_normal.mat']);
if ~isfile(file_channel)
    error('get_ftHeadmodel_and_ftLeadfield:anatomy', 'No cap %s on the head %s: make_anatomy(S, ''caps'', ''%s'') places it', ...
        montage, A.dir, A.dir);
end
elec = get_ftElec(montage);                         % labels and order of the montage (run_pipeline's channels)
on_head = BS_to_ft_channels(file_channel, []);
if ~isequal(lower(elec.label(:)), lower(on_head.label(:)))
    error('get_ftHeadmodel_and_ftLeadfield:labels', '%s: its channels differ from montage %s''s (labels or order)', ...
        file_channel, montage);
end
anatomy_key = file_fingerprint({A.head, A.outerskull, A.innerskull, A.cortex, file_channel});
if isfile(file_leadfield)
    w = whos('-file', file_leadfield);
    stale = ~ismember('anatomy_key', {w.name});
    if ~stale
        L = load(file_leadfield, 'anatomy_key');
        stale = ~strcmp(L.anatomy_key, anatomy_key);
    end
    if stale
        fprintf('%s: the head''s surfaces or cap changed, leadfield recomputed\n', file_leadfield);
    elseif free && ~ismember('G_free', {w.name})
        fprintf('%s: no three-orientation leadfield (S.orientation = ''free''): recomputed\n', file_leadfield);
        stale = true;
    end
    if stale, delete(file_leadfield); end
end
if ~isfile(file_leadfield)
    compute_leadfield(file_channel, file_leadfield, fullfile(A.dir, 'ftHeadmodel.mat'), A, S);
    save(file_leadfield, 'anatomy_key', '-append');
end
ftLeadfield = load_leadfield(file_leadfield, free);
if ~isequal(lower(elec.label(:)), lower(ftLeadfield.label(:)))
    error('get_ftHeadmodel_and_ftLeadfield:labels', '%s: the leadfield channels differ from the montage''s', file_leadfield);
end
end


function ftLeadfield = load_leadfield(file_leadfield, free)
w = whos('-file', file_leadfield);
if ismember('ftLeadfield', {w.name})
    load(file_leadfield, 'ftLeadfield');
elseif free && ismember('G_free', {w.name})              % compact file, three orientations
    L = load(file_leadfield, 'G_free', 'labels', 'src_pos');
    G = double(L.G_free);                                % x, y, z of each source
    ftLeadfield = compact_leadfield(L);
    ftLeadfield.leadfield = mat2cell(G, size(G, 1), 3 * ones(1, size(G, 2) / 3));
else                                                     % compact normal-orientation file
    L = load(file_leadfield, 'G_normal', 'labels', 'src_pos');
    G = double(L.G_normal);
    ftLeadfield = compact_leadfield(L);
    ftLeadfield.leadfield = num2cell(G, 1);
    ftLeadfield.ori = 'normal';
end
end


function lf = compact_leadfield(L)
% FieldTrip fields of a compact leadfield file (G_normal or G_free)
lf.pos = L.src_pos;
lf.inside = true(size(L.src_pos, 1), 1);
lf.unit = 'm';
lf.label = cellstr(L.labels(:));
lf.leadfielddimord = '{pos}_chan_ori';
end


function compute_leadfield(file_channel, file_leadfield, file_headmodel, A, S)
% FieldTrip + OpenMEEG programs (the original route), else Python OpenMEEG;
% A: the head (get_anatomy), its surfaces and cortex
cortex_file = A.cortex;
if ispc, [st, ~] = system('where om_assemble'); else, [st, ~] = system('which om_assemble'); end
if st == 0
    if isfile(file_headmodel)
        load(file_headmodel, 'ftHeadmodel');
    else
        ftGeometry = BS_to_ft_tess({A.head; A.outerskull; A.innerskull});
        cfg = [];
        cfg.method = 'openmeeg';
        cfg.conductivity = [1, 0.0125, 1];          % brain, skull, scalp (Brainstorm default)
        ftHeadmodel = ft_prepare_headmodel(cfg, ftGeometry);
        save(file_headmodel, 'ftHeadmodel');
    end
    ftElec = BS_to_ft_channels(file_channel, []);
    C = load(cortex_file, 'Vertices', 'VertNormals');
    cfg = [];
    cfg.elec = ftElec;
    cfg.sourcemodel = struct('pos', C.Vertices, 'inside', true(size(C.Vertices, 1), 1), 'unit', 'm');
    cfg.headmodel = ftHeadmodel;
    cfg.reducerank = 'no';
    lf = ft_prepare_leadfield(cfg);
    n = C.VertNormals ./ vecnorm(C.VertNormals, 2, 2);
    G_normal = zeros(numel(ftElec.label), numel(lf.leadfield), 'single');
    for i = 1:numel(lf.leadfield)
        G_normal(:, i) = lf.leadfield{i} * n(i, :)';
    end
    labels = ftElec.label; src_pos = C.Vertices; elec_pos = ftElec.chanpos;
    note = sprintf('FieldTrip ft_prepare_leadfield, OpenMEEG BEM, %s, conductivities 1/0.0125/1, normal orientation', A.name);
    save(file_leadfield, 'G_normal', 'labels', 'src_pos', 'elec_pos', 'note');
    if isfield(S, 'orientation') && strcmp(S.orientation, 'free')
        G_free = single(cell2mat(lf.leadfield));    % x, y, z of each source
        save(file_leadfield, 'G_free', '-append');
    end
    return
end
py = 'python3';
if ispc, py = 'python'; end                         % the name of the Windows installers
if isfield(S, 'openmeeg_python') && ~isempty(S.openmeeg_python), py = S.openmeeg_python; end
tool = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'tools', 'build_leadfield_openmeeg.py');
cmd = sprintf('"%s" "%s" --channel-file "%s" --out "%s"', py, tool, file_channel, file_leadfield);
if ~A.shipped, cmd = sprintf('%s --anatomy "%s" --cortex "%s"', cmd, A.dir, A.cortex); end
if isfield(S, 'orientation') && strcmp(S.orientation, 'free'), cmd = [cmd ' --free']; end
fprintf('Computing the leadfield of %s with OpenMEEG (a few minutes):\n  %s\n', file_channel, cmd);
if ispc, cmd = ['"' cmd '"']; end                  % cmd.exe /c removes the first and the last quote of a
[st, out] = system(cmd);                           % line that starts with one: these, not the path's
disp(out);
if st ~= 0 || ~isfile(file_leadfield)
    stopped = regexp(out, '(RuntimeError|AssertionError): [^\n]*', 'match', 'once');
    if ~isempty(stopped)                           % the builder ran and stopped at one of its checks
        error('get_ftHeadmodel_and_ftLeadfield:openmeeg', 'No leadfield for %s: %s (the builder''s output above)', ...
            file_channel, stopped);
    end
    error('get_ftHeadmodel_and_ftLeadfield:openmeeg', ['No leadfield for %s, and OpenMEEG is not available.\n' ...
        'Install the OpenMEEG programs (FieldTrip route) or the Python package (pip install openmeeg scipy)\n' ...
        'and set S.openmeeg_python to that Python (README, section 1: What you need).'], file_channel);
end
end
