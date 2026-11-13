function folders = make_anatomy(S, kind, source, out, varargin)
%MAKE_ANATOMY Head-model folders for stage 1 (S.anatomy), made by Brainstorm's own processes.
%   make_anatomy(S, 'template', name, out)
%       a Brainstorm template head in the folder out (then S.anatomy = out):
%       a name of Brainstorm's template list ('ICBM152' is part of
%       Brainstorm; 'Colin27_2016', 'ICBM152_2023b', 'FsAverage_2020', ...
%       are downloaded once, as Brainstorm's "Use template" does), or a
%       Brainstorm anatomy .zip or folder (e.g. an exported subject).
%       Brainstorm: db_set_template, and process_generate_bem when the
%       template has no skull surfaces.
%   make_anatomy(S, 't1', subjects, root)
%       individual heads from T1 MRIs: subjects = {ID, T1 file; ...} (one
%       row per subject; ID as stage 1 finds it in the file names,
%       subject_ids), one folder <root>/<ID> each (then S.anatomy = root,
%       S.anatomy_individual = true). Brainstorm: process_import_mri,
%       process_segment_cat12 (CAT12: cortex of 15,000 vertices, its
%       atlases Desikan-Killiany, Destrieux, DKT, ..., non-linear MNI
%       normalization; NAS/LPA/RPA from the MNI positions, Brainstorm's
%       default in scripts), process_generate_bem (scalp, outer and inner
%       skull of 1922 vertices, skull 4 mm). SPM12 and CAT12 are installed as
%       Brainstorm plugins the first time (downloaded). About 15-40 minutes
%       a subject. make_anatomy(..., 'cat12', folder): a CAT12 folder of
%       one's own (the folder of cat_version.m, e.g. when CAT12's server
%       cannot be reached: its releases are also on
%       github.com/ChristianGaser/cat12), given to Brainstorm as its CAT12
%       plugin (bst_plugin SetCustomPath; kept for the next runs).
%   make_anatomy(S, 'segmentation', subjects, root)
%       the same from segmentations already made: subjects = {ID, folder}
%       of CAT12 or FreeSurfer outputs (process_import_anatomy, then
%       process_mni_normalize when it gave no MNI normalization, and
%       process_generate_bem).
%   make_anatomy(S, 'caps', folders)
%       place the caps of S.montages that existing folders do not have yet
%       (a cap added to S.montages, or a folder written by hand); a cap
%       already there is kept (delete its channel_<montage>.mat to place it
%       again: a new file changes the head's key, so its outputs are
%       recomputed).
%   folders = make_anatomy(...) returns the folders made.
%
%   Each folder gets Brainstorm's files of that head (see get_anatomy):
%   tess_head.mat, tess_outerskull.mat, tess_innerskull.mat (BEM surfaces),
%   tess_cortex.mat (the subject's default cortex in Brainstorm, with its
%   atlases and sphere registration), subjectimage_T1.mat (the MRI and its
%   MNI normalization), a channel_<montage>.mat for each cap of S.montages,
%   and anatomy_info.txt (what was done, and the checks). The caps: the
%   pipeline's positions on the ICBM152 (get_ftElec) taken to this head
%   through the MNI normalizations (Brainstorm's cs_convert, as its "Add EEG
%   positions" with MNI coordinates, channel_add_loc) and projected on its
%   scalp (channel_project_scalp, as "Project electrodes on scalp"). Then
%   the leadfield of each cap is computed in the folder (stage 1's
%   get_ftHeadmodel_and_ftLeadfield, OpenMEEG; its checks: nested surfaces,
%   electrodes on the scalp, the sign of the gain), unless
%   make_anatomy(..., 'leadfield', false). A folder already complete is not
%   made again (delete it to redo it); its missing leadfields are computed.
%   Electrode positions digitised on a subject's head are not handled here:
%   a Brainstorm channel file with those positions, in the head's
%   coordinates (Brainstorm: import them on the subject and align them with
%   its scalp), saved in the folder as channel_<montage>.mat with the
%   montage's labels in its order, replaces the projected cap.
%
%   Brainstorm (the copy of S.toolboxes) runs headless in this MATLAB
%   (quit the Brainstorm interface first). Its user folder (plugins,
%   downloaded templates) is <S.toolboxes>/brainstorm_user/.brainstorm, kept
%   for the next runs, never ~/.brainstorm; its database is a protocol
%   deleted at the end.
%
%   Example (one other template for every recording):
%     S = pipeline_settings();
%     make_anatomy(S, 'template', 'Colin27_2016', fullfile(S.root, 'anatomy', 'Colin27_2016'));
%     S.anatomy = fullfile(S.root, 'anatomy', 'Colin27_2016');  run_all(S)
%   Example (individual heads):
%     make_anatomy(S, 't1', {'sub-01', '/mri/sub-01_T1w.nii.gz'; 'sub-02', '/mri/sub-02_T1w.nii.gz'}, ...
%         fullfile(S.root, 'anatomy', 'individual'));
%     S.anatomy = fullfile(S.root, 'anatomy', 'individual');  S.anatomy_individual = true;  run_all(S)
if nargin < 1 || isempty(S), S = pipeline_settings(); end
opts = struct('leadfield', true, 'nvertices', 15000, 'cat12', '');
for k = 1:2:numel(varargin), opts.(lower(varargin{k})) = varargin{k + 1}; end
if nargin < 4, out = ''; end
S.anatomy = '';  S.anatomy_individual = false;  S.anatomy_fallback = false;   % the folders may not exist yet
S = setup_toolboxes(S);
kind = lower(char(kind));
switch kind
    case 'template'
        if isempty(out), error('make_anatomy:out', 'Give the folder of the head: make_anatomy(S, ''template'', name, folder)'); end
        jobs = struct('folder', char(out), 'source', char(source), 'id', 'template');
    case {'t1', 'segmentation'}
        if ~iscell(source) || size(source, 2) ~= 2
            error('make_anatomy:subjects', 'subjects: {ID, T1 file or segmentation folder; ...}, one row per subject');
        end
        if isempty(out), error('make_anatomy:out', 'Give the folder of the subjects'' folders (S.anatomy)'); end
        jobs = struct('folder', fullfile(char(out), cellstr(source(:, 1))), 'source', cellstr(source(:, 2)), ...
            'id', cellstr(source(:, 1)));
    case 'caps'
        jobs = struct('folder', cellstr(source), 'source', '', 'id', 'caps');
    otherwise
        error('make_anatomy:kind', 'kind: ''template'', ''t1'', ''segmentation'' or ''caps''');
end
folders = {jobs.folder};
surfaces = {'tess_head.mat', 'tess_outerskull.mat', 'tess_innerskull.mat', 'tess_cortex.mat', 'subjectimage_T1.mat'};
done = arrayfun(@(j) all(isfile(fullfile(j.folder, [surfaces, {'anatomy_info.txt'}]))), jobs);
if ~strcmp(kind, 'caps') && any(done)
    fprintf('Already made (delete to redo): %s\n', strjoin(folders(done), ', '));
    jobs = jobs(~done);
end
if strcmp(kind, 'caps') && ~all(done)
    error('make_anatomy:folder', 'Not a folder made by make_anatomy: %s', strjoin(folders(~done), ', '));
end
for p = {jobs.source}
    if ~isempty(p{1}) && ~strcmp(kind, 'template') && ~isfile(p{1}) && ~isfolder(p{1})
        error('make_anatomy:source', 'Not found: %s', p{1});
    end
end

if ~isempty(jobs)
    make_folders(S, kind, jobs, opts);
end
if opts.leadfield                                   % the leadfields, with their checks (stage 1's code;
    for j = 1:numel(folders)                        % a leadfield already there is kept)
        S2 = S;  S2.anatomy = folders{j};
        A = get_anatomy(S2);
        for m = 1:numel(S.montages)
            get_ftHeadmodel_and_ftLeadfield(S.montages{m}, 'Desikan_Killiany_68', S2, A);
        end
    end
end
end


function make_folders(S, kind, jobs, opts)
% Brainstorm: each job's head, then its caps
c = start_brainstorm(S); %#ok<NASGU>              % stopped (its folders off the path) when this returns
sMni = load(fullfile(bst_get('BrainstormHomeDir'), 'defaults', 'anatomy', 'ICBM152', 'subjectimage_T1.mat'));
check_gifti(kind, jobs);
if strcmp(kind, 't1')
    install_plugin('spm12');
    if ~isempty(opts.cat12), cat12_folder(char(opts.cat12)); end
    install_plugin('cat12');
end
for j = 1:numel(jobs)
    folder = jobs(j).folder;
    t0 = tic;
    if strcmp(kind, 'caps')
        info = strsplit(fileread(fullfile(folder, 'anatomy_info.txt')), newline);
    else
        fprintf('\n=== %s: %s (%s) ===\n', folder, jobs(j).source, kind);
        info = make_head(kind, jobs(j), folder, opts);
    end
    placed = place_caps(S, folder, sMni, strcmp(kind, 'caps'));
    placed = placed(~cellfun(@isempty, placed));
    for c = placed                                    % a cap placed again: its line replaces the old one
        info = info(~startsWith(info, [extractBefore(c{1}, ':') ':']));
    end
    info = [info, placed]; %#ok<AGROW>
    write_info(folder, info);
    fprintf('%s done in %.1f min\n', folder, toc(t0) / 60);
end
end


function info = make_head(kind, job, folder, opts)
% the head in Brainstorm's database, then its files in the folder
subj = sprintf('S%s', char(datetime('now', 'Format', 'HHmmssSSS')));
[~, iSubject] = db_add_subject(subj, [], 0, 0);
switch kind
    case 'template'
        sTemplate = bst_get('AnatomyDefaults', job.source);
        if isempty(sTemplate)
            if ~isfile(job.source) && ~isfolder(job.source)
                list = bst_get('AnatomyDefaults');
                error('make_anatomy:template', 'Unknown template %s. Brainstorm''s: %s', job.source, strjoin({list.Name}, ', '));
            end
            [~, name] = fileparts(job.source);
            sTemplate = struct('Name', name, 'FilePath', job.source);
        end
        db_set_template(iSubject, sTemplate, 0);
        what = sprintf('Brainstorm template %s (%s)', sTemplate.Name, sTemplate.FilePath);
    case 't1'
        call_process('process_import_mri', 'subjectname', subj, 'mrifile', {job.source, 'ALL'});
        call_process('process_segment_cat12', 'subjectname', subj, 'nvertices', opts.nvertices, ...
            'tpmnii', {'', 'Nifti1'}, 'sphreg', 1, 'vol', 0, 'extramaps', 0, 'cerebellum', 0);
        what = sprintf('T1 %s, CAT12 segmentation (Brainstorm process_segment_cat12, %d vertices)', job.source, opts.nvertices);
    case 'segmentation'
        format = 'FreeSurfer-fast';
        if ~isempty(dir(fullfile(job.source, '**', 'cat_*.xml'))) || ~isempty(dir(fullfile(job.source, 'surf', '*.central.*.gii')))
            format = 'CAT12';
        end
        call_process('process_import_anatomy', 'subjectname', subj, ...
            'mrifile', {job.source, format}, 'nvertices', opts.nvertices);
        what = sprintf('%s segmentation %s (Brainstorm process_import_anatomy, %d vertices)', format, job.source, opts.nvertices);
end
sSubject = bst_get('Subject', iSubject);
if isempty(sSubject.Anatomy) || isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
    error('make_anatomy:brainstorm', '%s: Brainstorm gave no MRI, cortex or scalp (see its messages above)', job.source);
end
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
if ~isfield(sMri, 'NCS') || isempty(cs_convert(sMri, 'mni', 'scs', [0 0 0]))      % no MNI normalization yet
    call_process('process_mni_normalize', 'subjectname', subj, 'method', 'maff8');
    sSubject = bst_get('Subject', iSubject);
    sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
end
comments = {sSubject.Surface.Comment};
if ~any(startsWith(comments, 'bem_')) && (isempty(sSubject.iInnerSkull) || isempty(sSubject.iOuterSkull))
    call_process('process_generate_bem', 'subjectname', subj, ...
        'nscalp', 1922, 'nouter', 1922, 'ninner', 1922, 'thickness', 4, 'method', 'brainstorm');
    sSubject = bst_get('Subject', iSubject);
end
if ~isfolder(folder), mkdir(folder); end
files = {'tess_head.mat', 'tess_outerskull.mat', 'tess_innerskull.mat', 'tess_cortex.mat'};
kinds = {'iScalp', 'bem_head'; 'iOuterSkull', 'bem_outerskull'; 'iInnerSkull', 'bem_innerskull'; 'iCortex', ''};
nv = zeros(1, 4);
names = cell(1, 4);
for k = 1:4
    i = [];
    if ~isempty(kinds{k, 2}), i = find(startsWith({sSubject.Surface.Comment}, kinds{k, 2}), 1, 'last'); end
    if isempty(i), i = sSubject.(kinds{k, 1}); end
    if isempty(i), error('make_anatomy:brainstorm', '%s: no %s surface', job.source, files{k}); end
    sSurf = in_tess_bst(sSubject.Surface(i).FileName);
    sSurf.VertNormals = interface_normals(sSurf);
    if isfield(sSurf, 'NCS') && isfield(sSurf.NCS, 'iy')     % Brainstorm copies the MRI's deformation fields into the
        sSurf.NCS.y = [];  sSurf.NCS.iy = [];                % BEM head (75 MB with CAT12); subjectimage_T1.mat keeps them
    end
    save(fullfile(folder, files{k}), '-struct', 'sSurf');
    nv(k) = size(sSurf.Vertices, 1);
    names{k} = sSurf.Comment;
end
save(fullfile(folder, 'subjectimage_T1.mat'), '-struct', 'sMri');
mni = 'linear';
if isfield(sMri.NCS, 'iy') && ~isempty(sMri.NCS.iy), mni = 'non-linear'; end
info = {sprintf('Head model made by tools/make_anatomy.m, %s', char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm'))), ...
    sprintf('Source: %s', what), ...
    sprintf('Brainstorm %s', bst_version()), ...
    sprintf('MNI normalization of the MRI: %s', mni), ...
    sprintf('Surfaces: head %d V (%s), outer skull %d V (%s), inner skull %d V (%s), cortex %d V (%s)', ...
        nv(1), names{1}, nv(2), names{2}, nv(3), names{3}, nv(4), names{4})};
C = load(fullfile(folder, 'tess_cortex.mat'), 'Atlas');
atlases = {C.Atlas(~cellfun(@isempty, {C.Atlas.Scouts})).Name};
info{end + 1} = sprintf('Atlases of the cortex (S.atlases by name): %s', strjoin(atlases, ', '));
A = struct('name', job.id, 'dir', folder, 'shipped', false, 'fallback', false, 'cortex', fullfile(folder, 'tess_cortex.mat'));
try
    Sc = get_atlas('Desikan_Killiany_68', A);
    info{end + 1} = sprintf('Desikan_Killiany_68: %d regions, the labels of the shipped ICBM152', numel(Sc));
catch err
    info{end + 1} = sprintf('Desikan_Killiany_68 NOT AVAILABLE on this cortex: %s', err.message);
    warning('make_anatomy:atlas', '%s', info{end});
end
end


function check_gifti(kind, jobs)
% CAT12's surfaces are GIfTI files: Brainstorm before 13 October 2026 reads them with Java 8's
% sun.misc.BASE64Decoder, which Java 9 and later (MATLAB's since R2023b on Apple silicon) do not have
cat12 = strcmp(kind, 't1') || (strcmp(kind, 'segmentation') && ...
    any(arrayfun(@(J) ~isempty(dir(fullfile(J.source, '**', '*.gii'))), jobs)));
if ~cat12 || exist('bst_base64', 'file') == 2
    return;
end
try
    javaObject('sun.misc.BASE64Decoder');
catch
    error('make_anatomy:gifti', ['This Brainstorm cannot read CAT12''s surfaces with MATLAB''s Java (%s): ' ...
        'it needs Brainstorm''s fix of 13 October 2026 (commit 952475b, toolbox/misc/bst_base64.m), which ' ...
        'tools/fetch_toolboxes.sh applies to the pinned Brainstorm, or a Brainstorm of that date or later'], ...
        char(java.lang.System.getProperty('java.version')));
end
end

function install_plugin(name)
[ok, msg] = bst_plugin('Install', name, 0);
if ~ok, error('make_anatomy:plugin', 'Brainstorm could not install %s: %s', name, msg); end
end


function cat12_folder(folder)
% a CAT12 folder of one's own as Brainstorm's CAT12 plugin (its custom plugin
% path), linked into SPM12's toolbox folder as Brainstorm links the CAT12 it
% installs (bst_plugin LinkCatSpm, which skips plugins of a custom path)
if ~isfile(fullfile(folder, 'cat_version.m'))
    error('make_anatomy:cat12', '%s: no cat_version.m (give the CAT12 folder)', folder);
end
plugins = bst_get('UserPluginsDir');
if startsWith([folder filesep], [plugins filesep])
    error('make_anatomy:cat12', ['%s is inside Brainstorm''s plugins folder, where Brainstorm takes a folder ' ...
        'whose path starts with %s for the CAT12 it installs: put it elsewhere (e.g. next to the toolboxes)'], ...
        folder, fullfile(plugins, 'cat12'));
end
spm = bst_plugin('GetInstalled', 'spm12');
link = fullfile(spm.Path, spm.SubFolder, 'toolbox', 'cat12');
target = '';
if ~ispc, [~, target] = system(sprintf('readlink "%s"', link)); target = strtrim(target); end
if ~isfile(fullfile(link, 'cat12.m')) || (~isempty(target) && ~strcmp(target, folder))   % none, or a link elsewhere
    if ispc, cmd = sprintf('mklink /D "%s" "%s"', link, folder); else, cmd = sprintf('ln -sfn "%s" "%s"', folder, link); end
    [st, out] = system(cmd);
    if st ~= 0, error('make_anatomy:cat12', 'Could not link %s into SPM12: %s', folder, out); end
end
bst_plugin('SetCustomPath', 'cat12', folder);
end


function call_process(name, varargin)
% a Brainstorm process, as its scripts call it (bst_process CallProcess); the
% errors and warnings it leaves in Brainstorm's report are printed, an error stops
global GlobalData
R = GlobalData.ProcessReports.Reports;
n0 = 0;
if iscell(R) && size(R, 2) == 5, n0 = size(R, 1); end
bst_process('CallProcess', name, [], [], varargin{:});
R = GlobalData.ProcessReports.Reports;
if ~iscell(R) || size(R, 2) ~= 5 || size(R, 1) <= n0, return; end
R = R(n0 + 1:end, :);
for i = find(ismember(lower(R(:, 1)), {'error', 'warning'}))'
    fprintf('Brainstorm %s (%s): %s\n', R{i, 1}, name, strtrim(char(R{i, 4})));
end
err = strcmpi(R(:, 1), 'error');
if any(err)
    error('make_anatomy:brainstorm', '%s: %s', name, strjoin(cellfun(@(m) strtrim(char(m)), R(err, 4), 'UniformOutput', false), ' | '));
end
end


function N = interface_normals(sSurf)
% the vertex normals as Brainstorm's interface computes them (tess_normals:
% MATLAB's patch normals). Headless (server mode) it sums the faces' normals
% instead, which on the ICBM152 cortex differs from the patch normals by a
% median 5 degrees (up to 90, 23 vertices reversed); the patch normals of
% this headless MATLAB are those of the pipeline's cortex (0 degrees)
global GlobalData
gui = GlobalData.Program.GuiLevel;
GlobalData.Program.GuiLevel = 0;
c = onCleanup(@() restore_gui(gui));
N = tess_normals(sSurf.Vertices, sSurf.Faces, sSurf.VertConn);
end


function restore_gui(gui)
global GlobalData
GlobalData.Program.GuiLevel = gui;
end


function info = place_caps(S, folder, sMni, keep)
% the pipeline's caps (ICBM152 positions) on this head: MNI, then the scalp;
% keep: leave the caps already there
head = load(fullfile(folder, 'tess_head.mat'), 'Vertices', 'Faces');
sMri = load(fullfile(folder, 'subjectimage_T1.mat'));
% the linear part of the head's MNI normalization: a non-linear one (CAT12's
% deformation) is fitted to the brain, and outside it, at the scalp, moved
% electrodes near the vertex by up to 3 cm (Colin27's T1)
linear = 'linear';
if isfield(sMri.NCS, 'y') && ~isempty(sMri.NCS.y)
    if isfield(sMri.NCS, 'R') && ~isempty(sMri.NCS.R) && isfield(sMri.NCS, 'T') && ~isempty(sMri.NCS.T)
        sMri.NCS.y = [];  sMri.NCS.iy = [];
        linear = 'linear part of the';
    else
        linear = 'non-linear (no linear part)';
    end
end
info = {scalp_check(head, folder)};
for m = 1:numel(S.montages)
    montage = S.montages{m};
    if keep && isfile(fullfile(folder, ['channel_' file_tag(montage) '.mat']))
        fprintf('%s: cap %s already there, kept\n', folder, montage);
        continue
    end
    elec = get_ftElec(montage);
    P = cs_convert(sMni, 'scs', 'mni', elec.chanpos);
    P = cs_convert(sMri, 'mni', 'scs', P);
    if isempty(P) || any(~isfinite(P(:)))
        error('make_anatomy:mni', '%s: the MNI normalization of %s does not cover cap %s', folder, ...
            'subjectimage_T1.mat', montage);
    end
    P0 = P;
    P = channel_project_scalp(head.Vertices, P);
    d0 = sqrt(sum((P - P0) .^ 2, 2));               % how far the projection moved each electrode
    ChannelMat = db_template('channelmat');
    ChannelMat.Comment = sprintf('%s (%d) on %s: MNI positions, projected on the scalp', montage, numel(elec.label), folder);
    Channel = repmat(db_template('channeldesc'), 1, numel(elec.label));
    for i = 1:numel(elec.label)
        Channel(i).Name = elec.label{i};
        Channel(i).Type = 'EEG';
        Channel(i).Loc = P(i, :)';
        Channel(i).Weight = 1;
    end
    ChannelMat.Channel = Channel;
    save(fullfile(folder, ['channel_' file_tag(montage) '.mat']), '-struct', 'ChannelMat');
    info{end + 1} = sprintf(['Cap %s: %d electrodes through the MNI normalizations (%s MNI normalization of this ' ...
        'head), then projected on the scalp: moved by %.1f mm (median; max %.1f)'], montage, numel(elec.label), ...
        linear, 1e3 * median(d0), 1e3 * max(d0)); %#ok<AGROW>
    [dd, o] = sort(d0, 'descend');
    far = o(dd > 0.010);                              % where the head differs from the template's: worth a look
    if ~isempty(far)
        names = arrayfun(@(i) sprintf('%s %.0f', elec.label{i}, 1e3 * d0(i)), far(1:min(end, 12)), 'UniformOutput', false);
        info{end} = sprintf('%s; more than 10 mm: %d (%s%s mm)', info{end}, numel(far), strjoin(names, ', '), ...
            repmat(', ...', 1, numel(far) > 12));
    end
    fprintf('%s\n', info{end});
end
end


function line = scalp_check(head, folder)
% vertices of the scalp standing out from their neighbours (distance from the
% head's centre) by more than 8 mm: a spike or a fold of the surface, which
% moves the electrodes projected near it (e.g. CAT12's scalp of Colin27's T1:
% two bumps of up to 39 mm on top of the head)
V = head.Vertices;  F = double(head.Faces);
r = sqrt(sum((V - mean(V, 1)) .^ 2, 2));
nb = sparse(F(:, [1 2 3 2 3 1]), F(:, [2 3 1 1 2 3]), 1, size(V, 1), size(V, 1)) > 0;
excess = zeros(size(r));
for i = 1:numel(r)
    excess(i) = r(i) - median(r(nb(:, i)));
end
bad = find(abs(excess) > 0.008);
line = '';
if isempty(bad), return; end
line = sprintf(['Scalp: %d vertices stand out from their neighbours by more than 8 mm (up to %.0f mm): ' ...
    'check the head surface in Brainstorm'], numel(bad), 1e3 * max(abs(excess)));
warning('make_anatomy:scalp', '%s: %s', folder, line);
end


function write_info(folder, info)
fid = fopen(fullfile(folder, 'anatomy_info.txt'), 'w');
fprintf(fid, '%s\n', info{~cellfun(@isempty, info)});
fclose(fid);
end


function v = bst_version()
B = bst_get('Version');
v = sprintf('%s (%s)', B.Version, B.Date);
end


function c = start_brainstorm(S)
% Brainstorm headless; user folder <S.toolboxes>/brainstorm_user (kept),
% database: a protocol deleted at the end
bst = fileparts(fileparts(fileparts(which('bst_bandpass_hfilter'))));   % brainstorm3/toolbox/math
if isempty(bst), error('make_anatomy:brainstorm', 'Brainstorm is not on the path (setup_toolboxes)'); end
addpath(bst);
if brainstorm('status'), error('make_anatomy:running', 'Brainstorm is running: quit it first (brainstorm stop)'); end
user = fullfile(S.toolboxes, 'brainstorm_user');
if ~isfolder(user), mkdir(user); end
home = char(java.lang.System.getProperty('user.home'));
p0 = path;
java.lang.System.setProperty('user.home', user);    % Brainstorm's user folder: <user>/.brainstorm
c = onCleanup(@() stop_brainstorm(home, p0));
brainstorm('setpath');                              % its whole tree, until it stops
bst_startup(bst, -1, 'local', '');                  % server mode, database in <user>/.brainstorm/local_db
gui_brainstorm('CreateProtocol', sprintf('make_anatomy_%s', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))), 0, 0);
end


function stop_brainstorm(home, p0)
try, db_delete_protocol(0, 1); catch, end           % the protocol and its files
try, brainstorm('stop'); catch, end
java.lang.System.setProperty('user.home', home);
path(p0);                                           % only the two folders of setup_toolboxes again
end
