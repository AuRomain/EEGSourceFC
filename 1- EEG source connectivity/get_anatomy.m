function [A, problem] = get_anatomy(S, recording)
%GET_ANATOMY The head of a recording: the shipped ICBM152 or a folder of S.anatomy.
%   A = get_anatomy(S) the head given by the settings: the pipeline's
%   ICBM152 (inputs/, S.anatomy = '') or the folder S.anatomy, one head for
%   every recording (another template, made by tools/make_anatomy.m).
%   A = get_anatomy(S, recording) with S.anatomy_individual = true: the
%   folder <S.anatomy>/<ID> of the recording's subject (ID = subject_ids,
%   from S.id_regexp; recording = its file name). Without that folder the
%   recording is an error, unless S.anatomy_fallback = true: it then uses
%   the shipped ICBM152 (A.fallback = true; stage1_summary.csv says so).
%   [A, problem] = get_anatomy(...) returns the error message in problem
%   instead of stopping (run_all lists every recording without a head).
%
%   A folder holds Brainstorm files (tools/make_anatomy.m writes them):
%     tess_head.mat, tess_outerskull.mat, tess_innerskull.mat  BEM surfaces
%     tess_cortex.mat           sources (its vertices) and atlases (its
%                               Atlas field: Desikan-Killiany, Destrieux, ...)
%     channel_<montage>.mat     each cap of S.montages placed on this head
%                               (<montage> = file_tag: ASA_10-05_343, ...)
%     subjectimage_T1.mat       the MRI, with its MNI normalization (only
%                               for make_anatomy(S, 'caps', ...))
%   and the leadfields stage 1 computes there (ftLeadfield_<montage>_normal.mat).
%   Coordinates: Brainstorm's SCS of that head, in metres.
%
%   A.name      'ICBM152 (shipped)' or the folder's name
%   A.dir       the folder of the files (inputs/ for the shipped head)
%   A.shipped   true for the pipeline's own ICBM152
%   A.fallback  true when the subject has no folder and uses the shipped head
%   A.head, A.outerskull, A.innerskull, A.cortex  the files
%   A.key       '' for the shipped head (its files are part of the settings
%               key already); else the fingerprint of the folder's surfaces
%               and caps, part of each output's input key: outputs made with
%               another head are recomputed
problem = '';
in = fullfile(fileparts(mfilename('fullpath')), 'inputs');
A = struct('name', 'ICBM152 (shipped)', 'dir', in, 'shipped', true, 'fallback', false, ...
    'head', fullfile(in, 'tess_head.mat'), 'outerskull', fullfile(in, 'tess_outerskull.mat'), ...
    'innerskull', fullfile(in, 'tess_innerskull.mat'), 'cortex', fullfile(in, 'tess_cortex_pial_low.mat'), 'key', '');
if ~isfield(S, 'anatomy') || isempty(S.anatomy), return; end
folder = char(S.anatomy);
if isfield(S, 'anatomy_individual') && S.anatomy_individual
    if nargin < 2 || isempty(recording)
        error('get_anatomy:recording', 'S.anatomy_individual = true: give the recording (get_anatomy(S, file))');
    end
    [~, stem] = fileparts(char(recording));
    id = char(subject_ids(S, {stem}));
    folder = fullfile(folder, id);
    if ~isfolder(folder)
        problem = sprintf('%s: no head model folder %s (subject %s)', stem, folder, id);
        if isfield(S, 'anatomy_fallback') && S.anatomy_fallback
            A.fallback = true;
            A.name = 'ICBM152 (shipped, fallback)';
            problem = '';
            return
        end
        if nargout < 2
            error('get_anatomy:missing', ['%s. Make it with tools/make_anatomy.m, or set S.anatomy_fallback = true ' ...
                'to use the shipped ICBM152 for subjects without one'], problem);
        end
        return
    end
end
[~, name] = fileparts(folder);
A = struct('name', name, 'dir', folder, 'shipped', false, 'fallback', false, ...
    'head', fullfile(folder, 'tess_head.mat'), 'outerskull', fullfile(folder, 'tess_outerskull.mat'), ...
    'innerskull', fullfile(folder, 'tess_innerskull.mat'), 'cortex', fullfile(folder, 'tess_cortex.mat'), 'key', '');
caps = cellfun(@(m) fullfile(folder, ['channel_' file_tag(m) '.mat']), cellstr(S.montages), 'UniformOutput', false);
files = [{A.head, A.outerskull, A.innerskull, A.cortex}, caps(:)'];
missing = files(~cellfun(@isfile, files));
if ~isempty(missing)
    [~, n, e] = cellfun(@fileparts, missing, 'UniformOutput', false);
    problem = sprintf('%s: missing %s (tools/make_anatomy.m makes them; make_anatomy(S, ''caps'', folder) adds a cap)', ...
        folder, strjoin(strcat(n, e), ', '));
    if nargout < 2, error('get_anatomy:files', '%s', problem); end
    return
end
A.key = file_fingerprint(files);
end
