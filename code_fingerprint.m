function [h, code, toolboxes] = code_fingerprint(S, part)
%CODE_FINGERPRINT Fingerprint of the code behind an output, for cache keys.
%   h = code_fingerprint(S, part) hashes the .m files that shape the output:
%   'stage1'      "1- EEG source connectivity" (with its packages) and the
%                 shared helpers it calls (band_edges, file_tag, subject_ids)
%   'statistics'  "Analysis" and the helpers at the top of the pipeline
%                 folder (not the settings and demo scripts): the code behind
%                 the statistics, for their provenance (run_statistics;
%                 no toolbox part)
%   'clustering'  "Dynamic connectivity/fun/Clustering" (the clustering, the
%                 stability of the number of states and window_matrix, which
%                 scales the windows) and fun/ICA/load_and_concatenate.m, which
%                 reads them (the cached part of the states)
%   (the stage folders found by name, stage_folder: "1- EEG source
%   connectivity", "3- Dynamic connectivity", "4- Analysis")
%   plus the toolboxes they use (S.versions from setup_toolboxes: the
%   content fingerprints of toolbox_check, or the version strings when the
%   check is off): FieldTrip and Brainstorm for stage 1, modkmeans for the
%   clustering. A change of this code or of these toolboxes makes
%   cached outputs stale, so they are recomputed. The settings files,
%   the tests and the documentation are not part of it.
%   [h, code, toolboxes] = code_fingerprint(S, part) also gives the two
%   parts of h: the code's own fingerprint (its files, by name and content)
%   and the toolbox versions, h = hash_text([code toolboxes]); settings_key
%   puts a listed code fingerprint in place of the current one to accept
%   outputs of code that computes the same (equivalent_code.txt).
switch part
    case 'stage1'
        folders = {fullfile(S.root, '1- EEG source connectivity')};
        extra = fullfile(S.root, {'band_edges.m', 'file_tag.m', 'subject_ids.m'});
    case 'clustering'
        dyn = stage_folder(S.root, 'Dynamic connectivity');
        folders = {fullfile(dyn, 'fun', 'Clustering')};
        extra = {fullfile(dyn, 'fun', 'ICA', 'load_and_concatenate.m')};
    case 'statistics'
        folders = {stage_folder(S.root, 'Analysis')};
        d = dir(fullfile(S.root, '*.m'));
        top = setdiff({d.name}, {'pipeline_settings.m', 'run_demo.m', 'run_static_demo.m'});
        extra = fullfile(S.root, top);
    otherwise
        error('code_fingerprint:part', 'Unknown part %s', part);
end
files = extra;
for k = 1:numel(folders)
    d = dir(fullfile(folders{k}, '**', '*.m'));
    files = [files, fullfile({d.folder}, {d.name})]; %#ok<AGROW>
end
files = sort(strrep(strrep(files, [S.root filesep], ''), '\', '/'));   % relative names, with '/': the
                                                       % same in any folder and on Windows
code = [file_fingerprint(fullfile(S.root, files)) hash_text(strjoin(files, '|'))];
v = '';
tb = {'modkmeans'};
if strcmp(part, 'stage1'), tb = {'fieldtrip', 'brainstorm'}; end
if strcmp(part, 'statistics'), tb = {}; end
if isfield(S, 'versions')
    for f = tb
        if isfield(S.versions, 'fingerprints') && isfield(S.versions.fingerprints, f{1})
            v = [v f{1} '=' S.versions.fingerprints.(f{1}) ';']; %#ok<AGROW>
        elseif isfield(S.versions, f{1})
            v = [v f{1} '=' S.versions.(f{1}) ';']; %#ok<AGROW>
        end
    end
end
toolboxes = v;
h = hash_text([code v]);
end
