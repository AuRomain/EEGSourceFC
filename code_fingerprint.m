function [h, code, toolboxes] = code_fingerprint(S, part)
%CODE_FINGERPRINT Fingerprint of the code behind an output, for cache keys.
%   h = code_fingerprint(S, part) hashes the .m files that shape the output:
%   'stage1'      "1- EEG source connectivity" (with its packages) and the
%                 shared helpers it calls (band_edges, file_tag, subject_ids)
%   'clustering'  "2- Dynamic connectivity/fun/Clustering" (the clustering, the
%                 stability of the number of states and window_matrix, which
%                 scales the windows) and fun/ICA/load_and_concatenate.m, which
%                 reads them (the cached part of stage 2)
%   plus the toolboxes they use (S.versions from setup_toolboxes: the
%   content fingerprints of toolbox_check, or the version strings when the
%   check is off): FieldTrip and Brainstorm for stage 1, modkmeans for the
%   clustering. A change of this code or of these toolboxes makes
%   cached outputs stale, so they are recomputed (V3-6). The settings files,
%   the tests and the documentation are not part of it.
%   [h, code, toolboxes] = code_fingerprint(S, part) also gives the two
%   parts of h: the code's own fingerprint (its files, by name and content)
%   and the toolbox versions, h = hash_text([code toolboxes]); settings_key
%   puts an earlier code in place of the current one to accept outputs of
%   earlier code that computes the same (equivalent_code.txt).
switch part
    case 'stage1'
        folders = {fullfile(S.root, '1- EEG source connectivity')};
        extra = fullfile(S.root, {'band_edges.m', 'file_tag.m', 'subject_ids.m'});
    case 'clustering'
        folders = {fullfile(S.root, '2- Dynamic connectivity', 'fun', 'Clustering')};
        extra = {fullfile(S.root, '2- Dynamic connectivity', 'fun', 'ICA', 'load_and_concatenate.m')};
    otherwise
        error('code_fingerprint:part', 'Unknown part %s', part);
end
files = extra;
for k = 1:numel(folders)
    d = dir(fullfile(folders{k}, '**', '*.m'));
    files = [files, fullfile({d.folder}, {d.name})]; %#ok<AGROW>
end
files = sort(strrep(files, [S.root filesep], ''));     % relative names: the same in any folder
code = [file_fingerprint(fullfile(S.root, files)) hash_text(strjoin(files, '|'))];
v = '';
tb = {'modkmeans'};
if strcmp(part, 'stage1'), tb = {'fieldtrip', 'brainstorm'}; end
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
