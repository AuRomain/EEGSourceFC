function R = toolbox_check(S)
%TOOLBOX_CHECK Are the toolboxes on the path the versions validated with this pipeline?
%   R = toolbox_check(S) compares, for each toolbox of toolbox_manifest
%   found on the path, the content of its checked files with the
%   fingerprints of toolbox_reference.json (the versions validated with
%   the pipeline). R: one row per toolbox with its root folder, its
%   fingerprint, the reference and the status: 'validated', 'equivalent'
%   (a version listed under "equivalent" in toolbox_reference.json: a
%   published version checked to give the same results as the validated
%   one), 'differs' (another version, or edited files: results
%   may differ), 'not found' (not on the path: not needed by these
%   settings) or 'no reference'.
%   setup_toolboxes calls it and, by S.toolbox_check, warns ('warn', the
%   default), stops ('error') or skips the check ('off'); the fingerprints
%   of FieldTrip and Brainstorm enter the cache keys of stage 1, that of
%   modkmeans the clustering's (code_fingerprint).
M = toolbox_manifest();
ref = struct();
f = fullfile(S.root, 'toolbox_reference.json');
if isfile(f), ref = jsondecode(fileread(f)); end
n = numel(M);
root = strings(n, 1);  fp = strings(n, 1);  reference = strings(n, 1);  status = strings(n, 1);
for i = 1:n
    w = which(M(i).marker);
    if isempty(w), status(i) = "not found"; continue; end
    d = fileparts(w);
    for k = 1:M(i).up, d = fileparts(d); end
    root(i) = d;
    fp(i) = toolbox_fingerprint(d, M(i).files);
    if isfield(ref, M(i).name)
        reference(i) = ref.(M(i).name);
        eq = strings(0, 1);
        if isfield(ref, 'equivalent') && isfield(ref.equivalent, M(i).name)
            eq = string({ref.equivalent.(M(i).name).fingerprint});
        end
        if fp(i) == reference(i)
            status(i) = "validated";
        elseif ismember(fp(i), eq)
            status(i) = "equivalent";
        else
            status(i) = "differs";
        end
    else
        status(i) = "no reference";
    end
end
R = table(string({M.name})', root, fp, reference, status, string({M.version})', ...
    'VariableNames', {'toolbox', 'root', 'fingerprint', 'reference', 'status', 'validated_version'});
end

function h = toolbox_fingerprint(d, patterns)
% fingerprint of the checked files of a toolbox (names relative to its root and contents)
files = {};
for p = patterns
    l = dir(fullfile(d, p{1}));
    l = l(~[l.isdir]);
    rel = fileparts(p{1});
    files = [files, fullfile(rel, {l.name})]; %#ok<AGROW>
end
files = sort(strrep(files, '\', '/'));
if isempty(files), h = "none"; return; end
h = string(hash_text([file_fingerprint(fullfile(d, files)) strjoin(files, '|')]));
end
