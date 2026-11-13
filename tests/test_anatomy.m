function tests = test_anatomy
%TEST_ANATOMY Other heads in stage 1 (S.anatomy, get_anatomy, get_atlas, tools/make_anatomy.m).
%   From the pipeline folder:  results = runtests('tests/test_anatomy.m')
%   Two of the simulated recordings of the static tests (static_test_recordings),
%   static stage 1. A head folder holding the shipped ICBM152's files
%   (surfaces, cap, and its leadfield) must give exactly the shipped
%   results; individual heads, the fallback, atlases by name, the
%   settings and the stale leadfield are checked on it. make_anatomy with
%   Brainstorm's own ICBM152 template (no download; skipped without
%   Brainstorm's full tree) must give the shipped head.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'), fullfile(root, 'tools'));
S = static_test_settings();
S.parallel = false;
S.conn = {'wpli', 'aec_orth'};
S = setup_toolboxes(S);
all_data = static_test_recordings(S);
tmp = tempname;
mkdir(tmp);
S.data = fullfile(tmp, 'data');
for id = {'sub-control01', 'sub-patient01'}
    mkdir(fullfile(S.data, id{1}));
    copyfile(fullfile(all_data, id{1}), fullfile(S.data, id{1}));
end
copyfile(fullfile(all_data, 'participants.tsv'), S.data);
S.file_pattern = {'sub-*/eeg/*_eeg.set', 'sub-*/eeg/*_eeg.edf'};   % the patient's recording is in EDF
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.subs_info = fullfile(S.data, 'participants.tsv');
S.results = fullfile(tmp, 'shipped');
evalc('run_all(S);');
tc.TestData.S = S;
tc.TestData.tmp = tmp;
tc.TestData.head = shipped_copy(fullfile(tmp, 'heads', 'copy'));
end

function teardownOnce(tc)
if isfolder(tc.TestData.tmp), rmdir(tc.TestData.tmp, 's'); end
end

function folder = shipped_copy(folder)
% a head folder with the shipped ICBM152's files and leadfield (get_anatomy's layout)
in = fullfile(fileparts(fileparts(mfilename('fullpath'))), '1- EEG source connectivity', 'inputs');
mkdir(folder);
for f = {'tess_head', 'tess_outerskull', 'tess_innerskull'}
    copyfile(fullfile(in, [f{1} '.mat']), folder);
end
copyfile(fullfile(in, 'tess_cortex_pial_low.mat'), fullfile(folder, 'tess_cortex.mat'));
copyfile(fullfile(in, 'channel_ASA_10-05_343_bs.mat'), fullfile(folder, 'channel_ASA_10-05_343.mat'));
lf = fullfile(folder, 'ftLeadfield_ASA_10-05_343_normal.mat');
copyfile(fullfile(in, 'ftLeadfield_ASA_10-05_343_normal.mat'), lf);
anatomy_key = file_fingerprint(fullfile(folder, {'tess_head.mat', 'tess_outerskull.mat', 'tess_innerskull.mat', ...
    'tess_cortex.mat', 'channel_ASA_10-05_343.mat'})); %#ok<NASGU>   % as get_ftHeadmodel_and_ftLeadfield keeps it
save(lf, 'anatomy_key', '-append');
end

function w = output(results, rec, atlas)
if nargin < 3, atlas = 'Desikan_Killiany_68'; end
w = load(fullfile(results, 'Connectivity', 'cmats', rec, ...
    ['static_cmats_eloreta_ASA_10-05_343_' file_tag(atlas) '_alpha.mat']));
end

function quiet(f)
% f() with its printed output hidden
evalc('f()');
end

function err = caught(f)
% the error f() throws (verifyError returns f's outputs, not the error)
err = MException('test_anatomy:none', 'no error');
try
    f();
catch err
end
end

function Q = summary(results)
Q = readtable(fullfile(results, 'Connectivity', 'stage1_summary.csv'), 'TextType', 'string');
end

%% A folder holding the shipped head gives the shipped results
function test_folder_with_the_shipped_head(tc)
S = tc.TestData.S;
Q0 = summary(S.results);
verifyEqual(tc, Q0.anatomy, ["ICBM152 (shipped)"; "ICBM152 (shipped)"]);
S.anatomy = tc.TestData.head;
S.results = fullfile(tc.TestData.tmp, 'folder');
evalc('run_all(S);');
Q = summary(S.results);
verifyEqual(tc, Q.status, ["ok"; "ok"], strjoin(string(Q.message), ' | '));
verifyEqual(tc, Q.anatomy, ["copy"; "copy"]);
for rec = {'sub-control01_task-rest_eeg', 'sub-patient01_task-rest_eeg'}
    a = output(tc.TestData.S.results, rec{1});
    b = output(S.results, rec{1});
    verifyEqual(tc, b.static_full_cmats, a.static_full_cmats);   % exactly: the same leadfield and cortex
    verifyEqual(tc, b.static_power, a.static_power);
    verifyEqual(tc, b.info.roi_labels, a.info.roi_labels);
    verifyEqual(tc, b.info.settings_key, a.info.settings_key);   % the settings are the same,
    verifyNotEqual(tc, b.info.input_key, a.info.input_key);      % the head is part of the input key
    verifyEqual(tc, b.info.anatomy.name, 'copy');
    verifyEmpty(tc, a.info.anatomy.key);
end
% the states and the statistics accept these outputs with the same head, not with another
quiet(@() stage1_outputs(S, 'eloreta', 'alpha'));
S0 = S;  S0.anatomy = '';
verifyError(tc, @() quiet(@() stage1_outputs(S0, 'eloreta', 'alpha')), 'stage1_outputs:input');
% the same head again: skipped; the shipped head in this folder: recomputed
evalc('run_all(S);');
verifyEqual(tc, summary(S.results).status, ["skipped (up to date)"; "skipped (up to date)"]);
S.anatomy = '';
evalc('run_all(S);');
verifyEqual(tc, summary(S.results).status, ["ok"; "ok"]);
verifyEqual(tc, summary(S.results).anatomy, ["ICBM152 (shipped)"; "ICBM152 (shipped)"]);
end

%% One head per subject: missing ones stop the run, or use the shipped head
function test_individual_heads_and_fallback(tc)
S = tc.TestData.S;
root = fullfile(tc.TestData.tmp, 'individual');
mkdir(root);
copyfile(tc.TestData.head, fullfile(root, 'sub-control01'));   % only the control has a head
S.anatomy = root;
S.anatomy_individual = true;
err = caught(@() quiet(@() run_all(S)));
verifyEqual(tc, err.identifier, 'run_all:anatomy');
verifySubstring(tc, err.message, 'sub-patient01');
verifyFalse(tc, contains(err.message, 'sub-control01_task'));
verifyEqual(tc, get_anatomy(S, 'sub-control01_task-rest_eeg.set').name, 'sub-control01');
verifyError(tc, @() get_anatomy(S, 'sub-patient01_task-rest_eeg.set'), 'get_anatomy:missing');
S.anatomy_fallback = true;
S.results = fullfile(tc.TestData.tmp, 'fallback');             % a copy of the shipped results
copyfile(tc.TestData.S.results, S.results);
verifyWarning(tc, @() quiet(@() run_all(S)), 'run_all:anatomy_fallback');
Q = summary(S.results);
[~, i] = sort(Q.file);
Q = Q(i, :);
verifyEqual(tc, Q.anatomy, ["sub-control01"; "ICBM152 (shipped, fallback)"]);
verifyEqual(tc, Q.status, ["ok"; "skipped (up to date)"]);    % the fallback is the shipped head: reused
A = get_anatomy(S, 'sub-patient01_task-rest_eeg.set');
verifyTrue(tc, A.shipped && A.fallback);
% the folder exists but lacks a file: an error naming it
delete(fullfile(root, 'sub-control01', 'channel_ASA_10-05_343.mat'));
[~, problem] = get_anatomy(S, 'sub-control01_task-rest_eeg.set');
verifySubstring(tc, problem, 'channel_ASA_10-05_343.mat');
verifyError(tc, @() quiet(@() run_all(S)), 'run_all:anatomy');
end

%% Atlases by name, on any head, in the shipped head's order
function test_atlases_by_name(tc)
in = fullfile(tc.TestData.S.root, '1- EEG source connectivity', 'inputs');
dk = load(fullfile(in, 'scout_Desikan-Killiany_68.mat'));
Sc = get_atlas('Destrieux');
verifyEqual(tc, numel(Sc), 148);
verifyEqual(tc, numel(get_atlas('Schaefer_100_17net')), 100);
S = tc.TestData.S;
S.anatomy = tc.TestData.head;
A = get_anatomy(S);
Sc = get_atlas('Desikan_Killiany_68', A);                     % the cortex's 'Desikan-Killiany'
verifyEqual(tc, {Sc.Label}, {dk.Scouts.Label});
verifyEqual(tc, {Sc.Vertices}, {dk.Scouts.Vertices});
verifyError(tc, @() get_atlas('NoSuchAtlas', A), 'get_atlas:atlas');
% another head whose atlas lists the regions in another order, one of them
% missing, one extra: reordered like the shipped head's, the missing one an error
other = fullfile(tc.TestData.tmp, 'heads', 'shuffled');
copyfile(tc.TestData.head, other);
C = load(fullfile(other, 'tess_cortex.mat'));
i = find(strcmp({C.Atlas.Name}, 'Desikan-Killiany'));
C.Atlas(i).Scouts = C.Atlas(i).Scouts(end:-1:1);
C.Atlas(i).Scouts(end + 1) = C.Atlas(i).Scouts(1);
C.Atlas(i).Scouts(end).Label = 'unknown L';
save(fullfile(other, 'tess_cortex.mat'), '-struct', 'C');
S.anatomy = other;
B = get_anatomy(S);
Sc = get_atlas('Desikan_Killiany_68', B);
verifyEqual(tc, {Sc.Label}, {dk.Scouts.Label});
C.Atlas(i).Scouts(5) = [];
save(fullfile(other, 'tess_cortex.mat'), '-struct', 'C');
verifyError(tc, @() get_atlas('Desikan_Killiany_68', B), 'get_atlas:regions');
% through stage 1: Destrieux by name on the shipped head
S = tc.TestData.S;
S.atlases = {'Destrieux'};
S.results = fullfile(tc.TestData.tmp, 'destrieux');
S.data = fullfile(tc.TestData.tmp, 'one');
mkdir(fullfile(S.data, 'sub-control01'));
copyfile(fullfile(tc.TestData.S.data, 'sub-control01'), fullfile(S.data, 'sub-control01'));
evalc('run_all(S);');
w = output(S.results, 'sub-control01_task-rest_eeg', 'Destrieux');
verifyEqual(tc, size(w.static_full_cmats), [2, 148, 148]);
Sc = get_atlas('Destrieux');
verifyEqual(tc, w.info.roi_labels, {Sc.Label});
end

%% The settings
function test_anatomy_settings(tc)
S = tc.TestData.S;
S.anatomy = fullfile(tc.TestData.tmp, 'no_such_folder');
verifyError(tc, @() check_settings(S), 'check_settings:anatomy');
S.anatomy = '';
S.anatomy_individual = true;
verifyError(tc, @() check_settings(S), 'check_settings:anatomy');
S.anatomy = tc.TestData.tmp;
S.atlases = {fullfile(S.root, '1- EEG source connectivity', 'inputs', 'scout_Desikan-Killiany_68.mat')};
verifyError(tc, @() check_settings(S), 'check_settings:anatomy');      % a scout file is one cortex's
S = tc.TestData.S;
S.anatomy_fallback = true;
verifyWarning(tc, @() check_settings(S), 'check_settings:not_used');
S = tc.TestData.S;
S.anatomy_individual = 2;
verifyError(tc, @() check_settings(S), 'check_settings:anatomy_individual');
% the defaults: the shipped head, and the same settings key as without the anatomy fields
S = check_settings(rmfield(tc.TestData.S, {'anatomy', 'anatomy_individual', 'anatomy_fallback'}));
verifyEqual(tc, {S.anatomy, S.anatomy_individual, S.anatomy_fallback}, {'', false, false});
verifyEqual(tc, settings_key(S, '', 'eloreta', 'alpha'), settings_key(tc.TestData.S, '', 'eloreta', 'alpha'));
verifyTrue(tc, get_anatomy(S).shipped);
end

%% A leadfield is made again when the head's surfaces or cap change
function test_leadfield_of_a_changed_head(tc)
S = tc.TestData.S;
folder = fullfile(tc.TestData.tmp, 'heads', 'changed');
copyfile(tc.TestData.head, folder);
S.anatomy = folder;
A = get_anatomy(S);
[~, L] = get_ftHeadmodel_and_ftLeadfield('ASA_10-05_343', 'Desikan_Killiany_68', S, A);   % accepted: its key
verifyEqual(tc, numel(L.leadfield), 15002);
ch = load(fullfile(folder, 'channel_ASA_10-05_343.mat'));
ch.Channel(1).Loc = ch.Channel(1).Loc * 1.001;                  % the cap moved
save(fullfile(folder, 'channel_ASA_10-05_343.mat'), '-struct', 'ch');
A = get_anatomy(S);
[st, ~] = system('which om_assemble');
assumeTrue(tc, st ~= 0, 'OpenMEEG programs installed: the FieldTrip route would compute the leadfield');
S.openmeeg_python = fullfile(tc.TestData.tmp, 'no_python');     % so that it cannot be made here
verifyError(tc, @() quiet(@() get_ftHeadmodel_and_ftLeadfield('ASA_10-05_343', 'Desikan_Killiany_68', S, A)), ...
    'get_ftHeadmodel_and_ftLeadfield:openmeeg');
% the stale leadfield is kept until a new one is complete (runs sharing
% the cache never read a file being written), no partial file
% is left, and the stale one is never used: the next call tries again
verifyTrue(tc, isfile(fullfile(folder, 'ftLeadfield_ASA_10-05_343_normal.mat')));
verifyEmpty(tc, dir(fullfile(folder, '*_partial_*')));
verifyError(tc, @() quiet(@() get_ftHeadmodel_and_ftLeadfield('ASA_10-05_343', 'Desikan_Killiany_68', S, A)), ...
    'get_ftHeadmodel_and_ftLeadfield:openmeeg');
% a cap whose labels are not the montage's: an error
ch.Channel = ch.Channel([2 1 3:end]);
save(fullfile(folder, 'channel_ASA_10-05_343.mat'), '-struct', 'ch');
A = get_anatomy(S);
verifyError(tc, @() get_ftHeadmodel_and_ftLeadfield('ASA_10-05_343', 'Desikan_Killiany_68', S, A), ...
    'get_ftHeadmodel_and_ftLeadfield:labels');
end

%% make_anatomy with Brainstorm's ICBM152 template: the shipped head
function test_make_anatomy_template(tc)
S = tc.TestData.S;
bst = fileparts(fileparts(fileparts(which('bst_bandpass_hfilter'))));
assumeTrue(tc, isfile(fullfile(bst, 'defaults', 'anatomy', 'ICBM152', 'tess_head.mat')), ...
    'Brainstorm''s ICBM152 template not found');
folder = fullfile(tc.TestData.tmp, 'heads', 'ICBM152');
evalc('make_anatomy(S, ''template'', ''ICBM152'', folder, ''leadfield'', false);');
verifyEmpty(tc, which('db_add_subject'));                        % Brainstorm's tree off the path again
in = fullfile(S.root, '1- EEG source connectivity', 'inputs');
for f = {'tess_head.mat', 'tess_outerskull.mat', 'tess_innerskull.mat', 'tess_cortex.mat', ...
        'subjectimage_T1.mat', 'channel_ASA_10-05_343.mat', 'anatomy_info.txt'}
    verifyTrue(tc, isfile(fullfile(folder, f{1})), f{1});
end
C0 = load(fullfile(in, 'tess_cortex_pial_low.mat'), 'Vertices', 'VertNormals');
C = load(fullfile(folder, 'tess_cortex.mat'), 'Vertices', 'VertNormals');
verifyEqual(tc, C.Vertices, C0.Vertices);
verifyEqual(tc, C.VertNormals, C0.VertNormals, 'AbsTol', 1e-9);   % Brainstorm's interface normals
H0 = load(fullfile(in, 'tess_head.mat'), 'Vertices');
H = load(fullfile(folder, 'tess_head.mat'), 'Vertices');
verifyEqual(tc, H.Vertices, H0.Vertices);
S.anatomy = folder;
A = get_anatomy(S);
dk = load(fullfile(in, 'scout_Desikan-Killiany_68.mat'));
Sc = get_atlas('Desikan_Killiany_68', A);
verifyEqual(tc, {Sc.Label}, {dk.Scouts.Label});
e0 = get_ftElec('ASA_10-05_343');
e = BS_to_ft_channels(fullfile(folder, 'channel_ASA_10-05_343.mat'), []);
verifyEqual(tc, e.label, e0.label);
d = sqrt(sum((e.chanpos - e0.chanpos) .^ 2, 2));
verifyLessThan(tc, median(d), 0.002);                            % the MNI round trip on the same head: ~0
verifySubstring(tc, fileread(fullfile(folder, 'anatomy_info.txt')), 'Desikan_Killiany_68: 68 regions');
end
