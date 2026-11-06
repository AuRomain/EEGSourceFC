function tests = test_same_code_as_v_dynamic
%TEST_SAME_CODE_AS_V_DYNAMIC The files shared with the full pipeline are unchanged copies of it.
%   From v_static/:  results = runtests('tests/test_same_code_as_v_dynamic.m')
%   tools/v_dynamic_files.txt lists every file of this folder copied from
%   v_dynamic by tools/sync_from_v_dynamic.sh, with its SHA-256. Each
%   must still have it (an edit here would make this folder another
%   pipeline: change the file in v_dynamic, test it there, and sync); when
%   v_dynamic is next to this folder, each must also equal v_dynamic's
%   file (else: run the sync). Every other code file must be one of those
%   written for this folder.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
txt = splitlines(strtrim(fileread(fullfile(root, 'tools', 'v_dynamic_files.txt'))));
txt = txt(~startsWith(txt, '#'));
tok = regexp(txt, '^([0-9a-f]{64})  (.+)$', 'tokens', 'once');
verifyTrue(tc, all(~cellfun(@isempty, tok)), 'tools/v_dynamic_files.txt: a line is not "SHA-256  path"');
tok = vertcat(tok{:});
tc.TestData.root = root;
tc.TestData.sha = tok(:, 1);
tc.TestData.files = tok(:, 2);
end

function test_shared_files_are_unchanged(tc)
root = tc.TestData.root;
for i = 1:numel(tc.TestData.files)
    f = fullfile(root, tc.TestData.files{i});
    verifyTrue(tc, isfile(f), ['missing: ' tc.TestData.files{i}]);
    if isfile(f)
        verifyEqual(tc, sha256(f), tc.TestData.sha{i}, ['changed since the copy from v_dynamic: ' tc.TestData.files{i}]);
    end
end
verifyGreaterThan(tc, numel(tc.TestData.files), 100);
end

function test_shared_files_equal_v_dynamic(tc)
v_dynamic = fullfile(fileparts(tc.TestData.root), 'v_dynamic');
assumeTrue(tc, isfolder(v_dynamic), 'v_dynamic is not next to this folder');
for i = 1:numel(tc.TestData.files)
    f = tc.TestData.files{i};
    verifyTrue(tc, isfile(fullfile(v_dynamic, f)) && strcmp(sha256(fullfile(v_dynamic, f)), tc.TestData.sha{i}), ...
        ['differs from v_dynamic (bash tools/sync_from_v_dynamic.sh): ' f]);
end
end

function test_no_other_code(tc)
% every code file is shared with v_dynamic or written for this folder: no
% leftover of an earlier version. A study's settings function in the top
% folder (settings_<study>.m, README, section 4: Your study, step by step) is the user's.
root = tc.TestData.root;
own = {'pipeline_settings.m', 'tests/test_same_code_as_v_dynamic.m', 'tests/test_template.m', 'tools/sync_from_v_dynamic.sh'};
d = [dir(fullfile(root, '**', '*.m')); dir(fullfile(root, '**', '*.py')); dir(fullfile(root, '**', '*.sh'))];
rel = strrep(strrep(fullfile({d.folder}, {d.name}), [root filesep], ''), filesep, '/');
rel = rel(~startsWith(rel, {'results', 'simulated', 'examples/'}) & cellfun(@isempty, regexp(rel, '^settings_[^/]*\.m$')));
other = setdiff(rel, [tc.TestData.files; own(:)]);
verifyEmpty(tc, other, ['code files neither shared with v_dynamic nor listed as this folder''s own: ' strjoin(other, ', ')]);
end

function test_no_v_dynamic_file_missing(tc)
% every function called by a code file here is here too: with v_dynamic's code at
% the end of the path, no direct dependency resolves into v_dynamic (a file that
% the list of tools/sync_from_v_dynamic.sh would have missed), except those that
% run_statistics.m calls only with S.dynamic = true (the states and the
% windows), the stage 2 of the test that is skipped here, and the duration
% check of tests/test_reference_statistics.m, skipped here (no states)
dynamic_only = {'2- Dynamic connectivity/compute_modkmeans_states_meta_crit.m', ...
    '3- Analysis/fun/duration_rank_summary.m', ...
    '3- Analysis/compute_export_nodes_edges.m', '3- Analysis/compute_language_network_strength.m', ...
    '3- Analysis/duration_checks.m', '3- Analysis/node_strength_RSN.m', ...
    '3- Analysis/permu_ttest_permu_corr_GT_multcomp_BNS.m', ...
    '3- Analysis/permu_ttest_permu_corr_language_GT_multcomp_BNS.m', ...
    '3- Analysis/permu_ttest_permu_corr_microfeatures_multcomp_BNS.m', '3- Analysis/plot_state_networks.m', ...
    '3- Analysis/state_free_features.m', '3- Analysis/state_network_contrast.m', '3- Analysis/state_reliability.m', ...
    '3- Analysis/static_language_GT_metrics_stats_permu.m', '3- Analysis/surrogate_variability.m', ...
    '3- Analysis/window_measures.m'};
root = tc.TestData.root;
v_dynamic = fullfile(fileparts(root), 'v_dynamic');
assumeTrue(tc, isfolder(v_dynamic), 'v_dynamic is not next to this folder');
p = path;
c = onCleanup(@() path(p));
addpath(v_dynamic, genpath(fullfile(v_dynamic, '1- EEG source connectivity')), genpath(fullfile(v_dynamic, '2- Dynamic connectivity')), ...
    genpath(fullfile(v_dynamic, '3- Analysis')), fullfile(v_dynamic, 'simulation'), fullfile(v_dynamic, 'tools'), fullfile(v_dynamic, 'tests'), '-end');
addpath(root, genpath(fullfile(root, '1- EEG source connectivity')), genpath(fullfile(root, '3- Analysis')), ...
    fullfile(root, 'simulation'), fullfile(root, 'tools'), fullfile(root, 'examples'), fullfile(root, 'tests'), '-begin');
d = dir(fullfile(root, '**', '*.m'));
f = fullfile({d.folder}, {d.name});
f = f(~contains(f, [filesep 'results' filesep]));
dep = matlab.codetools.requiredFilesAndProducts(f, 'toponly');
in_v_dynamic = strrep(strrep(dep(startsWith(dep, [v_dynamic filesep])), [v_dynamic filesep], ''), filesep, '/');
missing = setdiff(in_v_dynamic, dynamic_only);
verifyEmpty(tc, missing, ['called here but only in v_dynamic (add them to tools/sync_from_v_dynamic.sh): ' strjoin(missing, ', ')]);
verifyGreaterThan(tc, numel(f), 100);
end

function h = sha256(file)
fid = fopen(file, 'r');
b = fread(fid, inf, '*uint8');
fclose(fid);
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(b);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
end
