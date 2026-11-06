function S = setup_toolboxes(S)
%SETUP_TOOLBOXES Put the pipeline folders and the toolboxes it uses on the path.
%   S = setup_toolboxes() or setup_toolboxes(S) (S from pipeline_settings).
%   A toolbox already on the MATLAB path is used as it is; otherwise it is
%   looked for in S.toolboxes (the Toolboxes folder next to this pipeline),
%   or in S.toolbox_paths.<name> when set (the toolbox's root folder, e.g.
%   .../brainstorm3, or the folder of its functions). Toolboxes (see README for links):
%     FieldTrip      reading EEG files, eLORETA and LCMV
%     Brainstorm 3   band-pass filter, scout values (two folders only: the
%                    whole tree would shadow MATLAB functions)
%     modkmeans      Poulsen's modified k-means and calc_fitmeas
%     CVIK           cluster validity indices of the criteria function
%     BCT            graph metrics of the language lists
%                    (these three only with states: not in a static-only
%                    run, S.dynamic = false)
%     NBS 1.2        network-based statistic: included (3- Analysis/fun/NBS,
%                    the two functions used, unchanged; GPL-3), not looked for
%     restingIAF     individual alpha frequency (Corcoran et al. 2018;
%                    github.com/corcorana/restingIAF), only when a band is
%                    anchored on it (S.bands or S.psd.bands, band_edges) or
%                    with the spectral profile of stage 3 (S.spectral_profile)
%     NMF_neurodevelopment  only with S.energy_entropy (off by default):
%                    signal_entropy for the energy-entropy features. Not put
%                    on the path: its histogram.m would shadow MATLAB's
%                    histogram; normalized_energy_entropy adds it only while it runs.
%   Also records the versions used in S.versions (saved with the results),
%   and compares the toolboxes with the versions validated with the
%   pipeline (toolbox_check, toolbox_reference.json; S.toolbox_check:
%   'warn' (default), 'error' or 'off'): a toolbox that differs gets a
%   warning naming it; tools/fetch_toolboxes.sh gets the validated versions.
%   On a cap file whose channels have other names than the 10-05 ones (an
%   ANT equidistant or EGI cap), S.iaf.channels gets the cap's electrodes
%   nearest to the positions of its names (cap_channels), with a message.
if nargin < 1, S = pipeline_settings(); end
S = check_settings(S);
tb = S.toolboxes;
% name, folder under S.toolboxes, function that must then be found
% and the folder of that function inside the toolbox's own root (S.toolbox_paths
% may give the root, e.g. brainstorm3, or that folder, e.g. brainstorm3/toolbox/math)
need = {
    'fieldtrip',  'fieldtrip-toolbox',                                         'ft_sourceanalysis',     ''
    'brainstorm', fullfile('brainstorm3', 'toolbox', 'math'),                  'bst_bandpass_hfilter',  fullfile('toolbox', 'math')
    'modkmeans',  'modkmeans',                                                 'calc_fitmeas',          ''
    'cvik',       fullfile('GitHub resources', 'cvik-toolbox-master', 'cvi'),  'chindex',               'cvi'
    'bct',        fullfile('BCT', '2019_03_03_BCT'),                           'clustering_coef_wu',    ''
    };
if ~S.dynamic                                      % static only: no states, no clustering
    need = need(1:2, :);
end
psd_iaf = any(arrayfun(@(b) anchored(b, S.iaf), S.psd.bands));   % PSD bands anchored on the IAF
if S.iaf_used || S.spectral_profile || psd_iaf     % individual bands, spectral profile, PSD
    need(end + 1, :) = {'restingiaf', 'restingIAF', 'restingIAF', ''};
end
missing = {};
for k = 1:size(need, 1)
    [name, sub, fun, inner] = need{k, :};
    if ~isempty(which(fun)) && ~strcmp(name, 'fieldtrip'), continue; end
    folder = fullfile(tb, sub);
    if isfield(S, 'toolbox_paths') && isfield(S.toolbox_paths, name)
        folder = S.toolbox_paths.(name);
        if ~isfile(fullfile(folder, [fun '.m'])) && isfile(fullfile(folder, inner, [fun '.m']))
            folder = fullfile(folder, inner);              % the toolbox's root was given
        end
    end
    if ~isfolder(folder) && isempty(which(fun))
        missing{end + 1} = sprintf('%s (%s)', name, folder); %#ok<AGROW>
        continue
    end
    switch name
        case 'fieldtrip'
            if isempty(which('ft_defaults')), addpath(folder); end
            ft_defaults;
        case 'brainstorm'
            addpath(folder);                                   % bst_bandpass_hfilter, bst_scout_value
            addpath(fullfile(fileparts(folder), 'core'));      % bst_get (prefixed names: no shadowing)
        case 'cvik'
            addpath(folder);
            addpath(fullfile(fileparts(folder), 'proximity'));  % eucdist, proxconfig
        otherwise
            addpath(folder);
    end
end
if ~isempty(missing)
    error('setup_toolboxes:missing', 'Toolboxes not found: %s.\nSee README.md, sections 1 (What you need) and 2 (Installation).', ...
        strjoin(missing, '; '));
end
global ft_default %#ok<GVMIS>
ft_default.checksize = 1e7;
ft_default.showcallinfo = 'no';
ft_default.trackcallinfo = 'no';

% signal_entropy (only with S.energy_entropy): found, not put on the path (see above)
S.signal_entropy_dir = '';
if S.energy_entropy
    S.signal_entropy_dir = fullfile(tb, 'GitHub resources', 'NMF_neurodevelopment-master', 'matlab_functions');
    if isfield(S, 'toolbox_paths') && isfield(S.toolbox_paths, 'signal_entropy')
        S.signal_entropy_dir = S.toolbox_paths.signal_entropy;
    end
    if ~isfile(fullfile(S.signal_entropy_dir, 'signal_entropy.m'))
        error('setup_toolboxes:entropy', ['S.energy_entropy is on but signal_entropy.m is not in %s ' ...
            '(set S.toolbox_paths.signal_entropy, or S.energy_entropy = false)'], S.signal_entropy_dir);
    end
end

% pipeline folders (those present: a static-only copy of the pipeline,
% v_static, has no "2- Dynamic connectivity")
folders = {S.root, ...
    fullfile(S.root, '1- EEG source connectivity'), fullfile(S.root, '1b- PSD analysis'), ...
    fullfile(S.root, '2- Dynamic connectivity'), fullfile(S.root, '2- Dynamic connectivity', 'fun'), ...
    fullfile(S.root, '2- Dynamic connectivity', 'fun', 'Clustering'), ...
    fullfile(S.root, '2- Dynamic connectivity', 'fun', 'Backfitting'), ...
    fullfile(S.root, '2- Dynamic connectivity', 'fun', 'ICA'), ...
    fullfile(S.root, '3- Analysis'), fullfile(S.root, '3- Analysis', 'fun'), ...
    fullfile(S.root, '3- Analysis', 'fun', 'NBS')};
folders = folders(cellfun(@isfolder, folders));
addpath(folders{:});

% the posterior channels of the IAF (S.iaf.channels, 10-05 names) on a cap
% file whose channels have other names (an ANT equidistant or EGI cap): its
% electrodes nearest to their 10-05 positions, added to the names
% (cap_channels). Caps with 10-20 names: S.iaf unchanged.
for m = 1:numel(S.montages)
    if ~isfile(S.montages{m}), continue; end         % 'ASA_10-05_343', 'EEGLAB_61': 10-20 names
    elec = get_ftElec(S.montages{m});
    [ch, how] = cap_channels(S.iaf.channels, elec.label, elec);
    if ~isempty(how)
        S.iaf.channels = [cellstr(S.iaf.channels(:)'), ch];
        fprintf('S.iaf.channels: %s has none of these names; by position: %s\n', S.montages{m}, how);
    elseif isempty(ch) && (S.iaf_used || S.spectral_profile || psd_iaf)
        warning('setup_toolboxes:iaf_channels', ['%s has none of the posterior channels of S.iaf.channels, ' ...
            'by name or within 20 mm of their 10-05 positions: set S.iaf.channels to its posterior channels'], ...
            S.montages{m});
    end
end

% versions, for the provenance files written by each stage
S.versions.matlab = version;
S.versions.computer = computer;
S.versions.fieldtrip = toolbox_version('ft_defaults');
S.versions.brainstorm = toolbox_version('bst_bandpass_hfilter');
S.versions.pipeline = pipeline_version(S.root);
if ~strcmp(S.toolbox_check, 'off')
    R = toolbox_check(S);
    S.versions.toolboxes = R;
    for i = find(R.status ~= "not found")'
        S.versions.fingerprints.(char(R.toolbox(i))) = char(R.fingerprint(i));
    end
    bad = R.status == "differs" | R.status == "no reference";
    if any(bad)
        msg = sprintf(['Toolbox(es) not the versions validated with this pipeline: %s. Results may differ ' ...
            '(README, sections 1 and 2; tools/fetch_toolboxes.sh gets the validated versions)'], ...
            strjoin(R.toolbox(bad) + " (" + R.status(bad) + ", " + R.root(bad) + ")", '; '));
        if strcmp(S.toolbox_check, 'error'), error('setup_toolboxes:version', '%s', msg); end
        warning('setup_toolboxes:version', '%s', msg);
    end
    eqv = R.status == "equivalent";
    if any(eqv)
        fprintf('Published versions equivalent to the validated copies (toolbox_reference.json): %s\n', ...
            strjoin(R.toolbox(eqv), ', '));
    end
end
end

function v = toolbox_version(fun)
% git commit, or Brainstorm's doc/version.txt, or the folder's date
v = 'unknown';
d = fileparts(which(fun));
for up = 1:4
    if isfolder(fullfile(d, '.git'))
        % --no-pager: with a terminal (matlab -batch from a shell) git would
        % otherwise open a pager and wait for a key
        [st, out] = system(sprintf('git --no-pager -C "%s" log -1 --format="%%h %%cs"', d));
        if st == 0, v = ['git ' strtrim(out)]; end
        return
    end
    if isfile(fullfile(d, 'doc', 'version.txt'))
        txt = splitlines(fileread(fullfile(d, 'doc', 'version.txt')));
        v = strtrim(erase(strjoin(txt(1:min(3, end)), ' '), '%'));
        return
    end
    d = fileparts(d);
end
info = dir(which(fun));
if ~isempty(info), v = ['file date ' char(datetime(info.datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd'))]; end
end

function v = pipeline_version(root)
v = 'not a git repository';
if isfolder(fullfile(root, '.git'))
    [st, out] = system(sprintf('git --no-pager -C "%s" describe --always --dirty', root));
    if st == 0, v = strtrim(out); end
end
end

function tf = anchored(band, opt)
[~, ~, tf] = band_edges(band, 10, opt);
end
