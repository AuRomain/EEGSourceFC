function subs = get_subs(path, pattern)
%GET_SUBS File names of the recordings in a directory.
%   subs = get_subs(path, pattern) lists the files matching PATTERN, a
%   pattern or a cell of patterns (default '*.set'), e.g. {'*.edf', '*.vhdr'}.
%   A pattern may name subfolders ('sub-*/eeg/*_eeg.set' for a BIDS
%   dataset, '**/*.vhdr' for any depth): the file is then given with its
%   path relative to PATH (callers read fullfile(path, subs{i}) and name
%   the recording by its file name without extension), and two files with
%   the same name in different folders are an error, since their outputs
%   would share a folder.
%   Only matching files are returned, so system files such as .DS_Store or
%   project_state.mat cannot be taken for subjects; the list is
%   sorted, so the subject order never depends on the file system. EEGLAB
%   .fdt and BrainVision .eeg/.vmrk companions are not listed separately.
%   An EGI recording is a folder <name>.mff (FieldTrip reads it with
%   mffmatlabio, which needs Java): such folders are listed like files
%   ('*.mff', or one named in full), and the files inside them never are;
%   other folders never are listed.
if nargin < 2, pattern = '*.set'; end
pattern = cellstr(pattern);
base = dir(path);
if isempty(base)
    error('get_subs:none', 'Folder %s not found or empty', path);
end
base = base(1).folder;                 % PATH as dir writes it, to make the paths relative
subs = {};
for p = 1:numel(pattern)
    f = fullfile(path, pattern{p});
    if endsWith(f, '.mff', 'IgnoreCase', true) && isfolder(f)
        % one .mff named in full: the folder itself (dir would list its files)
        [folder, name, ext] = fileparts(f);
        d = dir(fullfile(folder, [name ext '*']));
        d = d(strcmp({d.name}, [name ext]));
    else
        d = dir(f);
    end
    d = d((~[d.isdir] | endsWith({d.name}, '.mff', 'IgnoreCase', true)) & ~startsWith({d.name}, '.') & ...
        ~endsWith({d.folder}, '.mff', 'IgnoreCase', true));          % nor the files inside a .mff
    for i = 1:numel(d)
        if strcmp(d(i).folder, base)
            subs{end + 1, 1} = d(i).name; %#ok<AGROW>
        else                               % in a subfolder: its path relative to PATH
            subs{end + 1, 1} = fullfile(d(i).folder(numel(base) + 2:end), d(i).name); %#ok<AGROW>
        end
    end
end
subs = subs(~endsWith(subs, {'.fdt', '.vmrk', '.eeg'}, 'IgnoreCase', true) | ...
    any(endsWith(pattern, {'.eeg'}, 'IgnoreCase', true)));
if isempty(subs)
    error('get_subs:none', 'No files matching %s in %s', strjoin(pattern, ', '), path);
end
subs = unique(subs);                   % sorted, and a file matched by two patterns counted once
[~, stems] = cellfun(@fileparts, subs, 'UniformOutput', false);
[u, ~, j] = unique(stems);
twice = u(accumarray(j, 1) > 1);
if ~isempty(twice)
    error('get_subs:duplicate', ['Recordings with the same name in different folders (their outputs ' ...
        'would share a folder): %s'], strjoin(subs(ismember(stems, twice)), ', '));
end
end
