function elec = BS_to_ft_channels(bs_channel_file, eeg_id)
%BS_TO_FT_CHANNELS Brainstorm channel file -> FieldTrip electrode structure.
%   elec = BS_to_ft_channels(file, eeg_id) keeps the channels eeg_id
%   (all EEG channels when eeg_id is empty). A channel without a position
%   now stops the run instead of being placed at the origin (0, 0, 0).
%   FILE can also be a CSV with columns label,x,y,z (metres, ICBM152
%   coordinates), for caps that are not in Brainstorm's format.
[~, ~, ext] = fileparts(bs_channel_file);
if any(strcmpi(ext, {'.csv', '.txt', '.tsv'}))
    T = readtable(bs_channel_file, 'TextType', 'string');
    T.Properties.VariableNames = lower(T.Properties.VariableNames);
    if ~all(ismember({'label', 'x', 'y', 'z'}, T.Properties.VariableNames))
        error('BS_to_ft_channels:csv', '%s needs the columns label, x, y, z (metres)', bs_channel_file);
    end
    if nargin < 2 || isempty(eeg_id), eeg_id = 1:height(T); end
    elec = struct('label', {cellstr(T.label(eeg_id))}, 'unit', 'm', 'chanpos', [T.x(eeg_id), T.y(eeg_id), T.z(eeg_id)]);
else
    a = load(bs_channel_file);
    Channel = a.Channel;
    if nargin < 2 || isempty(eeg_id)
        eeg_id = find(strcmpi({Channel.Type}, 'EEG'));
    end
    elec = struct();
    elec.label = {Channel(eeg_id).Name}';
    elec.unit = 'm';
    elec.chanpos = zeros(numel(eeg_id), 3);
    for i = 1:numel(eeg_id)
        loc = Channel(eeg_id(i)).Loc;
        if numel(loc) < 3
            error('BS_to_ft_channels:position', 'Channel %s has no position', Channel(eeg_id(i)).Name);
        end
        elec.chanpos(i, :) = loc(:, 1)';
    end
end
if any(~isfinite(elec.chanpos(:))) || any(all(elec.chanpos == 0, 2))
    error('BS_to_ft_channels:position', '%s: a channel has no position', bs_channel_file);
end
if numel(unique(lower(elec.label))) < numel(elec.label)
    error('BS_to_ft_channels:labels', '%s: duplicate channel labels', bs_channel_file);
end
elec.elecpos = elec.chanpos;
elec.tra = eye(numel(elec.label));
end
