function nodes = load_network_table(S, roi_labels)
%LOAD_NETWORK_TABLE Coordinates and network assignment of each region, in the order of the matrices.
%   nodes = load_network_table(S, roi_labels) reads S.network_table (for the
%   DK-68 atlas: DK68_networks.csv with label, short_label, x_mni, y_mni,
%   z_mni, the lab's assignment rsn_lab and the Yeo-7 networks rsn_yeo7)
%   and returns its rows matched to roi_labels by label (the original read
%   cellfile_node.mat and assumed its rows were in the atlas order).
%   Empty (with a message) when there is no table or a region is missing,
%   so the network scripts can be skipped for another atlas.
nodes = table();
if isempty(S.network_table) || ~isfile(S.network_table)
    fprintf('No network table (S.network_table): network profiles and BrainNet files skipped\n');
    return
end
T = readtable(S.network_table, 'TextType', 'string');
[ok, row] = ismember(lower(strtrim(string(roi_labels(:)))), lower(strtrim(T.label)));
if ~all(ok)
    fprintf('Network table %s lacks regions (%s ...): network profiles and BrainNet files skipped\n', ...
        S.network_table, strjoin(string(roi_labels(find(~ok, 3))), ', '));
    return
end
nodes = T(row, :);
end
