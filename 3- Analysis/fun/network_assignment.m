function [lab, nets, column] = network_assignment(nodes)
%NETWORK_ASSIGNMENT The network of each region, for the network features and tests.
%   [lab, nets, column] = network_assignment(nodes) with nodes the network
%   table rows of the regions (load_network_table). The assignment is the
%   column rsn_yeo7 (Yeo et al. 2011, 7 networks) when the table has it,
%   otherwise its first column rsn_* (a table for another atlas):
%   lab   network of each region (string, nRegions x 1)
%   nets  the networks, "Other" left out, in the order of the column
%         <column>_code when the table has it (Yeo-7: VIS, SMN, DAN, SAL,
%         LIM, FPN, DMN), else of first appearance
%   Empty when nodes is empty or has no rsn_* column.
lab = strings(0, 1);  nets = strings(0, 1);  column = "";
if isempty(nodes), return; end
v = string(nodes.Properties.VariableNames);
cols = v(startsWith(v, 'rsn_') & ~endsWith(v, {'_code', '_fraction', '_coverage'}));
if isempty(cols), return; end
column = cols(1);
if ismember("rsn_yeo7", cols), column = "rsn_yeo7"; end
lab = string(nodes.(column));
lab = lab(:);
code = column + "_code";
if ismember(code, v)
    [~, order] = sort(nodes.(code));
    nets = unique(lab(order), 'stable');
else
    nets = unique(lab, 'stable');
end
nets = nets(nets ~= "Other");
end
