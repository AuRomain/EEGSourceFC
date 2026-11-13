function rgb = network_colours(names)
%NETWORK_COLOURS Colour of each network: the Yeo et al. (2011) colours, grey otherwise.
%   rgb = network_colours(names): one row (0-1) per name. The names of rsn_lab
%   take the colour of the matching Yeo-7 network (SAN = SAL, MOT = SMN);
%   limbic is darkened from Yeo's pale green so that it shows on white.
names = upper(string(names(:)));
tab = {
    'VIS', [120 18 134]
    'SMN', [70 130 180];   'MOT', [70 130 180]
    'DAN', [0 118 14]
    'SAL', [196 58 250];   'SAN', [196 58 250];   'VAN', [196 58 250]
    'LIM', [150 175 80]
    'FPN', [230 148 34];   'CON', [230 148 34]
    'DMN', [205 62 78]
    'AUD', [0 160 160]
    };
rgb = repmat([150 150 150], numel(names), 1);
for i = 1:size(tab, 1)
    rgb(names == tab{i, 1}, :) = repmat(tab{i, 2}, nnz(names == tab{i, 1}), 1);
end
rgb = rgb / 255;
end
