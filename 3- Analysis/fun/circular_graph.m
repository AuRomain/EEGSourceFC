function circular_graph(ax, A, nodes, net, mode)
%CIRCULAR_GRAPH Connectogram: regions on a circle, grouped by network, with the edges of A.
%   circular_graph(ax, A, nodes, net, mode) draws the nonzero edges of A
%   (regions x regions, symmetric) into the axes AX. Left-hemisphere
%   regions are on the left half, right-hemisphere regions on the right
%   half, mirrored; in each half the regions are grouped by network (column
%   NET of the network table NODES, in the order of first appearance) and,
%   inside a network, ordered front to back (y_mni). Each region's dot has
%   its network's colour (network_colours), on an arc of that colour, and a
%   size that grows with the summed absolute weight of its drawn edges.
%   mode 'weight': edges in grey, darker and wider when stronger;
%   'signed': positive edges red, negative blue, wider when larger.
n = height(nodes);
lab = string(nodes.(net));
nets = unique(lab, 'stable');
nets = [nets(nets ~= "Other"); nets(nets == "Other")];
[~, g] = ismember(lab, nets);
hemi = repmat('R', n, 1);
names = strtrim(string(nodes.label));
if ismember('x_mni', nodes.Properties.VariableNames), hemi(nodes.x_mni < 0) = 'L'; end
hemi(endsWith(names, " L")) = 'L';  hemi(endsWith(names, " R")) = 'R';
front = zeros(n, 1);
if ismember('y_mni', nodes.Properties.VariableNames), front = -nodes.y_mni; end
theta = zeros(n, 1);
gap = 6 * pi / 180;                                     % between the hemispheres, at top and bottom
for h = 'LR'
    m = find(hemi == h);
    [~, o] = sortrows([g(m), front(m)]);
    m = m(o);
    t = linspace(pi / 2 + gap, 3 * pi / 2 - gap, numel(m))';   % left: top to bottom, anticlockwise
    if h == 'R', t = pi - t; end                                % right: mirrored
    theta(m) = t;
end
xy = [cos(theta), sin(theta)];
col = network_colours(nets);
hold(ax, 'on');
% network arcs
for h = 'LR'
    for k = 1:numel(nets)
        m = find(hemi == h & g == k);
        if isempty(m), continue; end
        t = linspace(min(theta(m)) - 0.02, max(theta(m)) + 0.02, 30);
        plot(ax, 1.04 * cos(t), 1.04 * sin(t), '-', 'Color', col(k, :), 'LineWidth', 3);
    end
end
% edges, weakest first
A(1:n + 1:end) = 0;
[i, j] = find(triu(A ~= 0, 1));
w = A(sub2ind([n n], i, j));
[~, o] = sort(abs(w));
i = i(o); j = j(o); w = w(o);
wmax = max(abs(w));
t = linspace(0, 1, 40)';
for e = 1:numel(w)
    p0 = xy(i(e), :);  p2 = xy(j(e), :);
    p1 = (p0 + p2) / 2 * (1 - norm(p0 - p2) / 2);          % near regions: shallow; far: through the centre
    c = (1 - t) .^ 2 .* p0 + 2 * (1 - t) .* t .* p1 + t .^ 2 .* p2;
    s = abs(w(e)) / wmax;
    if strcmp(mode, 'signed')
        rgb = [0.80 0.12 0.15];
        if w(e) < 0, rgb = [0.10 0.35 0.75]; end
        plot(ax, c(:, 1), c(:, 2), '-', 'Color', [rgb, 0.25 + 0.65 * s], 'LineWidth', 0.4 + 2.2 * s);
    else
        plot(ax, c(:, 1), c(:, 2), '-', 'Color', [repmat(0.7 - 0.6 * s, 1, 3), 0.3 + 0.6 * s], 'LineWidth', 0.4 + 2.2 * s);
    end
end
% regions
strength = sum(abs(A), 2);
sz = 6 + 40 * strength / max([strength; eps]);
scatter(ax, xy(:, 1), xy(:, 2), sz, col(g, :), 'filled', 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.25);
short = names;
if ismember('short_label', nodes.Properties.VariableNames), short = string(nodes.short_label); end
for r = 1:n
    deg = theta(r) * 180 / pi;
    if hemi(r) == 'L'
        text(ax, 1.09 * xy(r, 1), 1.09 * xy(r, 2), short(r), 'Rotation', deg - 180, ...
            'HorizontalAlignment', 'right', 'FontSize', 4.5, 'Interpreter', 'none');
    else
        text(ax, 1.09 * xy(r, 1), 1.09 * xy(r, 2), short(r), 'Rotation', deg, ...
            'HorizontalAlignment', 'left', 'FontSize', 4.5, 'Interpreter', 'none');
    end
end
text(ax, -1.42, 1.32, 'L', 'FontWeight', 'bold');
text(ax, 1.35, 1.32, 'R', 'FontWeight', 'bold');
axis(ax, 'equal', 'off');
xlim(ax, [-1.45 1.45]);  ylim(ax, [-1.4 1.4]);
end
