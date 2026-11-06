function glass_brain(ax, G, xyz, A, rgb, view_name, mode)
%GLASS_BRAIN Edges and regions inside a transparent cortex (BrainNet Viewer style).
%   glass_brain(ax, G, xyz, A, rgb, view_name, mode): the cortex of G
%   (atlas_regions: ICBM152 in MNI mm), the regions at XYZ (MNI) with
%   colours RGB (regions x 3) and a size that grows with the summed absolute
%   weight of their edges, and the nonzero edges of A as straight lines;
%   mode 'weight' (grey) or 'signed' (red positive, blue negative), as in
%   circular_graph. view_name: 'left', 'right' or 'top'. Drawn with MATLAB
%   only: BrainNet Viewer 1.7 (2019) stops in R2026a (its whitebg call), so
%   it cannot make these images in batch; its .node/.edge files are still
%   written (compute_export_nodes_edges) for use on an older MATLAB.
n = size(xyz, 1);
hold(ax, 'on');
patch(ax, 'Faces', G.faces, 'Vertices', G.vertices, 'FaceColor', [0.82 0.82 0.82], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.10, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.6);
A(1:n + 1:end) = 0;
[i, j] = find(triu(A ~= 0, 1));
w = A(sub2ind([n n], i, j));
[~, o] = sort(abs(w));
i = i(o); j = j(o); w = w(o);
wmax = max(abs(w));
for e = 1:numel(w)
    s = abs(w(e)) / wmax;
    if strcmp(mode, 'signed')
        c = [0.80 0.12 0.15];
        if w(e) < 0, c = [0.10 0.35 0.75]; end
        c = [c, 0.3 + 0.6 * s];
    else
        c = [repmat(0.7 - 0.6 * s, 1, 3), 0.35 + 0.55 * s];
    end
    plot3(ax, xyz([i(e) j(e)], 1), xyz([i(e) j(e)], 2), xyz([i(e) j(e)], 3), '-', 'Color', c, 'LineWidth', 0.5 + 2.5 * s);
end
strength = sum(abs(A), 2);
scatter3(ax, xyz(:, 1), xyz(:, 2), xyz(:, 3), 8 + 60 * strength / max([strength; eps]), rgb, 'filled', ...
    'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.25);
switch view_name
    case 'left', view(ax, -90, 0);
    case 'right', view(ax, 90, 0);
    case 'top', view(ax, 0, 90);
    otherwise, error('glass_brain:view', 'view: left, right or top');
end
camlight(ax, 'headlight');
axis(ax, 'equal', 'off');
end
