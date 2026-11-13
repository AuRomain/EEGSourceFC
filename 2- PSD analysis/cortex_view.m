function cortex_view(ax, geo, P, clim, cmap, view_name)
%CORTEX_VIEW One view of the cortex coloured by region (or by vertex), with marks.
%   cortex_view(ax, geo, P, clim, cmap, view_name), geo from
%   cortex_geometry; view_name 'left lateral', 'left medial', 'right medial',
%   'right lateral' (that hemisphere), 'dorsal', 'ventral' (both, the front
%   at the top), or '3d' (both, lit from every side, to be turned).
%   P.val    one value per region (NaN: grey), coloured on clim with cmap;
%   P.sig, P.clu   regions marked by a filled dot (significant) or an open
%            circle (in a significant cluster); the other regions are paler
%   P.vertex (optional) one value per vertex of the cortex instead: the
%            colour follows the vertices (interpolated inside each
%            triangle), no region is paler or marked (P.val, P.sig and
%            P.clu are then not used)
hold(ax, 'on');
both = ismember(view_name, {'dorsal', 'ventral', '3d'});
left = startsWith(view_name, 'left');
lateral = endsWith(view_name, 'lateral');
h = 1 + ~left;
grey = [0.78 0.78 0.78];
col = repmat(grey, size(geo.V, 1), 1);
by_vertex = isfield(P, 'vertex') && ~isempty(P.vertex);
if by_vertex                                          % every vertex its own value
    v = double(P.vertex(:));
    col = colour(v, clim, cmap);
    col(isnan(v), :) = repmat(grey, nnz(isnan(v)), 1);
    P.sig = false(max([geo.region; 0]), 1);  P.clu = P.sig;    % no marks
else
    in = geo.region > 0;
    v = P.val(geo.region(in));
    c = colour(v, clim, cmap);
    c(isnan(v), :) = repmat(grey, nnz(isnan(v)), 1);
    pale = ~(P.sig(geo.region(in)) | P.clu(geo.region(in)));
    c(pale, :) = 0.6 * c(pale, :) + 0.4 * repmat(grey, nnz(pale), 1);
    col(in, :) = c;
end
if both, F = geo.F; else, F = geo.F(all(geo.vh(geo.F) == h, 2), :); end
patch(ax, 'Faces', F, 'Vertices', geo.V, 'FaceVertexCData', col, 'FaceColor', 'interp', 'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', 'AmbientStrength', 0.45, 'DiffuseStrength', 0.65, 'SpecularStrength', 0.05);
up = [0 0 1] - geo.lr * geo.lr(3);  up = up / norm(up);   % superior
ant = cross(geo.lr, up);                              % anterior (left x superior)
switch view_name
    case 'dorsal',  cam = up;  top = ant;
    case 'ventral', cam = -up;  top = ant;
    case '3d',      cam = geo.lr + 0.35 * up - 0.25 * ant;  top = up;   % from the left, a little above
    otherwise                                         % from the left (+lr) or the right
        if left == lateral, cam = geo.lr; else, cam = -geo.lr; end
        top = up;
end
view(ax, cam);
camup(ax, top);
axis(ax, 'equal', 'off');
camlight(ax, 'headlight');                            % in 3-D it follows the camera (cortex_3d)
% marks: in a fixed view on each region's vertex nearest the camera (a region
% that cannot be seen from there is hidden by the surface in front of it);
% in 3-D on its most outward vertex
if both, k = find(P.sig | P.clu); else, k = find(geo.hemi(:) == h & (P.sig | P.clu)); end
k = k(:);
pt = geo.point(k);
if ~strcmp(view_name, '3d')
    toward = cam(:) / norm(cam);
    for i = 1:numel(k)
        vk = find(geo.region == k(i));
        [~, j] = max(geo.V(vk, :) * toward);
        pt(i) = vk(j);
    end
    x = geo.V(pt, :) + 0.002 * toward';               % 2 mm towards the camera
else
    x = geo.V(pt, :) + 0.002 * geo.N(pt, :);
end
s = P.sig(k);
plot3(ax, x(~s, 1), x(~s, 2), x(~s, 3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot3(ax, x(s, 1), x(s, 2), x(s, 3), 'o', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
colormap(ax, cmap);
caxis(ax, clim); %#ok<CAXIS>
end

function c = colour(v, clim, cmap)
n = size(cmap, 1);
i = round((v - clim(1)) / diff(clim) * (n - 1)) + 1;
i = min(max(i, 1), n);
i(isnan(i)) = 1;
c = cmap(i, :);
end
