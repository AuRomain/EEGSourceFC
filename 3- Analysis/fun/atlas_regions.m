function G = atlas_regions(S, roi_labels, nodes)
%ATLAS_REGIONS Geometry of the atlas regions, for the spin test and the brain figures.
%   G = atlas_regions(S, roi_labels, nodes) for the regions of the run's
%   atlas (S.atlases{1}) in the order of roi_labels, from the cortex of
%   stage 1 (inputs/tess_cortex_pial_low.mat, ICBM152):
%     G.hemi      'L' or 'R' (the label's " L"/" R" suffix, else the side
%                 of the region's MNI x)
%     G.sphere    unit vector to the region's centroid on the FreeSurfer
%                 sphere of its hemisphere (the cortex's Reg.Sphere), for
%                 spin_permutations; [] when the cortex has no sphere
%     G.vertices, G.faces  the cortex in MNI coordinates (mm): the affine
%                 transform that best maps the regions' cortex centroids
%                 to the x_mni, y_mni, z_mni of NODES (the network table);
%                 G.fit_mm is the median distance left (about 2 mm for
%                 DK-68, whose table centroids were computed on another
%                 surface); [] without table centroids
%   Fields are empty, with a message, when the atlas or a region is
%   missing. The scouts are read as stage 1 reads them (get_atlas: an atlas
%   file, or an atlas of the cortex by name). Another head for every
%   recording (S.anatomy, a folder): its cortex; one head per subject
%   (S.anatomy_individual): the shipped ICBM152's, whose regions are the
%   same (get_atlas orders every head's regions like it).
G = struct('hemi', '', 'sphere', [], 'vertices', [], 'faces', [], 'fit_mm', NaN);
atlas = S.atlases{1};
A = [];                                            % the shipped ICBM152
if isfield(S, 'anatomy') && ~isempty(S.anatomy) && ~(isfield(S, 'anatomy_individual') && S.anatomy_individual)
    A = get_anatomy(S);
end
try
    [Sc, C] = get_atlas(atlas, A, {'Vertices', 'Faces', 'Reg'});
catch err
    fprintf('Atlas %s: %s: no spin test and no brain figures\n', atlas, err.message);
    return
end
[ok, row] = ismember(lower(strtrim(string(roi_labels(:)))), lower(strtrim(string({Sc.Label}))));
if ~all(ok)
    fprintf('Atlas %s lacks regions (%s): no spin test and no brain figures\n', atlas, ...
        strjoin(string(roi_labels(find(~ok, 3))), ', '));
    return
end
Sc = Sc(row);
n = numel(Sc);
cen = zeros(n, 3);
for i = 1:n, cen(i, :) = mean(C.Vertices(Sc(i).Vertices, :), 1); end

lab = strtrim(string(roi_labels(:)));
hemi = repmat(' ', n, 1);
hemi(endsWith(lab, [" L", "_L", ".L", "-L", "lh"])) = 'L';
hemi(endsWith(lab, [" R", "_R", ".R", "-R", "rh"])) = 'R';
has_xyz = all(ismember({'x_mni', 'y_mni', 'z_mni'}, nodes.Properties.VariableNames));
if has_xyz
    xyz = [nodes.x_mni, nodes.y_mni, nodes.z_mni];
    T = [cen, ones(n, 1)] \ xyz;                      % affine, cortex coordinates to MNI
    G.vertices = [C.Vertices, ones(size(C.Vertices, 1), 1)] * T;
    G.faces = C.Faces;
    G.fit_mm = median(sqrt(sum(([cen, ones(n, 1)] * T - xyz) .^ 2, 2)));
    hemi(hemi == ' ' & xyz(:, 1) < 0) = 'L';
    hemi(hemi == ' ') = 'R';
end
if any(hemi == ' ')
    fprintf('Hemisphere of some regions unknown (no L/R in the labels, no MNI centroids): no spin test\n');
    G.hemi = hemi;
    return
end
G.hemi = hemi;
if ~isfield(C, 'Reg') || ~isfield(C.Reg, 'Sphere') || isempty(C.Reg.Sphere.Vertices)
    fprintf('The cortex has no sphere registration: no spin test\n');
    return
end
P = C.Reg.Sphere.Vertices;
G.sphere = zeros(n, 3);
for h = 'LR'
    m = find(hemi == h);
    X = P(unique([Sc(m).Vertices]), :);
    s = [2 * X, ones(size(X, 1), 1)] \ sum(X .^ 2, 2);     % least-squares sphere: centre s(1:3)
    for i = m'
        u = mean(P(Sc(i).Vertices, :), 1) - s(1:3)';
        G.sphere(i, :) = u / norm(u);
    end
end
end
