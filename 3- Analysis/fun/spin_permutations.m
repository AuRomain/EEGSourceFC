function perm = spin_permutations(G, n_rot, rs)
%SPIN_PERMUTATIONS Reorderings of the regions that keep their spatial layout (spin test).
%   perm = spin_permutations(G, n_rot, rs) with G from atlas_regions:
%   n_rot x n_regions; in rotation r, region i takes the label (network) of
%   region perm(r, i). Each rotation turns the left hemisphere's sphere
%   centroids (G.sphere) by a uniformly random rotation and the right
%   hemisphere's by its mirror image (Alexander-Bloch et al. 2018,
%   NeuroImage 178:540), so homologous regions move alike. The rotated
%   regions are then matched one-to-one to the regions of their hemisphere
%   with the least total distance (Hungarian algorithm, MATLAB matchpairs;
%   Vasa et al. 2018, Cereb Cortex 28:281): every label is used once, so
%   the networks keep their sizes, and a network made of neighbouring
%   regions stays made of neighbours. Random labels that keep neighbours
%   together are the null for connectivity that is high between close
%   regions (source leakage, distance), which random sets of regions are not.
%   The rotations are made as BrainSpace's spin_permutations makes them
%   (Vos de Wael et al. 2020, Commun Biol 3:103; github.com/MICA-MNI/BrainSpace:
%   QR of a 3 x 3 Gaussian matrix, signs fixed, a reflection turned into a
%   rotation, the right hemisphere's mirrored), which then takes each
%   rotated vertex's nearest vertex; the one-to-one assignment of regions
%   is that of Vasa et al.'s rotate_parcellation
%   (github.com/frantisekvasa/rotate_parcellation).
%   rs: a RandStream. Empty when G has no sphere.
perm = [];
if isempty(G.sphere), return; end
n = size(G.sphere, 1);
perm = zeros(n_rot, n);
mirror = diag([-1 1 1]);                       % x is left-right on the FreeSurfer sphere
L = find(G.hemi == 'L');
R = find(G.hemi == 'R');
for r = 1:n_rot
    [Q, U] = qr(randn(rs, 3));
    Q = Q * diag(sign(diag(U)));               % uniform on the orthogonal group
    if det(Q) < 0, Q(:, 1) = -Q(:, 1); end     % a rotation, not a reflection
    perm(r, L) = match(G.sphere(L, :), G.sphere(L, :) * Q', L);
    perm(r, R) = match(G.sphere(R, :), G.sphere(R, :) * (mirror * Q * mirror)', R);
end
end

function p = match(orig, rotated, idx)
% region idx(i) is at the position nearest (jointly) to rotated region idx(j)
D = sqrt(max(0, 2 - 2 * (orig * rotated')));  % chord distances on the unit sphere
M = matchpairs(D, 1e3);                          % every pair matched (distances are at most 2)
p = zeros(1, numel(idx));
p(M(:, 1)) = idx(M(:, 2));
end
