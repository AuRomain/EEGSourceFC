function g = roi_patterns(fm, radius)
%ROI_PATTERNS Scalp pattern of a coherent cortical patch around each region's seed.
%   g = sim.roi_patterns(fm, radius), fm from sim.forward_model: nChan x
%   nRoi, the mean leadfield of the vertices within RADIUS (m) of the
%   region's seed, with the sign of each vertex aligned to the seed's
%   normal (all dipoles pointing the same way relative to the cortical
%   surface). (From v_crosscheck.)
nroi = numel(fm.scouts);
g = zeros(size(fm.G, 1), nroi);
for r = 1:nroi
    sc = fm.scouts(r);
    v = sc.Vertices(:);
    v = v(vecnorm(fm.src_pos(v, :) - fm.src_pos(sc.Seed, :), 2, 2) < radius);
    fl = sign(fm.src_nrm(v, :) * fm.src_nrm(sc.Seed, :)');
    fl(fl == 0) = 1;
    g(:, r) = fm.G(:, v) * fl / numel(v);
end
end
