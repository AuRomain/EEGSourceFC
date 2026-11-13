function M = to_matrices(E, ii, jj, n)
%TO_MATRICES Edge values -> symmetric matrices with a zero diagonal.
%   M = connectivity.to_matrices(E, ii, jj, n) with E nEdges x nWindows
%   (edges ii(k) - jj(k)) returns nWindows x n x n (the cmats layout).
nw = size(E, 2);
M = zeros(nw, n, n, 'like', E);
lin = sub2ind([n, n], ii, jj);
for w = 1:nw
    A = zeros(n, 'like', E);
    A(lin) = E(:, w);
    M(w, :, :) = A + A.';
end
end
