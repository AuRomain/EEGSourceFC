function C = cluster_perm_test(t_obs, T, df, A, alpha_c)
%CLUSTER_PERM_TEST Cluster-mass permutation test over neighbouring channels or regions (FieldTrip's clusterstat).
%   C = cluster_perm_test(t_obs, T, df, A, alpha_c)
%   t_obs    1 x K statistics of the data: the z of perm_ancova_maxT (its
%            Welch t or Welch-type G converted with their own degrees of
%            freedom, as PALM; df Inf) or of perm_partial_corr_maxT (its
%            slope over the HC3 standard error, as z on its degrees of
%            freedom; df Inf); NaN = not tested
%   T        n_perm x K the same statistics in each permutation: the joint
%            permutations of perm_ancova_maxT or perm_partial_corr_maxT (the
%            same permutation of the subjects in every column)
%   df       1 x K degrees of freedom of each column (subjects minus the
%            parameters of its model; Inf for z, the threshold then the
%            normal's, as PALM's cluster-forming threshold on the z scale)
%   A        K x K adjacency (true for neighbours; made symmetric, its
%            diagonal ignored)
%   alpha_c  cluster-forming threshold: a column is supra-threshold when
%            its statistic passes the two-sided p < alpha_c of Student's t
%            with its df (0.05: |t| >= 2.0 for 40 subjects, |z| >= 1.96)
%   The test of FieldTrip's ft_statistics_montecarlo with cfg.correctm =
%   'cluster', clusterstatistic 'maxsum', clusterthreshold 'parametric'
%   (one critical value per column: clustercritval = [tinv(alpha_c/2, df)
%   tinv(1 - alpha_c/2, df)]), tail 0, correcttail 'prob' (Maris &
%   Oostenveld 2007, J Neurosci Methods 164:177), adapted from FieldTrip's
%   private/clusterstat.m and findcluster.m (commit 0beabe6, the tested
%   version; GPL-3) to take the statistics of all permutations at once:
%   - supra-threshold: t >= the upper critical value (positive) or <= the
%     lower (negative) (clusterstat lines 191-196);
%   - clusters: the connected sets, in A, of the supra-threshold columns of
%     one sign (findcluster with combineClusters); a cluster's statistic is
%     the sum of its t (maxsum);
%   - null distributions: in each permutation the largest positive
%     cluster sum and the smallest (most negative) negative one, 0 when
%     there is none (clusterstat lines 273-370);
%   - p of a positive cluster: (1 + the permutations whose largest sum is
%     above its sum) / (n_perm + 1); of a negative cluster, below
%     (clusterstat lines 437-444 and 517-524);
%   - two-sided: p doubled, at most 1 (ft_statistics_montecarlo lines
%     406-419, correcttail 'prob').
%   So cluster_p controls the family-wise error over the clusters of the
%   map, both signs. tests/test_psd checks the clusters, sums and p against
%   FieldTrip's own clusterstat. A cluster says that the effect is present
%   somewhere in it, not at each of its channels or regions.
%   C (table, one row per column): cluster (its number, 1 = the largest
%   |sum|; 0 = in no cluster), cluster_sign (+1, -1; 0), cluster_size
%   (columns in it), cluster_mass (the sum) and cluster_p (two-sided; NaN
%   outside clusters).
t_obs = t_obs(:)';
K = numel(t_obs);
if size(T, 2) ~= K || numel(df) ~= K || ~isequal(size(A), [K K])
    error('cluster_perm_test:size', 't_obs (1 x %d), T (n_perm x K), df (1 x K) and A (K x K) must agree', K);
end
n_perm = size(T, 1);
df = double(df(:)');
thr = inf(1, K);                                  % the upper critical value; the lower is -thr
ok = df > 0;                                      % Inf: tinv gives the normal quantile
thr(ok) = tinv(1 - alpha_c / 2, df(ok));
A = logical(A);
A = A | A';
A(1:K + 1:end) = false;
deg = sum(A, 2);
Nb = repmat(K + 1, K, max([deg; 1]));             % neighbour lists, padded with an empty column K + 1
for k = 1:K
    j = find(A(k, :));
    Nb(k, 1:numel(j)) = j;
end
% the data's clusters, positive and negative numbered apart
lab = zeros(1, K);  sgn = zeros(1, K);
for s = [1 -1]
    L = components(s * t_obs >= thr, Nb);
    new = L > 0;
    lab(new) = L(new) + (s < 0) * K;
    sgn(new) = s;
end
ids = unique(lab(lab > 0));
mass = zeros(1, numel(ids));  csgn = zeros(1, numel(ids));
for i = 1:numel(ids)
    mass(i) = sum(t_obs(lab == ids(i)));
    csgn(i) = sgn(find(lab == ids(i), 1));
end
% the null distributions: the largest positive and the smallest negative sum of each permutation
pos = zeros(n_perm, 1);  neg = zeros(n_perm, 1);
for chunk = 1:2000:n_perm                         % in pieces, to bound the memory
    rows = chunk:min(chunk + 1999, n_perm);
    [pos(rows), neg(rows)] = extreme_sums(T(rows, :), thr, Nb);
end
[~, order] = sort(abs(mass), 'descend');
C = table(zeros(K, 1), zeros(K, 1), zeros(K, 1), nan(K, 1), nan(K, 1), ...
    'VariableNames', {'cluster', 'cluster_sign', 'cluster_size', 'cluster_mass', 'cluster_p'});
for r = 1:numel(order)
    i = order(r);
    in = lab == ids(i);
    if csgn(i) > 0
        p = (sum(pos > mass(i)) + 1) / (n_perm + 1);
    else
        p = (sum(neg < mass(i)) + 1) / (n_perm + 1);
    end
    C.cluster(in) = r;
    C.cluster_sign(in) = csgn(i);
    C.cluster_size(in) = nnz(in);
    C.cluster_mass(in) = mass(i);
    C.cluster_p(in) = min(1, 2 * p);
end
end

function L = components(M, Nb)
% connected components of the true entries of each row of M (n x K) in the
% graph of the neighbour lists Nb (K x d, padded with K + 1): each true
% entry gets the smallest column index of its component, false entries 0
% (what findcluster and combineClusters give for channels alone, for all
% rows at once: minimum-label propagation)
[n, K] = size(M);
L = repmat(1:K, n, 1);
L(~M) = Inf;
L = [L, inf(n, 1)];
while true
    Lnew = L(:, 1:K);
    for d = 1:size(Nb, 2)
        Lnew = min(Lnew, L(:, Nb(:, d)));
    end
    Lnew(~M) = Inf;
    if isequal(Lnew, L(:, 1:K)), break; end
    L(:, 1:K) = Lnew;
end
L = L(:, 1:K);
L(isinf(L)) = 0;
end

function [pos, neg] = extreme_sums(T, thr, Nb)
% the largest positive and the smallest negative cluster sum of each row
% of T (0 when a row has no cluster of that sign)
[n, K] = size(T);
pos = zeros(n, 1);  neg = zeros(n, 1);
for s = [1 -1]
    M = s * T >= thr;                             % NaN (not tested) is never in a cluster
    if ~any(M(:)), continue; end
    L = components(M, Nb);
    [r, c] = find(L);
    lin = sub2ind([n, K], r, c);
    ms = accumarray([r, L(lin)], s * T(lin), [n, K]);   % the |sum| of each cluster, at its label
    if s > 0
        pos = max(ms, [], 2);
    else
        neg = -max(ms, [], 2);
    end
end
end
