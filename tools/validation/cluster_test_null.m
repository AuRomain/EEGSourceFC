function cluster_test_null(n_sim, S)
%CLUSTER_TEST_NULL Family-wise error of the channel cluster test (cluster_perm_test) under the null.
%   cluster_test_null()          the documented size
%   cluster_test_null(n_sim)     fewer datasets, for a quick check
%   cluster_test_null(n_sim, S)  with settings S (pipeline_settings), e.g.
%                                S.toolbox_paths when the toolboxes are not
%                                in the Toolboxes folder next to the repository
%
%   Reproduces the paragraph "The cluster test with FieldTrip's rules,
%   under the null ..." after the table of the section "Where the methods
%   come from" of docs/validation.md (also quoted as "with FieldTrip's
%   rules: family-wise error 4.2%" in "Statistics of stage 4 and of the PSD
%   analysis"), with the current perm_ancova_maxT, cluster_perm_test and
%   psd_adjacency of this repository.
%   Data: 54 channels of the default cap (S.montages{1}, 'ASA_10-05_343'),
%   spatially correlated (Gaussian in the distance between the
%   electrodes, scale half the median distance), FieldTrip's
%   triangulation neighbours (psd_adjacency, ft_prepare_neighbours);
%   25 + 25 recordings, no difference; Welch t with 1000 permutations
%   (perm_ancova_maxT), clusters at the two-sided 0.05 threshold.
%   Needs FieldTrip (and the toolboxes setup_toolboxes requires for S).
%   Default (that of the documented paragraph): n_sim = 2000 datasets.
%   Seeds: the data from RandStream('mrg32k3a', 'Seed', 11), in the order
%   of the datasets; perm_ancova_maxT's seed is the dataset's number.
%   Prints the number of channels and of datasets, the share of datasets
%   with a cluster p <= 0.05 (family-wise error) with its 95%
%   Clopper-Pearson interval, and the share of datasets in which at least
%   one cluster formed.
%   Run time with the defaults: about 3 minutes (one core).
%   The cluster masses are sums of z (each permutation's Welch t as z on
%   its own degrees of freedom); with the defaults: 83 of 2,000 (4.2%, 95%
%   CI 3.3-5.1%), at least one cluster in 76%, as documented.
if nargin < 1 || isempty(n_sim), n_sim = 2000; end
validation_setup();
if nargin < 2 || isempty(S), S = pipeline_settings(); end
S.parallel = false;
S = setup_toolboxes(S);
el = get_ftElec(S.montages{1});
chans = {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8','CP5','CP1', ...
       'CP2','CP6','P7','P3','Pz','P4','P8','O1','Oz','O2','AF3','AF4','FC3','FC4','C1','C2','CP3','CP4', ...
       'P1','P2','PO3','PO4','F1','F2','FT7','FT8','TP7','TP8','PO7','PO8','F5','F6','C5','C6','P5','P6'};
chans = chans(ismember(lower(chans), lower(el.label)));
A = psd_adjacency('sensor', chans, S);
[~, row] = ismember(lower(chans), lower(el.label));
pos = el.chanpos(row, :);
K = numel(chans);
Dm = squareform(pdist(pos));
Cs = exp(-(Dm / (0.5 * median(Dm(:)))) .^ 2);           % spatial correlation
Lc = chol(Cs + 1e-6 * eye(K), 'lower');
rs = RandStream('mrg32k3a', 'Seed', 11);
g = [false(25, 1); true(25, 1)];
any05 = false(n_sim, 1);  n_cl = zeros(n_sim, 1);
for i = 1:n_sim
    Y = randn(rs, 50, K) * Lc';
    [~, D] = perm_ancova_maxT(Y, g, [], struct('n_perm', 1000, 'seed', i, 'n_boot', 0));
    C = cluster_perm_test(D.t_obs, D.T, D.df, A, 0.05);
    any05(i) = any(C.cluster_p <= 0.05);
    n_cl(i) = max(C.cluster);
end
[~, ci] = binofit(sum(any05), n_sim);
fprintf(['Cluster test under the null, %d channels, %d datasets: family-wise error (a cluster p <= 0.05) %.3f ' ...
    '(%d of %d; 95%% CI %.3f-%.3f); datasets with at least one cluster %.2f\n'], ...
    K, n_sim, mean(any05), sum(any05), n_sim, ci(1), ci(2), mean(n_cl > 0));
end
