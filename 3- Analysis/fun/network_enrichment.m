function [P, stat, G] = network_enrichment(W, X, grp, nodes, Gatlas, opt)
%NETWORK_ENRICHMENT Where a group difference of the edges concentrates: networks against random and spin nulls.
%   [P, stat, G] = network_enrichment(W, X, grp, nodes, Gatlas, opt)
%   W       nROI x nROI x nRecordings, one measure, the recordings of the
%           comparison
%   X, grp  design [1, A, covariates] and group A (logical) of these
%           recordings, as compute_NBS
%   nodes   network table rows of the regions (network_assignment)
%   Gatlas  atlas_regions (Gatlas.sphere for the spin test; may be empty)
%   opt     n_perm, seed
%   Every edge's group statistic A - B is the Welch-type G of the model
%   (welch_G; Welch t without covariates), G (nROI x nROI), as z on its own
%   Aspin-Welch degrees of freedom (t_to_z, as PALM; 5 November 2026), as
%   in the NBS: stat (nROI x nROI), on which the blocks are scored. A block is a network (the edges inside it) or a pair of
%   networks (the edges between them); its score is the mean statistic of
%   its edges. The score is compared with two nulls, two-sided:
%     random  the network labels shuffled over the regions (n_perm times)
%     spin    the labels rotated on the sphere (spin_permutations, n_perm
%             times), which keeps the networks' sizes AND spatial layout
%   z = (score - null mean) / null SD; p = share of the null at least as
%   far from its mean; q = Benjamini-Hochberg over the blocks, per null.
%   Neighbouring regions share signal (source leakage) and differ together,
%   so a compact network can stand out against random regions for that
%   reason alone: the spin null keeps this and is the one to read. Both
%   nulls keep the edge statistics fixed: they say where a group
%   difference lies, beyond regions laid out alike, not whether a block
%   differs between the groups (the family 'networks' tests that, with
%   permutations of the recordings).
P = table();
[lab, nets] = network_assignment(nodes);
n = size(W, 1);
low = find(tril(true(n), -1));
Y = reshape(W, n * n, []);
Y = Y(low, :)';                                    % recordings x edges
c = zeros(size(X, 2), 1);  c(2) = 1;
[g, ~, df] = welch_G(Y, X, c, grp);
zg = t_to_z(g, df);
g(isnan(g)) = 0;  zg(isnan(zg)) = 0;
stat = zeros(n);  G = stat;
stat(low) = zg;  G(low) = g;
stat = stat + stat';  G = G + G';
if numel(nets) < 2, return; end
rs = RandStream('mrg32k3a', 'Seed', opt.seed);
perm_rand = zeros(opt.n_perm, n);
for r = 1:opt.n_perm, perm_rand(r, :) = randperm(rs, n); end
perm_spin = spin_permutations(Gatlas, opt.n_perm, rs);
[bi, bj] = find(triu(true(numel(nets))));
E = double(lab == nets');
score = block_scores(stat, E, bi, bj);
[z_r, p_r] = null_compare(score, stat, lab, nets, perm_rand, bi, bj);
if isempty(perm_spin)
    z_s = nan(size(score));  p_s = z_s;
else
    [z_s, p_s] = null_compare(score, stat, lab, nets, perm_spin, bi, bj);
end
kind = repmat("between", numel(bi), 1);
kind(bi == bj) = "within";
q = sum(E, 1)';
edges = q(bi) .* q(bj);                            % edges between two networks
edges(bi == bj) = q(bi(bi == bj)) .* (q(bi(bi == bj)) - 1) / 2;   % inside one
P = table(nets(bi), nets(bj), kind, edges, score, z_r, p_r, z_s, p_s, ...
    'VariableNames', {'network_a', 'network_b', 'type', 'n_edges', 'mean_stat', 'z_random', 'p_random', ...
    'z_spin', 'p_spin'});
P = P(P.n_edges > 0, :);                           % a network of one region has no edge inside
P = addvars(P, fdr_bh(P.p_random), 'After', 'p_random', 'NewVariableNames', 'q_random');
P = addvars(P, fdr_bh(P.p_spin), 'After', 'p_spin', 'NewVariableNames', 'q_spin');
end

function s = block_scores(T, E, bi, bj)
% mean of T over the edges of each block (T symmetric, zero diagonal)
n = size(T, 1);
cnt = E' * (1 - eye(n)) * E;
S = (E' * T * E) ./ cnt;
s = S(sub2ind(size(S), bi, bj));
end

function [z, p] = null_compare(score, T, lab, nets, perms, bi, bj)
nperm = size(perms, 1);
null = zeros(numel(score), nperm);
for r = 1:nperm
    null(:, r) = block_scores(T, double(lab(perms(r, :)) == nets'), bi, bj);
end
mu = mean(null, 2, 'omitnan');
sd = std(null, 0, 2, 'omitnan');
z = (score - mu) ./ sd;
p = (1 + sum(abs(null - mu) >= abs(score - mu), 2)) / (nperm + 1);
end
