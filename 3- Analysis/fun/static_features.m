function F = static_features(D, k, nodes, families, part)
%STATIC_FEATURES The features of one static connectivity measure for every recording (long table).
%   F = static_features(D, k, nodes, families) with D from load_static, k
%   the measure (index into D.measures), nodes the network table rows of
%   the regions (load_network_table; empty: no network features) and
%   families a cell of S.static_families.
%   F = static_features(..., part): part 0 = D.conn (default), 1 or 2 =
%   the odd or even blocks of the whole recording (for the split-half
%   reliability, static_reliability).
%   F: table id, family, feature, state (0: a feature without states, for
%   the statistics of permu_ttest_permu_corr_BNS), value (an empty table
%   with these columns when no family gives a feature). Families:
%     global    Strength    mean node strength
%               ClustCoeff  mean weighted clustering coefficient
%               CPL         characteristic path length (lengths 1 / weight)
%               (gt_window_metrics: the graph metrics of the windows of
%               the states, BCT's strengths_und, clustering_coef_wu, distance_wei and
%               charpath, on the positive part of a signed measure such as AEC)
%               dwPLI: Strength only, from the signed weights (BCT's
%               strengths_und_sign: positive minus negative strength). dwPLI
%               estimates the squared wPLI without bias, so without coupling
%               it is centred on 0 and negative half the time: its positive
%               part is biased upwards by the noise, its
%               sum is not. ClustCoeff and CPL need non-negative weights and
%               are not computed for it
%     networks  within_<net>       mean of the edges inside a network
%               between_<a>_<b>    mean of the edges between two networks
%               (networks of network_assignment: Yeo-7 for DK-68; "Other"
%               left out; a network of one region has no within value;
%               signed measures keep their sign)
%     regions   strength_<region>  node strength of every region, the sum of
%                                  its positive weights (dwPLI: of all its
%                                  weights; the mean over the regions is
%                                  Strength)
%     power     logpower_all       mean over the regions of log10 band power
%               logpower_<net>     the same over each network's regions
%   A recording whose matrix has a NaN edge (a measure that could not be
%   computed for it, e.g. the symmetric orthogonalisation with fewer
%   independent signals than regions, or a split half too short) gets NaN
%   connectivity features, missing in the statistics.
if nargin < 5, part = 0; end
n = numel(D.roi_labels);
nr = numel(D.ids);
if part == 0
    W = reshape(D.conn(k, :, :, :), n, n, nr);
    P = D.power;
else
    W = reshape(D.halves(part, k, :, :, :), n, n, nr);
    P = reshape(D.power_halves(:, part, :), n, nr);
end
low = find(tril(true(n), -1));
E = reshape(W, n * n, nr);
E = E(low, :);                                     % edges x recordings
missing = any(isnan(E), 1);                        % not computed for this recording: NaN features
signed = strcmp(D.measures{k}, 'dwpli');           % centred on 0 without coupling: all its weights
rows = {};
add = @(fam, names, V) table(repmat(D.ids(:), numel(names), 1), repmat(string(fam), nr * numel(names), 1), ...
    repelem(string(names(:)), nr, 1), zeros(nr * numel(names), 1), reshape(V', [], 1), ...
    'VariableNames', {'id', 'family', 'feature', 'state', 'value'});   % V: features x recordings
[lab, nets] = network_assignment(nodes);
if ismember('global', families) && signed
    V = 2 * sum(E, 1) / n;                         % mean over the regions of strengths_und_sign's Spos - Sneg
    V(:, missing) = NaN;
    rows{end + 1} = add('global', "Strength", V);
elseif ismember('global', families)
    G = gt_window_metrics(W);
    V = [G.strength'; G.clustering'; G.path_length'];
    V(:, missing) = NaN;
    rows{end + 1} = add('global', ["Strength", "ClustCoeff", "CPL"], V);
end
if ismember('networks', families) && ~isempty(nets)
    names = strings(0, 1);  V = zeros(0, nr);
    for a = 1:numel(nets)
        for b = a:numel(nets)
            ia = lab == nets(a);  ib = lab == nets(b);
            if a == b
                if nnz(ia) < 2, continue; end
                M = false(n);  M(ia, ia) = true;
                names(end + 1, 1) = "within_" + nets(a); %#ok<AGROW>
            else
                M = false(n);  M(ia, ib) = true;  M(ib, ia) = true;
                names(end + 1, 1) = "between_" + nets(a) + "_" + nets(b); %#ok<AGROW>
            end
            M = M(low);                            % the block's edges, counted once
            V(end + 1, :) = mean(E(M, :), 1); %#ok<AGROW>
        end
    end
    rows{end + 1} = add('networks', names, V);
end
if ismember('regions', families)
    names = "strength_" + regexprep(string(D.roi_labels(:)), '[^A-Za-z0-9]+', '_');
    V = zeros(n, nr);
    for r = 1:nr
        if signed
            Wr = W(:, :, r);  Wr(1:n + 1:end) = 0;  % as strengths_und_sign (Spos - Sneg), without the diagonal
            V(:, r) = sum(Wr, 2);
        else
            V(:, r) = sum(max(W(:, :, r), 0), 2);  % as strengths_und
        end
    end
    V(:, missing) = NaN;
    rows{end + 1} = add('regions', names, V);
end
if ismember('power', families)
    L = log10(P);
    names = "logpower_all";  V = mean(L, 1);
    for a = 1:numel(nets)
        names(end + 1, 1) = "logpower_" + nets(a); %#ok<AGROW>
        V(end + 1, :) = mean(L(lab == nets(a), :), 1); %#ok<AGROW>
    end
    rows{end + 1} = add('power', names, V);
end
F = vertcat(rows{:});
if isempty(F)                                      % no family gave a feature: an empty table, its columns kept
    F = add('none', strings(0, 1), zeros(0, nr));
end
end
