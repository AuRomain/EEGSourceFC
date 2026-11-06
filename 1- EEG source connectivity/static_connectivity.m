function [M, power, notes] = static_connectivity(Z, seg, measures, srate, fmin, fmax)
%STATIC_CONNECTIVITY Connectivity between region signals over a whole recording.
%   [M, power, notes] = static_connectivity(Z, seg, measures, srate, fmin, fmax)
%   Z         nROI x nSamples: analytic signals of the band-passed region
%             signals (complex, computed by run_pipeline before the filter
%             edges were trimmed), or the band-passed signals (real: Hilbert
%             transform of each continuous piece)
%   seg       nPieces x 2, first and last sample of each continuous piece
%             used; other samples are ignored (the halves of the split-half
%             reliability are pieces of the same Z: split_halves.m)
%   measures  cell of names (S.conn):
%     'plv'          phase locking value |<u_i conj(u_j)>|, u = z / |z|
%                    (Lachaux et al. 1999)
%     'wpli'         weighted phase lag index |sum Im(s_ij)| / sum |Im(s_ij)|,
%                    s_ij = z_i conj(z_j) at every sample (Vinck et al. 2011)
%     'pli'          phase lag index |<sign Im(s_ij)>| (Stam et al. 2007)
%     'ciplv'        corrected imaginary PLV |Im P| / sqrt(1 - (Re P)^2),
%                    P = <u_i conj(u_j)> (Bruna et al. 2018)
%     'dwpli'        debiased squared wPLI (Vinck et al. 2011): the recording
%                    is cut into blocks of tau / 2, tau = 1 / bandwidth, each
%                    block's sum of Im(s_ij) is one observation, and pairs of
%                    blocks closer than tau are left out as well as each
%                    block with itself (neighbouring blocks are correlated);
%                    estimates wPLI^2: about 0 without coupling, can be negative;
%                    NaN when the samples hold fewer than 4 blocks (a split
%                    half too short: notes; it stopped the recording)
%     'aec'          Pearson correlation of the amplitude envelopes |z|,
%                    signed (Brookes et al. 2011)
%     'aec_orth'     AEC after pairwise orthogonalisation: the zero-lag
%                    part of x is regressed out of y, and the two directions
%                    averaged (Brookes et al. 2012; Hipp et al. 2012)
%     'aec_symorth'  AEC after symmetric orthogonalisation of all regions
%                    together (Colclough et al. 2015); needs as many
%                    independent region signals as regions (at most the data
%                    rank of the channels minus 1 for the average reference;
%                    connectivity.symmetric_orthogonalisation): otherwise NaN
%                    (notes; it stopped the whole recording)
%     'plv_symorth'  PLV after symmetric orthogonalisation (the same condition)
%   M      nMeasures x nROI x nROI: symmetric, zero diagonal (NaN edges for a
%          measure that could not be computed)
%   power  nROI x 1: mean squared envelope |z|^2 over the samples used, the
%          band power of each region (units of the inverse solution, squared).
%          For a band-limited signal x with analytic signal z = x + i H(x),
%          |z|^2 = x^2 + H(x)^2 and H(x) has the power of x: the mean of
%          |z|^2 is twice the variance of the band-passed signal
%   notes  1 x nMeasures cell: '' or why the measure is NaN
%   wPLI, PLI, ciPLV, dwPLI and the orthogonalised measures ignore
%   zero-lag coupling, which is what source leakage produces; PLV and AEC
%   do not. get_all_connectivity_mats gives these as static_full_cmats,
%   with or without windows (S.dynamic); the windowed functions of
%   +connectivity give the same values over all the samples, except dwPLI,
%   whose blocks they fit into the windows (tests/test_static_estimators.m).
n = size(Z, 1);
if nargin < 2 || isempty(seg), seg = [1, size(Z, 2)]; end
measures = cellstr(measures);
in = false(1, size(Z, 2));
for k = 1:size(seg, 1), in(seg(k, 1):seg(k, 2)) = true; end
if isreal(Z)                                       % Hilbert transform of each piece
    Zr = Z;
    Z = complex(zeros(size(Zr)));
    for k = 1:size(seg, 1)
        idx = seg(k, 1):seg(k, 2);
        Z(:, idx) = hilbert(Zr(:, idx).').';
    end
    clear Zr
end
Zs = Z;
if ~all(in), Zs = Z(:, in); end                    % the samples used
N = size(Zs, 2);
[ii, jj] = find(tril(true(n), -1));
lin = sub2ind([n, n], ii, jj);
M = zeros(numel(measures), n, n);
power = mean(abs(Zs) .^ 2, 2);
notes = repmat({''}, 1, numel(measures));
U = [];
for k = 1:numel(measures)
    switch measures{k}
        case 'plv'
            if isempty(U), U = unit(Zs); end
            e = plv(U, lin, N);
        case 'wpli'
            e = wpli(Zs, ii, jj, lin);
        case 'pli'
            e = pli(Zs, ii, jj, N);
        case 'ciplv'
            if isempty(U), U = unit(Zs); end
            P = U * U' / N;                        % <u_i conj(u_j)>
            e = abs(imag(P(lin))) ./ sqrt(max(1 - real(P(lin)) .^ 2, eps));
        case 'dwpli'
            [e, notes{k}] = dwpli(Z, seg, ii, jj, round(srate / (fmax - fmin) / 2));
        case 'aec'
            e = aec(abs(Zs), lin);
        case 'aec_orth'
            e = aec_orth(Zs, lin);
        case 'aec_symorth'
            [Zo, notes{k}] = connectivity.symmetric_orthogonalisation(Zs, '''aec_orth''');
            e = nan(numel(lin), 1);
            if isempty(notes{k}), e = aec(abs(Zo), lin); end
        case 'plv_symorth'
            [Zo, notes{k}] = connectivity.symmetric_orthogonalisation(Zs, '''wpli'' or ''ciplv''');
            e = nan(numel(lin), 1);
            if isempty(notes{k}), e = plv(unit(Zo), lin, N); end
        otherwise
            error('static_connectivity:measure', 'Unknown connectivity measure %s', measures{k});
    end
    A = zeros(n);
    A(lin) = e;
    M(k, :, :) = A + A.';
end
end

function U = unit(Z)
% unit phasors; a sample without amplitude has no phase and adds nothing
U = Z ./ abs(Z);
U(abs(Z) == 0) = 0;
end

function e = plv(U, lin, N)
P = U * U';                                        % sum of u_i conj(u_j): one product
e = abs(P(lin)) / N;
end

function e = wpli(Z, ii, jj, lin)
% numerator: one product; denominator, sum of |Im(s_ij)|: edges in chunks
C = Z * Z';
num = imag(C(lin));
den = zeros(numel(ii), 1);
X = real(Z);  Y = imag(Z);
chunk = max(1, min(256, floor(2e8 / (8 * size(Z, 2)))));   % about 200 MB per array
for c0 = 1:chunk:numel(ii)
    c = c0:min(c0 + chunk - 1, numel(ii));
    den(c) = sum(abs(Y(ii(c), :) .* X(jj(c), :) - X(ii(c), :) .* Y(jj(c), :)), 2);
end
e = min(abs(num) ./ den, 1);                       % at most 1 (sums in different orders)
e(den == 0) = 0;                                   % no imaginary part at all
end

function e = pli(Z, ii, jj, N)
e = zeros(numel(ii), 1);
chunk = max(1, min(256, floor(2e8 / (8 * size(Z, 2)))));   % edges per pass, as wpli
for c0 = 1:chunk:numel(ii)
    c = c0:min(c0 + chunk - 1, numel(ii));
    e(c) = abs(sum(sign(imag(Z(ii(c), :) .* conj(Z(jj(c), :)))), 2)) / N;
end
end

function [e, note] = dwpli(Z, seg, ii, jj, B)
% blocks of B samples inside the pieces; X_b = sum of Im(s_ij) over block b
first = cell(1, size(seg, 1));
for k = 1:size(seg, 1), first{k} = seg(k, 1):B:(seg(k, 2) - B + 1); end
first = [first{:}];
note = '';
if numel(first) < 4                                % too short (a split half): NaN, not an error
    note = sprintf('dwPLI needs at least 4 blocks of %d samples; the data have %d', B, numel(first));
    e = nan(numel(ii), 1);
    return
end
e = zeros(numel(ii), 1);
chunk = max(1, min(256, floor(2e8 / (8 * size(Z, 2)))));
for c0 = 1:chunk:numel(ii)
    c = c0:min(c0 + chunk - 1, numel(ii));
    cs = cumsum([zeros(numel(c), 1), imag(Z(ii(c), :) .* conj(Z(jj(c), :)))], 2);
    X = cs(:, first + B) - cs(:, first);           % edges x blocks
    A = abs(X);
    num = sum(X, 2) .^ 2 - sum(X .^ 2, 2);         % every pair of different blocks ...
    den = sum(A, 2) .^ 2 - sum(A .^ 2, 2);
    for lag = 1:2                                  % ... except those less than tau apart
        num = num - 2 * sum(X(:, 1:end - lag) .* X(:, 1 + lag:end), 2);
        den = den - 2 * sum(A(:, 1:end - lag) .* A(:, 1 + lag:end), 2);
    end
    d = num ./ den;
    d(den <= 0) = 0;
    e(c) = d;
end
end

function e = aec(env, lin)
R = corr(env.');
e = R(lin);
end

function e = aec_orth(Z, lin)
% for every x, the leakage of x into every y (y - b x, b by least squares
% on the real signals) is removed, and the envelope of x correlated with
% the corrected envelopes; the analytic signal of y - b x is z_y - b z_x.
% The signals y in chunks of about 200 MB (all at once up to about 180,000
% samples)
n = size(Z, 1);
X = real(Z);
env = abs(Z);
F = zeros(n);
chunk = max(1, min(n, floor(2e8 / (16 * size(Z, 2)))));
for i = 1:n
    b = (X * X(i, :)') / (X(i, :) * X(i, :)');
    for c0 = 1:chunk:n
        c = c0:min(c0 + chunk - 1, n);
        F(i, c) = corr(env(i, :)', abs(Z(c, :) - b(c) * Z(i, :))');
    end
end
F = (F + F') / 2;                                  % both directions
e = F(lin);
end

