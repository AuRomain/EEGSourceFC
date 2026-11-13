function [M, power, notes] = static_connectivity(Z, seg, measures, srate, fmin, fmax, opts)
%STATIC_CONNECTIVITY Connectivity between region signals over a whole recording.
%   [M, power, notes] = static_connectivity(Z, seg, measures, srate, fmin, fmax)
%   [M, power, notes] = static_connectivity(..., opts): opts.dwpli_epoch_cycles
%   or opts.dwpli_epoch_length, dwPLI's epochs (connectivity.dwpli_epoch;
%   default 6 cycles of the band centre): S, as get_all_connectivity_mats
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
%     'dwpli'        debiased wPLI (Vinck et al. 2011), computed by FieldTrip
%                    with epochs as trials: each piece cut into epochs
%                    (connectivity.dwpli_epoch: 6 cycles of the band centre
%                    by default), their Fourier spectra (ft_freqanalysis,
%                    mtmfft, Hanning taper), ft_connectivity_wpli with
%                    'debias' at each Fourier frequency of the band, and the
%                    mean over these frequencies (connectivity.dwpli_epochs);
%                    estimates wPLI^2: about 0 without coupling, can be
%                    negative; NaN with fewer than 10 epochs (a split half
%                    too short: notes)
%     'aec'          Pearson correlation of the amplitude envelopes |z|,
%                    signed (Brookes et al. 2011)
%     'aec_orth'     AEC after pairwise orthogonalisation: the zero-lag
%                    part of x is regressed out of y over the whole
%                    recording, and the two directions averaged (the
%                    regression form of Brookes et al. 2012 and O'Neill et
%                    al. 2017; not Hipp et al. 2012, who orthogonalise the
%                    analytic signals sample by sample)
%     'aec_symorth'  AEC after symmetric orthogonalisation of all regions
%                    together (Colclough et al. 2015); needs as many
%                    independent region signals as regions (at most the data
%                    rank of the channels minus 1 for the average reference;
%                    connectivity.symmetric_orthogonalisation): otherwise NaN
%                    (notes)
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
%   +connectivity give the same values over all the samples, dwPLI over the
%   same epochs (tests/test_static_estimators.m).
if nargin < 7, opts = struct(); end
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
        case 'dwpli'                                % epochs of every piece, from its start
            ep = connectivity.dwpli_epoch(srate, fmin, fmax, opts);
            [d, ~, notes{k}] = connectivity.dwpli_epochs(Z, srate, ...
                connectivity.window_starts(seg, ep.samples, ep.samples), ep);
            e = d(1, lin)';
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

