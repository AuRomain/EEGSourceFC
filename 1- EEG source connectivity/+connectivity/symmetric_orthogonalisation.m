function [Z, note, r] = symmetric_orthogonalisation(Z, alternative)
%SYMMETRIC_ORTHOGONALISATION Analytic signals of the closest mutually orthogonal signals (Colclough et al. 2015).
%   [Z, note, r] = connectivity.symmetric_orthogonalisation(Z, alternative)
%   Z            nSignals x nSamples analytic signals (complex) of
%                band-passed signals, over the samples to use
%   alternative  the measure to suggest instead, for the note
%   Z            the analytic signals of the orthogonalised signals, or []
%                when the signals are not independent enough (note says why)
%   note         '' or the reason the orthogonalisation was not done
%   r            the rank of the real signals
%   The closest orthogonal signals are a linear map of the real signals
%   (L = X * W, ROInets.closest_orthogonal_matrix) and the Hilbert
%   transform is linear, so their analytic signals are Z * W: those
%   computed before the filter edges were trimmed.
%   It needs as many independent signals as regions. The rank is counted
%   as run_pipeline's data_rank counts the channels': eigenvalues of the
%   covariance of the real signals above 1e-12 of the largest (singular
%   values above 1e-6 of the largest), so that the rounding of data stored
%   as 32-bit floats does not count as signal. MATLAB's rank() (tolerance
%   at machine precision) took that rounding for independent signals, and
%   the orthogonalisation then amplified it. Regions are mixtures of at most
%   the data rank of the channels (minus 1 for the average reference,
%   unless the data were average-referenced already): of the OpenNeuro
%   recordings (data ranks 9-49 of 60 channels, median 33), none can be
%   orthogonalised with 68 regions, and with 34 regions 73-81 of 138 (data
%   rank below 34, or 34) cannot. The callers then
%   give NaN for the measure (for that recording only: run_pipeline names
%   it, stage1_summary.csv lists it in measures_nan), where an error used
%   to stop the whole recording.
note = '';
X = double(real(Z));
Xc = X - mean(X, 2);
e = eig(Xc * Xc' / max(1, size(Xc, 2)), 'vector');
r = nnz(e > 1e-12 * max(e));
n = size(X, 1);
if r < n
    note = sprintf(['symmetric orthogonalisation of %d regions needs %d independent signals; the region ' ...
        'signals have rank %d (at most the data rank of the channels minus 1, for the average reference: ' ...
        'data_rank in stage1_summary.csv). Use an atlas with fewer regions than that, or %s'], n, n, r, alternative);
    Z = [];
    return
end
L = ROInets.closest_orthogonal_matrix(real(Z)')';
W = real(Z)' \ L';                                 % L' = X' * W (exact: L is in the span of X)
Z = (Z.' * W).';
end
