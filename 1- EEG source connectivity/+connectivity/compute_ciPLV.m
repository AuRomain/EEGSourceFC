function ciplv = compute_ciPLV(data)
%COMPUTE_CIPLV Corrected imaginary phase-locking value (Bruna et al. 2018).
%   ciplv = connectivity.compute_ciPLV(data) with data nSignals x nSamples
%   (x nTrials), analytic (complex) or real (phase taken by Hilbert).
%   ciPLV = |Im(mean(e^{i dphi}))| / sqrt(1 - Re(mean(e^{i dphi}))^2): the
%   imaginary part of the PLV, insensitive to zero-lag (volume conduction)
%   coupling. Zero diagonal.
%   Reference implementation (tests); ciplv_sliding_window computes the
%   same values in windows.
if isreal(data), data = hilbert(data.').'; end
[nc, ns, nt] = size(data);
ndat = data ./ abs(data);
ciplv = zeros(nc, nc, nt);
for t = 1:nt
    P = ndat(:, :, t) * ndat(:, :, t)' / ns;
    C = abs(imag(P)) ./ sqrt(1 - real(P) .^ 2);
    C(1:nc + 1:end) = 0;
    ciplv(:, :, t) = C;
end
end
