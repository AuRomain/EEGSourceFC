function x = peak_noise(n_sig, n, fs, f0, sd, rs)
%PEAK_NOISE Independent noise with a Gaussian spectral peak, unit RMS per row.
%   x = sim.peak_noise(n_sig, n, fs, f0, sd, rs): random phases and a
%   spectrum exp(-(f - f0)^2 / (2 sd^2)) (amplitude), so the power spectrum
%   peaks at f0 Hz with an SD of sd / sqrt(2) Hz: an alpha rhythm whose peak
%   frequency is known (the simulated individual alpha frequency).
nf = 2 ^ nextpow2(n);
f = (0:nf / 2)' / nf * fs;
amp = exp(-(f - f0) .^ 2 / (2 * sd ^ 2));
X = (randn(rs, nf / 2 + 1, n_sig) + 1i * randn(rs, nf / 2 + 1, n_sig)) .* amp;
X(1, :) = 0;
X(end, :) = real(X(end, :));
x = ifft([X; conj(X(end - 1:-1:2, :))], 'symmetric');
x = x(1:n, :).';
x = x ./ std(x, 0, 2);
end
