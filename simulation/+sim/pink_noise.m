function x = pink_noise(n_sig, n, fs, rs)
%PINK_NOISE Independent 1/f noise (above 0.5 Hz), unit RMS per row.
%   (From v_crosscheck.)
nf = 2 ^ nextpow2(n);
f = (0:nf / 2)' / nf * fs;
amp = zeros(size(f));
amp(f >= 0.5) = 1 ./ sqrt(f(f >= 0.5));
X = (randn(rs, nf / 2 + 1, n_sig) + 1i * randn(rs, nf / 2 + 1, n_sig)) .* amp;
X(1, :) = 0;
X(end, :) = real(X(end, :));
x = ifft([X; conj(X(end - 1:-1:2, :))], 'symmetric');
x = x(1:n, :).';
x = x ./ std(x, 0, 2);
end
