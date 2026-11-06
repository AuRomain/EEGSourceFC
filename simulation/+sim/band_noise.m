function x = band_noise(n_sig, n, fs, band, rs)
%BAND_NOISE Independent band-limited Gaussian noise, unit RMS per row.
%   x = sim.band_noise(n_sig, n, fs, band, rs): white noise from the
%   RandStream rs, band-passed with Brainstorm's filter (transition 25% of
%   the lower edge, as stage 1), one second discarded at each edge.
w = randn(rs, n_sig, n + 2 * fs);
tran = min([max(0.25 * band(1), 2), band(1), fs / 2 - band(2)]);
x = bst_bandpass_hfilter(w, fs, band(1), band(2), 1, 0, [], tran, 'bst-hfilter-2019');
x = x(:, fs + 1:fs + n);
x = x ./ std(x, 0, 2);
end
