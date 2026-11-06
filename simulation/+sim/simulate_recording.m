function [X, roi] = simulate_recording(lab, masks, lags, g, G, fs, p, rs, alpha_mask, iaf)
%SIMULATE_RECORDING Sensor data of one simulated recording (see sim.default_params).
%   [X, roi] = sim.simulate_recording(lab, masks, lags, g, G, fs, p, rs)
%   lab    nSamples x 1 state labels (sim.state_sequence)
%   masks  nRoi x K, ROIs of each state's network; lags nRoi x K phase lags
%   g      nChan x nRoi ROI patterns (sim.roi_patterns); G full leadfield
%   alpha_mask, iaf  optional: regions carrying an alpha rhythm around the
%          subject's individual alpha frequency iaf (p.alpha_rms > 0)
%   Returns X (nChan x nSamples, microvolts, average reference) and the ROI
%   source signals ROI (nRoi x nSamples), for checks at the source level.
%   (From v_crosscheck.)
n = numel(lab);
[nroi, K] = size(masks);
w = round(p.crossfade_s * fs);
c = zeros(K, n);
for k = 1:K
    c(k, :) = movmean(double(lab' == k), max(w, 1));      % 0..1, smooth on/off of state k
end
osc = sim.band_noise(nroi, n, fs, p.band, rs);
drv = hilbert(sim.band_noise(K, n, fs, p.band, rs).').';
jit = randn(rs, nroi, K) * p.lag_jitter_rad;
beta = osc;
for k = 1:K
    m = find(masks(:, k));
    shared = real(drv(k, :) .* exp(1i * (lags(m, k) + jit(m, k))));
    ck = p.coupling * c(k, :);
    gain = 1 + (p.state_gain - 1) * c(k, :);             % network amplitude while its state is on
    beta(m, :) = gain .* (sqrt(1 - ck .^ 2) .* osc(m, :) + ck .* shared);
end
roi = p.osc_rms * beta + p.pink_rms * sim.pink_noise(nroi, n, fs, rs);
if nargin >= 10 && p.alpha_rms > 0                       % alpha rhythm at the subject's IAF
    roi(alpha_mask, :) = roi(alpha_mask, :) + p.alpha_rms * ...
        sim.peak_noise(nnz(alpha_mask), n, fs, iaf, p.alpha_sd_hz, rs);
end
X = g * roi;
src = randperm(rs, size(G, 2), p.brain_noise_sources);
X = X + p.brain_noise_rms * G(:, src) * sim.pink_noise(numel(src), n, fs, rs) / sqrt(numel(src) / nroi);
X = X + p.sensor_noise * std(X(:)) * randn(rs, size(X));
X = X - mean(X, 1);
s = std(X(:));
X = X / s * p.scale_uv;
end
