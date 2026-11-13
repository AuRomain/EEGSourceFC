function [Pm, Pv] = psd_vertex_spectra(W, C, vert, free)
%PSD_VERTEX_SPECTRA The PSD of every source from the channels' cross-spectra, and each region's mean of them.
%   [Pm, Pv] = psd_vertex_spectra(W, C, vert, free)
%   W     the inverse filters (sources x channels; with free orientation
%         three rows per source, its x, y and z: psd_source_signals)
%   C     the channels' cross-spectra (channels x channels x frequencies,
%         real part: psd_cross_spectra)
%   vert  each region's sources (1 x regions cell of vertex indices)
%   free  true: three rows per source, whose PSDs are added
%   Pv    sources x frequencies: w' C(f) w, the PSD of the source's signal
%         w' X (psd_welch of it, without computing it)
%   Pm    regions x frequencies: the mean of Pv over each region's sources
%   The PSD of a source is that of its current density (LORETA-KEY's
%   spectra; FieldTrip's source power from a cross-spectral density; the
%   PSD of every source in Brainstorm); with free orientation the sum over
%   its three orientations (the squared norm of the current density
%   vector).
nf = size(C, 3);
Pv = zeros(size(W, 1), nf);
for k = 1:nf
    Pv(:, k) = sum((W * C(:, :, k)) .* W, 2);
end
if free, Pv = reshape(sum(reshape(Pv, 3, [], nf), 1), [], nf); end
Pm = cell2mat(cellfun(@(v) mean(Pv(v, :), 1), vert(:), 'UniformOutput', false));
end
