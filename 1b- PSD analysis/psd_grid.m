function f = psd_grid(o)
%PSD_GRID The frequencies of the PSD step: S.psd.freq_range in steps of 1 / S.psd.segment Hz.
%   f = psd_grid(S.psd) (1-45 Hz in 0.25 Hz steps by default), the grid
%   every recording's spectra are given on (psd_spectra).
f = o.freq_range(1):1 / o.segment:o.freq_range(2) + 1e-9;
f = round(f * 1e6) / 1e6;
end
