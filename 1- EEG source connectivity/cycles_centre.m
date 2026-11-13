function c = cycles_centre(fmin, fmax, opts)
%CYCLES_CENTRE The frequency (Hz) whose cycles set a band's windows, epochs and filter.
%   c = cycles_centre(fmin, fmax, opts): the band centre, (fmin + fmax) / 2;
%   for a band anchored on the individual alpha frequency (IAF), the centre
%   of the band at an IAF of 10 Hz (opts.band_ref, its edges then: set by
%   run_pipeline and kept in each output's info.band_ref), the same for every
%   recording. Windows and epochs given in cycles (S.window_cycles, dwPLI's
%   epochs and windows) and the filter's transition (bandpass_pieces) then
%   have the same length in seconds whatever the recording's IAF, so the
%   bias of short-window estimates, which falls with window length times
%   bandwidth, does not follow the IAF (windows of 6 cycles of each
%   recording's own centre would be shorter at a higher IAF, and the window
%   mean of the phase measures would rise with the IAF).
if nargin > 2 && isfield(opts, 'band_ref') && ~isempty(opts.band_ref)
    c = 0.5 * (opts.band_ref(1) + opts.band_ref(2));
else
    c = 0.5 * (fmin + fmax);
end
end
