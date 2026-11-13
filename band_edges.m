function [fmin, fmax, anchored] = band_edges(band, iaf, opt)
%BAND_EDGES Edges (Hz) of a frequency band: fixed, or anchored on the individual alpha frequency.
%   [fmin, fmax, anchored] = band_edges(band, iaf, opt)
%   band   an element of S.bands; its fmin and fmax are each either a
%          number (Hz) or text relative to the recording's individual alpha
%          frequency (IAF):
%            'iaf'           the IAF
%            'iaf+2', 'iaf-6'  the IAF plus or minus a number of Hz
%            'iaf+bins', 'iaf-bins'  plus or minus S.iaf.bins bins of
%                            S.iaf.bin_hz Hz (the individual alpha band)
%   iaf    the recording's IAF (Hz); NaN when it could not be estimated
%          (anchored edges are then NaN)
%   opt    S.iaf (bins, bin_hz)
%   anchored  true when an edge refers to the IAF (individual bands are
%          optional: S.individual_bands)
%   The edges are not checked here: check_settings checks them at a typical
%   IAF (10 Hz), run_pipeline for each recording (an anchored band narrower
%   than 1 Hz for a recording's IAF is left out for that recording).
%   Example (docs/individual_alpha_frequency.md): with 2 bins of 1 Hz,
%   {'iaf-bins', 'iaf+bins'} gives 7-11 Hz for an IAF of 9 Hz.
[fmin, a1] = one_edge(band.fmin, iaf, opt, band.name);
[fmax, a2] = one_edge(band.fmax, iaf, opt, band.name);
anchored = a1 || a2;
end

function [f, anchored] = one_edge(e, iaf, opt, name)
anchored = false;
if isnumeric(e) && isscalar(e)
    f = e;
    return
end
if ~(ischar(e) || isstring(e))
    error('band_edges:edge', 'Band %s: an edge must be a number (Hz) or text such as ''iaf-bins''', name);
end
tok = regexp(lower(strrep(char(e), ' ', '')), '^iaf([+-](?:bins|[\d.]+))?$', 'tokens', 'once');
offset = 0;
if ~isempty(tok) && ~isempty(tok{1})                % '+bins', '-6', ...
    if strcmp(tok{1}(2:end), 'bins'), offset = opt.bins * opt.bin_hz; else, offset = str2double(tok{1}(2:end)); end
    if tok{1}(1) == '-', offset = -offset; end
end
if isempty(tok) || isnan(offset)
    error('band_edges:edge', ['Band %s: edge ''%s'' not understood (a number, ''iaf'', ''iaf+2'', ' ...
        '''iaf-6'', ''iaf+bins'' or ''iaf-bins'')'], name, char(e));
end
anchored = true;
f = iaf + offset;
end
