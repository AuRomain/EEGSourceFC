function [chosen, how] = cap_channels(wanted, labels, elec, max_mm)
%CAP_CHANNELS The channels of a cap that stand for a set of 10-05 names.
%   [chosen, how] = cap_channels(wanted, labels, elec, max_mm)
%   wanted   10-05 names of a set of channels: the posterior channels of the
%            IAF (S.iaf.channels), the frontal, frontopolar or temporal
%            channels of the spectral profile and of the window measures
%   labels   the channels to choose from (a recording's, or a cap's)
%   elec     the cap's electrodes (get_ftElec: label, chanpos in metres, in
%            the ICBM152 coordinates of the 10-05 cap); [] = by name only
%   max_mm   largest distance from a name's 10-05 position (default 20 mm)
%   chosen   the names of WANTED among LABELS (any case; in the spelling of
%            LABELS, in the order of WANTED). Only when LABELS
%            has none of them (a cap with other names: an ANT equidistant
%            cap, Z1 ... RR14; an EGI HydroCel cap, E1 ... E256), each name
%            is matched by position: the channel of LABELS nearest to its
%            position on Brainstorm's ASA 10-05 cap, within max_mm, each
%            channel once (two names at the same position, T3 and T7, give
%            one channel). A cap whose names are 10-20 names is matched
%            by name only.
%   how      '' (by name, or none found), or the matches by position, e.g.
%            'Z9 for Oz (5 mm), L9 for O1 (12 mm)'
%   A cap with some of the names (e.g. Cz and Oz among E1 ... E128) keeps
%   only those, by name: the others are not matched by position. Check
%   the channels a run used (its log) on such caps, and give S.iaf.channels
%   (or the setting of the measure) in the cap's own names.
if nargin < 3, elec = []; end
if nargin < 4 || isempty(max_mm), max_mm = 20; end
wanted = cellstr(wanted);  wanted = wanted(:)';
labels = cellstr(labels);  labels = labels(:)';
how = '';
[found, at] = ismember(lower(wanted), lower(labels));
chosen = labels(at(found));
if ~isempty(chosen) || isempty(elec) || isempty(labels), return; end
T = template_cap();
[in_template, it] = ismember(lower(wanted), lower(T.label));
[on_cap, ic] = ismember(lower(labels), lower(elec.label));
cand = labels(on_cap);
P = elec.chanpos(ic(on_cap), :) * unit_to_mm(elec);
parts = {};
for k = find(in_template)
    [d, j] = min(sqrt(sum((P - T.chanpos(it(k), :)) .^ 2, 2)));
    if ~isempty(d) && d <= max_mm && ~ismember(cand{j}, chosen)
        chosen{end + 1} = cand{j}; %#ok<AGROW>
        parts{end + 1} = sprintf('%s for %s (%.0f mm)', cand{j}, T.label{it(k)}, d); %#ok<AGROW>
    end
end
how = strjoin(parts, ', ');
end

function T = template_cap()
% Brainstorm's ICBM152 ASA 10-05 cap (343 positions), in mm
persistent cap
if isempty(cap)
    cap = get_ftElec('ASA_10-05_343');
    cap.chanpos = cap.chanpos * unit_to_mm(cap);
end
T = cap;
end

function s = unit_to_mm(elec)
s = 1000;                                            % get_ftElec: metres
if isfield(elec, 'unit')
    switch elec.unit
        case 'cm', s = 10;
        case 'mm', s = 1;
    end
end
end
