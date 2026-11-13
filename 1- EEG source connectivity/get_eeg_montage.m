function [eeg, labels, idx_montage, report] = get_eeg_montage(eeg, eeg_labels, montage_labels, channel_subset)
%GET_EEG_MONTAGE Keep the channels of the recording that the montage knows, in the montage's order.
%   [eeg, labels, idx_montage, report] = get_eeg_montage(eeg, eeg_labels, montage_labels, channel_subset)
%   eeg             nChannels x nSamples, rows labelled by eeg_labels (from the file)
%   montage_labels  labels of the montage (get_ftElec)
%   channel_subset  labels to keep ({} = every matched channel), e.g. a
%                   32- or 19-channel subset, for any montage
%   eeg             rows of the matched channels, in montage order
%   labels          their montage labels; idx_montage their rows in the montage
%   report          channels of the file left out, channels of the subset
%                   missing, reference suffixes of the matched labels
%   Labels are matched, ignoring case, spaces, an
%   'EEG' prefix and a '-REF'/'-LE'/'-AR'/'-M1' suffix (EDF style), with the
%   old 10-20 names T3/T4/T5/T6 accepted for T7/T8/P7/P8 when the montage only
%   has the new ones (and the reverse). Matched channels with different
%   reference suffixes (e.g. C3-A2 and C4-A1: each against the other
%   mastoid) are not on one common reference, which the source analysis
%   needs: a warning (get_eeg_montage:mixed_reference).
if nargin < 4, channel_subset = {}; end
key_eeg = normalise(eeg_labels);
key_mon = normalise(montage_labels);
[found, loc] = ismember(key_mon, key_eeg);
alias = {'t3', 't7'; 't4', 't8'; 't5', 'p7'; 't6', 'p8'};
for a = 1:size(alias, 1)
    for d = 1:2
        m = strcmp(key_mon, alias{a, d}) & ~found;
        [f2, l2] = ismember(alias{a, 3 - d}, key_eeg);
        if any(m) && f2 && ~any(strcmp(key_mon, alias{a, 3 - d}))
            found(m) = true;  loc(m) = l2;
        end
    end
end
if ~isempty(channel_subset)
    in_subset = ismember(key_mon, normalise(channel_subset));
    report.subset_missing = setdiff(normalise(channel_subset), key_mon(found & in_subset));
    found = found & in_subset;
else
    report.subset_missing = {};
end
idx_montage = find(found);
labels = montage_labels(idx_montage);
eeg = eeg(loc(idx_montage), :);
report.left_out = eeg_labels(~ismember(1:numel(eeg_labels), loc(idx_montage)));
suf = regexp(lower(strtrim(cellstr(eeg_labels(loc(idx_montage))))), ...
    ['[\s_-]+(' ref_suffixes ')$'], 'tokens', 'once');
suf = unique(cellfun(@(t) [t{:}], suf(~cellfun(@isempty, suf)), 'UniformOutput', false));
report.reference_suffixes = suf(:)';
if numel(suf) > 1
    warning('get_eeg_montage:mixed_reference', ['Channels with different reference suffixes (%s): ' ...
        'not one common reference; re-reference the recording before the source analysis'], strjoin(upper(suf), ', '));
end
end

function r = ref_suffixes
r = 'ref|le|ar|avg|m1|m2|a1|a2|cz';
end

function k = normalise(l)
k = lower(strtrim(cellstr(l(:))));
k = regexprep(k, '^eeg[\s_-]*(?=[a-z])', '');              % 'EEG Fp1', 'EEGFp1' (EEG001 stays)
k = regexprep(k, ['[\s_-]+(' ref_suffixes ')$'], '');      % 'Fp1-REF' (Cz, A1 stay: a separator is required)
k = regexprep(k, '\s+', '');
end
