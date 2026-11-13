function [H, block] = window_mean_halves(cmats, windows, seg, onset, srate, block_s, sz)
%WINDOW_MEAN_HALVES The window mean of each measure on the odd and on the even blocks of a recording.
%   [H, block] = window_mean_halves(cmats, windows, seg, onset, srate, block_s, sz)
%   cmats    1 x nMeasures cell: nWindows x nROI x nROI windows of each measure
%   windows  1 x nMeasures cell: their info.windows (samples, first_sample)
%   seg      nPieces x 2 continuous pieces (samples of the analysed data)
%   onset    time of each piece's first sample in the recording (s)
%   block_s  shortest block (s; S.reliability_block)
%   sz       size(static_cmats)
%   H        2 x nMeasures x nROI x nROI: the mean of the windows that lie
%            wholly inside the pieces of each half (split_halves: the odd
%            and the even blocks); NaN for a half without samples or
%            without a whole window of the measure
%   block    1 x nMeasures block length of each measure (s)
%   Each measure has blocks of max(block_s, 3 x its window): a window lies
%   wholly inside one block only, so the halves share no data, and with
%   blocks of at least 3 windows at least about 2/3 of the windows are used
%   (a window of length w with its start anywhere in a block of length b
%   fits in it with probability (b - w) / b). With blocks of block_s (10 s)
%   for every measure, AEC's 6 s windows would use under half of their
%   windows, and windows longer than a block none (NaN halves).
H = nan([2, sz]);
block = nan(1, numel(cmats));
for k = 1:numel(cmats)
    block(k) = max(block_s, 3 * windows{k}.samples / srate);
    halves = split_halves(seg, onset, srate, block(k));
    first = windows{k}.first_sample(:);
    last = first + windows{k}.samples - 1;
    for h = find(~cellfun(@isempty, halves))
        p = halves{h};
        in = any(first >= p(:, 1).' & last <= p(:, 2).', 2);
        if any(in), H(h, k, :, :) = mean(double(cmats{k}(in, :, :)), 1); end
    end
end
end
