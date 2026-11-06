function halves = split_halves(seg, onset_s, srate, block_s)
%SPLIT_HALVES The samples of a recording in two halves: odd and even blocks of time.
%   halves = split_halves(seg, onset_s, srate, block_s)
%   seg      nPieces x 2 continuous pieces (samples of the analysed data)
%   onset_s  time of each piece's first sample in the recording (s)
%   block_s  block length (s; S.reliability_block, 10)
%   halves   1 x 2 cell: the pieces (first and last samples) of the odd and
%            of the even blocks, blocks counted from the start of the
%            recording, so a cut does not shift them
%   Connectivity computed on each half (static_connectivity with these
%   pieces) gives two estimates from interleaved stretches of the same
%   recording: their agreement across recordings is the split-half
%   reliability of a measure (static_reliability.m). Interleaved blocks, rather
%   than first and second half, keep slow drifts (drowsiness, electrode
%   drift) out of the difference between the halves.
halves = {zeros(0, 2), zeros(0, 2)};
for k = 1:size(seg, 1)
    s = seg(k, 1):seg(k, 2);
    t = onset_s(k) + (s - seg(k, 1)) / srate;
    h = mod(floor(t / block_s + 1e-9), 2);         % 0: odd blocks (the first), 1: even
    edges = [1, find(diff(h) ~= 0) + 1, numel(s) + 1];
    for r = 1:numel(edges) - 1
        halves{h(edges(r)) + 1}(end + 1, :) = s([edges(r), edges(r + 1) - 1]);
    end
end
end
