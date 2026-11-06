function starts = window_starts(seg, win_samples, step_samples)
%WINDOW_STARTS First sample of every window that lies inside one continuous piece.
%   starts = connectivity.window_starts(seg, win_samples, step_samples)
%   seg: nPieces x 2 (first and last sample of each piece; [1, nSamples]
%   for a continuous recording). Windows start at the beginning of each
%   piece and never cross a cut (boundary event, epoch limit), so no window
%   mixes data from both sides of a discontinuity. For one piece the
%   windows are the original ones: 1, 1 + step, ... up to nSamples - win + 1.
starts = cell(1, size(seg, 1));
for k = 1:size(seg, 1)
    starts{k} = seg(k, 1):step_samples:(seg(k, 2) - win_samples + 1);
end
starts = [starts{:}];
end
