function R = individual_alpha_frequency(X, labels, fs, segments, opt)
%INDIVIDUAL_ALPHA_FREQUENCY Individual alpha frequency (IAF) of one recording.
%   R = individual_alpha_frequency(X, labels, fs, segments, opt)
%   X         channels x samples: the recording's channels matched to the
%             cap, average-referenced (as stage 1 uses them)
%   labels    their labels; fs the sampling rate; segments the continuous
%             pieces (rows [first last] sample)
%   opt       S.iaf: channels (posterior), search ([7 13] Hz), cmin,
%             sgf_frame, sgf_order, measure ('paf' or 'cog'), no_peak
%             ('exclude' or 'cog'), checks ({'irasa', 'specparam'}),
%             max_disagreement (Hz), edge_check (true: a peak within 0.5 Hz
%             of the search range's edge, where the IRASA oscillatory
%             spectrum has no maximum inside the range, is not an alpha
%             peak; needs the 'irasa' check; see "Edge check" below)
%   R.paf, R.cog          peak alpha frequency and centre of gravity of
%                         restingIAF (Corcoran et al. 2018, Psychophysiology
%                         55:e13064; github.com/corcorana/restingIAF):
%                         Welch spectra (4 s), Savitzky-Golay smoothing,
%                         peak criteria, weighted over the channels with a
%                         peak (at least opt.cmin); NaN when not resolved
%   R.n_paf, R.n_cog      channels contributing
%   R.irasa, R.specparam  checks with FieldTrip's ft_freqanalysis, on 4 s
%                         pieces within continuous stretches: the peak of
%                         the oscillatory spectrum after IRASA (Wen & Liu
%                         2016), and the strongest specparam/FOOOF peak in the
%                         search range (Donoghue et al. 2020; FieldTrip calls
%                         Brainstorm's MATLAB implementation, process_fooof);
%                         median over the channels
%   R.irasa_edge          share of the posterior channels whose IRASA
%                         oscillatory spectrum has its largest value in the
%                         search range on the range's bin nearest the IAF's
%                         edge (edge check; NaN when not computed)
%   R.value               the IAF used: opt.measure; when it is NaN and
%                         opt.no_peak is 'cog', the CoG; NaN otherwise (the
%                         recording is then left out of the anchored bands)
%   R.source              'paf', 'cog' or 'none'
%   R.flag                text: why there is no IAF, or a warning (the IAF
%                         and a check more than opt.max_disagreement apart:
%                         kept, to be looked at)
%   R.channels            the channels used
%   Continuous pieces shorter than 4 s are not used; the others are joined
%   for restingIAF (which reads one continuous signal).
%   Edge check: the IRASA spectrum "has no maximum inside the range" in a
%   channel when its largest value in the search range is on the range's
%   first (last) frequency bin. The IAF is rejected when it lies within
%   0.5 Hz of that edge and at least half of the posterior channels have
%   their maximum on that edge's bin. Bins, not frequencies, are compared,
%   channel by channel: the bins of 4 s pieces padded to a power of 2 are
%   0.244 Hz apart at 250, 500 or 1000 Hz, so 7 and 13 Hz are not bins (the
%   check compared the channels' median peak frequency with 7 and 13 Hz
%   within 0.01 Hz: never true at these rates, and at 256 or 512 Hz not
%   when exactly half of an even number of channels were on the edge bin,
%   the median then lying halfway to the next bin).
R = struct('paf', NaN, 'cog', NaN, 'n_paf', 0, 'n_cog', 0, 'irasa', NaN, 'specparam', NaN, ...
    'irasa_edge', NaN, 'value', NaN, 'source', "none", 'flag', "", 'channels', {{}});
ch = find(ismember(lower(string(labels)), lower(string(opt.channels))));
R.channels = cellstr(labels(ch));
if numel(ch) < opt.cmin
    R.flag = sprintf('only %d of the posterior channels (S.iaf.channels) present; at least %d needed', ...
        numel(ch), opt.cmin);
    return
end
n4 = round(4 * fs);
keep = segments(segments(:, 2) - segments(:, 1) + 1 >= n4, :);
if isempty(keep)
    R.flag = "no continuous piece of 4 s";
    return
end
x = double(X(ch, :));
idx = cell2mat(arrayfun(@(k) keep(k, 1):keep(k, 2), (1:size(keep, 1))', 'UniformOutput', false)');

% restingIAF (Corcoran et al. 2018)
fsi = round(fs);
pSum = restingIAF(x(:, idx), numel(ch), opt.cmin, [1 40], fsi, opt.search, opt.sgf_frame, opt.sgf_order);
R.paf = pSum.paf;  R.cog = pSum.cog;
R.n_paf = pSum.pSel;  R.n_cog = pSum.gSel;

% checks with FieldTrip (IRASA, specparam), 4 s pieces
if ~isempty(opt.checks)
    data = struct('label', {cellstr(labels(ch))}, 'fsample', fs, 'trial', {{}}, 'time', {{}}, ...
        'sampleinfo', zeros(0, 2));
    for k = 1:size(keep, 1)
        for t0 = keep(k, 1):n4:(keep(k, 2) - n4 + 1)
            data.trial{end + 1} = x(:, t0:t0 + n4 - 1);
            data.time{end + 1} = (0:n4 - 1) / fs;
            data.sampleinfo(end + 1, :) = [t0, t0 + n4 - 1];
        end
    end
    a = @(f) f >= opt.search(1) & f <= opt.search(2);
    if ismember('irasa', opt.checks)
        cfg = struct('method', 'irasa', 'taper', 'hanning', 'foilim', [2 40], 'pad', 'nextpow2', ...
            'output', 'fractal', 'feedback', 'no');
        frac = ft_freqanalysis(cfg, data);
        cfg.output = 'original';
        orig = ft_freqanalysis(cfg, data);
        osc = orig.powspctrm - frac.powspctrm;
        in = find(a(orig.freq));
        fa = orig.freq(in);
        [~, k] = max(osc(:, in), [], 2);
        R.irasa = median(fa(k));
        at_edge = [k == 1, k == numel(in)];      % per channel: the in-range maximum on the first, last bin
    end
    if ismember('specparam', opt.checks)
        cfg = struct('method', 'mtmfft', 'taper', 'hanning', 'foilim', [2 40], 'pad', 4, ...
            'output', 'fooof_peaks', 'keeptrials', 'no', 'feedback', 'no');
        fr = ft_freqanalysis(cfg, data);
        cf = nan(numel(ch), 1);
        for c = 1:numel(ch)
            pk = fr.fooofparams(c).peak_params;           % [centre frequency, power, bandwidth]
            if ~isempty(pk), pk = pk(a(pk(:, 1)), :); end
            if ~isempty(pk), [~, j] = max(pk(:, 2)); cf(c) = pk(j, 1); end
        end
        R.specparam = median(cf, 'omitnan');
    end
end

% the IAF used
flags = strings(0, 1);
v = R.(lower(opt.measure));
% a "peak" near the edge of the search range while the oscillatory spectrum
% (IRASA) rises out of the range at that edge is the flank of a
% neighbouring rhythm (beta above 13 Hz, theta below 7 Hz), not an alpha
% peak: in simulated recordings without alpha, every method found one at
% 12.4-13.1 Hz. Decided channel by channel on the IRASA bins (help above)
if opt.edge_check && ~isnan(v) && ~isnan(R.irasa) && min(abs(v - opt.search)) < 0.5
    [~, e] = min(abs(v - opt.search));            % the edge the "peak" is near
    R.irasa_edge = mean(at_edge(:, e));
    if R.irasa_edge >= 0.5
        flags(end + 1) = sprintf(['%s %.2f Hz at the edge of the search range, where the IRASA spectrum ' ...
            'has its maximum in the range on the edge bin in %d of %d channels: not an alpha peak'], ...
            upper(opt.measure), v, nnz(at_edge(:, e)), numel(ch));
        v = NaN;
        R.cog = NaN;                              % on the same flank: not a fallback either
    end
end
if ~isnan(v)
    R.value = v;  R.source = string(lower(opt.measure));
elseif strcmp(opt.no_peak, 'cog') && ~isnan(R.cog)
    R.value = R.cog;  R.source = "cog";
    flags(end + 1) = sprintf('no %s: CoG used', upper(opt.measure));
else
    flags(end + 1) = sprintf('no %s on %d or more channels: left out of the anchored bands', ...
        upper(opt.measure), opt.cmin);
end
% the checks disagree: kept, but to be looked at (stage1_summary.csv). In
% simulation, a weak alpha next to strong beta misled restingIAF and IRASA
% (search range only) and not specparam (separate peaks); in the 8
% development recordings specparam differed by just over 0.5 Hz twice
for c = {'irasa', 'specparam'}
    if ~isnan(R.value) && ~isnan(R.(c{1})) && abs(R.value - R.(c{1})) > opt.max_disagreement
        flags(end + 1) = sprintf('IAF %.2f Hz and %s peak %.2f Hz differ by more than %.2g Hz', ...
            R.value, c{1}, R.(c{1}), opt.max_disagreement); %#ok<AGROW>
    end
end
R.flag = strjoin(flags, '; ');
end
