function run_psd_analysis(S)
%RUN_PSD_ANALYSIS  Stage 1b: power spectra at the scalp and in the regions, group statistics and maps.
%   run_psd_analysis() uses pipeline_settings(); run_psd_analysis(S) uses
%   the settings S (a variable S in the workspace is not used).
%   run_full_pipeline runs it after stage 1 when S.psd.run is true. It
%   reads the recordings itself (stage 1's outputs are not needed), with
%   the settings of S.psd (pipeline_settings.m, section "PSD analysis"):
%   1. Spectra (psd_spectra; kept with keys in results/PSD/spectra, only
%      what changed is recomputed): each recording read as stage 1 reads it
%      (analysed stretch, continuous pieces, the channels of the first cap,
%      flat channels left out, average reference); the PSD of each channel
%      ('sensor') and of each region of the first atlas for every inverse
%      method of S.inv ('source', the regions' signals of stage 1:
%      psd_source_signals), by Welch's method (MATLAB's pwelch on each
%      continuous piece, Hann segments of S.psd.segment s, S.psd.overlap;
%      psd_welch), and its specparam fit over S.psd.fit_range (psd_specparam).
%   2. Measures (psd_measures) of each channel or region, for the bands of
%      S.psd.bands (their own setting, not the connectivity bands of
%      S.bands; edges can be anchored on the IAF as there, band_edges):
%      'absolute' (log10 band power), 'relative' (% of S.psd.total),
%      'periodic' (log10 power above the 1/f fit), 'aperiodic' (exponent
%      and offset of the fit) and 'alpha_frequency' (specparam's alpha
%      peak, or the centre of gravity above the fit without one), as
%      listed in S.psd.measures. Channels: those of every recording (a
%      warning names those left out).
%   3. Statistics (psd_statistics), the methods of stage 3: each level x
%      measure x band is a family whose states are the channels or regions
%      (permu_ttest_permu_corr_BNS): every comparison of S.comparisons
%      (Welch t, or ANCOVA with S.covariates and Freedman-Lane
%      permutations; Hedges' g) and correlation set (S.correlations,
%      S.covariate_sets), S.psd.n_perm permutations (S.n_perm when empty),
%      max-T over the channels or regions (p_maxT) and over the family
%      (p_maxT_family, family_p), Benjamini-Hochberg q within the family
%      (q_fdr, a description), and the families of the PSD selected by
%      Benjamini-Bogomolov (family_adjust: p_bb, significant); with
%      S.psd.cluster, a cluster
%      test over neighbouring channels or regions on the same permutations
%      (cluster_perm_test, psd_adjacency; cluster_p, then Benjamini-
%      Hochberg over all clusters of the PSD: cluster_q, cluster_significant).
%   4. Figures (psd_figures, S.figures): group spectra (absolute and above
%      the 1/f fit), box plots of the measures averaged over channels or
%      regions, scalp and cortex maps of Hedges' g (or r) with the
%      significant channels, regions and clusters marked, and the
%      significant clusters' values.
%   Outputs in results/PSD/<S.run_name, or main>/:
%     psd_log.txt, provenance_psd.mat    the log, the settings and versions
%     recordings.csv                     each recording: what was analysed,
%                                        fit quality, IAF, failures
%     values_<level>_<measure>.csv       the measures: one row per
%                                        recording, a column per band x channel/region
%     psd_values.mat                     the same and the spectra (recording x
%                                        unit x frequency), for other plots
%     Statistical analysis/              summary.csv, families.csv,
%                                        clusters.csv, neighbours_<level>.csv,
%                                        <level>/<family>_group_comparison.csv,
%                                        <family>_correlations.csv (.xlsx too)
%     figures/                           the figures
%   <level> is 'sensor' or 'source_<inverse>' (e.g. source_eloreta).
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if nargin < 1, S = []; end
S = setup_toolboxes(stage_settings(S, 'run_psd_analysis'));
o = S.psd;
name = 'main';
if isfield(S, 'run_name') && ~isempty(S.run_name), name = char(S.run_name); end
out = fullfile(S.results, 'PSD', name);
if ~isfolder(out), mkdir(out); end
diary(fullfile(out, 'psd_log.txt'));
stop_diary = onCleanup(@() diary('off'));
fprintf(['PSD analysis, %s\nLevels: %s; bands: %s; measures: %s\nSpectra: Welch, %g s Hann segments, %g%% ' ...
    'overlap, %g-%g Hz; specparam %g-%g Hz; relative power: %% of %g-%g Hz\n'], char(datetime('now')), ...
    strjoin(o.levels, ', '), strjoin(arrayfun(@band_text, o.bands, 'UniformOutput', false), ', '), ...
    strjoin(o.measures, ', '), o.segment, 100 * o.overlap, o.freq_range, o.fit_range, o.total);
t0 = tic;
[R, stems] = psd_spectra(S);
L = assemble(S, R, stems);
write_values(S, L, R, stems, out);
stats = struct('families', table(), 'summary', table(), 'clusters', table(), 'res', struct([]));
has_tests = ~isempty(S.comparisons) || ~isempty(S.correlations);
if has_tests && ~isempty(S.subs_info)
    try
        stats = psd_statistics(S, L, fullfile(out, 'Statistical analysis'));
    catch err
        fprintf(2, 'PSD statistics failed: %s\n', err.getReport('extended', 'hyperlinks', 'off'));
    end
else
    fprintf('No comparisons or correlations (S.comparisons, S.correlations): spectra and measures only\n');
end
if S.figures
    try
        psd_figures(S, L, stats, fullfile(out, 'figures'));
    catch err
        fprintf(2, 'PSD figures failed: %s\n', err.getReport('extended', 'hyperlinks', 'off'));
    end
end
provenance = struct('stage', '1b', 'date', char(datetime('now')), 'settings', S, ...
    'recordings', {stems}, 'families', stats.families, 'summary', stats.summary); %#ok<NASGU>
save(fullfile(out, 'provenance_psd.mat'), 'provenance');
fprintf('PSD analysis finished in %.1f min: %s\n', toc(t0) / 60, out);
end


function L = assemble(S, R, stems)
% one entry per level (sensor, source_<inverse>): the recordings with
% spectra, the channels or regions, the spectra and the measures
o = S.psd;
n = numel(R);
ids = subject_ids(S, stems);
names = {};  lev = {};  inv = {};
if ismember('sensor', o.levels), names{end + 1} = 'sensor'; lev{end + 1} = 'sensor'; inv{end + 1} = ''; end
if ismember('source', o.levels)
    for iv = 1:numel(S.inv)
        names{end + 1} = ['source_' S.inv{iv}]; lev{end + 1} = 'source'; inv{end + 1} = S.inv{iv}; %#ok<AGROW>
    end
end
iaf = nan(n, 1);
for s = 1:n
    if ~isfield(R(s).sensor, 'error') && isfield(R(s).sensor, 'iaf'), iaf(s) = R(s).sensor.iaf.value; end
end
L = struct('name', {}, 'level', {}, 'inverse', {}, 'units', {}, 'rec', {}, 'stems', {}, 'ids', {}, ...
    'freq', {}, 'pow', {}, 'flat', {}, 'M', {}, 'iaf', {}, 'failed', {});
for k = 1:numel(names)
    spec = cell(n, 1);  failed = strings(0, 1);
    for s = 1:n
        if strcmp(lev{k}, 'sensor')
            x = R(s).sensor;
        else
            x = R(s).source;
            if ~isfield(x, 'error')
                x = x.inverses(strcmp({x.inverses.inverse}, inv{k}));
                if isempty(x), x = struct('error', "not computed for this inverse method"); end
            end
        end
        if isfield(x, 'error')
            failed(end + 1) = string(stems{s}) + ": " + x.error; %#ok<AGROW>
        else
            spec{s} = x;
        end
    end
    ok = find(~cellfun(@isempty, spec));
    if ~isempty(failed)
        warning('psd:failed', '%s: %d recording(s) left out: %s', names{k}, numel(failed), strjoin(failed, '; '));
    end
    if isempty(ok)
        warning('psd:none', '%s: no recording has spectra', names{k});
        continue
    end
    % the channels (regions) of every recording, in the order of the first
    all_labels = cellfun(@(x) lower(x.labels), spec(ok), 'UniformOutput', false);
    units = spec{ok(1)}.labels;
    present = cell2mat(cellfun(@(l) ismember(lower(units), l), all_labels(:), 'UniformOutput', false));
    common = units(all(present, 1));
    if strcmp(lev{k}, 'sensor')
        elec = get_ftElec(S.montages{1});                 % the cap's order
        [~, ord] = ismember(lower(common), lower(elec.label));
        [~, s_] = sort(ord);  common = common(s_);
        everyone = unique([all_labels{:}]);
        left = setdiff(everyone, lower(common));
        if ~isempty(left)
            cnt = cellfun(@(c) sum(cellfun(@(l) ismember(c, l), all_labels)), left);
            warning('psd:channels', ['Channels left out of the PSD maps and tests (not in every recording; ' ...
                'recordings with them): %s'], strjoin(compose('%s (%d)', string(left(:)), cnt(:)), ', '));
        end
    elseif numel(common) < numel(units)
        error('psd:regions', '%s: the recordings have different regions', names{k});
    end
    K = numel(common);
    freq = spec{ok(1)}.freq;
    F = numel(freq);
    nr = numel(ok);
    pow = nan(nr, K, F, 'single');  flat = pow;
    nb = numel(o.bands);
    M = struct('absolute', nan(nr, K, nb), 'relative', nan(nr, K, nb), 'periodic', nan(nr, K, nb), ...
        'exponent', nan(nr, K), 'offset', nan(nr, K), 'fit_r2', nan(nr, K), 'alpha_frequency', nan(nr, K), ...
        'alpha_peak', false(nr, K), 'edges', nan(nr, nb, 2));
    for i = 1:nr
        x = spec{ok(i)};
        [~, row] = ismember(lower(common), lower(x.labels));
        P = double(x.pow(row, :));
        fit = x.fit;
        for f = {'offset', 'exponent', 'r2', 'error'}, fit.(f{1}) = fit.(f{1})(row); end
        fit.peaks = fit.peaks(row);  fit.ap = fit.ap(row, :);
        m = psd_measures(freq, P, fit, o, iaf(ok(i)), S.iaf);
        pow(i, :, :) = P;
        flat(i, :, :) = log10(P) - fit.ap;
        for f = {'absolute', 'relative', 'periodic'}, M.(f{1})(i, :, :) = m.(f{1}); end
        for f = {'exponent', 'offset', 'alpha_frequency', 'alpha_peak'}, M.(f{1})(i, :) = m.(f{1}); end
        M.fit_r2(i, :) = m.fit_r2;
        M.edges(i, :, :) = m.edges;
    end
    L(end + 1) = struct('name', names{k}, 'level', lev{k}, 'inverse', inv{k}, 'units', {common(:)'}, ...
        'rec', ok(:), 'stems', {stems(ok)}, 'ids', ids(ok), 'freq', freq, 'pow', pow, 'flat', flat, ...
        'M', M, 'iaf', iaf(ok), 'failed', failed); %#ok<AGROW>
    fprintf('%s: %d recordings, %d %s; median specparam R^2 %.3f\n', names{k}, nr, K, ...
        unit_word(lev{k}), median(M.fit_r2(:), 'omitnan'));
end
if isempty(L), error('psd:none', 'No recording has spectra (see the messages above)'); end
end

function write_values(S, L, R, stems, out)
% recordings.csv (each recording), values_<level>_<measure>.csv and psd_values.mat
o = S.psd;
n = numel(stems);
Q = table(string(stems(:)), subject_ids(S, stems), 'VariableNames', {'recording', 'ID'});
Q.srate = nan(n, 1);  Q.analysed_s = nan(n, 1);  Q.n_pieces = nan(n, 1);  Q.n_segments = nan(n, 1);
Q.n_channels = nan(n, 1);  Q.iaf = nan(n, 1);  Q.iaf_flag = strings(n, 1);  Q.status = repmat("ok", n, 1);
for s = 1:n
    x = R(s).sensor;
    if isfield(x, 'error'), Q.status(s) = "failed: " + x.error; continue; end
    Q.srate(s) = x.info.srate;  Q.analysed_s(s) = x.info.analysed_s;  Q.n_pieces(s) = x.info.n_pieces;
    Q.n_segments(s) = x.n_segments;  Q.n_channels(s) = numel(x.labels);
    Q.iaf(s) = x.iaf.value;
    if isfield(x.iaf, 'flag'), Q.iaf_flag(s) = string(x.iaf.flag); end
    if isfield(R(s), 'source') && isfield(R(s).source, 'error')
        Q.status(s) = "regions failed: " + R(s).source.error;
    end
end
for k = 1:numel(L)                                   % fit quality and alpha peaks of each level
    v = nan(n, 1);  a = nan(n, 1);
    v(L(k).rec) = median(L(k).M.fit_r2, 2, 'omitnan');
    a(L(k).rec) = mean(L(k).M.alpha_peak, 2);
    Q.(['fit_r2_' L(k).name]) = v;
    Q.(['alpha_peak_share_' L(k).name]) = a;
end
writetable(Q, fullfile(out, 'recordings.csv'));
bands = {o.bands.name};
for k = 1:numel(L)
    base = table(string(L(k).stems(:)), L(k).ids(:), 'VariableNames', {'recording', 'ID'});
    for m = {'absolute', 'relative', 'periodic'}
        V = L(k).M.(m{1});
        cols = strings(1, 0);  X = zeros(size(V, 1), 0);
        for b = 1:numel(bands)
            cols = [cols, string(bands{b}) + "__" + string(L(k).units)]; %#ok<AGROW>
            X = [X, V(:, :, b)]; %#ok<AGROW>
        end
        T = [base, array2table(X, 'VariableNames', cellstr(matlab.lang.makeUniqueStrings(cols)))];
        writetable(T, fullfile(out, sprintf('values_%s_%s.csv', L(k).name, m{1})));
    end
    for m = {'exponent', 'offset', 'alpha_frequency', 'fit_r2'}
        T = [base, array2table(L(k).M.(m{1}), 'VariableNames', cellstr(matlab.lang.makeUniqueStrings(string(L(k).units))))];
        writetable(T, fullfile(out, sprintf('values_%s_%s.csv', L(k).name, m{1})));
    end
end
psd = rmfield(L, 'rec'); %#ok<NASGU>
save(fullfile(out, 'psd_values.mat'), 'psd', 'bands', '-v7.3');
end

function s = band_text(b)
s = sprintf('%s %s-%s Hz', b.name, edge_text(b.fmin), edge_text(b.fmax));
end

function s = edge_text(e)
if isnumeric(e), s = sprintf('%g', e); else, s = char(e); end
end

function w = unit_word(level)
w = 'channels';
if strcmp(level, 'source'), w = 'regions'; end
end
