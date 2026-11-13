function failures = run_psd_analysis(S)
%RUN_PSD_ANALYSIS  Stage 2: power spectra at the scalp and in the regions, group statistics and maps.
%   run_psd_analysis() uses pipeline_settings(); run_psd_analysis(S) uses
%   the settings S (a variable S in the workspace is not used).
%   failures = run_psd_analysis(...): the steps that failed (string array:
%   what, message; empty when none), as run_statistics. A failing step
%   (statistics, figures, region check, vertex maps) does not stop the
%   others, unless S.stop_on_error = true; they are listed in failures.csv
%   and run_full_pipeline then ends "Done with errors" (when only figures
%   failed: "every result written, but N figure step(s) failed",
%   failure_summary). Every run starts
%   from an empty results/PSD/<run>/ (the spectra, kept in results/PSD/spectra,
%   are reused), so no table, figure or log of an earlier run is left
%   beside the new ones.
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
%      peak; NaN without one), as listed in S.psd.measures; the centre of
%      gravity of the power above the fit in the alpha range is reported,
%      not tested. Channels: those of every recording (a warning names
%      those left out). Each recording's spectra are checked for narrow
%      lines (psd_line_check) and spectra that rise over the fit range
%      (a negative exponent): recordings.csv and a warning; so are bands
%      whose periodic power starts within 2 Hz of the start of the fit
%      (near_fit_start), where it depends on where the fit starts.
%   3. Statistics (psd_statistics), the methods of the statistics: each level x
%      measure x band is a family whose states are the channels or regions
%      (permu_ttest_permu_corr_BNS): every comparison of S.comparisons
%      (Welch t, or ANCOVA with S.covariates and Freedman-Lane
%      permutations; Hedges' g, adjusted with covariates) and correlation
%      set (S.correlations, S.covariate_sets; the slope over its HC3
%      standard error, S.correlation_type), S.psd.n_perm permutations (S.n_perm when empty),
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
%      regions, scalp and cortex maps of g (or r) with the
%      significant channels, regions and clusters marked, and the
%      significant clusters' values.
%   5. The regions' two spectra (psd_region_check): the power of each
%      region's signal against the mean of its sources' power, which no
%      cancellation lowers (their ratio, correlation and group g).
%      S.psd.region_spectra = 'vertex_mean' runs steps 2 to 4 on the
%      second (levels source_<inverse>_vertexmean).
%   6. The power of every source of the cortex (psd_vertex_maps, with
%      S.psd.vertex_maps): values, group means and g (r) per band on the
%      cortex, LORETA-style; descriptive, not tested.
%   Outputs in results/PSD/<S.run_name, or main>/:
%     psd_log.txt, provenance_psd.mat    the log, the settings and versions
%     failures.csv                       the steps that failed (none: no file)
%     recordings.csv                     each recording: what was analysed,
%                                        fit quality, IAF, lines, rising
%                                        spectra, bands near the fit's
%                                        start (near_fit_start), failures
%     values_<level>_<measure>.csv       the measures: one row per
%                                        recording, a column per band x channel/region
%     psd_values.mat                     the same and the spectra (recording x
%                                        unit x frequency), for other plots
%     region_power_check.csv             the regions' two spectra compared
%                                        (psd_region_check)
%     vertex_values_<inverse>.mat        every source's band power, group
%                                        means, g and r (psd_vertex_maps)
%     Statistical analysis/              summary.csv, families.csv,
%                                        clusters.csv, neighbours_<level>.csv,
%                                        <level>/<family>_group_comparison.csv,
%                                        <family>_correlations.csv (.xlsx too)
%     figures/                           the figures
%   <level> is 'sensor' or 'source_<inverse>' (e.g. source_eloreta;
%   source_<inverse>_vertexmean with S.psd.region_spectra = 'vertex_mean').
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
if nargin < 1, S = []; end
S = setup_toolboxes(stage_settings(S, 'run_psd_analysis'));
o = S.psd;
name = 'main';
if isfield(S, 'run_name') && ~isempty(S.run_name), name = char(S.run_name); end
out = fullfile(S.results, 'PSD', name);
if isfolder(out), rmdir(out, 's'); end              % a new run: nothing of the earlier one kept
mkdir(out);
diary(fullfile(out, 'psd_log.txt'));
stop_diary = onCleanup(@() diary('off'));
fprintf(['PSD analysis, %s\nLevels: %s; bands: %s; measures: %s\nSpectra: Welch, %g s Hann segments, %g%% ' ...
    'overlap, %g-%g Hz; specparam %g-%g Hz; relative power: %% of %g-%g Hz\n'], char(datetime('now')), ...
    strjoin(o.levels, ', '), strjoin(arrayfun(@band_text, o.bands, 'UniformOutput', false), ', '), ...
    strjoin(o.measures, ', '), o.segment, 100 * o.overlap, o.freq_range, o.fit_range, o.total);
t0 = tic;
failures = strings(0, 2);                           % what, message (failures.csv)
stop = isfield(S, 'stop_on_error') && S.stop_on_error;
[R, stems] = psd_spectra(S);
L = assemble(S, R, stems);
write_values(S, L, R, stems, out);
stats = struct('families', table(), 'summary', table(), 'clusters', table(), 'res', struct([]));
has_tests = ~isempty(S.comparisons) || ~isempty(S.correlations);
if has_tests && ~isempty(S.subs_info)
    try
        check_subject_table(S, false);           % the table against the recordings (also when run alone)
        stats = psd_statistics(S, L, fullfile(out, 'Statistical analysis'));
    catch err
        failures = note_failure(failures, 'PSD statistics', err, stop);
    end
else
    fprintf('No comparisons or correlations (S.comparisons, S.correlations): spectra and measures only\n');
end
if S.figures
    try
        psd_figures(S, L, stats, fullfile(out, 'figures'));
    catch err
        failures = note_failure(failures, 'PSD figures', err, stop);
    end
end
try                                                  % after psd_figures, which empties the figures folder
    psd_region_check(S, L, out);
catch err
    failures = note_failure(failures, 'region power check', err, stop);
end
try
    psd_vertex_maps(S, L, out);
catch err
    failures = note_failure(failures, 'vertex maps', err, stop);
end
if ~isempty(failures)
    writetable(table(failures(:, 1), failures(:, 2), 'VariableNames', {'what', 'message'}), ...
        fullfile(out, 'failures.csv'));
end
provenance = struct('stage', 2, 'date', char(datetime('now')), 'settings', S, ...
    'recordings', {stems}, 'families', stats.families, 'summary', stats.summary, ...
    'failures', failures);
save(fullfile(out, 'provenance_psd.mat'), 'provenance');
if isempty(failures)
    fprintf('PSD analysis finished in %.1f min: %s\n', toc(t0) / 60, out);
else
    fprintf(2, 'PSD analysis finished in %.1f min WITH %d FAILURE(S) (failures.csv):\n', toc(t0) / 60, ...
        size(failures, 1));
    fprintf(2, '  %s: %s\n', failures');
    fprintf('Outputs: %s\n', out);
end
end

function failures = note_failure(failures, what, err, stop)
% a step that failed: named on the console and in failures.csv; with
% S.stop_on_error the error stops the stage
fprintf(2, '%s failed: %s\n', what, err.getReport('extended', 'hyperlinks', 'off'));
if stop
    diary off
    rethrow(err);
end
failures(end + 1, :) = [string(what), string(err.message)];
end

function L = assemble(S, R, stems)
% one entry per level (sensor, source_<inverse>): the recordings with
% spectra, the channels or regions, the spectra and the measures
o = S.psd;
n = numel(R);
ids = subject_ids(S, stems);
names = {};  lev = {};  inv = {};
if ismember('sensor', o.levels), names{end + 1} = 'sensor'; lev{end + 1} = 'sensor'; inv{end + 1} = ''; end
vm = strcmp(o.region_spectra, 'vertex_mean');         % the regions' spectra: the mean of their sources' PSDs
if ismember('source', o.levels)
    for iv = 1:numel(S.inv)
        names{end + 1} = ['source_' S.inv{iv} repmat('_vertexmean', 1, vm)]; lev{end + 1} = 'source'; %#ok<AGROW>
        inv{end + 1} = S.inv{iv}; %#ok<AGROW>
    end
end
iaf = nan(n, 1);
for s = 1:n
    if ~isfield(R(s).sensor, 'error') && isfield(R(s).sensor, 'iaf'), iaf(s) = R(s).sensor.iaf.value; end
end
L = struct('name', {}, 'level', {}, 'inverse', {}, 'units', {}, 'rec', {}, 'stems', {}, 'ids', {}, ...
    'freq', {}, 'pow', {}, 'flat', {}, 'M', {}, 'iaf', {}, 'failed', {}, 'check', {});
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
        'alpha_peak', false(nr, K), 'alpha_cog', nan(nr, K), 'edges', nan(nr, nb, 2));
    chk = [];                                        % the regions: both spectra's measures (psd_region_check)
    if strcmp(lev{k}, 'source')
        e = struct('absolute', nan(nr, K, nb), 'relative', nan(nr, K, nb), 'periodic', nan(nr, K, nb));
        chk = struct('signal', e, 'vertex_mean', e);
    end
    for i = 1:nr
        x = spec{ok(i)};
        [~, row] = ismember(lower(common), lower(x.labels));
        if strcmp(lev{k}, 'source') && vm
            P = double(x.pow_vertex_mean(row, :));
            fit = fit_rows(x.fit_vertex_mean, row);
        else
            P = double(x.pow(row, :));
            fit = fit_rows(x.fit, row);
        end
        m = psd_measures(freq, P, fit, o, iaf(ok(i)), S.iaf);
        if ~isempty(chk) && isfield(x, 'pow_vertex_mean')
            kinds = {'signal', 'vertex_mean'};
            if vm, Po = x.pow;  fo = x.fit; else, Po = x.pow_vertex_mean;  fo = x.fit_vertex_mean; end   % the other one
            other = psd_measures(freq, double(Po(row, :)), fit_rows(fo, row), o, iaf(ok(i)), S.iaf);
            for f = {'absolute', 'relative', 'periodic'}
                chk.(kinds{1 + vm}).(f{1})(i, :, :) = m.(f{1});
                chk.(kinds{2 - vm}).(f{1})(i, :, :) = other.(f{1});
            end
        end
        pow(i, :, :) = P;
        flat(i, :, :) = log10(P) - fit.ap;
        for f = {'absolute', 'relative', 'periodic'}, M.(f{1})(i, :, :) = m.(f{1}); end
        for f = {'exponent', 'offset', 'alpha_frequency', 'alpha_peak', 'alpha_cog'}, M.(f{1})(i, :) = m.(f{1}); end
        M.fit_r2(i, :) = m.fit_r2;
        M.edges(i, :, :) = m.edges;
    end
    L(end + 1) = struct('name', names{k}, 'level', lev{k}, 'inverse', inv{k}, 'units', {common(:)'}, ...
        'rec', ok(:), 'stems', {stems(ok)}, 'ids', ids(ok), 'freq', freq, 'pow', pow, 'flat', flat, ...
        'M', M, 'iaf', iaf(ok), 'failed', failed, 'check', chk); %#ok<AGROW>
    fprintf('%s: %d recordings, %d %s; median specparam R^2 %.3f\n', names{k}, nr, K, ...
        unit_word(lev{k}), median(M.fit_r2(:), 'omitnan'));
end
if isempty(L), error('psd:none', 'No recording has spectra (see the messages above)'); end
end

function fit = fit_rows(fit, row)
% the specparam fit (psd_specparam) of the units row, in that order
for f = {'offset', 'exponent', 'r2', 'error'}, fit.(f{1}) = fit.(f{1})(row); end
fit.peaks = fit.peaks(row);  fit.ap = fit.ap(row, :);
end

function write_values(S, L, R, stems, out)
% recordings.csv (each recording), values_<level>_<measure>.csv and psd_values.mat
o = S.psd;
n = numel(stems);
Q = table(string(stems(:)), subject_ids(S, stems), 'VariableNames', {'recording', 'ID'});
Q.srate = nan(n, 1);  Q.analysed_s = nan(n, 1);  Q.n_pieces = nan(n, 1);  Q.n_segments = nan(n, 1);
Q.n_channels = nan(n, 1);  Q.iaf = nan(n, 1);  Q.iaf_flag = strings(n, 1);  Q.status = repmat("ok", n, 1);
Q.line_hz = nan(n, 1);  Q.line_ratio = nan(n, 1);  Q.line_channels = nan(n, 1);
for s = 1:n
    x = R(s).sensor;
    if isfield(x, 'error'), Q.status(s) = "failed: " + x.error; continue; end
    Q.srate(s) = x.info.srate;  Q.analysed_s(s) = x.info.analysed_s;  Q.n_pieces(s) = x.info.n_pieces;
    Q.n_segments(s) = x.n_segments;  Q.n_channels(s) = numel(x.labels);
    [Q.line_hz(s), Q.line_ratio(s), ~, Q.line_channels(s)] = psd_line_check(x.freq, x.pow);   % narrow lines
    Q.iaf(s) = x.iaf.value;
    if isfield(x.iaf, 'flag'), Q.iaf_flag(s) = string(x.iaf.flag); end
    if isfield(R(s), 'source') && isfield(R(s).source, 'error')
        Q.status(s) = "regions failed: " + R(s).source.error;
    end
end
fl = find(Q.line_ratio >= 2);
if ~isempty(fl)
    warning('psd:line', ['%d recording(s) with a narrow line in the mean of their channels'' spectra, at ' ...
        'least twice the power 1-2 Hz away (recordings.csv: line_hz, line_ratio, line_channels): %s. A line ' ...
        'enters the band powers, the total of relative power and the specparam fit; remove it before the ' ...
        'pipeline (docs/psd.md, "Line noise")'], numel(fl), ...
        strjoin(compose('%s (%.2f Hz, %.1f)', Q.recording(fl), Q.line_hz(fl), Q.line_ratio(fl)), ', '));
end
for k = 1:numel(L)                                   % fit quality, alpha peaks, rising spectra of each level
    v = nan(n, 1);  a = nan(n, 1);  u = nan(n, 1);
    v(L(k).rec) = median(L(k).M.fit_r2, 2, 'omitnan');
    a(L(k).rec) = mean(L(k).M.alpha_peak, 2);
    u(L(k).rec) = sum(L(k).M.exponent < 0, 2);
    Q.(['fit_r2_' L(k).name]) = v;
    Q.(['alpha_peak_share_' L(k).name]) = a;
    Q.(['rising_' L(k).name]) = u;
    if any(u > 0)
        warning('psd:rising', ['%s: %d %s in %d recording(s) have a spectrum that rises over the fit range ' ...
            '(a negative exponent; recordings.csv: rising_%s), usually an artefact such as muscle activity'], ...
            L(k).name, sum(u, 'omitnan'), unit_word(L(k).level), nnz(u > 0), L(k).name);
    end
end
if ismember('periodic', o.measures)                  % bands whose periodic power starts near the fit's start
    [Q.near_fit_start, msg] = near_fit_start(o, L, n);
    if ~isempty(msg)
        warning('psd:fit_start', ['Periodic power of %s starts within 2 Hz of the start of the specparam fit ' ...
            '(S.psd.fit_range(1) = %g Hz), where it depends on where the fit starts: the aperiodic fit is ' ...
            'least constrained at its edge. Repeat a periodic result of these bands with a later fit start ' ...
            '(e.g. S.psd.fit_range = [%g %g]), and check that the groups do not differ in the recordings ' ...
            'concerned (recordings.csv: near_fit_start; docs/psd.md, "Reading the measures")'], msg, ...
            o.fit_range(1), o.fit_range(1) + 1, o.fit_range(2));
    end
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
    for m = {'exponent', 'offset', 'alpha_frequency', 'alpha_cog', 'fit_r2'}
        T = [base, array2table(L(k).M.(m{1}), 'VariableNames', cellstr(matlab.lang.makeUniqueStrings(string(L(k).units))))];
        writetable(T, fullfile(out, sprintf('values_%s_%s.csv', L(k).name, m{1})));
    end
end
psd = rmfield(L, 'rec');
save(fullfile(out, 'psd_values.mat'), 'psd', 'bands', '-v7.3');
end

function [col, msg] = near_fit_start(o, L, n)
% per recording, the bands whose part inside the fit range, over which psd_measures takes periodic power,
% starts less than 2 Hz above the fit's first frequency (the edges of the first level: every level has the
% same, from the recording's IAF); msg: each such band with its number of recordings ('' if none)
col = strings(n, 1);  msg = '';
f = L(1).freq;
fit_lo = min(f(f >= o.fit_range(1) - 1e-9));  fit_hi = max(f(f <= o.fit_range(2) + 1e-9));
lo = L(1).M.edges(:, :, 1);  hi = L(1).M.edges(:, :, 2);   % recordings x bands; NaN without an IAF
a = max(lo, fit_lo);  z = min(hi, fit_hi);
measured = hi > lo & lo >= f(1) - 1e-9 & hi <= f(end) + 1e-9 & z - a >= 1;   % as psd_measures
hit = measured & a < fit_lo + 2 - 1e-9;
names = string({o.bands.name});
for i = 1:size(hit, 1), col(L(1).rec(i)) = strjoin(names(hit(i, :)), ';'); end
cnt = sum(hit, 1);
b = find(cnt);
if ~isempty(b)
    msg = strjoin(compose('%s (%d of %d recordings)', names(b)', cnt(b)', repmat(size(hit, 1), numel(b), 1)), ', ');
end
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
