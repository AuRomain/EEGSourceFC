function S = check_settings(S)
%CHECK_SETTINGS Validate the settings and fill the optional ones.
%   S = check_settings(S) with S from pipeline_settings (or a copy edited
%   for another dataset). Optional fields that are missing get the
%   defaults below, so a settings file for a new dataset only needs the
%   fields that differ. Errors name the field to fix. A field that the
%   pipeline does not know (e.g. a misspelled name) is an error that
%   suggests the closest known name.

% optional fields and their defaults (S.dynamic: true where the folder
% of the brain network states is, false in a static-only pipeline)
has_states = ~isempty(stage_folder(S.root, 'Dynamic connectivity'));
def = struct( ...
    'toolboxes', fullfile(fileparts(S.root), 'Toolboxes'), ...
    'id_regexp', '', 'srate', [], 'resample', [], 'core_start', 0, 'core_duration', [], ...
    'edge_trim', 1, 'filter_transition', 'window', 'filter_span', 2, 'boundary_events', {{'boundary', 'New Segment'}}, 'bad_events', {{'bad*'}}, 'min_segment', 4, ...
    'channel_subset', {{}}, 'min_channels', 16, 'eloreta_lambda', 0.05, 'lcmv_lambda', '5%', 'lcmv_kappa', [], ...
    'wmne', struct('weightExp', 0.5, 'weightLimit', 10, 'SNR', 3), ...
    'window_cycles', 6, 'window_length', [], 'window_overlap', 0.5, 'aec_window', 6, 'aec_step', 0.5, 'measure_windows', struct(), ...
    'dwpli_epoch_cycles', 6, 'dwpli_epoch_length', [], ...
    'anatomy', '', 'anatomy_individual', false, 'anatomy_fallback', false, ...
    'save_scouts', true, 'orientation', 'normal', 'scout_function', 'mean', 'toolbox_check', 'warn', 'drop_flat_channels', true, 'parallel', true, 'max_workers', 8, ...
    'subs_info', '', 'id_column', 'ID', 'group_column', 'group', 'group_files', {{}}, ...
    'state_conn', 'wpli', 'NCs', 3:7, 'k_selection', 'stability', 'k_fixed', [], 'k_alternatives', 0.1, 'stability', struct('splits', 20, 'reps', 50), 'reps', 200, 'max_iterations', 1000, 'clustering', 'modkmeans', 'clustering_engine', 'fast', 'null_test', struct('n_surrogates', 19, 'reps', 20), 'continuum_test', struct('n_null', 19, 'nulls', {{'recording', 'pooled'}}, 'n_windows', 20000, 'reps', 20, 'splits', 10, 'stab_reps', 5), 'minTime', 5, 'minTime_s', [], 'minTime_dwpli', 0, ...
    'energy_entropy', false, 'network_table', '', 'language_file', '', 'seed', 20261015, ...
    'n_perm', 10000, 'alpha', 0.05, 'correlation_type', 'pearson', 'correlation_test', 'permutation', 'label_sets', {{'raw', 'smoothed'}}, ...
    'comparisons', struct('name', {}, 'A', {}, 'B', {}), 'covariates', {{}}, 'categorical_covariates', {{}}, ...
    'correlations', struct('name', {}, 'columns', {}, 'subjects', {}, 'covariates', {}, 'zero_as_missing', {}), ...
    'covariate_sets', struct('name', {}, 'correlations', {}, 'covariates', {}, 'columns', {}), ...
    'nbs', struct('thresh', 3.1, 'perms', 5000, 'size', 'extent'), ...
    'network_tests', struct('n_perm', 10000, 'n_boot', 1000, 'min_windows', 10), 'figures', true, 'state_free', true, 'spectral_profile', true, ...
    'window_measures', true, 'surrogate_variability', true, 'surrogate_variability_n', 5, 'sequence_features', true, 'sequence_windows', 'shortest', 'sequence_surrogates', 200, 'reliability', true, 'reliability_block', 10, 'duration_checks', true, 'duration_strata', '', ...
    'dynamic', has_states, 'static_measures', {{}}, 'static_matrices', {{}}, 'static_families', {{'global'}}, ...
    'enrichment', struct('n_perm', 10000));
required = {'root', 'results', 'data', 'file_pattern', 'montages', 'atlases', 'inv', 'bands', 'conn', 'groups'};
% other fields the pipeline reads: optional ones without a default here,
% a former name, and fields that check_settings and setup_toolboxes add
other = {'iaf', 'psd', 'individual_bands', 'bands_individual', 'run_name', 'stop_on_error', 'toolbox_paths', ...
         'openmeeg_python', 'modkmeans_engine', 'co', 'iaf_used', 'versions', 'signal_entropy_dir', 'n_null'};
known = [fieldnames(def)', required, other];
unknown = setdiff(fieldnames(S)', known);
if ~isempty(unknown)
    msg = strings(1, numel(unknown));
    for i = 1:numel(unknown)
        d = cellfun(@(k) edit_distance(lower(unknown{i}), lower(k)), known);
        [dm, j] = min(d);
        msg(i) = "S." + unknown{i};
        if dm <= max(2, floor(numel(unknown{i}) / 3)), msg(i) = msg(i) + " (did you mean S." + known{j} + "?)"; end
    end
    error('check_settings:unknown', 'Unknown setting(s): %s. See pipeline_settings.m for the fields', strjoin(msg, ', '));
end
if isfield(S, 'n_null')                             % small-worldness's null graphs: not computed
    warning('check_settings:n_null', 'S.n_null is no longer used: small-worldness is not computed');
    S = rmfield(S, 'n_null');
end
given = fieldnames(S);
for f = fieldnames(def)'
    if ~isfield(S, f{1}), S.(f{1}) = def.(f{1}); end
end
if ~S.dynamic                                       % static only: the defaults of the static-only template
    if ~ismember('save_scouts', given), S.save_scouts = false; end     % nothing reads the region time series
    if ~ismember('static_families', given), S.static_families = {'global', 'networks', 'regions', 'power'}; end
    if ~ismember('filter_transition', given), S.filter_transition = 'brainstorm'; end   % Brainstorm's filter
    if ~ismember('min_segment', given), S.min_segment = 6; end
end
S.max_workers = local_workers(S);
for f = required
    if ~isfield(S, f{1}) || isempty(S.(f{1}))
        error('check_settings:missing', 'Setting S.%s is required (see pipeline_settings.m)', f{1});
    end
end

% types
S.file_pattern = cellstr(S.file_pattern);
S.montages = cellstr(S.montages);
S.atlases = cellstr(S.atlases);
S.inv = lower(cellstr(S.inv));
S.conn = lower(cellstr(S.conn));
if ~(ischar(S.groups) || isstring(S.groups) || iscellstr(S.groups))   % e.g. [1 2] or {1, 2}
    error('check_settings:groups', ['S.groups: the group names as text, as they are written in the column ' ...
        'S.group_column (compared as text, case ignored), e.g. {''1'', ''2''} for a column holding 1 and 2']);
end
S.groups = cellstr(S.groups);
S.boundary_events = cellstr(S.boundary_events);
S.bad_events = cellstr(S.bad_events);
S.channel_subset = cellstr(S.channel_subset);
S.label_sets = cellstr(S.label_sets);
S.covariates = cellstr(S.covariates);
if isempty(S.group_files)                         % 'ctrl', 'psych' for two groups (fixed names),
    S.group_files = default_group_files(S.groups); % whatever the groups are; the group names otherwise
end
S.group_files = cellstr(S.group_files);
if numel(S.group_files) ~= numel(S.groups)
    error('check_settings:groups', 'S.group_files needs one name per group in S.groups');
end
if ~all(cellfun(@isvarname, S.group_files)) || numel(unique(S.group_files)) < numel(S.group_files)
    error('check_settings:groups', ['S.group_files: one distinct name per group, of letters, digits and _ ' ...
        '(they are part of the names of files and variables of the states), e.g. {''ctrl'', ''psych''}']);
end

% values
ok_inv = {'eloreta', 'lcmv', 'wmne'};
bad = setdiff(S.inv, ok_inv);
if ~isempty(bad), error('check_settings:inv', 'Unknown inverse method %s (use %s)', strjoin(bad, ', '), strjoin(ok_inv, ', ')); end
ok_conn = {'plv', 'wpli', 'aec', 'pli', 'ciplv', 'aec_orth', 'aec_symorth', 'plv_symorth', 'dwpli'};
bad = setdiff(S.conn, ok_conn);
if ~isempty(bad), error('check_settings:conn', 'Unknown connectivity measure %s (use %s)', strjoin(bad, ', '), strjoin(ok_conn, ', ')); end
if numel(unique(S.conn)) < numel(S.conn), error('check_settings:conn', 'S.conn lists a measure twice'); end
if ~isscalar(S.dynamic) || ~(islogical(S.dynamic) || isequal(S.dynamic, 0) || isequal(S.dynamic, 1))
    error('check_settings:dynamic', 'S.dynamic is true (windows and brain network states) or false (static only)');
end
S.dynamic = logical(S.dynamic);
if S.dynamic && ~has_states
    error('check_settings:dynamic', ['S.dynamic = true needs the brain network states, not part ' ...
        'of this copy of the pipeline (%s): set S.dynamic = false'], S.root);
end
if ~S.dynamic                                      % settings of the windows and states, set but not used
    dyn = {'NCs', 'k_selection', 'k_fixed', 'k_alternatives', ...
        'stability', 'reps', 'max_iterations', 'clustering', 'clustering_engine', 'null_test', 'continuum_test', ...
        'minTime', 'minTime_s', 'minTime_dwpli', 'energy_entropy', 'label_sets', 'network_tests', 'state_free', 'window_measures', ...
        'surrogate_variability', 'surrogate_variability_n', 'duration_checks', 'duration_strata', ...
        'sequence_features', 'sequence_windows', 'sequence_surrogates'};
    window_mean = isempty(S.static_matrices) || ...  % the default, or asked for: the windows are computed then
        ((iscellstr(S.static_matrices) || ischar(S.static_matrices) || isstring(S.static_matrices)) && ...
        any(strcmp(cellstr(S.static_matrices), 'window_mean')));
    if ~window_mean, dyn = [{'window_overlap', 'aec_window', 'aec_step', 'measure_windows'}, dyn]; end
    if ~window_mean && ~strcmpi(char(string(S.filter_transition)), 'window')    % the 'window' filter rule reads these
        dyn = [dyn, {'window_cycles', 'window_length'}];
    end
    if ~strcmpi(char(string(S.filter_transition)), 'window'), dyn = [dyn, {'filter_span'}]; end
    unused = dyn(cellfun(@(f) ~isequal(S.(f), def.(f)), dyn));
    if ~isempty(unused)
        warning('check_settings:not_used', ['S.dynamic = false (static only): %s not used (windows and ' ...
            'brain network states only)'], strjoin("S." + string(unused), ', '));
    end
end
S.co = find(strcmp(S.conn, lower(S.state_conn)));      % index into cmats{}
if isempty(S.co) && ~S.dynamic                     % static only: no states; S.state_conn is only the
    S.state_conn = S.conn{1};                      % measure of the language networks (load_static_matrices)
    S.co = 1;
end
if S.dynamic && ismember(lower(S.state_conn), {'aec', 'aec_orth', 'aec_symorth', 'dwpli'})
    warning('check_settings:signed', ['S.state_conn = ''%s'' is signed: modkmeans ignores the sign of a ' ...
        'pattern, and graph metrics use its positive part only'], S.state_conn);
end
if isempty(S.co)
    error('check_settings:state_conn', 'S.state_conn = ''%s'' is not in S.conn', S.state_conn);
end
% the static analyses of the statistics (static_analysis): measures, matrices, families
S.static_measures = lower(cellstr(S.static_measures));
if isempty(S.static_measures)                      % the state measure; static only: every measure
    S.static_measures = S.conn(S.co);
    if ~S.dynamic, S.static_measures = S.conn; end
end
bad = setdiff(S.static_measures, S.conn);
if ~isempty(bad)
    error('check_settings:static_measures', ['S.static_measures: %s not in S.conn (settings already ' ...
        'checked keep the measures of the S.conn they had: set S.static_measures = {} after changing S.conn)'], ...
        strjoin(bad, ', '));
end
S.static_matrices = cellstr(S.static_matrices);
if isempty(S.static_matrices)                     % the default: the mean of the windows,
    S.static_matrices = {'window_mean'};           % with or without the states
end
if ~all(ismember(S.static_matrices, {'full', 'window_mean'})) || numel(unique(S.static_matrices)) < numel(S.static_matrices)
    error('check_settings:static_matrices', ['S.static_matrices: ''full'' (each measure over the whole recording) ' ...
        'and/or ''window_mean'' (the mean of its windows: S.window_cycles or S.window_length, S.window_overlap; ' ...
        'AEC: S.aec_window, S.aec_step)']);
end
S.static_families = lower(cellstr(S.static_families));
bad = setdiff(S.static_families, {'global', 'networks', 'regions', 'power'});
if ~isempty(bad)
    error('check_settings:static_families', 'Unknown static feature family %s (use global, networks, regions, power)', ...
        strjoin(bad, ', '));
end
if isempty(S.static_families)                      % allowed: the edges only
    warning('check_settings:static_families', ['S.static_families is empty: no feature family is tested; the ' ...
        'static analyses test the edges only (NBS, network enrichment) and their reliability']);
end
% individual alpha frequency (used only by bands anchored on it, band_edges)
iaf_def = struct('bins', 2, 'bin_hz', 1, 'measure', 'paf', 'no_peak', 'exclude', ...
    'channels', {{'O1', 'Oz', 'O2', 'PO3', 'POz', 'PO4', 'PO7', 'PO8', 'P3', 'Pz', 'P4'}}, ...
    'search', [7 13], 'cmin', 3, 'sgf_frame', 11, 'sgf_order', 5, 'checks', {{'irasa', 'specparam'}}, ...
    'max_disagreement', 0.5, 'edge_check', true);
if ~isfield(S, 'iaf') || ~isstruct(S.iaf), S.iaf = struct(); end
for f = fieldnames(iaf_def)'
    if ~isfield(S.iaf, f{1}), S.iaf.(f{1}) = iaf_def.(f{1}); end
end
S.iaf.checks = cellstr(S.iaf.checks);
if ~ismember(S.iaf.measure, {'paf', 'cog'}) || ~ismember(S.iaf.no_peak, {'exclude', 'cog'}) || ...
        ~all(ismember(S.iaf.checks, {'irasa', 'specparam'})) || ~(S.iaf.bins >= 0 && S.iaf.bin_hz > 0)
    error('check_settings:iaf', ['S.iaf: measure ''paf''/''cog'', no_peak ''exclude''/''cog'', ' ...
        'checks among {''irasa'', ''specparam''}, bins >= 0, bin_hz > 0']);
end
if ~isfield(S, 'individual_bands'), S.individual_bands = false; end
if S.individual_bands                              % the switch of pipeline_settings.m
    if ~isfield(S, 'bands_individual') || isempty(S.bands_individual)
        error('check_settings:bands', 'S.individual_bands is true but S.bands_individual is not set');
    end
    S.bands = S.bands_individual;
end
S.iaf_used = false;
for b = 1:numel(S.bands)
    B = S.bands(b);
    if ~(isfield(B, 'name') && isfield(B, 'fmin') && isfield(B, 'fmax'))
        error('check_settings:bands', 'S.bands(%d) needs name, fmin and fmax', b);
    end
    [fmin, fmax, anchored] = band_edges(B, 10, S.iaf);   % anchored edges: at a typical IAF
    if ~(fmin > 0 && fmax > fmin)
        error('check_settings:bands', 'S.bands(%d) (%s) needs fmin > 0 and fmax > fmin (for an IAF of 10 Hz: %g-%g Hz)', ...
            b, B.name, fmin, fmax);
    end
    S.iaf_used = S.iaf_used || anchored;
end
S = check_psd(S);
if isempty(S.window_length) && ~(S.window_cycles > 0)
    error('check_settings:window', 'Set S.window_cycles > 0 or a fixed S.window_length');
end
% the epochs (trials) of the debiased wPLI (connectivity.dwpli_epoch)
if ~(is_num(S.dwpli_epoch_cycles) && S.dwpli_epoch_cycles > 0) || ...
        ~(isempty(S.dwpli_epoch_length) || (is_num(S.dwpli_epoch_length) && S.dwpli_epoch_length > 0))
    error('check_settings:dwpli', ['S.dwpli_epoch_cycles: cycles of the band centre (> 0, default 6); ' ...
        'S.dwpli_epoch_length: [] (cycles) or seconds (> 0)']);
end
if ismember('dwpli', S.conn)
    windows = S.dynamic || isempty(S.static_matrices) || any(strcmp(cellstr(S.static_matrices), 'window_mean'));
    % bands anchored on the IAF: epochs and windows in cycles of the band's
    % centre at IAF 10 Hz, the same in every recording (cycles_centre)
    for b = 1:numel(S.bands)
        [fmin, fmax] = band_edges(S.bands(b), 10, S.iaf);
        ep = connectivity_epoch(S, fmin, fmax);
        if isempty(ep.freqs)
            error('check_settings:dwpli', ['dwPLI in band %s (%g-%g Hz): epochs of %.3g s have no Fourier frequency ' ...
                'inside the band (one every %.3g Hz): longer epochs (S.dwpli_epoch_cycles or S.dwpli_epoch_length)'], ...
                S.bands(b).name, fmin, fmax, ep.seconds, 1 / ep.seconds);
        end
        o = measure_settings(S, 'dwpli', fmin, fmax);
        w = o.window_length;
        if isempty(w), w = o.window_cycles / (0.5 * (fmin + fmax)); end
        if windows && floor(round(w * 1e6) / round(ep.seconds * 1e6)) < ep.min_epochs
            error('check_settings:dwpli', ['dwPLI in band %s: its windows of %.3g s hold %d epochs of %.3g s; ' ...
                'the debiased wPLI needs at least %d per window: longer windows (S.measure_windows.dwpli) or ' ...
                'shorter epochs'], S.bands(b).name, w, floor(w / ep.seconds), ep.seconds, ep.min_epochs);
        end
    end
elseif ~isequal(S.dwpli_epoch_cycles, def.dwpli_epoch_cycles) || ~isequal(S.dwpli_epoch_length, def.dwpli_epoch_length)
    warning('check_settings:not_used', 'S.dwpli_epoch_cycles and S.dwpli_epoch_length are not used: no ''dwpli'' in S.conn');
end
% windows of single measures (S.measure_windows): a field
% named after a measure of S.conn holds the window settings it uses
% instead of the common ones
if isempty(S.measure_windows), S.measure_windows = struct(); end
if ~isstruct(S.measure_windows) || ~isscalar(S.measure_windows)
    error('check_settings:measure_windows', ['S.measure_windows: a struct with a field for each measure that has ' ...
        'windows of its own, e.g. struct(''wpli'', struct(''window_length'', 2))']);
end
w_rules = {
    'window_cycles',  @(v) is_num(v) && v > 0,                           '> 0'
    'window_length',  @(v) isempty(v) || (is_num(v) && v > 0),          '[] or > 0 (s)'
    'window_overlap', @(v) is_num(v) && v >= 0 && v < 1,                 'in [0, 1)'
    'aec_window',     @(v) is_num(v) && v > 0,                           '> 0 (s)'
    'aec_step',       @(v) is_num(v) && v > 0,                           '> 0 (s)'
    };
for m = fieldnames(S.measure_windows)'
    if ~ismember(m{1}, S.conn)
        error('check_settings:measure_windows', 'S.measure_windows.%s: %s is not a measure of S.conn', m{1}, m{1});
    end
    ok = {'window_cycles', 'window_length', 'window_overlap'};       % phase measures
    if startsWith(m{1}, 'aec'), ok = {'aec_window', 'aec_step'}; end
    mw = S.measure_windows.(m{1});
    if ~isstruct(mw) || ~isscalar(mw) || isempty(fieldnames(mw)) || ~all(ismember(fieldnames(mw), ok))
        error('check_settings:measure_windows', 'S.measure_windows.%s: a struct of %s', m{1}, strjoin(ok, ', '));
    end
    for i = find(isfield(mw, w_rules(:, 1)'))
        if ~w_rules{i, 2}(mw.(w_rules{i, 1}))
            error('check_settings:measure_windows', 'S.measure_windows.%s.%s must be %s', m{1}, w_rules{i, 1}, w_rules{i, 3});
        end
    end
end
% the windows of the clustered measure must fit in the shortest piece kept
% (S.min_segment), or recordings cut into short pieces have none
if S.dynamic
    o = measure_settings(S, S.state_conn);
    w_state = o.window_length;
    if startsWith(lower(char(S.state_conn)), 'aec'), w_state = o.aec_window; end
    if strcmpi(S.state_conn, 'dwpli')                  % its own windows: the longest over the bands
        w_state = 0;
        for b = 1:numel(S.bands)
            [fmin, fmax] = band_edges(S.bands(b), 10, S.iaf);
            o = measure_settings(S, 'dwpli', fmin, fmax);
            wb = o.window_length;
            if isempty(wb), wb = o.window_cycles / (0.5 * (fmin + fmax)); end
            w_state = max(w_state, wb);
        end
    end
    if ~isempty(w_state) && w_state > S.min_segment
        warning('check_settings:window_vs_segment', ['The windows of %s (%g s) are longer than S.min_segment ' ...
            '(%g s): pieces of the recordings between them give no window, and a recording made only of such ' ...
            'pieces is left out of the states and their statistics. Set S.min_segment >= %g if your recordings are cut often'], ...
            S.state_conn, w_state, S.min_segment, w_state);
    end
end
% static only with the window mean: the windows of every measure should
% fit in the shortest piece kept, or such pieces add nothing to the mean
% and a recording made only of them has none (NaN; run_pipeline warns)
if ~S.dynamic && ismember('window_mean', S.static_matrices)
    w_max = zeros(1, numel(S.conn));
    for k = 1:numel(S.conn)
        o = measure_settings(S, S.conn{k});
        if startsWith(S.conn{k}, 'aec')
            w_max(k) = o.aec_window;
        elseif strcmp(S.conn{k}, 'dwpli')          % its own windows, band by band
            for b = 1:numel(S.bands)
                [fmin, fmax] = band_edges(S.bands(b), 10, S.iaf);
                ob = measure_settings(S, 'dwpli', fmin, fmax);
                wb = ob.window_length;
                if isempty(wb), wb = ob.window_cycles / (0.5 * (fmin + fmax)); end
                w_max(k) = max(w_max(k), wb);
            end
        elseif ~isempty(o.window_length)
            w_max(k) = o.window_length;
        else                                      % cycles: the longest window is that of the lowest band
            for b = 1:numel(S.bands)
                [fmin, fmax] = band_edges(S.bands(b), 10, S.iaf);
                w_max(k) = max(w_max(k), o.window_cycles / (0.5 * (fmin + fmax)));
            end
        end
    end
    long = w_max > S.min_segment;
    if any(long)
        warning('check_settings:window_vs_segment', ['The windows of %s (up to %g s) are longer than ' ...
            'S.min_segment (%g s): pieces of the recordings between them add nothing to the window mean, and a ' ...
            'recording made only of such pieces has none (NaN). Set S.min_segment >= %g if your recordings are cut often'], ...
            strjoin(S.conn(long), ', '), max(w_max(long)), S.min_segment, max(w_max(long)));
    end
end
if ~ismember(S.toolbox_check, {'warn', 'error', 'off'})
    error('check_settings:toolbox_check', 'S.toolbox_check: ''warn'' (default), ''error'' or ''off''');
end
if ~ismember(S.orientation, {'normal', 'free'})
    error('check_settings:orientation', 'S.orientation: ''normal'' (default) or ''free''');
end
S.scout_function = lower(char(string(S.scout_function)));
if ~ismember(S.scout_function, {'mean', 'pca', 'max', 'power', 'rms'})
    error('check_settings:scout_function', ['S.scout_function: ''mean'' (default), ''pca'', ''max'', ''power'' ' ...
        'or ''rms'' (Brainstorm''s scout functions; get_scouts_timeseries)']);
end
if strcmp(S.orientation, 'free') && ~strcmp(S.scout_function, 'pca')
    error('check_settings:scout_function', ['S.orientation = ''free'' takes the first principal component of ' ...
        'each region''s sources: set S.scout_function = ''pca''']);
end
phase = intersect(S.conn, {'plv', 'wpli', 'pli', 'ciplv', 'dwpli', 'plv_symorth'}, 'stable');
if ismember(S.scout_function, {'power', 'rms'}) && ~isempty(phase)
    what = struct('power', 'mean of its squared sources', 'rms', 'root mean square of its sources');
    warning('check_settings:scout_function', ['S.scout_function = ''%s'': each region''s signal is the %s, ' ...
        'without sign, and its band-passed part is the band''s share of the fluctuations of the region''s power, ' ...
        'not the region''s oscillation; %s do not measure the phase of the band''s oscillation'], ...
        S.scout_function, what.(S.scout_function), strjoin(upper(phase), ', '));
end
if ~(S.window_overlap >= 0 && S.window_overlap < 1)
    error('check_settings:window', 'S.window_overlap must be in [0, 1)');
end
if ~(ischar(S.clustering) || isstring(S.clustering)) || ~isvarname(char(S.clustering))
    error('check_settings:clustering', 'S.clustering: the name of a method of clustering_method.m (''modkmeans'')');
end
S.clustering = char(S.clustering);
if ~(isnumeric(S.filter_span) && isscalar(S.filter_span) && S.filter_span > 0)
    error('check_settings:filter_span', 'S.filter_span: the number of windows a filter may spread over (e.g. 2)');
end
if ~isstruct(S.null_test) || ~all(isfield(S.null_test, {'n_surrogates', 'reps'}))
    error('check_settings:null_test', 'S.null_test: struct(''n_surrogates'', 19, ''reps'', 20)');
end
if S.null_test.n_surrogates > 0 && S.null_test.n_surrogates < 19
    warning('check_settings:null_test', ['With %d surrogates the p-value cannot go below %.2f: ' ...
        '19 is the minimum for p = 0.05'], S.null_test.n_surrogates, 1 / (S.null_test.n_surrogates + 1));
end
if S.dynamic && S.null_test.n_surrogates > 0 && ~S.save_scouts
    error('check_settings:null_test', 'The surrogate test (S.null_test) needs the region time series: S.save_scouts = true');
end
if isfield(S, 'modkmeans_engine'), S.clustering_engine = S.modkmeans_engine; end   % a former name, still accepted
if ~ismember(S.clustering_engine, {'fast', 'toolbox'})
    error('check_settings:clustering_engine', 'S.clustering_engine: ''fast'' or ''toolbox''');
end
if ~ismember(S.k_selection, {'stability', 'median_vote', 'meta_criterion'})
    error('check_settings:k_selection', 'S.k_selection: ''stability'', ''median_vote'' or ''meta_criterion''');
end
if ~isstruct(S.stability) || ~all(isfield(S.stability, {'splits', 'reps'})) || ~is_int(S.stability.splits) ...
        || S.stability.splits < 1 || ~is_int(S.stability.reps) || S.stability.reps < 1
    error('check_settings:stability', 'S.stability: struct(''splits'', 20, ''reps'', 50), integers >= 1');
end
if ~isnumeric(S.NCs) || any(S.NCs < 2) || numel(S.NCs) < 3 || any(S.NCs ~= round(S.NCs)) || any(diff(S.NCs(:)) ~= 1)
    error('check_settings:NCs', ['S.NCs needs at least three consecutive numbers of states in increasing order, ' ...
        'all >= 2 (e.g. 3:7; the criteria compare neighbours)']);
end
if ~isempty(S.k_fixed) && ~(is_int(S.k_fixed) && ismember(S.k_fixed, S.NCs))
    error('check_settings:k_fixed', 'S.k_fixed: [] (chosen by S.k_selection) or one of the numbers of states tested (S.NCs = %s)', ...
        mat2str(S.NCs));
end
if ~is_num(S.k_alternatives) || S.k_alternatives < 0 || S.k_alternatives >= 1
    error('check_settings:k_alternatives', 'S.k_alternatives: a stability margin from 0 (no flag) to below 1 (0.1)');
end
% the head model (get_anatomy): the shipped ICBM152, a folder, or one folder per subject
if ~(ischar(S.anatomy) || (isstring(S.anatomy) && isscalar(S.anatomy)))
    error('check_settings:anatomy', 'S.anatomy is '''' (the shipped ICBM152) or a folder (tools/make_anatomy.m)');
end
S.anatomy = char(S.anatomy);
for f = {'anatomy_individual', 'anatomy_fallback'}
    if ~isscalar(S.(f{1})) || ~(islogical(S.(f{1})) || isequal(S.(f{1}), 0) || isequal(S.(f{1}), 1))
        error(['check_settings:' f{1}], 'S.%s is true or false', f{1});
    end
    S.(f{1}) = logical(S.(f{1}));
end
if ~isempty(S.anatomy) && ~isfolder(S.anatomy)
    error('check_settings:anatomy', 'S.anatomy: folder %s not found (tools/make_anatomy.m makes it)', S.anatomy);
end
if S.anatomy_individual && isempty(S.anatomy)
    error('check_settings:anatomy', 'S.anatomy_individual = true needs S.anatomy, the folder of the subjects'' folders');
end
if S.anatomy_individual && any(cellfun(@isfile, S.atlases))
    error('check_settings:anatomy', ['S.anatomy_individual = true: a scout file is defined on one cortex; give atlases ' ...
        'by name (''Desikan_Killiany_68'', ''Destrieux'', ...: the atlases of each subject''s cortex)']);
end
if S.anatomy_fallback && ~S.anatomy_individual
    warning('check_settings:not_used', 'S.anatomy_fallback is used only with S.anatomy_individual = true');
end
% numbers in their ranges: name, test, what is expected
rules = {
    'core_start',     @(v) is_num(v) && v >= 0,                          '>= 0 (s)'
    'core_duration',  @(v) isempty(v) || (is_num(v) && v > 0),          '[] or > 0 (s)'
    'edge_trim',      @(v) is_num(v) && v >= 0,                          '>= 0 (s)'
    'min_segment',    @(v) is_num(v) && v > 0,                           '> 0 (s)'
    'srate',          @(v) isempty(v) || (is_num(v) && v > 0),          '[] or > 0 (Hz)'
    'resample',       @(v) isempty(v) || (is_num(v) && v > 0),          '[] or > 0 (Hz)'
    'min_channels',   @(v) is_int(v) && v >= 1,                          'an integer >= 1'
    'window_cycles',  @(v) is_num(v) && v > 0,                           '> 0'
    'window_length',  @(v) isempty(v) || (is_num(v) && v > 0),          '[] or > 0 (s)'
    'aec_window',     @(v) is_num(v) && v > 0,                           '> 0 (s)'
    'aec_step',       @(v) is_num(v) && v > 0,                           '> 0 (s)'
    'reps',           @(v) is_int(v) && v >= 1,                          'an integer >= 1'
    'max_iterations', @(v) is_int(v) && v >= 1,                          'an integer >= 1'
    'minTime',        @(v) is_int(v) && v >= 0,                          'an integer >= 0 (windows)'
    'minTime_s',      @(v) isempty(v) || (is_num(v) && v >= 0),         '[] or >= 0 (s)'
    'minTime_dwpli',  @(v) is_int(v) && v >= 0,                          'an integer >= 0 (windows)'
    'seed',           @(v) is_int(v) && v >= 0,                          'an integer >= 0'
    'n_perm',         @(v) is_int(v) && v >= 1,                          'an integer >= 1'
    'alpha',          @(v) is_num(v) && v > 0 && v < 1,                  'in (0, 1)'
    'max_workers',    @(v) is_int(v) && v >= 1,                          'an integer >= 1'
    'filter_span',    @(v) is_num(v) && v > 0,                           '> 0'
    'lcmv_kappa',     @(v) isempty(v) || strcmp(string(v), "rank") || (is_int(v) && v >= 1), '[], ''rank'' or an integer >= 1'
    };
for i = 1:size(rules, 1)
    if ~rules{i, 2}(S.(rules{i, 1}))
        error('check_settings:range', 'S.%s must be %s', rules{i, 1}, rules{i, 3});
    end
end
if ~isempty(S.lcmv_kappa) && ~isnumeric(S.lcmv_kappa), S.lcmv_kappa = 'rank'; end   % one spelling in the key
ct = S.continuum_test;
ct_fields = {'n_null', 'nulls', 'n_windows', 'reps', 'splits', 'stab_reps'};
if ~isstruct(ct) || ~all(isfield(ct, ct_fields)) || ~is_int(ct.n_null) || ct.n_null < 0 || ~is_int(ct.n_windows) || ...
        ct.n_windows < 100 || ~is_int(ct.reps) || ct.reps < 1 || ~is_int(ct.splits) || ct.splits < 1 || ...
        ~is_int(ct.stab_reps) || ct.stab_reps < 1 || isempty(ct.nulls) || ~all(ismember(cellstr(ct.nulls), {'recording', 'pooled'}))
    error('check_settings:continuum_test', ['S.continuum_test: struct(''n_null'', 19, ''nulls'', {{''recording'', ''pooled''}}, ' ...
        '''n_windows'', 20000, ''reps'', 20, ''splits'', 10, ''stab_reps'', 5); n_null = 0: no test']);
end
S.continuum_test.nulls = cellstr(ct.nulls);
if ~islogical(S.state_free) && ~isequal(S.state_free, 0) && ~isequal(S.state_free, 1)
    error('check_settings:state_free', 'S.state_free is true or false');
end
if ~islogical(S.spectral_profile) && ~isequal(S.spectral_profile, 0) && ~isequal(S.spectral_profile, 1)
    error('check_settings:spectral_profile', 'S.spectral_profile is true or false');
end
if ~((ischar(S.correlation_type) || (isstring(S.correlation_type) && isscalar(S.correlation_type))) && ...
        ismember(lower(char(S.correlation_type)), {'pearson', 'spearman'}))
    error('check_settings:correlation_type', 'S.correlation_type is ''pearson'' or ''spearman''');
end
S.correlation_type = lower(char(S.correlation_type));
if ~((ischar(S.correlation_test) || (isstring(S.correlation_test) && isscalar(S.correlation_test))) && ...
        ismember(lower(char(S.correlation_test)), {'permutation', 'wild_bootstrap'}))
    error('check_settings:correlation_test', 'S.correlation_test is ''permutation'' or ''wild_bootstrap''');
end
S.correlation_test = lower(char(S.correlation_test));
for f = {'window_measures', 'surrogate_variability', 'reliability', 'duration_checks', 'sequence_features'}
    if ~islogical(S.(f{1})) && ~isequal(S.(f{1}), 0) && ~isequal(S.(f{1}), 1)
        error(['check_settings:' f{1}], 'S.%s is true or false', f{1});
    end
end
sw = S.sequence_windows;
if ~((is_int(sw) && sw >= 10) || ((ischar(sw) || (isstring(sw) && isscalar(sw))) && ismember(lower(char(sw)), {'shortest', 'all'})))
    error('check_settings:sequence_windows', ['S.sequence_windows: ''shortest'' (the shortest recording''s ' ...
        'number of windows), ''all'', or a number of windows (10 or more)']);
end
if ~is_int(sw), S.sequence_windows = lower(char(sw)); end
if ~is_int(S.sequence_surrogates) || S.sequence_surrogates < 19
    error('check_settings:sequence_surrogates', ['S.sequence_surrogates is the number of first-order Markov surrogates of each recording ' ...
        '(the syntax test, the null G of memory and non-stationarity), 19 or more']);
end
if ~is_int(S.surrogate_variability_n) || S.surrogate_variability_n < 1
    error('check_settings:surrogate_variability_n', 'S.surrogate_variability_n is the number of surrogates per recording, 1 or more');
end
if ~isnumeric(S.reliability_block) || ~isscalar(S.reliability_block) || ~(S.reliability_block > 0)
    error('check_settings:reliability_block', 'S.reliability_block is a block length in seconds, above 0');
end
if ~(ischar(S.duration_strata) || (isstring(S.duration_strata) && isscalar(S.duration_strata)))
    error('check_settings:duration_strata', 'S.duration_strata is a column of the subject table, or '''' for none');
end
if ~is_int(S.null_test.n_surrogates) || S.null_test.n_surrogates < 0 || ~is_int(S.null_test.reps) || S.null_test.reps < 1
    error('check_settings:null_test', 'S.null_test: n_surrogates an integer >= 0, reps an integer >= 1');
end
bad = setdiff(S.label_sets, {'raw', 'smoothed'});
if isempty(S.label_sets) || ~isempty(bad)
    error('check_settings:label_sets', 'S.label_sets: ''raw'' and/or ''smoothed''');
end
% the smoothing of dwPLI states: its windows (60 cycles,
% each over at least 10 epochs) are ten times those of the other phase
% measures, so S.minTime's 5 windows would reject segments ten times
% longer (up to 180 cycles); its own
% threshold S.minTime_dwpli (default 0: none), and S.minTime / S.minTime_s
% for the other measures
if S.dynamic && strcmpi(S.state_conn, 'dwpli')
    if ~isequal(S.minTime, def.minTime) || ~isempty(S.minTime_s)
        warning('check_settings:not_used', ['S.minTime and S.minTime_s are not used with dwPLI states: ' ...
            'S.minTime_dwpli (%d windows) instead'], S.minTime_dwpli);
    end
    if S.minTime_dwpli == 0 && all(ismember({'raw', 'smoothed'}, S.label_sets))
        S.label_sets = {'raw'};
        warning('check_settings:smoothing', ['dwPLI states are not smoothed (S.minTime_dwpli = 0): their ' ...
            'smoothed labels are the raw ones, so the statistics test them once, as ''raw''']);
    end
elseif S.dynamic && ~isequal(S.minTime_dwpli, def.minTime_dwpli)
    warning('check_settings:not_used', 'S.minTime_dwpli is not used: the states are of %s, not dwPLI', S.state_conn);
end
names = {S.bands.name};                            % band names are parts of file names
if ~all(cellfun(@(n) (ischar(n) || isstring(n)) && ~isempty(regexp(char(n), '^[A-Za-z][A-Za-z0-9_]*$', 'once')), names)) ...
        || numel(unique(cellfun(@char, names, 'UniformOutput', false))) < numel(names)
    error('check_settings:bands', 'S.bands: names must be unique, start with a letter and use only letters, digits and _');
end
if ~isstruct(S.nbs) || ~all(isfield(S.nbs, {'thresh', 'perms', 'size'})) || ~(isnumeric(S.nbs.thresh) && isscalar(S.nbs.thresh) ...
        && S.nbs.thresh > 0) || ~(isnumeric(S.nbs.perms) && isscalar(S.nbs.perms) && S.nbs.perms >= 1) ...
        || ~ismember(lower(char(S.nbs.size)), {'extent', 'intensity'})
    error('check_settings:nbs', ['S.nbs: struct(''thresh'', 3.1, ''perms'', 5000, ''size'', ''extent''): ' ...
        'threshold > 0, at least 1 permutation, size ''extent'' or ''intensity''']);
end
S.nbs.size = lower(char(S.nbs.size));             % NBS 1.2 compares 'Intensity' exactly
if ~isstruct(S.network_tests) || ~all(isfield(S.network_tests, {'n_perm', 'n_boot', 'min_windows'})) ...
        || ~is_int(S.network_tests.n_perm) || S.network_tests.n_perm < 1 || ~is_int(S.network_tests.n_boot) ...
        || S.network_tests.n_boot < 1 || ~is_int(S.network_tests.min_windows) || S.network_tests.min_windows < 1
    error('check_settings:network_tests', ['S.network_tests: struct(''n_perm'', 10000, ''n_boot'', 1000, ' ...
        '''min_windows'', 10), integers >= 1']);
end
if ~isstruct(S.enrichment) || ~isfield(S.enrichment, 'n_perm') || ~is_int(S.enrichment.n_perm) || S.enrichment.n_perm < 1
    error('check_settings:enrichment', 'S.enrichment: struct(''n_perm'', 10000), an integer >= 1');
end
if isfield(S, 'run_name') && ~isempty(S.run_name) && (~(ischar(S.run_name) || isstring(S.run_name)) || ...
        isempty(regexp(char(S.run_name), '^[A-Za-z0-9][A-Za-z0-9_.-]*$', 'once')))
    error('check_settings:run_name', 'S.run_name: a folder name (letters, digits, _ . -)');
end
if isfield(S, 'run_name') && strcmpi(char(S.run_name), 'spectra')   % results/PSD/spectra keeps the spectra
    error('check_settings:run_name', ['S.run_name ''spectra'': the PSD analysis keeps its spectra in ' ...
        'results/PSD/spectra, and each run empties results/PSD/<run name>: choose another name']);
end
for k = 1:numel(S.comparisons)
    if ~all(isfield(S.comparisons(k), {'name', 'A', 'B'}))
        error('check_settings:comparisons', 'S.comparisons(%d) needs name, A and B', k);
    end
end
for k = 1:numel(S.correlations)
    if ~isfield(S.correlations(k), 'zero_as_missing') || isempty(S.correlations(k).zero_as_missing)
        S.correlations(k).zero_as_missing = false;
    end
    S.correlations(k).covariates = cellstr(S.correlations(k).covariates);
end
for k = 1:numel(S.covariate_sets)
    S.covariate_sets(k).covariates = cellstr(S.covariate_sets(k).covariates);
end
S.categorical_covariates = cellstr(S.categorical_covariates);
named = S.covariates;                              % the covariates named anywhere
for k = 1:numel(S.correlations), named = [named, S.correlations(k).covariates]; end %#ok<AGROW>
for k = 1:numel(S.covariate_sets), named = [named, S.covariate_sets(k).covariates]; end %#ok<AGROW>
bad = setdiff(S.categorical_covariates, named);
if ~isempty(bad)
    error('check_settings:categorical_covariates', ['S.categorical_covariates names covariates of S.covariates, ' ...
        'S.correlations or S.covariate_sets; not among them: %s'], strjoin(bad, ', '));
end
end


function S = check_psd(S)
% the PSD analysis (run_psd_analysis, stage 2): S.psd, its missing fields
% filled with the defaults below (pipeline_settings.m, section "PSD analysis")
def = struct('run', true, 'levels', {{'sensor', 'source'}}, ...
    'bands', struct('name', {'delta', 'theta', 'alpha', 'alpha1', 'alpha2', 'beta', 'beta1', 'beta2'}, ...
                    'fmin', {1, 4, 7, 7, 10, 13, 13, 20}, 'fmax', {4, 7, 13, 10, 13, 30, 20, 30}), ...
    'total', [1 30], 'measures', {{'absolute', 'relative', 'periodic', 'aperiodic', 'alpha_frequency'}}, ...
    'segment', 4, 'overlap', 0.5, 'freq_range', [1 45], 'fit_range', [2 40], ...
    'cluster', struct('run', true, 'alpha', 0.05), 'n_perm', [], 'region_spectra', 'signal', ...
    'vertex_maps', true, 'power_line', '50');  % (a struct array value is one field: not split)
if ~isfield(S, 'psd') || isempty(S.psd), S.psd = struct(); end
if ~isstruct(S.psd) || ~isscalar(S.psd)
    error('check_settings:psd', 'S.psd is a struct (pipeline_settings.m, section "PSD analysis")');
end
bad = setdiff(fieldnames(S.psd), fieldnames(def));
if ~isempty(bad)
    error('check_settings:psd', 'Unknown field(s) of S.psd: %s (fields: %s)', strjoin(bad, ', '), ...
        strjoin(fieldnames(def), ', '));
end
for f = fieldnames(def)'
    if ~isfield(S.psd, f{1}), S.psd.(f{1}) = def.(f{1}); end
end
o = S.psd;
if ~isscalar(o.run) || ~(islogical(o.run) || isequal(o.run, 0) || isequal(o.run, 1))
    error('check_settings:psd', 'S.psd.run is true (run_full_pipeline runs the PSD analysis after stage 1) or false');
end
o.run = logical(o.run);
o.levels = lower(cellstr(o.levels));
if isempty(o.levels) || ~all(ismember(o.levels, {'sensor', 'source'}))
    error('check_settings:psd', 'S.psd.levels: ''sensor'' (channels) and/or ''source'' (regions)');
end
o.measures = lower(cellstr(o.measures));
known = {'absolute', 'relative', 'periodic', 'aperiodic', 'alpha_frequency'};
if isempty(o.measures) || ~all(ismember(o.measures, known))
    error('check_settings:psd', 'S.psd.measures: among %s', strjoin(known, ', '));
end
range_ok = @(v) isnumeric(v) && numel(v) == 2 && all(isfinite(v)) && v(1) > 0 && v(2) > v(1);
if ~(is_num(o.segment) && o.segment > 0) || ~(is_num(o.overlap) && o.overlap >= 0 && o.overlap < 1)
    error('check_settings:psd', 'S.psd.segment > 0 (s) and S.psd.overlap in [0, 1)');
end
if ~range_ok(o.freq_range)
    error('check_settings:psd', 'S.psd.freq_range: [fmin fmax] Hz, 0 < fmin < fmax');
end
if ~range_ok(o.fit_range) || o.fit_range(1) < o.freq_range(1) || o.fit_range(2) > o.freq_range(2) ...
        || diff(o.fit_range) < 4
    error('check_settings:psd', 'S.psd.fit_range: [fmin fmax] Hz inside S.psd.freq_range, at least 4 Hz wide');
end
if ~range_ok(o.total) || o.total(1) < o.freq_range(1) - 1e-9 || o.total(2) > o.freq_range(2) + 1e-9
    error('check_settings:psd', 'S.psd.total (relative power: %% of this range): [fmin fmax] Hz inside S.psd.freq_range');
end
B = o.bands;
if ~isstruct(B) || isempty(B) || ~all(isfield(B, {'name', 'fmin', 'fmax'}))
    error('check_settings:psd', 'S.psd.bands: struct(''name'', ..., ''fmin'', ..., ''fmax'', ...), at least one band');
end
names = {B.name};
if ~all(cellfun(@(n) (ischar(n) || isstring(n)) && ~isempty(regexp(char(n), '^[A-Za-z][A-Za-z0-9_]*$', 'once')), names)) ...
        || numel(unique(cellfun(@char, names, 'UniformOutput', false))) < numel(names)
    error('check_settings:psd', 'S.psd.bands: names must be unique, start with a letter and use only letters, digits and _');
end
for b = 1:numel(B)
    [lo, hi] = band_edges(B(b), 10, S.iaf);          % anchored edges: at a typical IAF
    if ~(lo > 0 && hi > lo) || lo < o.freq_range(1) - 1e-9 || hi > o.freq_range(2) + 1e-9
        error('check_settings:psd', ['S.psd.bands(%d) (%s): fmin > 0, fmax > fmin, inside S.psd.freq_range ' ...
            '(%g-%g Hz; for an IAF of 10 Hz: %g-%g Hz)'], b, char(B(b).name), o.freq_range(1), o.freq_range(2), lo, hi);
    end
end
c = o.cluster;
if ~isstruct(c) || ~all(isfield(c, {'run', 'alpha'})) || ~(islogical(c.run) || isequal(c.run, 0) || ...
        isequal(c.run, 1)) || ~(is_num(c.alpha) && c.alpha > 0 && c.alpha < 1)
    error('check_settings:psd', 'S.psd.cluster: struct(''run'', true, ''alpha'', 0.05) (cluster-forming p, in (0, 1))');
end
o.cluster.run = logical(c.run);
if ~(isempty(o.n_perm) || (is_int(o.n_perm) && o.n_perm >= 1))
    error('check_settings:psd', 'S.psd.n_perm: [] (S.n_perm) or an integer >= 1, the permutations of the PSD statistics');
end
if ~(ischar(o.region_spectra) || (isstring(o.region_spectra) && isscalar(o.region_spectra))) || ...
        ~ismember(char(o.region_spectra), {'signal', 'vertex_mean'})
    error('check_settings:psd', ['S.psd.region_spectra: ''signal'' (the PSD of each region''s signal, as stage 1 ' ...
        'makes it) or ''vertex_mean'' (the mean of its sources'' PSDs)']);
end
o.region_spectra = char(o.region_spectra);
if ~(islogical(o.vertex_maps) || isequal(o.vertex_maps, 0) || isequal(o.vertex_maps, 1)) || ~isscalar(o.vertex_maps)
    error('check_settings:psd', 'S.psd.vertex_maps is true (maps of every source''s band power) or false');
end
o.vertex_maps = logical(o.vertex_maps);
pl = o.power_line;
if isnumeric(pl) && isscalar(pl) && ismember(pl, [50 60]), pl = sprintf('%d', pl); end
if ~(ischar(pl) || (isstring(pl) && isscalar(pl))) || ~ismember(lower(char(pl)), {'50', '60', 'none'})
    error('check_settings:psd', ['S.psd.power_line: ''50'' or ''60'' (Hz of the mains, whose harmonics ' ...
        'the specparam fit leaves out) or ''none''']);
end
o.power_line = lower(char(pl));
if ismember('relative', o.measures)                  % relative power outside the total: a ratio, not a share
    out = {};
    for b = 1:numel(B)
        [lo, hi] = band_edges(B(b), 10, S.iaf);
        if lo < o.total(1) - 1e-9 || hi > o.total(2) + 1e-9, out{end + 1} = char(B(b).name); end %#ok<AGROW>
    end
    if ~isempty(out)
        warning('check_settings:psd_total', ['S.psd.bands outside S.psd.total (%g-%g Hz): %s. Their relative ' ...
            'power is a ratio to the power over S.psd.total, not a share of it (it can exceed 100%%); ' ...
            'widen S.psd.total to cover them'], o.total(1), o.total(2), strjoin(out, ', '));
    end
end
S.psd = o;
end

function tf = is_num(v)
tf = isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v);
end

function o = measure_settings(S, measure, fmin, fmax)
% the window settings of one measure: the common ones, with those of
% S.measure_windows.(measure) instead (window_cycles alone: cycles), as
% get_all_connectivity_mats; dwPLI (with its band) starts from its own:
% 60 cycles of the band centre, half overlapping, or 10 epochs when that
% holds fewer
o = S;
if strcmp(measure, 'dwpli') && nargin > 2
    ep = connectivity_epoch(S, fmin, fmax);
    o.window_length = max(60 / (0.5 * (fmin + fmax)), ep.min_epochs * ep.seconds);
    o.window_overlap = 0.5;
end
if ~isfield(S.measure_windows, measure), return; end
mw = S.measure_windows.(measure);
if isfield(mw, 'window_cycles') && ~isfield(mw, 'window_length'), o.window_length = []; end
for f = fieldnames(mw)'
    o.(f{1}) = mw.(f{1});
end
end

function tf = is_int(v)
tf = is_num(v) && v == round(v);
end

function e = connectivity_epoch(S, fmin, fmax)
% dwPLI's epochs in a band (connectivity.dwpli_epoch, at 1000 Hz: the
% sampling rate is not known yet; the length in seconds and the Fourier
% frequencies differ by less than a sample)
if isempty(which('connectivity.dwpli_epoch'))       % stage 1's folder, not yet on the path
    addpath(fullfile(S.root, '1- EEG source connectivity'));
end
e = connectivity.dwpli_epoch(1000, fmin, fmax, S);
end

function d = edit_distance(a, b)
% Levenshtein distance between two texts
m = numel(a);  n = numel(b);
D = zeros(m + 1, n + 1);
D(:, 1) = 0:m;  D(1, :) = 0:n;
for i = 1:m
    for j = 1:n
        D(i + 1, j + 1) = min([D(i, j + 1) + 1, D(i + 1, j) + 1, D(i, j) + (a(i) ~= b(j))]);
    end
end
d = D(end, end);
end

function n = local_workers(S)
% S.max_workers, at most the workers of this computer's local cluster (parpool
% errors when asked for more; 8 is the default). Looked up once per MATLAB session.
% Without the Parallel Computing Toolbox (getCurrentTask and parcluster are its
% functions) parfor runs as a for loop: S.max_workers as set.
persistent n_local has_pct
n = S.max_workers;
if isempty(has_pct), has_pct = ~isempty(ver('parallel')); end
if ~S.parallel || ~has_pct || ~isempty(getCurrentTask()), return; end
if isempty(n_local)
    n_local = Inf;
    try
        n_local = parcluster('Processes').NumWorkers;
    catch
    end
    if n > n_local
        fprintf('S.max_workers = %d, but this computer runs at most %d workers: %d used\n', n, n_local, n_local);
    end
end
n = min(n, n_local);
end
