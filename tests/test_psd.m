function tests = test_psd
%TEST_PSD The PSD analysis (stage 2, run_psd_analysis): spectra, measures, clusters, a planted difference.
%   From the pipeline folder:  results = runtests('tests/test_psd.m')
%   - psd_welch against signals of known power (a sine, white noise) and
%     pwelch;
%   - psd_specparam and psd_measures against spectra of known exponent,
%     offset and alpha peak; bands anchored on the IAF;
%   - cluster_perm_test: a planted cluster found, none under the null, the
%     signs apart, columns not tested left out; the same clusters, sums
%     and p as FieldTrip's own clusterstat (the installed FieldTrip's
%     private/clusterstat.m, findcluster.m and combineClusters.m, copied
%     to a temporary folder), on the same permutations;
%   - the neighbours of channels and regions (Brainstorm's VertConn and the
%     triangles of the cortex give the same region pairs);
%   - psd_source_signals equal to stage 1's region time series (src_data of
%     Scouts_*.mat) for a recording of static_test_recordings;
%   - psd_cross_spectra: its diagonal equal to psd_welch, the whole matrix
%     to MATLAB's cpsd of each piece;
%   - psd_vertex_spectra: every source's PSD equal to psd_welch of its
%     signal (eLORETA, wMNE, LCMV; normal and free orientation), eLORETA's
%     equal to FieldTrip's ft_inverse_eloreta on the cross-spectra and to
%     Brainstorm's PSD of a source map and of a scout ('after'), each
%     region's signal power at most the mean of its sources' power;
%   - the settings S.psd (defaults, errors);
%   - the whole step on simulated recordings (8 controls, 8 patients, the
%     64 channels of static_test_recordings' cap, 60 s at 256 Hz) whose
%     patients have a stronger alpha rhythm in the occipital regions and a
%     score that follows it: the posterior channels and occipital regions
%     must differ (tests or clusters), more than the frontal ones; the
%     score must correlate there; beta2 (20-30 Hz), not planted, must not;
%     the region check (region_power_check.csv) and the sources' maps
%     (vertex_values_<inverse>.mat: the largest g in the occipital
%     sources, g as hedges_g); S.psd.region_spectra = 'vertex_mean' from
%     the same kept spectra.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'simulation'));
addpath(stage_folder(root, 'PSD analysis'));
S = static_test_settings();
S.parallel = false;
tc.TestData.S = setup_toolboxes(S);
tc.TestData.root = root;
end

function test_welch_scaling(tc)
fs = 250;  t = (0:fs * 120 - 1) / fs;
rs = RandStream('mrg32k3a', 'Seed', 1);
x = [2 * sin(2 * pi * 10 * t) + randn(rs, size(t)); randn(rs, size(t))];
[P, f, n] = psd_welch(x, fs, [1 numel(t)], 4, 0.5, [0 fs / 2]);
verifyEqual(tc, n, 59);                                   % 120 s in 4 s segments, 50% overlap
verifyEqual(tc, f(2) - f(1), 0.25, 'AbsTol', 1e-12);
verifyEqual(tc, trapz(f, P(2, :)), 1, 'RelTol', 0.02);    % white noise: its variance
verifyEqual(tc, trapz(f, P(1, :)), 3, 'RelTol', 0.02);    % + a sine of amplitude 2: variance 2
in = f >= 8 & f <= 12;
verifyEqual(tc, trapz(f(in), P(1, in)) - trapz(f(in), P(2, in)), 2, 'RelTol', 0.03);
[~, ~, n2] = psd_welch(x, fs, [1 2000; 2001 numel(t)], 4, 0.5, [1 45]);   % no segment across a cut
verifyEqual(tc, n2, 3 + 55);
Pr = pwelch(x' - mean(x', 1), hann(1000), 500, 1000, fs);   % one piece: pwelch itself
verifyEqual(tc, P, Pr', 'RelTol', 1e-12);
verifyError(tc, @() psd_welch(x, fs, [1 500], 4, 0.5, [1 45]), 'psd_welch:pieces');
end

function test_specparam_and_measures(tc)
S = tc.TestData.S;
f = 1:0.25:45;
ap = 2 - 1.5 * log10(f);
P = [10 .^ (ap + 0.6 * exp(-(f - 10) .^ 2 / 2)); 10 .^ ap; 10 .^ (ap + 0.1 * (f >= 8 & f <= 11))];
fit = psd_specparam(f, P, [2 40]);
verifyEqual(tc, fit.exponent(1:2), [1.5; 1.5], 'AbsTol', 0.01);
verifyEqual(tc, fit.offset(1:2), [2; 2], 'AbsTol', 0.01);
verifyEqual(tc, fit.peaks{1}(1, 1:2), [10, 0.6], 'AbsTol', 0.02);
verifyEmpty(tc, fit.peaks{2});
verifyTrue(tc, all(isnan(fit.ap(:, f < 2 | f > 40)), 'all'));
M = psd_measures(f, P, fit, S.psd, NaN, S.iaf);
b = @(name) find(strcmp({S.psd.bands.name}, name));
verifyEqual(tc, M.alpha_frequency(1), 10, 'AbsTol', 0.02);    % the peak
verifyTrue(tc, M.alpha_peak(1));
verifyTrue(tc, isnan(M.alpha_frequency(2)));                  % pure 1/f: no alpha frequency
verifyTrue(tc, isnan(M.alpha_cog(2)));                        % nor a centre of gravity
verifyEqual(tc, M.alpha_cog(1), 10, 'AbsTol', 0.25);         % below the peak: power above the 1/f fit,
                                                              % larger below 10 Hz (9.89 here)
verifyEqual(tc, M.alpha_cog(3), 9.5, 'AbsTol', 0.5);         % a plateau: its centre of gravity
if M.alpha_peak(3)                                            % the frequency only when a peak is fitted
    verifyEqual(tc, M.alpha_frequency(3), 9.5, 'AbsTol', 0.5);
else
    verifyTrue(tc, isnan(M.alpha_frequency(3)));
end
verifyGreaterThan(tc, M.periodic(1, b('alpha')), 0.2);
verifyEqual(tc, M.periodic(2, :), zeros(1, numel(S.psd.bands)), 'AbsTol', 0.01);
% adjacent bands add up; relative power is a share of 1-30 Hz
pw = 10 .^ M.absolute;
verifyEqual(tc, pw(:, b('alpha1')) + pw(:, b('alpha2')), pw(:, b('alpha')), 'RelTol', 1e-10);
verifyEqual(tc, sum(M.relative(:, [b('delta'), b('theta'), b('alpha'), b('beta')]), 2), [100; 100; 100], 'AbsTol', 1e-8);
in = f >= 7 & f <= 13;
verifyEqual(tc, pw(2, b('alpha')), trapz(f(in), P(2, in)), 'RelTol', 1e-12);
% a band anchored on the IAF: this recording's edges, or NaN without an IAF
o = S.psd;
o.bands(end + 1) = struct('name', 'alpha_iaf', 'fmin', 'iaf-bins', 'fmax', 'iaf+bins');
M0 = psd_measures(f, P, fit, o, NaN, S.iaf);
verifyTrue(tc, all(isnan(M0.absolute(:, end))));
M1 = psd_measures(f, P, fit, o, 9.5, S.iaf);
verifyEqual(tc, M1.edges(end, :), [7.5 11.5]);
in = f >= 7.5 & f <= 11.5;
verifyEqual(tc, 10 .^ M1.absolute(2, end), trapz(f(in), P(2, in)), 'RelTol', 1e-12);
end

function test_periodic_lines_and_rising_spectra(tc)
% the periodic power is the mean from edge to edge, so a band is the
% width-weighted mean of its parts, also when anchored edges fall between
% grid points (the plain mean of the grid points inside is not); a narrow
% line is found and flagged, an alpha peak is not; a spectrum rising over
% the fit range gets a negative exponent and its own aperiodic curve
% (process_fooof returns |exponent|)
S = tc.TestData.S;
f = 1:0.25:45;
ap = 2 - 1.5 * log10(f);
P = [10 .^ (ap + 0.6 * exp(-(f - 9.8) .^ 2 / 2)); 10 .^ (ap + 0.4 * exp(-(f - 6) .^ 2 / 4))];
fit = psd_specparam(f, P, [2 40]);
o = S.psd;
o.bands = struct('name', {'a', 'a1', 'a2', 'w', 'w1', 'w2'}, 'fmin', {'iaf-2', 'iaf-2', 'iaf', 7, 7, 8.6}, ...
    'fmax', {'iaf+2', 'iaf', 'iaf+2', 13, 8.6, 13});
M = psd_measures(f, P, fit, o, 9.63, S.iaf);           % edges 7.63, 9.63, 11.63: between grid points
verifyEqual(tc, M.periodic(:, 1), (M.periodic(:, 2) + M.periodic(:, 3)) / 2, 'AbsTol', 1e-12);
verifyEqual(tc, M.periodic(:, 4), (1.6 * M.periodic(:, 5) + 4.4 * M.periodic(:, 6)) / 6, 'AbsTol', 1e-12);
% lines: a 25 Hz line 10 times its bin, an alpha peak and a narrow one (SD 0.35 Hz)
Q = [10 .^ ap; 10 .^ ap; 10 .^ (ap + 1.2 * exp(-(f - 10) .^ 2 / (2 * 0.35 ^ 2)))];
Q(1, f == 25) = 10 * Q(1, f == 25);
[hz, ratio, flagged, nu] = psd_line_check(f, Q(1, :));
verifyEqual(tc, hz, 25);
verifyEqual(tc, ratio, 10, 'RelTol', 0.05);
verifyTrue(tc, flagged);
verifyEqual(tc, nu, 1);
for k = 2:3
    [hz, ratio, flagged, nu] = psd_line_check(f, Q(k, :));
    verifyTrue(tc, isnan(hz) && ratio == 1 && ~flagged && nu == 0, sprintf('spectrum %d', k));
end
[~, ~, flagged, nu] = psd_line_check(f, Q);            % on the mean: diluted to about 4 times
verifyTrue(tc, flagged);
verifyEqual(tc, nu, 1);
% a spectrum rising over 2-40 Hz
up = 10 .^ (-1 + 0.8 * log10(f));
fr = psd_specparam(f, [up; 10 .^ ap], [2 40]);
verifyEqual(tc, fr.exponent, [-0.8; 1.5], 'AbsTol', 0.01);
in = f >= 2 & f <= 40;
verifyEqual(tc, fr.ap(1, in), log10(up(in)), 'AbsTol', 0.01);
Mr = psd_measures(f, up, fit_rows(fr, 1), S.psd, NaN, S.iaf);
verifyEqual(tc, Mr.periodic, zeros(1, numel(S.psd.bands)), 'AbsTol', 0.01);
% the mains: 'none' and '60' fit as '50' when no harmonic lies in 2-40 Hz
f2 = psd_specparam(f, P, [2 40], 'none');
f3 = psd_specparam(f, P, [2 40], '60');
verifyEqual(tc, [f2.exponent, f3.exponent], [fit.exponent, fit.exponent]);
end

function fit = fit_rows(fit, row)
for f = {'offset', 'exponent', 'r2', 'error'}, fit.(f{1}) = fit.(f{1})(row); end
fit.peaks = fit.peaks(row);  fit.ap = fit.ap(row, :);
end

function test_cluster_test(tc)
K = 20;  n = 30;
A = diag(true(K - 1, 1), 1);  A = A | A';                 % a line of 20 neighbours
g = [false(n, 1); true(n, 1)];
rs = RandStream('mrg32k3a', 'Seed', 2);
Y = randn(rs, 2 * n, K);
Y(g, 5:9) = Y(g, 5:9) + 1.2;                              % planted: 5 neighbours higher in B
Y(g, 14:15) = Y(g, 14:15) - 1.2;                          % and 2 lower
o = struct('n_perm', 2000, 'seed', 1, 'n_boot', 0);
[R, D] = perm_ancova_maxT(Y, g, [], o);
C = cluster_perm_test(D.t_obs, D.T, D.df, A, 0.05);
verifyEqual(tc, C.cluster(5:9), ones(5, 1));               % the largest cluster
verifyEqual(tc, C.cluster_sign(5:9), ones(5, 1));
verifyEqual(tc, C.cluster_size(5), 5);
verifyLessThan(tc, C.cluster_p(5), 0.01);
verifyEqual(tc, C.cluster_mass(5), sum(D.t_obs(5:9)), 'RelTol', 1e-12);   % a sum of z (perm_ancova_maxT gives z, as PALM)
verifyEqual(tc, D.t_obs, t_to_z(R.t', R.df'), 'RelTol', 1e-12);
verifyEqual(tc, unique(C.cluster_sign(14:15)), -1);       % the negative one, apart
verifyEqual(tc, D.df, inf(1, K));                         % z: the normal threshold
Y0 = randn(rs, 2 * n, K);                                 % the null: no small cluster p
[~, D0] = perm_ancova_maxT(Y0, g, [], o);
C0 = cluster_perm_test(D0.t_obs, D0.T, D0.df, A, 0.05);
verifyTrue(tc, all(C0.cluster_p(C0.cluster > 0) > 0.05));
Y(:, 7) = NaN;                                            % a column not tested splits the cluster
[~, D1] = perm_ancova_maxT(Y, g, [], o);
C1 = cluster_perm_test(D1.t_obs, D1.T, D1.df, A, 0.05);
verifyEqual(tc, C1.cluster(7), 0);
verifyNotEqual(tc, C1.cluster(6), C1.cluster(8));
% correlations: the z of t_hc3 on N - 2 - covariates degrees of freedom
x = randn(rs, 2 * n, 1);
Yc = randn(rs, 2 * n, K) + 0.8 * x .* ((1:K) >= 3 & (1:K) <= 6);
[Rc, Dc] = perm_partial_corr_maxT(Yc, x, [], struct('n_perm', 2000, 'seed', 1, 'type', 'Pearson'));
verifyEqual(tc, Dc.t_obs, t_to_z(Rc.t_hc3', 2 * n - 2), 'RelTol', 1e-10);
verifyEqual(tc, Dc.df, inf(1, K));
Cc = cluster_perm_test(Dc.t_obs, Dc.T, Dc.df, A, 0.05);
verifyEqual(tc, Cc.cluster(3:6), ones(4, 1));
verifyLessThan(tc, Cc.cluster_p(3), 0.01);
end

function test_cluster_test_equals_fieldtrip(tc)
% cluster_perm_test against FieldTrip's clusterstat (what
% ft_statistics_montecarlo calls with correctm 'cluster'), given the same
% statistics of the data and of every permutation
ft = fileparts(which('ft_defaults'));
assumeNotEmpty(tc, ft, 'FieldTrip not on the path');
tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
for f = {'clusterstat.m', 'findcluster.m', 'combineClusters.m'}
    copyfile(fullfile(ft, 'private', f{1}), tmp);
end
tc.applyFixture(matlab.unittest.fixtures.PathFixture(tmp));
rs = RandStream('mrg32k3a', 'Seed', 3);
K = 40;  n_perm = 3000;
xy = rand(rs, K, 2);
A = squareform(pdist(xy)) < 0.3;  A(1:K + 1:end) = false;
df = 30 + randi(rs, 10, 1, K);                        % a df per column, as with missing values
T = tinv(rand(rs, n_perm, K), 30);  T(:, 3) = NaN;   % a column not tested
for rep = 1:3
    t_obs = tinv(rand(rs, 1, K), 30);
    t_obs(1:10) = t_obs(1:10) + (rep - 1.5) * 3;      % a planted cluster: none, negative, positive
    t_obs(3) = NaN;
    t_obs(20:23) = t_obs(20:23) - 2.5;                % and a negative one
    C = cluster_perm_test(t_obs, T, df, A, 0.05);
    cfg = struct('dim', [K 1], 'numrandomization', n_perm, 'feedback', 'no', 'tail', 0, ...
        'clustertail', 0, 'clusterstatistic', 'maxsum', 'clusterthreshold', 'parametric', ...
        'connectivity', A, 'clustercritval', [tinv(0.025, df(:)), tinv(0.975, df(:))]);
    [~, st] = evalc('clusterstat(cfg, T'', t_obs'')');
    for s = [1 -1]
        if s > 0, f = 'pos'; else, f = 'neg'; end
        mine = C(C.cluster_sign == s, :);
        if ~isfield(st, [f 'clusters']) || isempty(st.([f 'clusters']))
            verifyEmpty(tc, mine, sprintf('rep %d: FieldTrip found no %s cluster', rep, f));
            continue
        end
        lab = st.([f 'clusterslabelmat'])(:)';
        verifyEqual(tc, C.cluster_sign(lab > 0)', repmat(s, 1, nnz(lab > 0)), ...
            sprintf('rep %d: the %s supra-threshold columns', rep, f));
        verifyEqual(tc, nnz(C.cluster_sign == s), nnz(lab > 0));
        cl = st.([f 'clusters']);
        for j = 1:numel(cl)
            in = lab == j;
            verifyEqual(tc, numel(unique(C.cluster(in))), 1, 'one FieldTrip cluster is one cluster here');
            verifyEqual(tc, C.cluster_size(find(in, 1)), nnz(in));
            verifyEqual(tc, C.cluster_mass(find(in, 1)), cl(j).clusterstat, 'RelTol', 1e-10);
            verifyEqual(tc, C.cluster_p(find(in, 1)), min(1, 2 * cl(j).prob), 'AbsTol', 1e-12, ...
                'p = FieldTrip''s, doubled (correcttail ''prob'')');
        end
    end
end
end

function test_neighbours(tc)
S = tc.TestData.S;
ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'C1', 'C2', 'C3', 'C4', 'O1', 'Oz', 'O2', 'Fp1', 'Fp2'};
[A, info] = psd_adjacency('sensor', ch, S);
verifyEqual(tc, A, A');
verifyFalse(tc, any(diag(A)));
c = @(x) strcmp(ch, x);
verifyTrue(tc, all(A(c('Cz'), c('C1') | c('C2') | c('FCz') | c('CPz'))));
verifyFalse(tc, A(c('Fp1'), c('Oz')));
verifyEqual(tc, info.n_neighbours, sum(A, 2));
L = load(fullfile(S.root, '1- EEG source connectivity', 'inputs', 'scout_Desikan-Killiany_68.mat'));
labels = {L.Scouts.Label};
[B, info] = psd_adjacency('source', labels, S);
r = @(x) strcmp(labels, x);
verifyTrue(tc, B(r('precentral L'), r('postcentral L')));
verifyTrue(tc, B(r('lateraloccipital R'), r('inferiorparietal R')));
verifyFalse(tc, B(r('precentral L'), r('precentral R')));         % never across the hemispheres
verifyFalse(tc, B(r('frontalpole L'), r('lateraloccipital L')));
verifyEqual(tc, info.hemi(r('cuneus L')), 1);
verifyEqual(tc, info.hemi(r('cuneus R')), 2);
verifyTrue(tc, all(info.n_neighbours >= 1));
% Brainstorm's vertex adjacency (VertConn) and the triangles' edges give the same region pairs
[Sc, C] = get_atlas(S.atlases{1}, [], {'Faces', 'VertConn'});
verifyNotEmpty(tc, C.VertConn, 'the shipped cortex has its VertConn');
region = zeros(size(C.VertConn, 1), 1);
for k = 1:numel(Sc), region(Sc(k).Vertices) = k; end
[vi, vj] = find(C.VertConn);
R = region(C.Faces);
pf = [R(:, [1 2]); R(:, [2 3]); R(:, [1 3])];
pv = [region(vi), region(vj)];
keep = @(q) unique(sort(q(all(q > 0, 2) & q(:, 1) ~= q(:, 2), :), 2), 'rows');
verifyEqual(tc, keep(pv), keep(pf));
end

function test_source_signals_as_stage1(tc)
% the regions' signals of the PSD step are those stage 1 saves
S = tc.TestData.S;
S.data = static_test_recordings(S);
S.file_pattern = {'sub-*/eeg/*_eeg.set'};
S.results = fullfile(tempdir, 'psd_test_stage1');
S.save_scouts = true;
c = onCleanup(@() rmdir_if(S.results));
subs = get_subs(S.data, S.file_pattern);
sub = subs{1};
out = fullfile(S.results, 'Connectivity');
mkdir(out);
evalc('qc = run_pipeline(sub, S.data, out, S);');
verifyEqual(tc, qc.status, "ok", qc.message);
[~, stem] = fileparts(sub);
w = load(fullfile(out, 'Scouts time-series', stem, 'Scouts_eloreta_ASA_10-05_343_Desikan_Killiany_68.mat'));
eeg = get_eeg_data(fullfile(S.data, sub), S);
ftElec = get_ftElec(S.montages{1});
[X, labels, idx] = get_eeg_montage(eeg.data, eeg.label, ftElec.label, S.channel_subset);
verifyEqual(tc, labels(:), w.channels(:));
X = X - mean(X, 1);
[src, roi] = psd_source_signals(X, idx, eeg.srate, sub, S);
verifyEqual(tc, roi(:), w.roi_labels(:));
verifyEqual(tc, double(src{1}), double(w.src_data), 'AbsTol', 1e-9 * max(abs(double(w.src_data(:)))));
end

function test_settings(tc)
S = tc.TestData.S;
verifyEqual(tc, {S.psd.bands.name}, {'delta', 'theta', 'alpha', 'alpha1', 'alpha2', 'beta', 'beta1', 'beta2'});
verifyEqual(tc, [S.psd.bands.fmin; S.psd.bands.fmax], [1 4 7 7 10 13 13 20; 4 7 13 10 13 30 20 30]);
verifyEqual(tc, S.psd.total, [1 30]);
verifyTrue(tc, S.psd.run && S.psd.cluster.run);
verifyEmpty(tc, S.psd.n_perm);                            % [] = S.n_perm
verifyEqual(tc, S.psd.region_spectra, 'signal');
verifyTrue(tc, S.psd.vertex_maps);
S2 = S;  S2.psd = struct('segment', 2);                   % partial: the rest from the defaults
S2 = check_settings(S2);
verifyEqual(tc, S2.psd.segment, 2);
verifyEqual(tc, S2.psd.fit_range, [2 40]);
S2 = S;  S2.psd = struct('bands', struct('name', 'alpha_iaf', 'fmin', 'iaf-bins', 'fmax', 'iaf+bins'));
S2 = check_settings(S2);
verifyEqual(tc, S2.psd.bands.name, 'alpha_iaf');
verifyEqual(tc, S.psd.power_line, '50');
S2 = S;  S2.psd = struct('power_line', 60);
verifyEqual(tc, getfield(check_settings(S2), 'psd', 'power_line'), '60'); %#ok<GFLD>
S2 = S;  S2.psd = struct('bands', struct('name', 'gamma', 'fmin', 30, 'fmax', 45));   % relative: not a share
verifyWarning(tc, @() check_settings(S2), 'check_settings:psd_total');
S2.psd.measures = {'absolute'};
verifyWarningFree(tc, @() check_settings(S2));
bad = {struct('measure', {{'absolute'}}), struct('measures', {{'abs'}}), struct('levels', {{'scalp'}}), ...
    struct('fit_range', [2 60]), struct('total', [0.5 30]), struct('overlap', 1), ...
    struct('bands', struct('name', 'gamma', 'fmin', 30, 'fmax', 60)), ...
    struct('bands', struct('name', {'a', 'a'}, 'fmin', {1, 2}, 'fmax', {2, 3})), ...
    struct('cluster', struct('run', true, 'alpha', 2)), struct('run', 'yes'), struct('n_perm', 0), struct('n_perm', 2.5), ...
    struct('region_spectra', 'mean'), struct('region_spectra', 1), struct('vertex_maps', 'yes'), ...
    struct('vertex_maps', [true false]), struct('power_line', '55'), struct('power_line', true)};
for i = 1:numel(bad)
    S2 = S;  S2.psd = bad{i};
    verifyError(tc, @() check_settings(S2), 'check_settings:psd', jsonencode(bad{i}));
end
end

function test_planted_alpha_difference(tc)
S = tc.TestData.S;
[folder, truth] = alpha_dataset(S);
S.data = folder;
S.file_pattern = {'sub-*/eeg/*_eeg.set'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.subs_info = fullfile(folder, 'participants.tsv');
S.id_column = 'participant_id';
S.results = fullfile(tempdir, 'psd_test_results');
S.run_name = 'test';
S.n_perm = 1000;
S.figures = true;
S.psd.bands = struct('name', {'alpha', 'beta2'}, 'fmin', {7, 20}, 'fmax', {13, 30});   % planted, and not
S.correlations = struct('name', 'Score', 'columns', '^score$', 'subjects', {{'group', '==', 'patient'}}, ...
    'covariates', {{}}, 'zero_as_missing', false);
out = fullfile(S.results, 'PSD', 'test');
c = onCleanup(@() rmdir_if(out));
log = evalc('run_psd_analysis(S);');
verifyFalse(tc, contains(log, 'failed'), log);
Q = readtable(fullfile(out, 'recordings.csv'), 'TextType', 'string');
verifyEqual(tc, height(Q), 16);
verifyEqual(tc, Q.status, repmat("ok", 16, 1));
verifyGreaterThan(tc, min(Q.fit_r2_sensor), 0.9);
nf = string(Q.near_fit_start);                              % an empty column is read as numbers
verifyTrue(tc, all(ismissing(nf) | nf == ""));               % alpha, beta2: far from 2 Hz
verifyFalse(tc, contains(regexprep(log, '\s+', ' '), 'start of the specparam fit'));
sa = fullfile(out, 'Statistical analysis');
for f = {'summary.csv', 'families.csv', 'clusters.csv', 'neighbours_sensor.csv', 'neighbours_source_eloreta.csv'}
    verifyTrue(tc, isfile(fullfile(sa, f{1})), f{1});
end
for f = {'spectra_sensor_patient_control.png', 'map_sensor_relative_patient_control.png', ...
        'map_source_eloreta_periodic_patient_control.png', 'boxplots_source_eloreta_absolute_patient_control.png', ...
        'map_sensor_aperiodic_alpha_patient_control.png'}
    verifyTrue(tc, isfile(fullfile(out, 'figures', f{1})), f{1});
end
f3 = fullfile(out, 'figures', 'map_source_eloreta_periodic_patient_control.fig');   % the 3-D cortex maps
verifyTrue(tc, isfile(f3));
verifyFalse(tc, isfile(fullfile(out, 'figures', 'map_sensor_periodic_patient_control.fig')));
fig = openfig(f3, 'new', 'invisible');
verifyEqual(tc, numel(findobj(fig, 'Type', 'axes')), 2);                  % a brain per band
verifySubstring(tc, fig.CreateFcn, 'rotate3d');                         % opens visible, to be turned
close(fig);
fams = readtable(fullfile(sa, 'summary.csv'), 'TextType', 'string');
verifyEqual(tc, height(fams), 2 * 9);                     % 2 levels x (3 measures x 2 bands + 3)
% the sensor level: posterior channels differ, more than frontal ones
G = readtable(fullfile(sa, 'sensor', 'psd_sensor_relative_alpha_group_comparison.csv'), 'TextType', 'string');
post = ismember(G.Unit, ["O1", "Oz", "O2", "PO7", "PO8"]);
front = ismember(G.Unit, ["Fp1", "Fpz", "Fp2", "AF7", "AF8"]);
verifyTrue(tc, all(G.g(post) > 0) && all(G.significant(post) | G.cluster_significant(post)));
verifyGreaterThan(tc, mean(G.g(post)) - mean(G.g(front)), 0.5);
% the regions: the occipital ones, more than the frontal ones
G = readtable(fullfile(sa, 'source_eloreta', 'psd_source_eloreta_periodic_alpha_group_comparison.csv'), 'TextType', 'string');
occ = startsWith(G.Unit, ["lateraloccipital", "cuneus", "pericalcarine", "lingual"]);
fr = startsWith(G.Unit, ["frontalpole", "rostralmiddlefrontal", "superiorfrontal", "medialorbitofrontal"]);
verifyTrue(tc, all(G.g(occ) > 0) && mean(G.significant(occ) | G.cluster_significant(occ)) >= 0.75);
verifyGreaterThan(tc, mean(G.g(occ)) - mean(G.g(fr)), 0.5);
% the score follows the patients' occipital alpha
C = readtable(fullfile(sa, 'source_eloreta', 'psd_source_eloreta_periodic_alpha_correlations.csv'), 'TextType', 'string');
verifyGreaterThan(tc, mean(C.r(startsWith(C.Unit, ["lateraloccipital", "cuneus"]))), 0.5);
% not planted: beta2 power
G = readtable(fullfile(sa, 'sensor', 'psd_sensor_absolute_beta2_group_comparison.csv'), 'TextType', 'string');
verifyFalse(tc, any(G.significant) || any(G.cluster_significant));
Cl = readtable(fullfile(sa, 'clusters.csv'), 'TextType', 'string');
verifyTrue(tc, any(Cl.cluster_significant & Cl.family == "psd_sensor_relative_alpha" & Cl.sign == 1));
verifyFalse(tc, any(Cl.cluster_significant & Cl.band == "beta2" & Cl.measure ~= "relative"));   % relative beta2
                                        % is a share of 1-30 Hz: more alpha lowers it (a true effect)
% the regions' power from their signal against the mean of their sources'
% power: never more (Cauchy-Schwarz), and the planted difference in both
R = readtable(fullfile(out, 'region_power_check.csv'), 'TextType', 'string');
verifyEqual(tc, height(R), 68 * 2);                   % regions x bands, one comparison
verifyLessThanOrEqual(tc, max(R.ratio_q3), 1 + 1e-9);
verifyGreaterThan(tc, min(R.ratio_q1), 0);
occ_r = R.band == "alpha" & startsWith(R.region, ["lateraloccipital", "cuneus", "pericalcarine", "lingual"]);
verifyTrue(tc, all(R.g_absolute_signal(occ_r) > 0) && all(R.g_absolute_vertex_mean(occ_r) > 0));
verifyTrue(tc, isfile(fullfile(out, 'figures', 'region_power_check_eloreta.png')));
% the sources' maps: values of every source, the largest g in the occipital
% sources, g as hedges_g, the figures
w = load(fullfile(out, 'vertex_values_eloreta.mat'));
Vk = w.Vk;
[Sc, Cx] = get_atlas(S.atlases{1}, [], {'Vertices'});
nv = size(Cx.Vertices, 1);
verifySize(tc, Vk.values.absolute, [16 nv 2]);
verifySize(tc, Vk.values.relative, [16 nv 2]);
verifyEqual(tc, numel(Vk.groups), 1);
g = double(Vk.groups.absolute.g(:, 1));                % alpha, patient - control
isocc = startsWith({Sc.Label}, {'lateraloccipital', 'cuneus', 'pericalcarine', 'lingual'});
isfr = startsWith({Sc.Label}, {'frontalpole', 'rostralmiddlefrontal', 'superiorfrontal', 'medialorbitofrontal'});
vo = [Sc(isocc).Vertices];  vf = [Sc(isfr).Vertices];
verifyGreaterThan(tc, median(g(vo)) - median(g(vf)), 0.5);
verifyGreaterThan(tc, median(g(vo)), 1);
[~, v] = max(g);
[~, row] = ismember(string(Vk.ids), truth.ids);
pat = truth.group(row) == "patient";
x = double(Vk.values.absolute(:, v, 1));
verifyEqual(tc, g(v), hedges_g(x(pat), x(~pat)), 'RelTol', 1e-5);
verifyEqual(tc, numel(Vk.correlations), 1);
for f = {'vertex_absolute_eloreta_patient_control.png', 'vertex_relative_eloreta_patient_control.png', ...
        'vertex_absolute_eloreta_patient_control.fig', 'vertex_absolute_eloreta_r_Score_score.png'}
    verifyTrue(tc, isfile(fullfile(out, 'figures', f{1})), f{1});
end
% S.psd.region_spectra = 'vertex_mean': the regions' statistics on the mean
% of their sources' power, from the same kept spectra
S2 = S;
S2.run_name = 'test_vm';
S2.figures = false;
S2.correlations = S2.correlations([]);
S2.psd.region_spectra = 'vertex_mean';
out2 = fullfile(S.results, 'PSD', 'test_vm');
c3 = onCleanup(@() rmdir_if(out2));
log3 = evalc('run_psd_analysis(S2);');
verifyFalse(tc, contains(log3, 'failed'), log3);
verifySubstring(tc, log3, 'to compute: 0 (sensor), 0 (source)');
fam2 = readtable(fullfile(out2, 'Statistical analysis', 'summary.csv'), 'TextType', 'string');
verifyTrue(tc, any(startsWith(fam2.family, "psd_source_eloreta_vertexmean_")));
verifyFalse(tc, any(startsWith(fam2.family, "psd_source_eloreta_") & ~contains(fam2.family, "vertexmean")));
G2 = readtable(fullfile(out2, 'Statistical analysis', 'source_eloreta_vertexmean', ...
    'psd_source_eloreta_vertexmean_absolute_alpha_group_comparison.csv'), 'TextType', 'string');
occ2 = startsWith(G2.Unit, ["lateraloccipital", "cuneus", "pericalcarine", "lingual"]);
verifyTrue(tc, all(G2.g(occ2) > 0));
R2 = readtable(fullfile(out2, 'region_power_check.csv'), 'TextType', 'string');
verifyEqual(tc, R2.ratio_median, R.ratio_median, 'AbsTol', 1e-12);    % the same check either way
% the spectra are kept: a second run computes none, also in a session where
% Brainstorm's specparam is not yet on the path (its file is in the key)
p0 = path;  c2 = onCleanup(@() path(p0));
bst = fileparts(which('process_fooof'));
rmpath(bst);
S.figures = false;  S.comparisons = S.comparisons([]);  S.correlations = S.correlations([]);   % spectra only
S.psd.bands = struct('name', {'delta', 'theta'}, 'fmin', {1, 4}, 'fmax', {4, 7});   % delta: 2-4 Hz in the fit
log2 = evalc('run_psd_analysis(S);');
verifySubstring(tc, log2, 'to compute: 0 (sensor), 0 (source)');
w2 = regexprep(log2, '\s+', ' ');                          % warnings may be wrapped
verifySubstring(tc, w2, 'Periodic power of delta (16 of 16 recordings) starts within 2 Hz of the start');
verifyFalse(tc, contains(w2, 'theta ('));                    % theta starts at 4 Hz, 2 Hz into the fit
Q2 = readtable(fullfile(out, 'recordings.csv'), 'TextType', 'string');
verifyEqual(tc, Q2.near_fit_start, repmat("delta", 16, 1));
verifyEqual(tc, truth.n_per_group, [8 8]);
end

function test_spectral_covariate_on_a_fresh_run(tc)
% a covariate from the spectral profile (the statistics compute it, after the PSD): on
% a fresh run it is not yet in the subject table, so the PSD statistics
% compute the profile first. With S.spectral_profile
% false: an error that names it. The log gives the sample sizes and the
% permutation resolution; its lines count q within each family
S = tc.TestData.S;
folder = alpha_dataset(S);
S.data = folder;
S.file_pattern = {'sub-*/eeg/*_eeg.set'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.subs_info = fullfile(folder, 'participants.tsv');
S.id_column = 'participant_id';
S.results = fullfile(tempdir, 'psd_test_results');
S.run_name = 'spectral_covariate';
S.n_perm = 200;  S.figures = false;  S.spectral_profile = true;
S.psd.n_perm = 150;                                                % its own
S.psd.levels = {'sensor'};  S.psd.measures = {'relative'};  S.psd.cluster.run = false;
S.psd.bands = struct('name', 'alpha', 'fmin', 7, 'fmax', 13);
S.covariates = {'age', 'exponent_all'};
S.correlations = S.correlations([]);  S.covariate_sets = S.covariate_sets([]);
csv = fullfile(S.results, 'Connectivity', 'spectral_profile.csv');
if isfile(csv), delete(csv); end                                   % as on a fresh run
out = fullfile(S.results, 'PSD', S.run_name);
c = onCleanup(@() rmdir_if(out));
log = evalc('run_psd_analysis(S);');
verifySubstring(tc, log, 'exponent_all come from the spectral profile: computed first');
verifySubstring(tc, log, 'Sample sizes');
verifySubstring(tc, log, 'Permutation resolution');
verifySubstring(tc, log, 'PSD statistics, 150 permutations');
verifySubstring(tc, log, 'within-family q < 0.05');
verifyTrue(tc, isfile(csv));
G = readtable(fullfile(out, 'Statistical analysis', 'sensor', 'psd_sensor_relative_alpha_group_comparison.csv'), ...
    'TextType', 'string');
verifyTrue(tc, all(G.StatisticalTest == "ANCOVA (age, exponent_all), Welch-type G, Freedman-Lane permutations"));
verifyTrue(tc, all(isfinite(G.p_unc)));
verifyEqual(tc, G.p_unc * 151, round(G.p_unc * 151), 'AbsTol', 1e-9);   % S.psd.n_perm, not S.n_perm
verifyTrue(tc, all(ismember({'p_maxT_family', 'family_p', 'p_bb', 'significant'}, G.Properties.VariableNames)));
S.spectral_profile = false;
verifyError(tc, @() prepare_spectral_covariates(S), 'prepare_spectral_covariates:off');
end

function test_failures_and_fresh_folder(tc)
% a failing step of the PSD analysis is returned, listed in failures.csv
% and stops it with S.stop_on_error; each run starts from an empty
% results/PSD/<run> (no tables of an earlier run are left); a recording that
% fails is saved without a key, so the next run tries it again
S = tc.TestData.S;
folder = alpha_dataset(S);
S.data = folder;
S.file_pattern = {'sub-*/eeg/*_eeg.set'};
S.id_regexp = '^(sub-[A-Za-z0-9]+)_';
S.subs_info = fullfile(folder, 'participants.tsv');
S.id_column = 'participant_id';
S.results = fullfile(tempdir, 'psd_test_failures');
c = onCleanup(@() rmdir_if(S.results));
S.run_name = 'failing';
S.n_perm = 100;  S.figures = false;
S.psd.levels = {'sensor'};  S.psd.measures = {'relative'};  S.psd.cluster.run = false;
S.psd.bands = struct('name', 'alpha', 'fmin', 7, 'fmax', 13);
S.correlations = S.correlations([]);  S.covariate_sets = S.covariate_sets([]);
S.covariates = {'exponent_all'};  S.spectral_profile = false;   % its statistics fail (no spectral profile)
out = fullfile(S.results, 'PSD', S.run_name);
stale = fullfile(out, 'Statistical analysis', 'summary.csv');
mkdir(fileparts(stale));  fclose(fopen(stale, 'w'));
failures = [];
log = evalc('failures = run_psd_analysis(S);');
verifySize(tc, failures, [1 2]);
verifyEqual(tc, failures(1, 1), "PSD statistics");
verifySubstring(tc, log, 'WITH 1 FAILURE(S)');
F = readtable(fullfile(out, 'failures.csv'), 'TextType', 'string', 'Delimiter', ',');
verifyEqual(tc, F.what, "PSD statistics");
verifyFalse(tc, isfile(stale));
S.stop_on_error = true;
verifyError(tc, @() run_psd_quietly(S), 'prepare_spectral_covariates:off');
% a recording that fails: kept without a key, computed again next time
S.stop_on_error = false;  S.covariates = {};
S.min_channels = 500;
evalc('psd_spectra(S);');
d = dir(fullfile(S.results, 'PSD', 'spectra', 'sensor', '*.mat'));
verifyNumElements(tc, d, 16);
w = load(fullfile(d(1).folder, d(1).name));
verifyTrue(tc, isfield(w.r, 'error'));
verifyEqual(tc, w.key, '');
log2 = evalc('psd_spectra(S);');
verifySubstring(tc, log2, 'to compute: 16 (sensor)');
end

function test_region_key_follows_lcmv_kappa(tc)
% LCMV's truncation (S.lcmv_kappa) is in the key of the regions' spectra
% (changing it does not reuse spectra of the other value)
S = tc.TestData.S;
folder = alpha_dataset(S);
S.data = folder;
S.file_pattern = {'sub-control01/eeg/*_eeg.set', 'sub-patient01/eeg/*_eeg.set'};
S.results = fullfile(tempdir, 'psd_test_kappa');
c = onCleanup(@() rmdir_if(S.results));
S.inv = {'lcmv'};  S.psd.vertex_maps = false;
log1 = evalc('psd_spectra(S);');
verifySubstring(tc, log1, 'to compute: 2 (sensor), 2 (source)');
S.lcmv_kappa = 40;
log2 = evalc('psd_spectra(S);');
verifySubstring(tc, log2, 'to compute: 0 (sensor), 2 (source)');
log3 = evalc('psd_spectra(S);');
verifySubstring(tc, log3, 'to compute: 0 (sensor), 0 (source)');
end

function run_psd_quietly(S) %#ok<INUSD> (used by evalc)
evalc('run_psd_analysis(S);');
end

function test_cross_spectra(tc)
% psd_cross_spectra: its diagonal is psd_welch's (pwelch's) spectra, the whole
% matrix MATLAB's cpsd ('mimo') of each piece weighted as psd_welch weights
% them; 0 Hz to Nyquist and an odd segment length too
fs = 256;  rs = RandStream('mrg32k3a', 'Seed', 3);
X = randn(rs, 6, fs * 100);  X(2, :) = X(2, :) + 0.7 * X(1, :);
seg = [1 9000; 9001 15000; 15001 size(X, 2)];
[C, f, n] = psd_cross_spectra(X, fs, seg, 4, 0.5, [1 45]);
[P, f2, n2] = psd_welch(X, fs, seg, 4, 0.5, [1 45]);
verifyEqual(tc, f, f2);
verifyEqual(tc, n, n2);
verifyEqual(tc, diagonals(C), P, 'RelTol', 1e-12);
L = 4 * fs;  acc = 0;  ns = 0;
for p = 1:size(seg, 1)
    x = X(:, seg(p, 1):seg(p, 2))';  x = x - mean(x, 1);
    k = fix((size(x, 1) - L / 2) / (L / 2));
    [Cp, fc] = cpsd(x, x, hann(L), L / 2, L, fs, 'mimo');
    acc = acc + k * real(Cp);  ns = ns + k;
end
keep = fc >= 1 - 1e-9 & fc <= 45 + 1e-9;
verifyEqual(tc, C, permute(acc(keep, :, :) / ns, [2 3 1]), 'AbsTol', 1e-12 * max(abs(C(:))));
verifyEqual(tc, C, permute(C, [2 1 3]));
for fs2 = [250 255]                                       % segments of 750 and 765 samples
    Z = randn(rs, 3, fs2 * 30);
    [C2, f3] = psd_cross_spectra(Z, fs2, [1 size(Z, 2)], 3, 0.5, [0 fs2 / 2]);
    [P2, f4] = psd_welch(Z, fs2, [1 size(Z, 2)], 3, 0.5, [0 fs2 / 2]);
    verifyEqual(tc, f3, f4, 'AbsTol', 1e-9);
    verifyEqual(tc, diagonals(C2), P2, 'RelTol', 1e-12);
end
verifyError(tc, @() psd_cross_spectra(X, fs, [1 500], 4, 0.5, [1 45]), 'psd_cross_spectra:pieces');
end

function test_vertex_spectra(tc)
% psd_vertex_spectra: the PSD of every source's signal w' X (psd_welch of it)
% without computing it, for eLORETA, wMNE and LCMV, normal and free
% orientation (free: the sum of the three); each region's mean of them; at
% least the power of the region's sign-flipped mean signal (no cancellation);
% for eLORETA equal to FieldTrip's source power from the cross-spectra
% (ft_inverse_eloreta with these filters; free orientation: FieldTrip keeps
% the dominant orientation's power, the largest eigenvalue, here the sum)
% and to Brainstorm's PSD of the source map and of a scout (brainstorm_psd)
S = tc.TestData.S;
[~, lf] = get_ftHeadmodel_and_ftLeadfield('EEGLAB_61', 'Desikan_Killiany_68', S);
el = get_ftElec('EEGLAB_61');
Cx = load(fullfile(S.root, '1- EEG source connectivity', 'inputs', 'tess_cortex_pial_low.mat'), 'VertNormals');
A = load(fullfile(S.root, '1- EEG source connectivity', 'inputs', 'scout_Desikan-Killiany_68.mat'));
lf.leadfield = cellfun(@(g) g - mean(g, 1), lf.leadfield, 'UniformOutput', false);
rs = RandStream('mrg32k3a', 'Seed', 5);
fs = 128;
X = randn(rs, 61, 61) * randn(rs, 61, fs * 60);           % correlated channels
X = X - mean(X, 1);
seg = [1 size(X, 2)];
C = psd_cross_spectra(X, fs, seg, 2, 0.5, [1 40]);
vert = arrayfun(@(sc) sc.Vertices(:)', A.Scouts(:)', 'UniformOutput', false);
r = 5;
pick = unique([vert{r}(1:3), 1, 7000, 15002]);
for o = {'normal', 'free'}
    free = strcmp(o{1}, 'free');
    F = get_inverse_solution(X, fs, eye(61), [], lf, el, Cx.VertNormals, {'eloreta', 'wmne', 'lcmv'}, ...
        struct('orientation', o{1}));
    for iv = 1:3
        W = reshape(F(iv, :, :), size(F, 2), size(F, 3));
        [Pm, Pv] = psd_vertex_spectra(W, C, vert, free);
        verifyEqual(tc, size(Pv, 1), numel(Cx.VertNormals) / 3);
        verifyEqual(tc, Pv(pick, :), brute(W, X, pick, free, fs, seg), 'RelTol', 1e-9);
        verifyEqual(tc, Pm(r, :), mean(brute(W, X, vert{r}, free, fs, seg), 1), 'RelTol', 1e-9);
        if ~free                                          % the sign-flipped mean: never more
            op = get_scouts_timeseries(W, 'Desikan_Killiany_68');
            Ps = psd_welch(op * X, fs, seg, 2, 0.5, [1 40]);
            verifyLessThanOrEqual(tc, Ps, Pm * (1 + 1e-9));
            verifyLessThan(tc, median(Ps(:) ./ Pm(:)), 0.9);   % and less in most regions
        end
        if iv == 1                                        % FieldTrip's eLORETA power
            k = 1 + 2 * free;                              % rows per source
            sm = struct('pos', zeros(numel(pick), 3), 'inside', true(numel(pick), 1));
            sm.filter = arrayfun(@(v) W(k * v + (1 - k:0), :), pick(:), 'UniformOutput', false);
            est = ft_inverse_eloreta(sm, [], [], [], C, 'keepmom', 'no');
            if free
                verifyLessThanOrEqual(tc, est.pow, Pv(pick, :) * (1 + 1e-9));
                verifyGreaterThanOrEqual(tc, est.pow, Pv(pick, :) / 3 * (1 - 1e-9));   % the largest of three
            else
                verifyEqual(tc, est.pow, Pv(pick, :), 'RelTol', 1e-10);
            end
            brainstorm_psd(tc, W, X, fs, free, vert{r});
        end
    end
end
end

function brainstorm_psd(tc, W, X, fs, free, v)
% Brainstorm's PSD of a source map (process_psd on an imaging kernel calls
% bst_psd with the kernel; free orientation: bst_source_orient 'sum'; a
% scout with ClusterFuncTime 'after': bst_scout_value 'mean' of its sources'
% PSDs) equal to psd_vertex_spectra on cross-spectra of Brainstorm's own
% segments (Hamming, the mean of each removed, 'physical' units)
bst = fileparts(fileparts(which('bst_bandpass_hfilter')));       % brainstorm3/toolbox
p0 = path;  c = onCleanup(@() path(p0));                         % its folders for this check only
addpath(fullfile(bst, 'timefreq'), fullfile(bst, 'inverse'), fullfile(bst, 'process', 'functions'));
Pb = bst_psd(X, fs, 2, 50, [], W, 'mean', 'physical');            % sources x 1 x frequencies
if free, Pb = bst_source_orient([], 3, [], Pb, 'sum', 'results'); end
Pb = permute(Pb, [1 3 2]);
L = round(2 * fs);  L = L - mod(L, 2);  ov = round(L * 50 / 100);
n = floor((size(X, 2) - ov) / (L - ov));
w = bst_window('hamming', L)';
nf = L / 2 + 1;
C = zeros(size(X, 1), size(X, 1), nf);
for k = 1:n
    x = X(:, (1:L) + (k - 1) * (L - ov));
    a = fft((x - mean(x, 2)) .* w, L, 2) * sqrt(2 / (fs * sum(w .^ 2)));
    a = a(:, 1:nf);
    a(:, [1 end]) = a(:, [1 end]) / sqrt(2);
    for f = 1:nf, C(:, :, f) = C(:, :, f) + real(a(:, f) * a(:, f)'); end
end
C = C / n;
Pc = permute(bst_psd(X, fs, 2, 50, [], [], 'mean', 'physical'), [1 3 2]);
verifyEqual(tc, cell2mat(arrayfun(@(f) diag(C(:, :, f)), 1:nf, 'UniformOutput', false)), Pc, 'RelTol', 1e-10);
[Pm, Pv] = psd_vertex_spectra(W, C, {v}, free);
verifyEqual(tc, Pv, Pb, 'RelTol', 1e-9);
Ps = cell2mat(arrayfun(@(f) bst_scout_value(Pb(v, f), 'mean'), 1:nf, 'UniformOutput', false));
verifyEqual(tc, Pm, Ps, 'RelTol', 1e-9);
end

function P = brute(W, X, v, free, fs, seg)
% the PSD of each source v from its time series (three added when free)
rows = v(:)';
if free, rows = reshape(3 * rows + (-2:0)', 1, []); end
P = psd_welch(W(rows, :) * X, fs, seg, 2, 0.5, [1 40]);
if free, P = reshape(sum(reshape(P, 3, [], size(P, 2)), 1), [], size(P, 2)); end
end

function D = diagonals(C)
D = cell2mat(arrayfun(@(k) diag(C(:, :, k)), 1:size(C, 3), 'UniformOutput', false));
end

function [folder, truth] = alpha_dataset(S)
% 8 controls and 8 patients: every region carries 1/f activity and an
% alpha rhythm (8-12 Hz); in the patients the alpha of the occipital
% regions (lateral occipital, cuneus, pericalcarine, lingual) is 1.5 to
% 2.5 times stronger, and their score follows that gain. The 64 channels
% of sim.static_params, simulated once per computer in a temporary folder.
folder = fullfile(tempdir, 'psd_test_recordings');
version = 2;                                          % of this simulation: older copies are made again
if isfile(fullfile(folder, 'truth.mat'))
    w = load(fullfile(folder, 'truth.mat'));
    truth = w.truth;
    if isfield(truth, 'version') && truth.version == version, return; end
    rmdir(folder, 's');
end
p = sim.static_params(S);
fm = sim.forward_model(S);
[~, row] = ismember(lower(p.channels), lower(fm.labels));
fm.G = fm.G(row, :) - mean(fm.G(row, :), 1);
fm.labels = fm.labels(row);
g = sim.roi_patterns(fm, 0.01);
labels = {fm.scouts.Label};
occ = startsWith(labels, {'lateraloccipital', 'cuneus', 'pericalcarine', 'lingual'})';
rs = RandStream('mrg32k3a', 'Seed', 11);
fs = 256;  n = 60 * fs;
ids = strings(16, 1);  grp = strings(16, 1);  gain = ones(16, 1);
for k = 1:16
    patient = k > 8;
    grp(k) = "control";  if patient, grp(k) = "patient"; end
    ids(k) = sprintf('sub-%s%02d', grp(k), mod(k - 1, 8) + 1);
    if patient, gain(k) = 1.5 + (k - 9) / 7; end         % 1.5 to 2.5
    a = 0.6 * ones(numel(labels), 1);
    a(occ) = 1.2 * gain(k);
    background = g * sim.pink_noise(numel(labels), n, fs, rs);
    X = background + g * (a .* sim.band_noise(numel(labels), n, fs, [8 12], rs));
    X = X + 0.05 * std(background(:)) * randn(rs, size(X));   % sensor noise: not scaled by the alpha
    X = X - mean(X, 1);
    if k == 1, scale = 10 / std(background(:)); end     % microvolts: the same scale for everyone
    X = scale * X;
    sim.write_eeg(fullfile(folder, char(ids(k)), 'eeg', char(ids(k) + "_task-rest_eeg")), X, fs, fm.labels, 'set');
end
score = nan(16, 1);
score(9:16) = round(20 + 8 * (gain(9:16) - 2) + randn(rs, 8, 1), 1);
T = table(ids, grp, round(25 + 3 * randn(rs, 16, 1), 1), score, 'VariableNames', {'participant_id', 'group', 'age', 'score'});
writetable(T, fullfile(folder, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');
truth = struct('ids', ids, 'group', grp, 'gain', gain, 'occipital', {labels(occ)}, 'n_per_group', [8 8], ...
    'version', version);
save(fullfile(folder, 'truth.mat'), 'truth');
end

function rmdir_if(d)
if isfolder(d), rmdir(d, 's'); end
end
