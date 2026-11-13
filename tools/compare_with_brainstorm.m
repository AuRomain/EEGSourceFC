function T = compare_with_brainstorm(X, srate, labels, opts)
%COMPARE_WITH_BRAINSTORM The phase measures of stage 1 against Brainstorm's bst_connectivity.
%   T = compare_with_brainstorm(X, srate, labels) compares PLV, ciPLV,
%   wPLI, PLI and dwPLI of this pipeline (+connectivity) with Brainstorm's
%   on broadband region time series: X nROI x nSamples (src_data of a
%   Scouts_*.mat file, which stage 1 saves with S.save_scouts = true), srate
%   in Hz, labels the region names.
%   - Both start from one analytic signal: the band-pass of Brainstorm's
%     connectivity processes (process_bandpass, 'bst-hfilter-2019', no
%     mirror) and MATLAB's hilbert over the whole stretch.
%   - Brainstorm: its NxN phase-synchronisation process (process_plv1n ->
%     bst_connectivity; Hilbert transform; magnitude) over the whole stretch
%     (time resolution 'none') and in windows of opts.window s with 50%
%     overlap ('windowed'). This pipeline: static_full_cmats
%     (static_connectivity) and the cmats of the same windows (+connectivity).
%     Expected: equal (cmats are single precision: 3e-8).
%   - dwPLI: Brainstorm's Hilbert dwPLI (not offered in its interface) is
%     above 1 on continuous data: its sums are divided by the number of
%     samples before the debiasing, which then subtracts N times too much.
%     The reference is Brainstorm's Fourier dwPLI over epochs of opts.epoch s
%     taken as trials (one window each; the frequencies of the band
%     averaged), the estimator of Vinck et al. (2011). This pipeline's dwPLI
%     is computed here on the same epochs (S.dwpli_epoch_length =
%     opts.epoch; FieldTrip's ft_freqanalysis and ft_connectivity_wpli on
%     the band-passed signals, connectivity.dwpli_epochs), Brainstorm's on
%     the broadband epochs with its own Fourier transform: close, not
%     equal (tests/test_reference_connectivity compares it with FieldTrip
%     itself: equal).
%   opts (struct, every field optional): band [8 13] Hz; seconds 120 (the
%   first seconds of X); regions 24 (evenly spaced rows of X: Brainstorm's
%   Hilbert path holds an nROI x nROI x nSamples array); window 0.5 s;
%   epoch 1 s.
%   Brainstorm (the copy on the path, added by setup_toolboxes) runs headless
%   with a temporary user folder and database, deleted at the end: the
%   user's Brainstorm folder (~/.brainstorm) and databases are not used.
%   Brainstorm must not be running in this MATLAB session; its folders stay
%   on the path afterwards.
%   T: measure, pipeline, brainstorm (what is compared), max_abs_diff, r
%   (over the edges, and the windows), mean_pipeline, mean_brainstorm,
%   min_brainstorm, max_brainstorm.
%   Example, from the pipeline folder with the toolboxes installed:
%       setup_toolboxes(pipeline_settings());  addpath('tools');
%       L = load(['results/Connectivity/Scouts time-series/<recording>/' ...
%           'Scouts_eloreta_<cap>_Desikan_Killiany_68.mat']);
%       T = compare_with_brainstorm(L.src_data, L.srate, L.roi_labels)
%   On three recordings of the OpenNeuro example: PLV, ciPLV, wPLI and PLI
%   equal to 1e-14 (whole stretch) and 3e-8 (479 windows); Brainstorm's
%   Hilbert dwPLI 1.2 to 2.0 (windows: 1 to 14). This pipeline's dwPLI
%   (FieldTrip's estimator) is compared with FieldTrip itself in
%   tests/test_reference_connectivity.
if nargin < 4, opts = struct(); end
def = struct('band', [8 13], 'seconds', 120, 'regions', 24, 'window', 0.5, 'epoch', 1);
for f = fieldnames(def)', if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end, end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, '1- EEG source connectivity'));             % static_connectivity, +connectivity
fmin = opts.band(1);  fmax = opts.band(2);
rows = round(linspace(1, size(X, 1), min(opts.regions, size(X, 1))));
X = double(X(rows, 1:min(end, round(opts.seconds * srate))));
labels = cellstr(labels(rows));
n = numel(rows);
nt = size(X, 2);
low = find(tril(true(n), -1));
c = start_brainstorm(); %#ok<NASGU>                                % stopped and deleted with c

%% one analytic signal, this pipeline's measures
Z = hilbert(process_bandpass('Compute', X, srate, fmin, fmax, 'bst-hfilter-2019', 0).').';
seg = [1, nt];
measures = {'plv', 'ciplv', 'wpli', 'pli', 'dwpli'};
M = static_connectivity(Z, seg, measures, srate, fmin, fmax, ...  % static_full_cmats; dwPLI on
    struct('dwpli_epoch_length', opts.epoch));                     % Brainstorm's epochs
[p, ~, w] = connectivity.plv_wpli_sliding_window(Z, srate, opts.window, opts.window / 2, seg);
ci = connectivity.ciplv_sliding_window(Z, srate, opts.window, opts.window / 2, seg);
q = connectivity.pli_sliding_window(Z, srate, opts.window, opts.window / 2, seg);
windows = struct('plv', p, 'ciplv', ci, 'wpli', w, 'pli', q);

%% Brainstorm
iStudy = db_add_condition('S1', 'stretch');
file = add_matrix(iStudy, X, srate, labels, 'stretch');
iEp = db_add_condition('S1', 'epochs');
ep = round(opts.epoch * srate);
ep_files = arrayfun(@(e) add_matrix(iEp, X(:, (e - 1) * ep + (1:ep)), srate, labels, sprintf('epoch %d', e)), ...
    1:floor(nt / ep), 'UniformOutput', false);
out = {};
for m = {'plv', 'ciplv', 'wpli', 'pli'}
    A = squeeze(M(strcmp(measures, m{1}), :, :));
    B = bst_measure(file, m{1}, 'hilbert', 'none', 2, n, opts);
    out(end + 1, :) = {m{1}, 'whole stretch (static_full_cmats)', 'Hilbert, whole stretch', ...
        compare(A(low), B(low))}; %#ok<AGROW>
    B = bst_measure(file, m{1}, 'hilbert', 'windowed', 2, n, opts);
    V = permute(double(windows.(m{1})), [2 3 1]);
    k = min(size(V, 3), size(B, 3));
    V = reshape(V(:, :, 1:k), n * n, k);  B = reshape(B(:, :, 1:k), n * n, k);
    out(end + 1, :) = {m{1}, sprintf('%d windows (cmats)', size(V, 2)), ...
        sprintf('Hilbert, %d windows of %g s', size(B, 2), opts.window), compare(V(low, :), B(low, :))}; %#ok<AGROW>
end
A = squeeze(M(5, :, :));
[B, freqs] = bst_measure(ep_files, 'dwpli', 'stft', 'none', 1, n, opts);
B = mean(B(:, :, 1, freqs >= fmin & freqs <= fmax), 4);
out(end + 1, :) = {'dwpli', sprintf('whole stretch, epochs of %g s (static_full_cmats)', opts.epoch), ...
    sprintf('Fourier, %d epochs of %g s as trials', numel(ep_files), opts.epoch), compare(A(low), B(low))};
B = bst_measure(file, 'dwpli', 'hilbert', 'none', 1, n, opts);
out(end + 1, :) = {'dwpli', '(Brainstorm''s values)', 'Hilbert, whole stretch', compare(B(low), B(low))};
T = [cell2table(out(:, 1:3), 'VariableNames', {'measure', 'pipeline', 'brainstorm'}), struct2table(vertcat(out{:, 4}))];
end

function c = start_brainstorm()
% Brainstorm headless, with its user folder and database in a temporary folder
bst = fileparts(fileparts(fileparts(which('bst_bandpass_hfilter'))));   % brainstorm3/toolbox/math
if isempty(bst), error('compare_with_brainstorm:brainstorm', 'Brainstorm is not on the path (setup_toolboxes)'); end
addpath(bst);
if brainstorm('status'), error('compare_with_brainstorm:running', 'Brainstorm is running: quit it first (brainstorm stop)'); end
tmp = tempname;
mkdir(tmp);
home = char(java.lang.System.getProperty('user.home'));
java.lang.System.setProperty('user.home', tmp);    % Brainstorm's user folder: <tmp>/.brainstorm
c = onCleanup(@() stop_brainstorm(home, tmp));
brainstorm('setpath');
bst_startup(bst, -1, 'local', '');                 % server mode, database in <tmp>/.brainstorm/local_db
gui_brainstorm('CreateProtocol', 'CompareWithBrainstorm', 0, 0);
db_add_subject('S1');
end

function stop_brainstorm(home, tmp)
try, brainstorm('stop'); catch, end
java.lang.System.setProperty('user.home', home);
[~, name] = fileparts(tmp);
p = strsplit(path, pathsep);
p = p(contains(p, name));                          % Brainstorm's folders in <tmp> (mex, process)
if ~isempty(p), rmpath(p{:}); end
if isfolder(tmp), rmdir(tmp, 's'); end
end

function file = add_matrix(iStudy, X, srate, labels, comment)
sMat = db_template('matrixmat');
sMat.Value = X;
sMat.Time = (0:size(X, 2) - 1) / srate;
sMat.Description = labels(:);
sMat.Comment = comment;
file = db_add(iStudy, sMat);
end

function [R, freqs] = bst_measure(files, method, tf, timeres, measure, n, opts)
% process_plv1n (-> bst_connectivity); R: n x n x time x frequency; measure 1 = signed, 2 = magnitude
if ~iscell(files), files = {files}; end
args = {'plvmethod', method, 'plvmeasure', measure, 'tfmeasure', tf, 'timeres', timeres};
if strcmp(tf, 'hilbert')
    args = [args, {'tfedit', struct('Freqs', {{'band', sprintf('%g, %g', opts.band), 'mean'}})}];
else
    args = [args, {'tfedit', struct('StftWinLen', opts.epoch, 'StftWinOvr', 0, 'StftFrqMax', 2 * opts.band(2))}];
end
if strcmp(timeres, 'windowed'), args = [args, {'avgwinlength', opts.window, 'avgwinoverlap', 50}]; end
mode = 'input';
if numel(files) > 1, mode = 'avg'; end                % the epochs as trials
out = bst_process('CallProcess', 'process_plv1n', files, [], args{:}, 'outputmode', mode);
if isempty(out), error('compare_with_brainstorm:bst', 'Brainstorm returned no result (%s, %s, %s)', method, tf, timeres); end
F = in_bst_timefreq(out(1).FileName, 0);
TF = F.TF;
if size(TF, 1) == n * (n + 1) / 2, TF = process_compress_sym('Expand', TF, n); end
R = reshape(TF, n, n, size(TF, 2), size(TF, 3));
for k = 1:n, R(k, k, :, :) = 0; end
freqs = F.Freqs;
end

function s = compare(a, b)
a = double(a(:));  b = double(b(:));
s = struct('max_abs_diff', max(abs(a - b)), 'r', corr(a, b), 'mean_pipeline', mean(a), ...
    'mean_brainstorm', mean(b), 'min_brainstorm', min(b), 'max_brainstorm', max(b));
end
