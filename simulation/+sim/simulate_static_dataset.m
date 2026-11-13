function truth = simulate_static_dataset(S, p)
%SIMULATE_STATIC_DATASET Simulated resting-state recordings with a planted static coupling difference (BIDS layout).
%   truth = sim.simulate_static_dataset(S, p) writes, in p.out_dir:
%     sub-<group><nn>/eeg/sub-<group><nn>_task-rest_eeg.<ext>  one recording
%                       per subject (p.format), in a BIDS-like layout
%     participants.tsv  participant_id, group, age, sex, score (the last
%                       group only: correlated with the planted coupling)
%     truth.mat         what was planted: each recording's coupling, the
%                       network's regions and phase lags
%   S: settings (pipeline_settings: the cap, atlas and network table of
%   stage 1; the recordings have the channels p.channels of the cap);
%   p: sim.static_params (the model is described there). The source model
%   is sim.simulate_recording with one "state" always on: the planted
%   network's regions share a driver with weight equal to the recording's
%   coupling. Nothing here comes from real participants.
rs = RandStream('mrg32k3a', 'Seed', p.seed);
fm = sim.forward_model(S);
if ~isempty(p.channels)                          % the cap's channels, average-referenced
    [ok, row] = ismember(lower(p.channels), lower(fm.labels));
    if ~all(ok)
        error('simulate_static_dataset:channels', 'Not in %s: %s', S.montages{1}, strjoin(p.channels(~ok), ', '));
    end
    fm.G = fm.G(row, :) - mean(fm.G(row, :), 1);
    fm.labels = fm.labels(row);
end
nroi = numel(fm.scouts);
mask = fm.network == p.network;
if ~any(mask), error('simulate_static_dataset:network', 'No region in network %s', p.network); end
lags = (rand(rs, nroi, 1) - 0.5) * pi;          % phase lags, shared by all recordings
g = sim.roi_patterns(fm, p.patch_radius_m);
fs = p.srate;
n = round(p.duration_s * fs);
lab = ones(n, 1);                                % one state, always on
fmts = cellstr(p.format);
rows = {};
coupling = [];
k = 0;
for gi = 1:numel(p.groups)
    for si = 1:p.n_per_group(gi)
        k = k + 1;
        id = sprintf('sub-%s%02d', p.groups{gi}, si);
        c = min(max(p.coupling(gi) + p.coupling_sd * randn(rs), 0), 0.99);
        q = p;
        q.coupling = c;
        X = sim.simulate_recording(lab, mask, lags, g, fm.G, fs, q, rs);
        fmt = fmts{mod(k - 1, numel(fmts)) + 1};
        file = sim.write_eeg(fullfile(p.out_dir, id, 'eeg', [id '_task-rest_eeg']), X, fs, fm.labels, fmt);
        rows(end + 1, :) = {string(id), string(p.groups{gi}), string(file)}; %#ok<AGROW>
        coupling(end + 1, 1) = c; %#ok<AGROW>
        fprintf('  simulated %s (%s, coupling %.2f)\n', id, fmt, c);
    end
end
ids = [rows{:, 1}]';
grp = [rows{:, 2}]';
% the subject table: age and sex unrelated to the EEG; in the last group a
% score correlated with the coupling (planted_r), missing in the others
nsub = numel(ids);
age = round(25 + 4 * randn(rs, nsub, 1), 1);
sex = repmat("F", nsub, 1);  sex(rand(rs, nsub, 1) < 0.5) = "M";
score = nan(nsub, 1);
last = grp == string(p.groups{end});
z = (coupling(last) - mean(coupling(last))) / std(coupling(last));
e = randn(rs, nnz(last), 1);
score(last) = round(20 + 5 * (p.planted_r * z + sqrt(1 - p.planted_r ^ 2) * e), 1);
T = table(ids, grp, age, sex, score, 'VariableNames', {'participant_id', 'group', 'age', 'sex', 'score'});
writetable(T, fullfile(p.out_dir, 'participants.tsv'), 'FileType', 'text', 'Delimiter', '\t');
truth = struct('params', p, 'ids', ids, 'group', grp, 'files', [rows{:, 3}]', 'coupling', coupling, ...
    'network', p.network, 'network_regions', {{fm.scouts(mask).Label}}, 'mask', mask, 'lags', lags, ...
    'roi_labels', {{fm.scouts.Label}}, 'score', score, 'created', char(datetime('now')));
save(fullfile(p.out_dir, 'truth.mat'), 'truth');
end
