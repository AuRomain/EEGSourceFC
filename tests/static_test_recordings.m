function [folder, truth] = static_test_recordings(S)
%STATIC_TEST_RECORDINGS Simulated recordings for the static-only tests (BIDS layout, three file formats).
%   [folder, truth] = static_test_recordings(S): 4 + 4 recordings of 60 s
%   at 256 Hz on a 64-channel 10-10 cap (sim.simulate_static_dataset), in
%   EEGLAB, EDF and BrainVision formats in turn, with participants.tsv,
%   simulated once per computer in a temporary folder. No private data.
folder = fullfile(tempdir, 'static_test_recordings');
if ~isfile(fullfile(folder, 'truth.mat'))
    addpath(fullfile(S.root, 'simulation'));
    p = sim.static_params(S);
    p.out_dir = folder;
    p.n_per_group = [4, 4];
    p.duration_s = 60;
    p.format = {'set', 'edf', 'vhdr'};
    p.seed = 5;
    evalc('sim.simulate_static_dataset(S, p);');
end
L = load(fullfile(folder, 'truth.mat'));
truth = L.truth;
end
