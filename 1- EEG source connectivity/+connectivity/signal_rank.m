function r = signal_rank(X)
%SIGNAL_RANK Independent signals in the rows of X, the rounding of 32-bit floats not counted (MNE-Python's rule).
%   r = connectivity.signal_rank(X), X nSignals x nSamples (real)
%   The singular values of the centred signals above n * eps('single')
%   times the largest, n the number of signals: the tolerance of
%   MNE-Python's estimate_rank with tol = 'float32' (mne/rank.py,
%   _estimate_rank_from_s: tol = s.shape[-1] * max_s * eps with eps =
%   np.finfo(np.float32).eps; rank = sum(s > tol)), the one MNE uses for
%   data saved as 32-bit floats, as EEGLAB saves them. Computed from the
%   eigenvalues e of the covariance (e proportional to s^2):
%   e > (n * eps('single'))^2 * max(e), 5.1e-11 of the largest for 60
%   signals, 2.4e-10 for 129. Signals of real recordings, at least the
%   amplifier's noise, lie far above it; the rounding of 32-bit floats far
%   below.
%   Used by run_pipeline (data_rank, the channels) and
%   connectivity.symmetric_orthogonalisation (the regions). A fixed limit
%   of 1e-12 of the largest eigenvalue (singular values 1e-6 of the
%   largest) could count the rounding of recordings saved as 32-bit floats
%   after cleaning as one more signal (tests/test_static_estimators.m,
%   test_signal_rank_rule).
X = double(X);
n = size(X, 1);
Xc = X - mean(X, 2);
e = eig(Xc * Xc' / max(1, size(Xc, 2)), 'vector');
r = nnz(e > (n * double(eps('single'))) ^ 2 * max(e));
end
