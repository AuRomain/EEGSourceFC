# Methods of EEGSourceFC

What each step computes, in more detail than the README.

## Stage 1: from the recordings to connectivity

1. **Reading** (FieldTrip). Channels are matched by label to the cap's
   template positions; other channels, flat channels and missing ones are
   left out (the leadfield follows the channels kept). The data rank, the
   number of independent signals in the channels kept (eigenvalues of
   their covariance, by MNE-Python's rule for data stored as 32-bit
   floats, `connectivity.signal_rank`; each ICA component removed or
   channel interpolated lowers it, the rounding of 32-bit floats does
   not), is reported in `stage1_summary.csv`,
   with a warning below `S.min_channels`. The data are
   re-referenced to the average. Boundary events and stretches of NaN cut
   the recording into continuous pieces. `S.core_start` and
   `S.core_duration` choose the analysed stretch.
2. **Inverse solution** on the ICBM152 template (OpenMEEG boundary-element
   head model, 15,002 cortical sources oriented normal to the surface):
   eLORETA (`S.eloreta_lambda`), LCMV (`S.lcmv_lambda`; `S.lcmv_kappa`:
   `[]` inverts the regularised covariance as FieldTrip does, `'rank'`
   only within the data rank, as MNE-Python's beamformers, for data with
   fewer independent signals than channels) or wMNE (`S.wmne`). The PCA
   weights below come from the broadband covariance: one set per
   recording, the same in every band. Each region's time series comes from its sources by
   Brainstorm's scout function `S.scout_function`: `'mean'` (the default:
   the mean, each source's sign flipped to the region's dominant
   orientation) or `'pca'` (the first principal component of the
   region's sources, from their covariance over the whole analysed
   recording); `'mean'` and `'pca'` are linear and folded into the
   inverse operator. Also `'max'` (the source of largest absolute value
   at each sample, with its sign), `'power'` (the mean of the squared
   sources) and `'rms'` (its square root), computed from the source time
   series. As in Brainstorm's connectivity processes, the function is
   applied to the broadband signals and the result is band-passed (step
   3): for power and RMS, which have no sign, that is the band's part of
   the fluctuations of the region's power, not the region's oscillation,
   so phase measures on them do not measure the band's phase (a
   warning).
3. **Band-pass** of each piece in each band (below), then its analytic
   signal (MATLAB's `hilbert`), then the edges of the piece are dropped.
4. **Connectivity over the whole recording** (`static_connectivity.m`),
   from the analytic signals of every piece kept:
   - PLV, PLI, wPLI, ciPLV from the phase differences;
   - debiased wPLI (`dwpli`, Vinck et al. 2011), computed by FieldTrip
     with epochs as trials: each piece cut into consecutive epochs of
     `S.dwpli_epoch_cycles` (6) cycles of the band centre (alpha 7–13 Hz:
     0.6 s) or `S.dwpli_epoch_length` seconds, `ft_freqanalysis` (mtmfft,
     Hanning taper) on the band-passed signals at the Fourier frequencies
     of the epochs inside the band, `ft_connectivity_wpli` with `debias`,
     then the mean over those frequencies (`+connectivity/dwpli_epochs.m`).
     An estimate needs at least 10 epochs, otherwise it is NaN with a
     note. It estimates the squared wPLI: about 0 without coupling, and it
     can be negative;
   - AEC (correlation of the amplitude envelopes), `aec_orth` (each pair
     orthogonalised in both directions, then averaged), `aec_symorth` and
     `plv_symorth` (symmetric orthogonalisation of all regions first).
   The same is computed on the odd and on the even blocks of
   `S.reliability_block` seconds (`split_halves.m`) for the reliability,
   with band power of each region.
5. **Window mean** (the default; `S.static_matrices` including
   `'window_mean'`): the same measures in windows, then averaged
   (`get_all_connectivity_mats.m`, the windowed functions of
   `+connectivity`): phase measures in windows of
   `S.window_cycles` cycles of the band centre or `S.window_length`
   seconds, overlapping by `S.window_overlap`; AEC measures in windows of
   `S.aec_window` seconds every `S.aec_step` seconds; dwPLI in windows of
   60 cycles of the band centre every 30 cycles (6 s every 3 s in
   alpha 7–13 Hz), each cut into consecutive epochs of 6 cycles, the
   trials of its estimate (10 per window, or 10 epochs when 60 cycles hold
   fewer, with epochs in seconds; counted from their length in seconds, so
   rounding each to whole samples does not lose one). An estimate needs at least 10 epochs: `S.measure_windows.dwpli`
   sets other windows, and a window of fewer is an error before stage 1). Windows never cross
   a cut; the analytic signal is that of step 3 (one Hilbert transform per
   piece). The mean of the windows inside each half gives the reliability
   of the window mean. A measure can have windows of its own
   (`S.measure_windows`, e.g. `struct('wpli', struct('window_length', 2))`);
   the others keep the common ones.

## Whole recording or windows

`S.static_matrices` chooses how each measure is estimated:
`'window_mean'` (the default, step 5), `'full'`
(step 4) or both, each analysed and reported on its own
(`features_window_mean.csv`, `features_full.csv`,
`static_reliability_window_mean.csv`, `static_reliability.csv`). Several
window settings for the same measure are separate runs, each with its own
`S.results`.

- `'full'` is how Brainstorm (time resolution "none"), FieldTrip and MNE
  estimate the static connectivity of continuous data, and the best
  estimate of a stable coupling: every sample adds to one estimate.
- `'window_mean'` averages the measure over short windows. The
  estimate of a phase measure from a short window is biased upwards: in
  0.6 s windows of independent alpha-band noise wPLI is about 0.5 (over
  minutes, near 0), and the bias grows with the strength of the rhythm,
  so it follows band power: in real recordings the window mean of wPLI is
  several times its whole-recording value and follows band power more
  closely than that value. Longer windows (`S.window_length = 2` or more) reduce the
  bias; `dwpli` has none (about 0 without coupling in every window:
  simulated alpha noise, 20 s windows, mean over pairs and windows below
  0.02). For AEC the two agree closely in real recordings, and windows of 6 s every 0.5 s
  are those of O'Neill et al. (2017).
- The window mean captures a coupling whose phase lead changes over time
  (the sign of the phase difference flips, or it drifts round the
  circle), which the whole recording averages out.
- A short stretch of much larger amplitude weighs more in the whole
  recording than in the window mean (in simulated noise, a 5 s burst ten
  times larger and coupled in two regions raised their whole-recording
  wPLI from 0.04 to 0.77 and AEC from 0 to 0.86, their window means
  hardly): mark such stretches as bad (`S.bad_events`) before either.
- On the demo's simulated recordings (`run_static_demo`: 16 + 16, 120 s,
  stronger lagged alpha coupling within the default mode network in the
  patients), wPLI over the whole recording found the planted difference
  (within-DMN g 2.17, family-wise p 0.0005). Its window mean did not with
  0.6 s windows (g 0.42) or 2 s windows (g 0.93), and did with 6 s windows
  (g 1.58, p 0.004). In a short window a weak coupling is small against
  the noise of the window's sum, whose absolute value wPLI takes, so it
  raises the window's expected value only to second order. The results of
  one matrix also depend on the other analysed with it, through the
  adjustment over families: 12 of the 28 network blocks of the whole
  recording were significant next to the 0.6 s window mean, 15 next to the
  2 s or 6 s one.
- Windows longer than `S.min_segment` add nothing from the pieces between
  the two lengths, and a recording made only of such pieces has no window
  mean (NaN, with a warning naming it); `check_settings` warns
  (`check_settings:window_vs_segment`). The template's `S.min_segment` is
  6 s so that the 6 s AEC windows fit in every piece kept (4 s is enough
  for `'full'`); for recordings that are each a single piece of minutes,
  the two give the same data.

**Which matrix for phase measures.** For wPLI, ciPLV, PLV and PLI the
whole recording (`'full'`) is recommended: `S.static_matrices =
{'full'}`, or both. Their window mean follows band power more than
coupling: in real recordings the window mean of wPLI can correlate only
weakly with its whole-recording value across recordings, be several
times larger, and be less reliable between split halves. The window mean
stays the default; `run_statistics` warns once per
run when it analyses the window mean of a phase measure
(`run_statistics:window_mean`, kept in `analysis_log.txt`). With it,
give band power (stage 2, `docs/psd.md`) as a covariate. AEC (6 s
windows) and dwPLI (debiased, windows of 60 cycles) are not part of the warning.

**Size of the bias.** A window of T seconds of a band B Hz wide holds
about T·B independent phase differences, and without coupling the phase
measures of N independent phase differences are of the order of 1/√N
(PLV: E[PLV²] = 1/N; the sample-size bias of Vinck et al. 2011), so the
bias of the window mean falls about as 1/√(T·B): 6 cycles of alpha 7–13
Hz (0.6 s, T·B = 3.6) gave about 0.5 without coupling, a band of 4 Hz in
the same window (T·B = 2.4) more, 2 s windows (T·B = 12) or 6 s (36)
much less, the whole recording (minutes) almost none. A narrow band in a
short window is the worst case. With bands anchored on the IAF
(`S.individual_bands`) the windows, dwPLI's epochs and the filter of
every recording are those of the band at IAF 10 Hz, not of each
recording's band centre (6 cycles of a higher frequency are shorter), so
the bias does not change with the IAF.

**Reliability is not validity.** A stable bias is reliable too: the
window mean of a phase measure follows band power, which is stable within
a recording, so its split-half reliability can be high without coupling,
and a high reliability does not show that a matrix measures coupling or
choose between matrices (`r_within_groups`, "Reliability").

## Filters

Each continuous piece is band-passed with Brainstorm's filter (even-order
linear-phase FIR, 60 dB stop band, mirrored edges, delay compensated),
turned into its analytic signal, and then its edges are discarded: at
least `S.edge_trim` s and at least the filter's transient (the
half-length holding 99.9% of the impulse response's energy). With
`S.filter_transition = 'brainstorm'` (Brainstorm's default: 0.5 Hz
transition bands when a band's lower edge is at or below 5 Hz, 1 Hz
above) the bands are sharp and the transients short enough:

| Band | Transition | Filter length | Transient | Attenuation 2 Hz below the band |
|---|---|---|---|---|
| delta 1–4 Hz | 0.5 Hz | 7.3 s | 1.8 s | −64 dB at 0 Hz |
| theta 4–7 Hz | 0.5 Hz | 7.3 s | 1.8 s | −90 dB |
| alpha 7–13 Hz | 1 Hz | 3.6 s | 0.9 s | −72 dB |
| beta1 13–20 Hz | 1 Hz | 3.6 s | 0.8 s | −68 dB |
| beta2 20–30 Hz | 1 Hz | 3.6 s | 0.8 s | −70 dB |
| gamma 30–45 Hz | 1 Hz | 3.6 s | 0.7 s | −69 dB |

(the same at 256, 512 and 1000 Hz). So at most 1.8 s are discarded at each
end of a piece, 1 s from alpha up, and a piece must be at least 8 s long
(9.6 s in delta and theta) to keep `S.min_segment` (6) s. `'eeglab'` (25%
of the lower edge, at least 2 Hz) gives shorter filters with wider
transition bands; `'window'` the narrowest transition whose
filter stays within two windows of 6 cycles of the band; a number sets
the width in Hz. A band anchored on the IAF gets the transition of the
band at IAF 10 Hz in every recording (at most its lower edge), so with
`'brainstorm'` the filter does not change in a step where a recording's
lower edge crosses 5 Hz. The edges of a band are the
edges of its pass band (full gain): the gain is half (−6 dB) half a
transition outside each edge and reaches the stop band a whole
transition outside (Kaiser window), so neighbouring bands (alpha 7–13,
beta1 13–20 Hz) both pass the frequencies around their shared edge, over
half the sum of their transitions between their −6 dB points.

## Statistics

For every band, inverse method and measure of `S.static_measures`
(`static_analysis.m`, `compute_NBS_analysis.m`):

- **Feature families** (`S.static_families`):
  - `global`: `Strength` (mean node strength), `ClustCoeff` (mean weighted
    clustering coefficient), `CPL` (characteristic path length, lengths
    1/weight); the Brain Connectivity Toolbox's definitions on the
    positive part of a signed measure such as AEC (`gt_window_metrics.m`;
    a test compares them with BCT when it is installed). dwPLI, centred
    on 0 without coupling, has `Strength` only, from its signed weights
    (BCT's `strengths_und_sign`), since its positive part is biased
    upwards by the noise;
  - `networks`: mean connectivity inside each network and between each
    pair of networks of `S.network_table` (Yeo-7 for DK-68);
  - `regions`: node strength of every region (dwPLI: of its signed
    weights);
  - `power`: log10 band power over all regions and each network (tested
    once per band and inverse method, since it does not depend on the
    measure).
  A recording for which a measure could not be computed (NaN matrix,
  `measures_nan` of `stage1_summary.csv`) has missing connectivity
  features and is left out of the NBS, with a warning naming it.
  `S.static_families = {}` tests the edges only (NBS, enrichment, edge
  tables) and their reliability, with a warning.
- **Sample sizes.** Before the families, the recordings of every
  comparison and correlation set are counted (those of stage 1, made or
  reused, with a row in the subject table and every covariate): a comparison needs at least 2
  per group and the covariates + 4 in all (2 residual degrees of freedom),
  a correlation the covariates + 4; a warning names those that cannot be
  tested, whose rows would be NaN (comparisons) or absent (correlations).
  Warnings also name a score with a recording of high leverage (Pearson
  with permutations: leverage above 3 p / n, `run_statistics:leverage`),
  and the recordings left out of a band anchored on the IAF, per group of
  every comparison (`run_statistics:bands_skipped`). The subject table is
  checked against the recordings when the statistics and the PSD
  statistics start, also when run alone; a covariate whose values are
  mostly numbers but not all is an error naming the others, and usual
  missing-value codes (−999, 99, ...) far from their column's values give
  a warning.
- **Tests.** Group comparisons: Welch t, or with covariates ANCOVA with
  the Welch-type G statistic and Freedman–Lane permutations
  (`S.n_perm`), the statistic of the data and of every permutation
  turned into z on its own Welch–Satterthwaite or Aspin–Welch degrees of
  freedom before the p-values are counted, as PALM does (Winkler et al.
  2014). Each comparison gives both groups' means and SDs (`mean_A`,
  `sd_A`, `mean_B`, `sd_B`: raw values on the subjects of the test), the
  difference (adjusted with covariates), and its effect size `g` with a
  bootstrap CI: Hedges' g, or with covariates the adjusted difference
  over the pooled SD of the raw values (times Hedges' J); `g_unadjusted`,
  the raw values' g, with its CI.
  Correlations: partial correlations (`S.correlation_type`: Pearson, or
  Spearman on ranks), in the subjects and with the covariates of each
  correlation set, tested by the slope of the feature on the score over
  its HC3 standard error (`t_hc3`), with permutations of the score
  residualised on the covariates (Smith's method, Winkler et al. 2014)
  or (`S.correlation_test =
  'wild_bootstrap'`) a restricted wild bootstrap: more robust than r to a
  spread of the feature that changes with the score, but with fewer than
  about 100 recordings and a few extreme scores the permutation test can
  still reject too often; the wild bootstrap or Spearman are then safer
  (`docs/validation.md`, "Correlations: the studentised slope"). Max-statistic over the features of a family.
- **Multiple testing.** A family is one band × inverse method × measure ×
  feature family × comparison or correlation set. Its tests run on the
  same permutations, and the family's largest statistic is compared with
  the largest of each permutation (max-statistic, Westfall–Young):
  `family_p` (the family as a whole) and `p_maxT_family` (each test).
  Families are selected by Benjamini–Hochberg over `family_p`, all the
  families of the run together, every band included (`families.csv`); a
  test of a selected family is significant when `p_bb` = `p_maxT_family`
  × F / R ≤ `S.alpha` (Benjamini–Bogomolov). `q_fdr`, the FDR q within the
  family, describes. A family can be selected only when 1 / (`S.n_perm` +
  1) ≤ `S.alpha` / F: with the six bands and several comparisons and
  correlation sets, 10,000 permutations can be too few (the OpenNeuro
  example uses 50,000; `docs/validation.md`, "Permutation resolution");
  a run where no family can be selected warns
  (`run_statistics:resolution`, `psd:resolution`). `families.csv` gives
  the 95% Clopper–Pearson interval of each max-statistic `family_p` from
  the number of permutations (`family_p_mc_low`, `family_p_mc_high`): a
  family whose interval holds `S.alpha` × R / F may be selected or not
  with another seed.
- **Edges.** The network-based statistic of every comparison (NBS 1.2's
  `NBSstats` on the z of every edge's Welch-type G, with Freedman–Lane
  permutations, primary threshold `S.nbs.thresh` on z, both directions; p = (1 + b)/(1 +
  permutations), and `p_both` = 2p, the one to compare with alpha);
  where the edge differences concentrate (network enrichment against
  random and spin-rotated network labels, `S.enrichment`); every edge's
  group means, difference, G and z; figures of the edge statistics and of
  each component with `p_both` < alpha (connectogram and brain view). The
  NBS p of a component is corrected over the edges by its permutation
  test (largest component), not over measures and bands; `q_run`
  (Benjamini–Hochberg over the components of every NBS of the run, in
  each `NBS_<matrix>.csv` and `NBS/components.csv`) is the exploratory
  view over all of them.
- **Failures.** A step that fails (a family, the features, the
  reliability, the NBS, a figure) does not stop the others: it is named
  on the console and in `Statistical analysis/failures.csv` (a family
  also in the column `status` of `summary.csv`), and the run ends "Done
  with errors"; when only figures failed, it ends "every result written,
  but N figure step(s) failed" (`S.stop_on_error = true` stops at the
  first one instead). The adjustment over families covers the families
  that ran.
- **Spectral profile** (`S.spectral_profile`): each recording at the
  scalp (IAF and alpha peak, 1/f exponent and offset, alpha power,
  muscle (temporal channels, T3/T4 included) and eye indices), tested as the family `spectral` and joined to
  the subject table, so any measure can be a covariate. On a cap without
  the 10-05 names of a set (an ANT equidistant or EGI cap), the set takes
  the channels nearest to their positions, within 20 mm (`cap_channels.m`).

## Reliability

Stage 1 computes every matrix also on the odd and on the even blocks of
`S.reliability_block` s of each recording (blocks counted from the start
of the recording, so a cut does not shift them; interleaved, so slow
drifts do not separate the halves; dwPLI over the epochs of each half,
NaN in a half of fewer than 10). `static_reliability.csv` gives, for
every feature, the correlation of the two halves across recordings (r)
and the Spearman–Brown reliability of the whole-recording value (`r_sb`);
for the edges, the median and quartiles of their `r_sb`; and how
reproducible each recording's pattern is (the median within-recording
correlation of the halves, and the share of recordings identified from
their other half). A comparison or correlation of a feature with low
reliability (`r_sb` below about 0.5) is attenuated: read it with care.

The window mean's halves are the mean of the windows that lie wholly
inside each half's blocks, each measure in blocks of
max(`S.reliability_block`, 3 × its window) (`static_halves_block`, s): the
halves share no data and about 2/3 of the windows are used (AEC's 6 s
windows in 18 s blocks, dwPLI's windows of 60 cycles in blocks of three (18 s in alpha), the phase
measures' short windows in 10 s blocks). The statistics leave out the
halves of outputs made before 1.2.0 (10 s blocks for every measure) and
say to run stage 1 again. A reliability is computed only when at least
max(3, half of the recordings) have both halves (NaN otherwise, with the
count in `note` and a warning); the log gives the block length, how many
recordings have both halves and why the others have none. The halves are
interleaved blocks of one session: within-session consistency, an upper
bound of test-retest reliability, with the groups pooled (a group
difference raises r): `r_within_groups` is r with each group's mean
removed from both halves. A stable bias is reliable too ("Which matrix
for phase measures"): reliability does not show that a feature measures
coupling.

## PSD analysis

The power spectra of every channel and region, their band power (absolute,
relative, above the 1/f background), 1/f exponent and offset and alpha
peak frequency (the centre of gravity of the alpha power above the 1/f
fit apart, untested), tested per channel and region with the statistics
above and a cluster test: `docs/psd.md`. Narrow lines (mains or
equipment) are flagged per recording (`recordings.csv`, `line_hz`,
`line_ratio`), the mains frequency is `S.psd.power_line`, and a spectrum
rising over the fit range keeps a negative 1/f exponent.

The power of every source point (vertex) comes from the channels'
cross-spectra and each source's filter, as LORETA's current density
(Frei et al. 2001): P_v(f) = w_v' C(f) w_v, the sum over the three
orientations of a free-orientation source. A region's signal (stage 1's
sign-flipped mean of its sources) cancels sources of opposite phase, so
its power is at most the mean of its sources' powers;
`region_power_check.csv` gives their ratio, region by region, and
`S.psd.region_spectra = 'vertex_mean'` runs the statistics on the mean of
the sources' powers instead. Maps of every source's band power (group
means, each comparison's g and each correlation's r) are descriptive, not
tested (`S.psd.vertex_maps`). Checked against FieldTrip's eLORETA source
power and Brainstorm's PSD of the source maps (`tests/test_psd.m`).
