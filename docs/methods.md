# Methods of v_static

What each step computes, in more detail than the README. The code is v_dynamic's
(`../v_dynamic`) run with `S.dynamic = false`.

## Stage 1: from the recordings to connectivity

1. **Reading** (FieldTrip). Channels are matched by label to the cap's
   template positions; other channels, flat channels and missing ones are
   left out (the leadfield follows the channels kept). The data rank, the
   number of independent signals in the channels kept (eigenvalues of
   their covariance above 1e-12 of the largest; each ICA component removed
   or channel interpolated lowers it), is reported in `stage1_summary.csv`,
   with a warning below `S.min_channels`. The data are
   re-referenced to the average. Boundary events and stretches of NaN cut
   the recording into continuous pieces. `S.core_start` and
   `S.core_duration` choose the analysed stretch.
2. **Inverse solution** on the ICBM152 template (OpenMEEG boundary-element
   head model, 15,002 cortical sources oriented normal to the surface):
   eLORETA (`S.eloreta_lambda`), LCMV (`S.lcmv_lambda`) or wMNE
   (`S.wmne`). Each region's time series comes from its sources by
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
   - debiased wPLI (`dwpli`): observations are blocks of half the band's
     decorrelation time, products of nearby blocks left out;
   - AEC (correlation of the amplitude envelopes), `aec_orth` (each pair
     orthogonalised in both directions, then averaged), `aec_symorth` and
     `plv_symorth` (symmetric orthogonalisation of all regions first).
   The same is computed on the odd and on the even blocks of
   `S.reliability_block` seconds (`split_halves.m`) for the reliability,
   with band power of each region.
5. **Window mean** (the default; `S.static_matrices` including
   `'window_mean'`): the same measures in windows, then averaged
   (`get_all_connectivity_mats.m`, the windowed functions of
   `+connectivity`, v_dynamic's windows): phase measures in windows of
   `S.window_cycles` cycles of the band centre or `S.window_length`
   seconds, overlapping by `S.window_overlap`; AEC measures in windows of
   `S.aec_window` seconds every `S.aec_step` seconds. Windows never cross
   a cut; the analytic signal is that of step 3 (one Hilbert transform per
   piece). The mean of the windows inside each half gives the reliability
   of the window mean. A measure can have windows of its own
   (`S.measure_windows`, e.g. `struct('wpli', struct('window_length', 2))`);
   the others keep the common ones.

## Whole recording or windows

`S.static_matrices` chooses how each measure is estimated:
`'window_mean'` (the default since 6 November 2026, step 5), `'full'`
(step 4) or both, each analysed and reported on its own
(`features_window_mean.csv`, `features_full.csv`,
`static_reliability_window_mean.csv`, `static_reliability.csv`). Several
window settings for the same measure are separate runs, each with its own
`S.results`.

- `'full'` is how Brainstorm (time resolution "none"), FieldTrip and MNE
  estimate the static connectivity of continuous data, and the best
  estimate of a stable coupling: every sample adds to one estimate.
- `'window_mean'` is the original pipeline's "static" matrix. The
  estimate of a phase measure from a short window is biased upwards: in
  0.6 s windows of independent alpha-band noise wPLI is about 0.5 (over
  minutes, near 0), and the bias grows with the strength of the rhythm,
  so it follows band power. On the OpenNeuro recordings (138, alpha) the
  mean wPLI was 0.15 over the whole recording and 0.55 as the mean of
  0.6 s windows, and its rank correlation with band power across
  recordings rose from 0.48 to 0.67. Longer windows (`S.window_length =
  2` or more) reduce the bias; `dwpli` removes most of it. For AEC the
  two agree closely (edge by edge, rank correlation 0.91 across the
  OpenNeuro recordings), and windows of 6 s every 0.5 s are those of
  O'Neill et al. (2017).
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
  6 s so that the 6 s AEC windows fit in every piece kept (4 s, as in
  v_dynamic, is enough for `'full'`); on the OpenNeuro recordings, each a
  single piece of about 5 minutes, the two give the same data.

## Filters

Each continuous piece is band-passed with Brainstorm's filter (even-order
linear-phase FIR, 60 dB stop band, mirrored edges, delay compensated),
turned into its analytic signal, and then its edges are discarded: at
least `S.edge_trim` s and at least the filter's transient (the
half-length holding 99.9% of the impulse response's energy). With
`S.filter_transition = 'brainstorm'` (Brainstorm's default: 0.5 Hz
transition bands below 5 Hz, 1 Hz above) the bands are sharp and the
transients short enough:

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
transition bands; `'window'` (v_dynamic's default) the narrowest transition whose
filter stays within two windows of 6 cycles of the band; a number sets
the width in Hz.

## Statistics

For every band, inverse method and measure of `S.static_measures`
(`static_analysis.m`, `compute_NBS_analysis.m`):

- **Feature families** (`S.static_families`):
  - `global`: `Strength` (mean node strength), `ClustCoeff` (mean weighted
    clustering coefficient), `CPL` (characteristic path length, lengths
    1/weight); the Brain Connectivity Toolbox's definitions on the
    positive part of a signed measure such as AEC (`gt_window_metrics.m`;
    a test compares them with BCT when it is installed);
  - `networks`: mean connectivity inside each network and between each
    pair of networks of `S.network_table` (Yeo-7 for DK-68);
  - `regions`: node strength of every region;
  - `power`: log10 band power over all regions and each network (tested
    once per band and inverse method, since it does not depend on the
    measure).
- **Sample sizes.** Before the families, the recordings of every
  comparison and correlation set are counted (those of stage 1, made or
  reused, with a row in the subject table and every covariate): a comparison needs at least 2
  per group and the covariates + 4 in all (2 residual degrees of freedom),
  a correlation the covariates + 4; a warning names those that cannot be
  tested, whose rows would be NaN (comparisons) or absent (correlations).
- **Tests.** Group comparisons: Welch t, or with covariates ANCOVA with
  the Welch-type G statistic and Freedman–Lane permutations
  (`S.n_perm`), the statistic of the data and of every permutation
  turned into z on its own Welch–Satterthwaite or Aspin–Welch degrees of
  freedom before the p-values are counted, as PALM does (Winkler et al.
  2014); Hedges' g with its CI. Correlations: partial correlations
  with permutations, in the subjects and with the covariates of each
  correlation set. Max-statistic over the features of a family.
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
  1) ≤ `S.alpha` / F: with the six bands' 258 families of the OpenNeuro
  example, 10,000 permutations are not enough (its settings use 50,000).
- **Edges.** The network-based statistic of every comparison (NBS 1.2's
  `NBSstats` on the z of every edge's Welch-type G, with Freedman–Lane
  permutations, primary threshold `S.nbs.thresh` on z, both directions; p = (1 + b)/(1 +
  permutations), and `p_both` = 2p, the one to compare with alpha);
  where the edge differences concentrate (network enrichment against
  random and spin-rotated network labels, `S.enrichment`); every edge's
  group means, difference, G and z; figures of the edge statistics and of
  each component with `p_both` < alpha (connectogram and brain view). The
  NBS p of a component is corrected over the edges by its permutation
  test (largest component), not over measures and bands.
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
drifts do not separate the halves). `static_reliability.csv` gives, for
every feature, the correlation of the two halves across recordings (r)
and the Spearman–Brown reliability of the whole-recording value (`r_sb`);
for the edges, the median and quartiles of their `r_sb`; and how
reproducible each recording's pattern is (the median within-recording
correlation of the halves, and the share of recordings identified from
their other half). A comparison or correlation of a feature with low
reliability (`r_sb` below about 0.5) is attenuated: read it with care.

## PSD analysis

The power spectra of every channel and region, their band power (absolute,
relative, above the 1/f background), 1/f exponent and offset and alpha
frequency, tested per channel and region with the statistics above and a
cluster test: `docs/psd.md`.
