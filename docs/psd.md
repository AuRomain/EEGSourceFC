# PSD analysis (stage 1b)

Band power of every channel and every cortical region, compared between
the groups and correlated with clinical scores, with scalp and cortex
maps. Code: `1b- PSD analysis/` (main function `run_psd_analysis.m`);
settings: `S.psd` (`pipeline_settings.m`, section "PSD analysis").
Written on 5 November 2026 from the lab's original scripts
(`v_original/1b- PSD analysis`, not distributed; what changed and why:
the last section).

```matlab
run_full_pipeline(S);        % stage 1, the PSD analysis, then the states (v_dynamic) and the statistics
run_psd_analysis(S);         % the PSD analysis alone (it does not need stage 1's outputs)
run_full_pipeline(S, 'from', '1b', 'to', '1b');   % the same, through the runner
```

`S.psd.run = false` switches it off in `run_full_pipeline`.

## What it computes

**1. Spectra** (`psd_spectra.m`). Each recording is read as stage 1 reads
it (`get_eeg_data`: the analysed stretch `S.core_start`/`S.core_duration`,
its continuous pieces; the channels of the first cap `S.montages{1}`,
`S.channel_subset`, flat channels left out, at least `S.min_channels`;
average reference). Then:

- `'sensor'`: the PSD of each channel by Welch's method: MATLAB's
  `pwelch` on each continuous piece, its mean removed (`psd_welch.m`),
  with Hann segments of `S.psd.segment` s (4 s: 0.25 Hz resolution) and
  `S.psd.overlap` (50%), so no segment crosses a cut; the pieces averaged
  with their number of segments as weight (the mean over all segments);
  one-sided, in µV²/Hz, over `S.psd.freq_range` (1–45 Hz) on a common
  grid. Always computed: it also gives the IAF.
- `'source'`: the same for each region of the first atlas
  (`S.atlases{1}`), for every inverse method of `S.inv`. The region
  signals are made exactly as in stage 1 (`psd_source_signals.m`: the
  same head model, the leadfield restricted to the channels used and
  average-referenced over them, the same inverse operator and scout
  function); the test `test_source_signals_as_stage1` checks that they
  equal stage 1's saved scouts.
- the **specparam fit** of every spectrum over `S.psd.fit_range` (2–40 Hz;
  `psd_specparam.m`: Brainstorm's `process_fooof('FOOOF_matlab')`, the
  'fixed' aperiodic mode, Brainstorm's default peak settings, 50 Hz line
  ignored), as the spectral profile of the statistics (`S.spectral_profile`).

The spectra are kept in `results/PSD/spectra/<sensor|source>/` with a key
(the recording's file content, the settings each level depends on, the
code of stage 1 and of this step): a level is recomputed only when its
key changes; changing a band, a measure or a test recomputes nothing. A
recording that fails (too few channels, no piece as long as a segment) is
reported in `recordings.csv` and left out with a warning; a failure of
the regions keeps the channels.

**2. Measures** (`psd_measures.m`), for each channel or region:

| Measure | What | Unit |
|---|---|---|
| `absolute` (per band) | log10 of the band power: the integral of the PSD from the band's lower to its upper edge (trapezoids with interpolated edges, so that adjacent bands add up exactly) | log10 µV² at the scalp |
| `relative` (per band) | the band power in % of the power over `S.psd.total` (1–30 Hz) | % |
| `periodic` (per band) | the power above the 1/f background: mean of log10 PSD minus the specparam aperiodic fit, over the band's part inside the fit range (NaN when less than 1 Hz of the band is inside) | log10 (0.3 = twice the background) |
| `aperiodic` | exponent and offset of the aperiodic fit (log10 P = offset − exponent · log10 f) | |
| `alpha_frequency` | the centre of the strongest specparam peak inside `S.iaf.search` (7–13 Hz); without one, the centre of gravity of the power above the fit in that range when that power is at least 5% of the range's power; NaN otherwise (a centre of gravity of noise means nothing) | Hz |

The bands are the PSD's own (`S.psd.bands`, not the connectivity bands of
`S.bands`): delta 1–4, theta 4–7, alpha 7–13, alpha1 7–10, alpha2 10–13,
beta 13–30, beta1 13–20, beta2 20–30 Hz, as in the original scripts. An
edge can be anchored on each recording's IAF as in `S.bands_individual`
(e.g. `'fmin', 'iaf-2'`; `band_edges.m`, `S.iaf`); the IAF is then
computed at the scalp as in stage 1 and listed in `recordings.csv`.
Channels: those present in every recording (a warning names those left
out). Regions: the atlas's.

**3. Statistics** (`psd_statistics.m`): the methods of stage 3
(`permu_ttest_permu_corr_BNS.m`). Each level × measure × band is a family
(`psd_<level>_<measure>_<band>`, e.g. `psd_sensor_relative_alpha`), whose
tests are the channels or regions:

- every comparison of `S.comparisons`: Welch t, or ANCOVA with
  `S.covariates` and Freedman-Lane permutations, the statistic of the data
  and of every permutation turned into z on its own Welch-Satterthwaite
  (Aspin-Welch with covariates) degrees of freedom, as PALM does (its
  `fastv` and `palm_gtoz`; Winkler et al. 2014); Hedges' g and its CI;
- every correlation set of `S.correlations` and `S.covariate_sets`:
  partial Pearson r with permutations;
- `S.psd.n_perm` permutations (`S.n_perm` when empty); `p_unc`, `p_maxT`
  (max-T over the channels or regions of the measure), `p_maxT_family`
  and `family_p` (max-statistic over the family: the comparison, or the
  correlation set with all its scores), `q_fdr` (Benjamini-Hochberg within
  the family, a description), and `p_bb`, `significant`: the families of
  the PSD selected by Benjamini-Hochberg over their `family_p`, and
  `p_bb` = `p_maxT_family` × F / R ≤ `S.alpha` inside them, over the PSD
  families only (`family_adjust.m`; the PSD is its own question, not
  adjusted with the connectivity families of stage 3);
- with `S.psd.cluster.run` (default), a **cluster test** on the same
  permutations: FieldTrip's (`ft_statistics_montecarlo` with `correctm =
  'cluster'`, `clusterstatistic = 'maxsum'`, a parametric threshold, two
  tails with `correcttail = 'prob'`; Maris & Oostenveld 2007), adapted in
  `cluster_perm_test.m` from FieldTrip's `clusterstat` to take every
  permutation at once. Neighbouring channels or regions whose statistic
  passes the two-sided `S.psd.cluster.alpha` (0.05) critical value form
  clusters of one sign (group comparisons: their z, |z| ≥ 1.96, PALM's
  cluster-forming threshold on the z scale; correlations: r turned to t,
  the critical t of its degrees of freedom); a cluster's mass is the sum
  of its statistics; a positive cluster's p is the share
  of permutations whose largest positive cluster is larger (a negative
  one's, whose most negative one is smaller), doubled for the two tails.
  `tests/test_psd.m` checks the clusters, masses and p against FieldTrip's
  own `clusterstat`. Neighbours: channels by FieldTrip's triangulation of
  the cap (`ft_prepare_neighbours`), regions with adjacent vertices on the
  cortex (Brainstorm's `VertConn` of the cortex surface), never across
  hemispheres (`psd_adjacency.m`; listed in `neighbours_<level>.csv`). The p of all clusters of the PSD are adjusted
  together (Benjamini-Hochberg: `cluster_q`, `cluster_significant`).
  A cluster says that the effect is spread over neighbouring channels or
  regions; it does not locate it at any single one of them.

The subjects are those of the subject table with a usable recording; a
subject with two recordings stops the statistics (as stage 3). Without
comparisons and correlations only the spectra and measures are written.
Before the tests the log gives the sample size of every comparison and
correlation set with the covariates each can use (`check_sample_sizes.m`,
as stage 3) and the permutation resolution (`resolution_note.m`): a
family can be selected only when 1 / (n_perm + 1) ≤ alpha / families,
whatever its size. The PSD has many families (levels × measures × bands ×
comparisons and correlation sets: 324 in the OpenNeuro example, which
needs at least 6,479 permutations and uses 50,000), hence its own
`S.psd.n_perm` when stage 3 needs fewer; the line says how many permutations may reach a family's
largest statistic for it to be selected, and how many permutations would
allow 9 (`docs/validation.md`, "Permutation resolution"). The counts of
each family's line are `q_fdr` < alpha within the family; the last line
gives the tests significant after the adjustment over the PSD families,
and the figures mark only those. A
covariate that is a measure of the spectral profile (e.g.
`alpha_pw_post`), which stage 3 computes later, is computed first
(`prepare_spectral_covariates.m`, the same cache as stage 3); with
`S.spectral_profile = false` the statistics stop with a message naming it
(before, the covariate was missing on a fresh run and every test failed).

**4. Figures** (`psd_figures.m`, `S.figures`), per level and comparison:

- `spectra_<level>_<cmp>.png`: each group's mean spectrum (± 95% CI) over
  the channels or regions, absolute and above the 1/f fit, the bands
  marked;
- `boxplots_<level>_<measure>_<cmp>.png`: each band's measure averaged
  over the channels or regions, per group, with Hedges' g (and
  `boxplots_<level>_aperiodic_alpha_<cmp>.png`: exponent, offset, alpha
  frequency);
- `map_<level>_<measure>_<cmp>.png`: maps of Hedges' g for every band,
  on the scalp (FieldTrip topography of the cap) or on the cortex (six
  views: left lateral and medial, right medial and lateral, the whole
  cortex from above and from below, where the left is on the right; each
  region coloured by its g, those neither significant nor in a
  significant cluster paler), significant channels or regions
  (`significant`) marked with a filled dot, those in a significant cluster
  with an open circle; the same for the exponent, offset and alpha frequency
  (`map_<level>_aperiodic_alpha_<cmp>.png`); maps of r for a correlation
  set and score when something is significant;
- `map_source_<inverse>_<measure>_<cmp>.fig` (and for the correlation
  maps): the same cortex maps in 3-D, one brain per band, to turn with the
  mouse in MATLAB (double-click the file, or `openfig`); the light follows
  the view. About 5 MB each;
- `cluster_<level>_<n>.png`: each significant cluster's values, per group.

## Outputs (`results/PSD/<S.run_name, or main>/`)

| File | Content |
|---|---|
| `psd_log.txt`, `provenance_psd.mat` | the log (channels left out, failures, the tests), the settings, recordings, families and summary |
| `recordings.csv` | one row per recording: status, rate, seconds analysed, pieces, Welch segments, channels, IAF (and its flag), median specparam R² and share of channels/regions with an alpha peak per level |
| `values_<level>_<measure>.csv` | the measures: one row per recording, a column per band and channel or region (`alpha__O1`) or per channel or region (exponent, offset, alpha frequency, fit R²) |
| `psd_values.mat` | the same and the spectra (recording × channel/region × frequency), for other plots |
| `Statistical analysis/summary.csv`, `families.csv`, `clusters.csv` | tests and significant results per family; the families selected; every cluster (family, test, sign, size, mass, p, q, mean effect, its channels or regions) |
| `Statistical analysis/<level>/<family>_group_comparison.csv`, `_correlations.csv` (.xlsx too) | one row per channel or region and test, with `Level`, `Measure`, `Band`, `Unit`, the statistics above and the cluster columns |
| `Statistical analysis/neighbours_<level>.csv` | the neighbours used by the cluster test |
| `figures/` | the figures above (`.png`; the cortex maps also as rotatable `.fig`) |

`<level>` is `sensor` or `source_<inverse>` (e.g. `source_eloreta`).

## Reading the measures

- `relative` power is a share: a band can change in relative power
  because another band changed (more alpha lowers every other band's
  share). Read it with `absolute` and `periodic`.
- `absolute` power depends on the 1/f background too: a change of the
  exponent or offset moves every band. `periodic` is the part above the
  background.

## Checks before reading the results

- `recordings.csv`: every recording `ok`; `fit_r2_*` above about 0.9 (a
  poor fit makes `periodic`, the exponent and the offset unreliable);
  `alpha_peak_share_*` (eyes-closed rest: most channels; eyes open: fewer).
- the log: channels left out (missing in some recordings), recordings
  left out.
- `summary.csv`: how many tests each family had (e.g. `alpha_frequency`
  has fewer when some channels have no alpha).

## Validation (5 November 2026)

- `tests/test_psd.m` (8 tests): Welch scaling (the integral of white noise
  is its variance, a sine adds its power, segments never cross a cut,
  one piece gives `pwelch` itself);
  specparam on synthetic spectra (exponent, offset and peak recovered, no
  peak in a pure 1/f spectrum); the band measures (adjacent bands add up,
  relative powers sum to 100%); the cluster test (a planted cluster
  found, p ≈ 0.001; none in noise; the same clusters, masses and p as
  FieldTrip's `clusterstat` on the same permutations); the neighbours
  (sensor: symmetric, every channel with neighbours; regions: never
  across hemispheres, Brainstorm's `VertConn` and the cortex's triangles
  giving the same pairs);
  the region signals equal stage 1's; the settings' checks; and the whole
  step on 8 + 8 simulated recordings with a planted occipital alpha
  increase in the patients that a score follows: the occipital channels
  and regions found, the score correlation found, nothing in the
  absolute and periodic beta2 power (relative beta2 falls, correctly: it
  is a share of 1–30 Hz, which holds more alpha); a second run, even in a
  session that has not yet loaded Brainstorm's specparam, recomputes no
  spectrum. That test takes about 4 minutes.
- `tests/test_static_integration.m`: the step inside `run_full_pipeline`
  on the three file formats of the BIDS test dataset.
- Two datasets never seen before: NIRP (4 recordings, 27 channels with
  the old names; too few for tests: spectra, measures and figures only)
  and Geneva (6 recordings, a 129-channel ANT cap without 10-20 names;
  tests, clusters and maps on the scalp and the cortex) ran without an
  error.

## From the original scripts

The lab's original scripts (`v_original/1b- PSD analysis`, Geneva
project, not distributed) had two computation
scripts (sensor, cortical) and two statistics-and-figures scripts. Kept:
the eight bands, absolute and relative power (1–30 Hz total), the sensor
and cortical levels, the group comparison per channel or region with FDR,
maps of the group difference with the significant channels or regions
marked, the blue-white-red colour scale. Changed:

| Original | Here | Why |
|---|---|---|
| a computer's paths (`D:\...`), rate 500 Hz and 240 s hard-coded, files by name pattern | the settings of the pipeline (`S.data`, `S.core_duration`, any rate) | general |
| `pwelch(x, fs, [], f, fs)`: 1 s Hamming segments (1 Hz resolution) evaluated on a 1/fs grid, across cuts | 4 s Hann segments (0.25 Hz) inside the continuous pieces | with 1 Hz resolution the power leaks over about 2 Hz, across the edges of the 3 Hz sub-bands (alpha1, alpha2); segments across a boundary mix discontinuous data |
| a recording with a NaN skipped, but the previous recording's values saved under its name | NaN stretches cut the recording (get_eeg_data); a failure is reported and the recording left out | silent wrong values |
| cortical PSD from the scout files saved by stage 1 (`save_scouts`) | the region signals made again as stage 1 makes them, kept with their key | no need to save the scouts; consistent with the settings in use |
| sensor data read as saved, no reference | average reference, flat channels left out, channels common to all recordings | as stage 1; the same channels for every subject |
| Wilcoxon rank-sum per channel, FDR within band × measure, no covariates, no correlations in the run | Welch t or ANCOVA with covariates and permutations, Hedges' g, max-T, FDR within family, Benjamini-Bogomolov over the PSD families, correlations, cluster test | the methods of stage 3; effect sizes; covariates (age); the many bands × measures × levels adjusted |
| mean difference masked by FDR on electrode or region-centroid dots | Hedges' g maps on the scalp topography and on the cortex surface (six views, and a rotatable 3-D figure), every value shown, significant ones marked | non-significant effects visible; units comparable across measures |
| absolute power in µV² | log10 µV², plus the periodic power and the 1/f fit | band power is skewed; the 1/f background shifts all bands together (a difference of absolute power may be a difference of the exponent) |
