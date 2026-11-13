# PSD analysis (stage 2)

Band power of every channel and every cortical region, compared between
the groups and correlated with clinical scores, with scalp and cortex
maps. Code: `2- PSD analysis/` (main function `run_psd_analysis.m`);
settings: `S.psd` (`pipeline_settings.m`, section "PSD analysis").

```matlab
run_full_pipeline(S);        % stage 1, the PSD analysis, then the states (if any) and the statistics
run_psd_analysis(S);         % the PSD analysis alone (it does not need stage 1's outputs)
run_full_pipeline(S, 'from', 2, 'to', 2);         % the same, through the runner
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
  with Hann segments of `S.psd.segment` s (4 s: 0.25 Hz resolution, fine
  enough for 3 Hz bands such as alpha1 and alpha2) and
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
- with the regions, the PSD of **every source** (vertex) of the cortex:
  P(f) = w′ C(f) w, from the channels' cross-spectra C
  (`psd_cross_spectra.m`: the same Welch segments and weights as the
  channels' PSD, which is its diagonal; MATLAB's `cpsd` of each piece) and
  the source's inverse filter w (`psd_vertex_spectra.m`; with
  `S.orientation = 'free'` the sum over its three orientations, the power
  of the current density). This is how LORETA makes its band images, from
  the cross-spectral matrix (Frei et al. 2001), and what Brainstorm's PSD
  of a full source map and FieldTrip's eLORETA on cross-spectra
  (`ft_inverse_eloreta`) give. Each region's **mean of its sources'
  PSDs** is kept beside the PSD of its signal; with `S.psd.vertex_maps`
  (default) every source's PSD is kept too, about 10 MB per recording and
  inverse (15,002 sources, 1–45 Hz at 0.25 Hz).
- the **specparam fit** of every spectrum over `S.psd.fit_range` (2–40 Hz;
  `psd_specparam.m`: Brainstorm's `process_fooof('FOOOF_matlab')`, the
  'fixed' aperiodic mode, Brainstorm's default peak settings; the mains
  frequency `S.psd.power_line` (50 Hz) and its harmonics ± 2 Hz left out
  where they lie in the fit range, which with 2–40 Hz they do not), as the
  spectral profile of the statistics (`S.spectral_profile`). A spectrum
  that rises over the fit range gets a negative exponent (specparam's
  MATLAB code returns the exponent's absolute value; the sign is taken
  from its aperiodic curve).

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
| `absolute` (per band) | log10 of the band power (band power is skewed): the integral of the PSD from the band's lower to its upper edge (trapezoids with interpolated edges, so that adjacent bands add up exactly) | log10 µV² at the scalp |
| `relative` (per band) | the band power in % of the power over `S.psd.total` (1–30 Hz) | % |
| `periodic` (per band) | the power above the 1/f background: mean of log10 PSD minus the specparam aperiodic fit over the band's part inside the fit range, from edge to edge (its trapezoid integral over its width, so a band is the width-weighted mean of its halves, whatever grid points an anchored band's edges fall between); NaN when less than 1 Hz of the band is inside | log10 (0.3 = twice the background) |
| `aperiodic` | exponent and offset of the aperiodic fit (log10 P = offset − exponent · log10 f) | |
| `alpha_frequency` | the centre of the strongest specparam peak inside `S.iaf.search` (7–13 Hz); NaN without one. The centre of gravity of the power above the fit in that range (when that power is at least 5% of the range's power) is reported in `values_<level>_alpha_cog.csv`, not tested: on spectra without alpha it lies near the middle of the range, so a group with weaker alpha would move toward it | Hz |

The bands are the PSD's own (`S.psd.bands`, not the connectivity bands of
`S.bands`): delta 1–4, theta 4–7, alpha 7–13, alpha1 7–10, alpha2 10–13,
beta 13–30, beta1 13–20, beta2 20–30 Hz. An
edge can be anchored on each recording's IAF as in `S.bands_individual`
(e.g. `'fmin', 'iaf-2'`; `band_edges.m`, `S.iaf`); the IAF is then
computed at the scalp as in stage 1 and listed in `recordings.csv`.
Channels: those present in every recording (a warning names those left
out). Regions: the atlas's.

**3. Statistics** (`psd_statistics.m`): the methods of the statistics
(`permu_ttest_permu_corr_BNS.m`). Each level × measure × band is a family
(`psd_<level>_<measure>_<band>`, e.g. `psd_sensor_relative_alpha`), whose
tests are the channels or regions:

- every comparison of `S.comparisons`: Welch t, or ANCOVA with
  `S.covariates` and Freedman-Lane permutations, the statistic of the data
  and of every permutation turned into z on its own Welch-Satterthwaite
  (Aspin-Welch with covariates) degrees of freedom, as PALM does (its
  `fastv` and `palm_gtoz`; Winkler et al. 2014); the effect size `g` and
  its CI (Hedges' g, adjusted with covariates) and `g_unadjusted`;
- every correlation set of `S.correlations` and `S.covariate_sets`:
  (partial) Pearson or Spearman r (`S.correlation_type`), tested with the
  slope over its HC3 standard error and permutations or a wild bootstrap
  (`S.correlation_test`);
- `S.psd.n_perm` permutations (`S.n_perm` when empty); `p_unc`, `p_maxT`
  (max-T over the channels or regions of the measure), `p_maxT_family`
  and `family_p` (max-statistic over the family: the comparison, or the
  correlation set with all its scores), `q_fdr` (Benjamini-Hochberg within
  the family, a description), and `p_bb`, `significant`: the families of
  the PSD selected by Benjamini-Hochberg over their `family_p`, and
  `p_bb` = `p_maxT_family` × F / R ≤ `S.alpha` inside them, over the PSD
  families only (`family_adjust.m`; the PSD is its own question, not
  adjusted with the connectivity families of the statistics);
- with `S.psd.cluster.run` (default), a **cluster test** on the same
  permutations: FieldTrip's (`ft_statistics_montecarlo` with `correctm =
  'cluster'`, `clusterstatistic = 'maxsum'`, a parametric threshold, two
  tails with `correcttail = 'prob'`; Maris & Oostenveld 2007), adapted in
  `cluster_perm_test.m` from FieldTrip's `clusterstat` to take every
  permutation at once. Neighbouring channels or regions whose statistic
  passes the two-sided `S.psd.cluster.alpha` (0.05) critical value form
  clusters of one sign (group comparisons: their z, |z| ≥ 1.96, PALM's
  cluster-forming threshold on the z scale; correlations: the slope over
  its HC3 standard error, as z on its degrees of freedom, |z| ≥ 1.96 too);
  a cluster's mass is the sum
  of its statistics; a positive cluster's p is the share
  of permutations whose largest positive cluster is larger (a negative
  one's, whose most negative one is smaller), doubled for the two tails.
  `tests/test_psd.m` checks the clusters, masses and p against FieldTrip's
  own `clusterstat`. Neighbours: channels by FieldTrip's triangulation of
  the cap (`ft_prepare_neighbours`), regions with adjacent vertices on the
  cortex (Brainstorm's `VertConn` of the cortex surface), never across
  hemispheres (`psd_adjacency.m`; listed in `neighbours_<level>.csv`).
  On a sparse or irregular cap the triangulation can still join distant
  channels at its edge (links up to the mean + 3 SD of the link lengths
  are kept): check `neighbours_sensor.csv`, or use the regions. The p of all clusters of the PSD are adjusted
  together: Benjamini-Hochberg over every cluster of every level, measure,
  band and test (`cluster_q`; `cluster_significant` = `cluster_q` ≤
  `S.alpha`), outside the adjustment over families of the other tests
  (exploratory, as the NBS).
  A cluster says that the effect is spread over neighbouring channels or
  regions; it does not locate it at any single one of them.

The subjects are those of the subject table with a usable recording; a
subject with two recordings stops the statistics (as it stops run_statistics). Without
comparisons and correlations only the spectra and measures are written.
Before the tests the log gives the sample size of every comparison and
correlation set with the covariates each can use (`check_sample_sizes.m`,
as run_statistics) and the permutation resolution (`resolution_note.m`): a
family can be selected only when 1 / (n_perm + 1) ≤ alpha / families,
whatever its size. The PSD has many families (levels × measures × bands ×
comparisons and correlation sets: 324 families, for instance, need at
least 6,479 permutations at `S.alpha` = 0.05), hence its own
`S.psd.n_perm` when run_statistics needs fewer; the line says how many permutations may reach a family's
largest statistic for it to be selected, and how many permutations would
allow 9 (`docs/validation.md`, "Permutation resolution"). The counts of
each family's line are `q_fdr` < alpha within the family; the last line
gives the tests significant after the adjustment over the PSD families,
and the figures mark only those. A
covariate that is a measure of the spectral profile (e.g.
`alpha_pw_post`), which run_statistics computes later, is computed first
(`prepare_spectral_covariates.m`, the same cache as run_statistics); with
`S.spectral_profile = false` the statistics stop with a message naming it.

**4. Figures** (`psd_figures.m`, `S.figures`), per level and comparison:

- `spectra_<level>_<cmp>.png`: each group's mean spectrum (± 95% CI) over
  the channels or regions, absolute and above the 1/f fit, the bands
  marked;
- `boxplots_<level>_<measure>_<cmp>.png`: each band's measure averaged
  over the channels or regions, per group, with Hedges' g (and
  `boxplots_<level>_aperiodic_alpha_<cmp>.png`: exponent, offset, alpha
  frequency);
- `map_<level>_<measure>_<cmp>.png`: maps of the tests' g (Hedges' g;
  the adjusted g with `S.covariates`) for every band,
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

**5. A region's signal against its sources' power** (`psd_region_check.m`,
always with the regions). A region's signal is the sign-flipped mean of
its sources (stage 1's scout function `'mean'`): sources whose activity is
out of phase cancel in it, and its power is at most the mean of its
sources' power (Cauchy-Schwarz; equal when they act together).
`region_power_check.csv` gives, per inverse, region and band, the ratio of
the two band powers over the recordings (`ratio_median`, `ratio_q1`,
`ratio_q3`), their correlation over the recordings for the absolute,
relative and periodic power (`r_absolute`, ...), and for each comparison
the unadjusted Hedges' g with each (`g_<measure>_signal`,
`g_<measure>_vertex_mean`; descriptive). `figures/region_power_check_<inverse>.png`
maps the median ratio; the log names the five lowest regions. A low ratio
is not an error: the region's signal is what the connectivity of stage 1
uses, and it keeps the part of the activity common to the region. A
region whose g differs between the two holds sources that do not act
together; its power is then better read from the mean of its sources'
power. `S.psd.region_spectra = 'vertex_mean'` runs the regions' measures,
statistics, clusters and figures on that mean instead (levels
`source_<inverse>_vertexmean`); the default `'signal'` keeps the PSD of
the signal stage 1 uses. The kept spectra hold both, so switching
recomputes no spectrum.

**6. Maps of every source's power** (`psd_vertex_maps.m`,
`S.psd.vertex_maps`, default true): from every source's PSD, its
`absolute` and `relative` power per band (as in 2, each recording's band
edges), in `vertex_values_<inverse>.mat` (recordings × sources × bands),
with each comparison's group means and Hedges' g (unadjusted) and each
correlation's r (`S.correlation_type`, no covariates), as LORETA's maps.
They are **descriptive, not tested**: the tests are on the channels and
regions. Figures: `vertex_<measure>_<inverse>_<cmp>.png` (each group's
mean and g per band; left, top and right views), the g maps in 3-D
(`.fig`), and `vertex_<measure>_<inverse>_r_<set>_<score>.png` (up to 10
scores). They need one cortex for everyone (the shipped head, or
`S.anatomy` without `S.anatomy_individual`); otherwise they are skipped
with a message, and steps 1 to 5 are unchanged.

## Outputs (`results/PSD/<S.run_name, or main>/`)

| File | Content |
|---|---|
| `psd_log.txt`, `provenance_psd.mat` | the log (channels left out, failures, the tests), the settings, recordings, families and summary |
| `recordings.csv` | one row per recording: status, rate, seconds analysed, pieces, Welch segments, channels, IAF (and its flag), the strongest narrow line (`line_hz`, `line_ratio`, `line_channels`: see Line noise), the bands whose periodic power starts less than 2 Hz above the start of the fit (`near_fit_start`: see Reading the measures), and per level the median specparam R², the share of channels/regions with an alpha peak and the number whose spectrum rises over the fit range (`rising_<level>`) |
| `values_<level>_<measure>.csv` | the measures: one row per recording, a column per band and channel or region (`alpha__O1`) or per channel or region (exponent, offset, alpha frequency, `alpha_cog` (not tested), fit R²) |
| `psd_values.mat` | the same and the spectra (recording × channel/region × frequency), for other plots |
| `Statistical analysis/summary.csv`, `families.csv`, `clusters.csv` | tests and significant results per family; the families selected; every cluster (family, test, sign, size, mass, p, q, mean effect, its channels or regions) |
| `Statistical analysis/<level>/<family>_group_comparison.csv`, `_correlations.csv` (.xlsx too) | one row per channel or region and test, with `Level`, `Measure`, `Band`, `Unit`, the statistics above and the cluster columns |
| `Statistical analysis/neighbours_<level>.csv` | the neighbours used by the cluster test |
| `region_power_check.csv` | per inverse, region and band: the region signal's band power over the mean of its sources' (median, quartiles), their correlations, g with each (step 5) |
| `vertex_values_<inverse>.mat` | every source's absolute and relative band power per recording, group means, g and r (step 6) |
| `figures/` | the figures above (`.png`; the cortex maps also as rotatable `.fig`) |

`<level>` is `sensor` or `source_<inverse>` (e.g. `source_eloreta`;
`source_<inverse>_vertexmean` with `S.psd.region_spectra = 'vertex_mean'`).

## Reading the measures

- `relative` power is a share: a band can change in relative power
  because another band changed (more alpha lowers every other band's
  share). Read it with `absolute` and `periodic`. A band outside
  `S.psd.total` (gamma with the default 1–30 Hz) gets a ratio to that
  total, not a share of it (a warning, `check_settings:psd_total`).
- `absolute` power depends on the 1/f background too: a change of the
  exponent or offset moves every band. `periodic` is the part above the
  background.
- `periodic` power in the first hertz or two of the fit range depends on
  where the fit starts: the aperiodic fit ('fixed', no knee) is least
  constrained at its edge, and the data are often high-pass filtered not
  far below it. A band whose part inside the fit starts less than 2 Hz
  above `S.psd.fit_range(1)` (delta with the defaults; an anchored band
  such as IAF − 6 to IAF − 4 Hz in recordings with a low IAF) is named
  in a warning (`psd:fit_start`, with its number of recordings) and in
  `recordings.csv` (`near_fit_start`). Before reading a periodic result
  of such a band, repeat it with a later fit start (e.g.
  `S.psd.fit_range = [3 40]`) and check that the groups do not differ in
  the recordings concerned: an anchored band reaches the fit start only
  in recordings with a low IAF, which may be more frequent in one group.

## Checks before reading the results

- `recordings.csv`: every recording `ok`; `fit_r2_*` above about 0.9 (a
  poor fit makes `periodic`, the exponent and the offset unreliable);
  `alpha_peak_share_*` (eyes-closed rest: most channels; eyes open: fewer).
- the log: channels left out (missing in some recordings), recordings
  left out.
- `summary.csv`: how many tests each family had (e.g. `alpha_frequency`
  has fewer when some channels have no alpha).
- `line_ratio` and `rising_<level>` in `recordings.csv`, and the warnings
  `psd:line` and `psd:rising` in the log (next section).

## Line noise and other narrow lines

Removing lines belongs to preprocessing, before the pipeline: for example
spectrum interpolation (FieldTrip's `ft_preproc_dftfilter` with
`'neighbour_fft'`) or a notch on each continuous piece. The PSD step
does not remove them, but looks for them (`psd_line_check.m`): on the
mean of each recording's channel spectra, a narrow peak is a local
maximum whose excess over the background (the median power 4–8 bins, 1–2
Hz, away) falls to at most 10% of itself 2 bins (0.5 Hz) away, where a
sinusoid's spectrum in Hann segments has its first zero. `line_hz` is its
frequency and `line_ratio` its power over the background (1: none); a
ratio of 2 or more (the line at least doubles its bin) gives the warning
`psd:line`. `line_channels` counts the channels whose own spectrum has
such a line. Rhythms such as alpha are broader than a line: they keep
much of their excess 0.5 Hz away, so they are not taken for lines.

A line in the analysed range enters the band powers, the total of
relative power, the periodic power and the specparam fit, where it can
take one of the three peaks. Remove it, or compare the results without the
flagged recordings. `S.psd.power_line` (`'50'`, `'60'` or `'none'`) is the
mains frequency whose first three harmonics ± 2 Hz the fit leaves out when
they lie in the fit range (with the default 2–40 Hz none does); it does
not remove a line elsewhere.

## Validation

- `tests/test_psd.m` (11 tests): Welch scaling (the integral of white noise
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
- Steps 5 and 6: the cross-spectra's diagonal
  equals the channels' PSD (1e-15) and the matrix MATLAB's `cpsd` of each
  piece; every source's PSD equals the PSD of its signal w′X (eLORETA,
  wMNE, LCMV; normal and free orientation; 1e-9), eLORETA's equals
  FieldTrip's `ft_inverse_eloreta` on the same cross-spectra (normal
  orientation; free: between its largest orientation's power and the sum
  of the three); each region's signal power is at most its sources' mean;
  in the simulated recordings the occipital sources have the largest g, g
  equals `hedges_g` at a source, and `'vertex_mean'` runs from the kept
  spectra.
- `tests/test_static_integration.m`: the step inside `run_full_pipeline`
  on the three file formats of the BIDS test dataset.
