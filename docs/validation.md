# Validation and tests of EEGSourceFC

## Tests

```matlab
results = runtests('tests');  table(results)
```

No data are needed: the recordings are simulated on the pipeline's forward
model.

- `test_static_estimators`: every connectivity measure against its direct
  formula, on one piece and on several; the known couplings (lagged,
  zero-lag, independent); the split halves, and the window mean's halves
  in blocks of each measure's windows; the measures equal the
  whole-recording outputs of the sliding-window functions.
- `test_static_settings`: the static settings and their checks (the
  defaults, the warning for settings of the windows, the scout function,
  `S.lcmv_kappa`, the warning for dwPLI in bands anchored on the alpha
  frequency), the outputs folder, what changes the cache keys.
  `test_template`: this folder's template and the examples.
- `test_static_stage1`: stage 1 on simulated recordings in a BIDS layout:
  the outputs, the cache (skipped when nothing changed, redone when a
  setting changes), the same recording in EEGLAB, BrainVision, EDF and BDF
  formats, the individual bands of a recording without alpha, flat
  channels and channel subsets, the data rank (recordings with 20 and 50
  of their 64 dimensions removed; reused outputs without a rank get it,
  nothing recomputed), a boundary event, LCMV's filters with
  `S.lcmv_kappa` against the beamformer formula on data of rank 40 on 64
  channels, seven measures through stage 1 (the stored matrices equal
  `static_connectivity` of the stored region time series), and the five
  scout functions (mean, PCA, max, power, RMS) against their definitions
  and through stage 1.
- `test_static_statistics`: the features against their definitions, the
  graph metrics against the Brain Connectivity Toolbox (when installed),
  the reliability (computed only with both halves in at least half of the
  recordings), and the enrichment, NBS and its figures on synthetic
  matrices with a planted network; temporal channels named T3/T4, groups
  from the IDs, the sample-size report, the channel sets of a cap without
  10-20 names.
- `test_static_integration`: both steps and the PSD analysis on a BIDS
  dataset (three file formats, `participants.tsv`), run by
  `run_full_pipeline` (the order of the stages, numbered 1–3, its `'from'`
  and `'to'`, its errors, the PSD analysis switched off, its outputs); the
  warning for the window mean of phase measures, once per run; the
  provenance of the statistics. Its test of the brain network states, a
  stage this pipeline does not have, is skipped here by design.
- `test_static_ground_truth`: 16 + 16 recordings with stronger DMN
  coupling in one group, and a null dataset without it; what both steps
  recover (below).
- `test_psd`: the PSD analysis (`docs/psd.md`): Welch's scaling and cuts
  (one piece equals MATLAB's `pwelch`), specparam and the band measures on
  synthetic spectra, the cluster test on a planted cluster, on noise and
  against FieldTrip's own `clusterstat`, the neighbours of channels and
  regions, the region signals against stage 1's saved scouts, the
  settings' checks, and the whole step on 8 + 8 simulated recordings with
  a planted occipital alpha increase that a score follows. The channels'
  cross-spectra (their diagonal is `pwelch`'s, to 1e-12) give the power
  of every source point: equal to the spectrum of each source's own
  signal (to 1e-9) for eLORETA, wMNE and LCMV, normal and free
  orientation (the sum of the three); equal to FieldTrip's eLORETA source
  power from the same cross-spectra (`ft_inverse_eloreta`, to 1e-10;
  with free orientation FieldTrip keeps the dominant orientation, between
  a third of the sum and the sum) and to Brainstorm's PSD of the source
  map and of a scout (`process_psd` on the imaging kernel); a region's
  signal never has more power than the mean of its sources'.
- `test_anatomy`: the other head models (`docs/anatomy.md`) on two
  simulated recordings: a head folder holding the shipped ICBM152's files
  gives exactly the shipped results; individual heads, the fallback, a
  missing folder or file; atlases by name and their order; the settings;
  a leadfield made again when its head changes; `make_anatomy` with
  Brainstorm's ICBM152 template gives the shipped cortex, scalp, atlas and
  cap (skipped without Brainstorm's full tree).
- `test_reference_connectivity`: PLV, ciPLV, wPLI, PLI (whole recording
  and windows) against Brainstorm's `process_plv1n`, AEC against
  Brainstorm's `bst_henv` and FieldTrip's `amplcorr`, the leakage-corrected
  AEC against ROInets (equal to 7 × 10⁻¹⁴). dwPLI (computed by FieldTrip,
  epochs of 6 cycles of the band centre as trials: `ft_freqanalysis`
  mtmfft with a Hanning taper, `ft_connectivity_wpli` with `debias`, the
  mean over the band's Fourier frequencies) against FieldTrip's own entry
  points on the same epochs (`ft_freqanalysis`,
  `ft_connectivityanalysis` `wpli_debiased`): equal over the whole
  recording (largest difference 0) and in each of 8 windows of 20 s
  (3 × 10⁻⁸, single precision); against FieldTrip's usual dwPLI of the
  broadband signals, r 1.000 over the edges, largest difference 0.027.
  `test_static_estimators` also checks it against Vinck et al.'s (2011)
  formula written without FieldTrip (10⁻¹⁰, whole recording, 2 s or
  10-cycle epochs, windows), and its `test_debiased_wpli` that it is about 0
  without coupling (200 s of simulated 8–12 Hz noise: 0.000, SD 0.005).
- `test_static_estimators`, `test_static_settings`, `test_static_stage1`,
  `test_static_integration`: the window mean (`S.static_matrices`) equals
  that of a run with windows; window sizes from cycles or a fixed length,
  the overlaps, the AEC windows; no window fits: NaN and the reason; the
  settings and the keys; end to end, saved with its split halves, the
  windows not saved, the whole recording unchanged, both matrices analysed
  by the statistics. Windows of one measure (`S.measure_windows`): equal
  to the sliding-window functions at those sizes, the other measures
  unchanged, through stage 1 too; the settings checked, the warnings, the
  keys.
- `test_reference_statistics`: the group tests against `ttest2`, `fitlm`,
  PERMUTOOLS and PALM (statistic, degrees of freedom and permutation p;
  the z of each permutation on its own degrees of freedom is PALM's), the
  correlations against `partialcorr` and PERMUTOOLS, the studentised slope
  of the correlations against statsmodels' HC3 t and the wild bootstrap
  against wildboottest ("Correlations: the studentised slope", below),
  Hedges' g against PERMUTOOLS, the adjusted g against `fitlm` and
  `hedges_g`, `t_to_z` against PALM's `palm_gtoz` (its tests of the brain
  network states are skipped here by design).

The references each test needs (a test whose reference is missing is
skipped, reported as incomplete, not failed):

| Test | Reference |
|---|---|
| `test_reference_connectivity` | Brainstorm's whole tree, able to start headless (with the `brainstorm.jar` that `tools/fetch_toolboxes.sh` fetches); FieldTrip; ROInets (the copy in var-EEG-FC, in `Toolboxes/GitHub resources`) |
| `test_reference_statistics` | PERMUTOOLS and PALM (commit 7c00e18), in `Toolboxes/GitHub resources` or on the path; statsmodels and wildboottest, in the Python virtual environment `Toolboxes/python-references` or in `python3` (stored statsmodels 0.14.5 values are compared without Python) |
| `test_static_statistics` | the Brain Connectivity Toolbox (`Toolboxes/BCT`), for the comparison of the graph metrics only |
| `test_anatomy` | Brainstorm's whole tree with its ICBM152 template, for `make_anatomy` |
| `test_psd` | FieldTrip and Brainstorm as fetched by `tools/fetch_toolboxes.sh` (`clusterstat`, `ft_inverse_eloreta`, Brainstorm's PSD functions) |

The suite runs on R2024b and the latest MATLAB release on Linux, Windows
and macOS (`.github/workflows/tests.yml`: the toolboxes fetched with
`tools/fetch_toolboxes.sh --static`, no Parallel Computing Toolbox, and on
Windows the Visual C++ 2010 runtime that FieldTrip's MEX files need).
Three tests are skipped by design, since they need the brain network
states: the windowed stage 2 of `test_static_integration`, and the
duration and adjusted Rand index checks of `test_reference_statistics`.
Run one MATLAB release at a time on one machine: `test_static_integration`
writes to `tempdir/static_test_results`, which the run that finishes first
deletes.

## Validation

- **Stage 1.** The measures are checked against their formulas and
  against reference software ("Tests" above; Brainstorm below). The
  whole-recording matrices equal those of a run with sliding windows on
  the same recording (`test_static_estimators`, 1e-10; dwPLI included).
- **Same as Brainstorm.** `tools/compare_with_brainstorm.m` gives
  Brainstorm's connectivity process (`process_plv1n`, `bst_connectivity`,
  Hilbert transform) and `static_connectivity` the same analytic signal:
  on three recordings of the OpenNeuro example (120 s, 24 regions, alpha
  8–13 Hz) PLV, ciPLV, wPLI and PLI are equal (largest difference 1e-14).
  dwPLI has no usable Brainstorm counterpart on a continuous recording:
  Brainstorm's Hilbert dwPLI (not offered in its interface) gives 1.2 to
  2.0 on these recordings, above the largest possible value, 1, because
  its sums are averaged before the debiasing. The tool also computes this
  pipeline's dwPLI on the epochs of Brainstorm's Fourier dwPLI (1 s epochs
  as trials, Vinck et al. 2011); the dwPLI itself is validated against
  FieldTrip (`test_reference_connectivity`).
- **Simulations** (`test_static_ground_truth`,
  `simulation/run_static_simulation.m`; z statistics on the Welch df,
  max-statistic family test). 16 + 16 recordings of 60 s, DMN coupling
  0.30 against 0.45: the measured within-DMN wPLI follows each
  recording's planted coupling (r = 0.77); the group difference is found
  (g = 1.78, family-wise p = 0.0025, `p_bb` = 0.004), but not confined to
  the DMN: the salience–DMN block differs more (g = 2.29), 8 other blocks
  are significant, and the NBS component (37 edges, in the planted
  direction, none in the other) has 31% of its regions in the DMN (source
  leakage). The planted score correlates r = 0.67 with within-DMN wPLI
  (0.60 planted) but is not significant after the correction over the
  family (family-wise p = 0.14). The null dataset (equal coupling, no
  correlation): no significant test in the run (0 of 214), no NBS
  component.
- **Usability.** Twelve simulated recordings made as unlike each other as
  possible, in one folder without a BIDS layout: 19 channels with the old
  10-20 names (T3–T6), 32 with labels such as `EEG Fp1-REF`, 64 and 128
  (10-05 positions), upper-case labels; 128, 200, 250, 256, 500, 512,
  1000 and 1024 Hz; EEGLAB, EDF, BDF, BrainVision and FIF files; ECG and
  VEOG channels, a flat channel, a boundary event, a 5 s stretch of NaN,
  an EEGLAB file of 2 s epochs; a semicolon-separated subject table with
  decimal commas and IDs 001–012. Both steps with the six bands, eLORETA
  and LCMV, wPLI and orthogonalised AEC, every family and NBS (stage 1 13
  minutes, statistics 3 minutes, on R2024b): 11 recordings analysed, the
  epoched one stopped with the message explaining the 8-s minimum. With
  the individual bands and wMNE, the recordings in which restingIAF finds
  no alpha peak (10 of 12: the simulation has no posterior alpha rhythm)
  keep the fixed bands; with the Destrieux atlas (148 regions) and no
  network table, the network features are skipped.
- **Band power across regions.** The band power of a region is not
  comparable between regions: with eLORETA and the mean of the region's
  sources, small medial temporal regions (parahippocampal, entorhinal,
  temporal pole) have the largest power in real recordings and large
  frontal ones the smallest. Compare a region between recordings or
  groups.
- **Data rank** (`stage1_summary.csv`, `data_rank`): the number of
  independent signals in the channels used, from the eigenvalues of their
  covariance, by MNE-Python's rule for data stored as 32-bit floats
  (`connectivity.signal_rank`: eigenvalues above 5 × 10⁻¹¹ of the largest
  for 60 channels, 2.4 × 10⁻¹⁰ for 129), so that the rounding of 32-bit
  floats does not count as one more signal. In ICA-cleaned recordings the
  gap between kept and removed eigenvalues is wide, so the limit does not
  decide the rank.
- **Other head models** (`docs/anatomy.md`). Three heads made by
  `tools/make_anatomy.m` with Brainstorm 3.250107 (with its fix 952475b)
  and compared with the shipped head through their leadfields (the 343
  electrodes of the ASA 10-05 cap, the 68 Desikan-Killiany regions; for
  each region the topography of its sources' summed gain, whose overall
  sign is arbitrary: |r| between heads, and the ratio of its norms).
  - **ICBM152 from Brainstorm** (its current template; 9 minutes with the
    leadfield): the cortex, scalp, atlas and cap of the shipped head
    exactly (`test_anatomy`); its skulls have 1922 vertices instead of 642
    and are another segmentation (inner skull 2007 against 2198 cm³, skull
    4.9 against 4.1 mm, scalp 9.8 against 8.9 mm), hence a gain about 9%
    lower (sources: r median 0.999; amplitude ratio 0.91, 5–95% 0.83–1.00).
  - **Colin27_2016** (a Brainstorm template, downloaded; 1.4 minutes): the
    cap moved by 1.3 mm (median; max 6.0) onto its scalp.
  - **CAT12 on Colin27's T1** (`'t1'`: an individual head as made for a
    subject, from the template's own MRI, so that it can be set against
    the Colin27_2016 head; CAT12 12.9; 29 minutes, the caps and leadfield
    5 more): cortex = CAT12's central surface (15,002 vertices), atlases
    from CAT12 (Desikan-Killiany, Destrieux, Schaefer, HCP_MMP1, ...),
    non-linear MNI normalization, BEM of 1922 vertices per surface. Its
    scalp has two bumps on top of the head (57 vertices standing out by
    more than 8 mm, up to 39 mm; none on the other heads: at most 7 mm),
    which `make_anatomy` reports: the four electrodes on them moved by
    12 to 34 mm (FCC3h, FCC4h, FCC2, FCC1); the cap otherwise by 2.0 mm
    (median). The cap goes through the linear part of the MNI
    normalization: through CAT12's non-linear deformation (fitted to the
    brain) the electrodes would move by 3.7 mm (median), nine by more than
    10 mm.

  | Heads | Regions: \|r\| median (min) | Regions with \|r\| < 0.9 | Norm ratio median (range) |
  |---|---|---|---|
  | ICBM152 (Brainstorm) vs shipped | 0.999 (0.926) | 0 of 68 | 0.91 |
  | Colin27_2016 vs shipped | 0.924 | 33 | 0.95 |
  | CAT12 on Colin27 vs shipped | 0.910 | 32 | 0.83 |
  | CAT12 on Colin27 vs Colin27_2016 | 0.969 (0.411) | 22 | 0.87 (0.68–1.14) |
  | CAT12 on Colin27 vs ICBM152 (Brainstorm) | 0.909 | 32 | 0.92 |

  Another template changes the topographies of about half of the regions
  (|r| < 0.9) more than another segmentation of the same template does;
  CAT12 on the same MRI comes closer to that template than to the
  ICBM152, as it should. Edge-level results therefore depend on the head
  model, and a study should report which one it used.
- **Not tested here**: computing a leadfield for a new cap on Windows
  (the tests themselves pass on Linux, Windows and macOS on GitHub:
  "Tests" above).

## Statistics

These checks were made with the pipeline's own functions (`perm_ancova_maxT`,
`perm_partial_corr_maxT`, `cluster_perm_test`, `family_adjust`), in
MATLAB R2026a. The simulations behind the tables of this section are in
`tools/validation/`, one MATLAB function per table, which reproduces it
with the documented seeds and sizes (its README lists them; the first
argument of each is the number of simulations).

### Unequal groups

With one group much larger than the other and the smaller group more
variable, the groups are not exchangeable under the null, and permutation
p-values of the raw Welch t or Welch-type G are somewhat too small in the
far tail. The statistic of the data and of every permutation is therefore
turned into z on its own degrees of freedom before counting, as PALM does
(Winkler et al. 2014), which removes the excess. 8 against 30
recordings, the smaller group with 3 times the SD, no difference; 8,000
simulations, 1,999 permutations each; share of p below each level:

| Test | p < 0.05 | p < 0.01 | p < 0.005 | p < 0.001 |
|---|---|---|---|---|
| Welch t, group labels permuted, z on its df | 0.041 | 0.0070 | 0.0045 | 0.0006 |
| ANCOVA (Welch-type G, Freedman–Lane, age), z on its df | 0.044 | 0.0070 | 0.0045 | 0.0008 |
| Welch t, raw | 0.059 | 0.018 | 0.0099 | 0.0026 |
| ANCOVA, raw G | 0.062 | 0.020 | 0.011 | 0.0026 |

The family-wise error of a few columns is 0.07 at 0.05 (10 against 25
recordings), and Freedman–Lane permutations with covariates are slightly
liberal at strong imbalance (about 0.055 at 0.05). The check of sample
sizes notes every comparison whose larger group is at least 3 times the
smaller, with the group SDs in the tables.

### Correlations: the studentised slope

A permutation test of r (or partial r) tests that feature and score are
independent, not that they are uncorrelated (DiCiccio & Romano 2017):
when the spread of the feature changes with the score but r is 0, it
rejects too often. `perm_partial_corr_maxT` therefore permutes, in place
of r, the slope of the feature on the score (covariates and intercept in
the model) over its HC3 standard error (MacKinnon & White 1985; Long &
Ervin 2000), each permutation studentised in the same way, then z on
N − 2 − q degrees of freedom. r is still reported, and `t_hc3` beside it.
`S.correlation_test = 'wild_bootstrap'` takes the null distribution of
the same statistic from a restricted wild bootstrap (Rademacher signs,
null imposed; Davidson & Flachaire 2008; MacKinnon 2013), which keeps the
link between the score and the spread that a permutation breaks.

Rejection at p < 0.05 under the null, age and sex as covariates, 499
permutations or draws, 2,000 datasets per row:

| n | SD of the feature | studentised slope, permutations (default) | studentised slope, wild bootstrap | r, permutations (for comparison) |
|---|---|---|---|---|
| 20 | 0.2 + 2\|x\| | 0.102 | 0.057 | 0.207 |
| 40 | 0.2 + 2\|x\| | 0.084 | 0.054 | 0.238 |
| 60 | 0.2 + 2\|x\| | 0.078 | 0.055 | 0.215 |
| 120 | 0.2 + 2\|x\| | 0.051 | 0.042 | 0.235 |
| 20 | 1 + 0.5\|x\| | 0.073 | 0.057 | 0.100 |
| 40 | 1 + 0.5\|x\| | 0.065 | 0.059 | 0.120 |
| 60 | 1 + 0.5\|x\| | 0.064 | 0.052 | 0.121 |
| 120 | 1 + 0.5\|x\| | 0.048 | 0.046 | 0.107 |
| 20–120 | constant | 0.045–0.059 | 0.047–0.060 | 0.045–0.055 |

Below about 100 recordings a strongly heteroscedastic feature still
gives 6–10% at 0.05 with permutations; the wild bootstrap stays within
0.042–0.060. With a constant spread both lose a little power (a slope of
0.4 SD found in 0.30, 0.60, 0.80 and 0.98 of the datasets at n = 20, 40,
60 and 120 with permutations, 0.28, 0.58, 0.79 and 0.98 with the wild
bootstrap, 0.34, 0.64, 0.83 and 0.99 for r). With fewer than about 100
recordings, or a feature whose spread visibly changes with the score,
use the wild bootstrap; plot the feature against the score when a
correlation matters.

Checks (`tests/test_reference_statistics.m`): `t_hc3` equals statsmodels'
`OLS(...).fit(cov_type='HC3')` t (1e-10; statsmodels 0.14.5, with and
without covariates, Pearson and Spearman); the permutation statistics
equal those redone one at a time with the textbook sandwich of the full
model; the wild bootstrap's p equals that of wildboottest 0.3.2's
`WildboottestHC` (Roodman et al. 2019; Rademacher, null imposed) within
Monte Carlo error; null datasets within 0.03 of 5%.

### Effect size with covariates

With covariates the group comparison tests the adjusted difference b1.
`g` is its effect size: b1 over the pooled within-group SD of the raw
values, times Hedges' exact J on nA + nB − 2 degrees of freedom (What
Works Clearinghouse handbook, version 5.0, 2022, appendix E), with a
percentile bootstrap CI (2,000 resamples within each group, the model
refitted on each); `g_unadjusted` is the raw g with its CI. Without
covariates both are Hedges' g. Checks: `g` equals `hedges_g` times
`fitlm`'s b1 over the raw difference of the means (1e-10). Simulation
(age as covariate, group B 5 years older on average, age raising the
feature by 0.05 SD per year, 1,000 datasets):

| Groups | true δ (target) | adjusted g, mean (SD) | 95% CI covers the target | raw g, mean |
|---|---|---|---|---|
| 20 + 20 | 0.5 (0.459) | 0.464 (0.316) | 0.932 | 0.684 |
| 10 + 30 | 0.8 (0.734) | 0.742 (0.374) | 0.923 | 0.963 |
| 30 + 30 | 0 (0) | 0.006 (0.257) | 0.941 | 0.232 |

The raw g mixes the group difference with the age difference between the
groups; the adjusted g does not. The percentile CI covers slightly less
than 95% in groups of 10–30.

### Multiple testing: the max-statistic family test

A family is one feature family and one comparison or correlation set
(for one band, inverse and measure). Every test of a family runs on the
same permutations, and each enters as a |z| (a Welch t or G turned into z
on its own degrees of freedom, the t of a partial r likewise), so that
tests of different N weigh the same. The family's largest |z| is compared
with the largest of each permutation (Westfall & Young 1993): `family_p`,
the family as a whole, and `p_maxT_family`, each test against the
family's maximum. The families are selected by Benjamini–Hochberg over
`family_p`, and a test of a selected family is significant when `p_bb` =
`p_maxT_family` × F / R ≤ alpha (Benjamini & Bogomolov 2014). A selected
family has at least one significant test, and both p-values reach
1 / (n_perm + 1) whatever the size of the family. `q_fdr` stays in the
tables as a description. The decision inside a family is family-wise: an
effect spread thinly over many tests, none strong, is found less often
than with FDR inside the family; the PSD's cluster tests cover spread
effects.

Checks: against PERMUTOOLS (Crosse et al. 2024, `permuttest2`, max
correction; 5 datasets, 19,999 permutations: family p equal within Monte
Carlo error, the tests' p-values correlated r ≥ 0.9999). Under the null
(1,000 simulated studies each): group comparison of 20 against 30
subjects, the smaller group twice as variable, age as covariate, 10%
missing values in one feature: family p ≤ 0.05 in 4.4% (≤ 0.01 in 1.0%);
partial correlations, N 40, 10% missing values: 4.5% (≤ 0.01 in 0.9%)
with the studentised slope (4.2% with permutations of r). The whole
adjustment over 10 null families of 5 tests: a significant test in 6.0%
of 400 runs (`tools/validation/family_maxT_checks.m`). Benjamini–Bogomolov
controls the error averaged over the selected families, under
independence of the families; the families of a run share recordings and
measures, so report F and R (`families.csv`).

### Permutation resolution

A permutation p is at least 1 / (n_perm + 1), and so are `family_p` and
`p_maxT_family`, whatever the size of the family. With F families of
which one has signal, Benjamini–Hochberg selects it when its `family_p` ≤
alpha / F, and its strongest test is then significant: a family is
selected alone when at most alpha (n_perm + 1) / F − 1 permutations reach
its largest statistic. So

    n_perm ≥ F / alpha − 1

at the least, and about 10 F / alpha leaves room (9 permutations). The
statistics state the rule at their start and the number for their
families at their end, the PSD statistics at their start
(`resolution_note.m`); the PSD has its own `S.psd.n_perm`. The template's
six bands and two measures, with its four feature families and the
spectral profile, give 43 families per comparison or correlation set, so
258 with six sets: 10,000 permutations allow none (the family's largest
statistic must beat every permutation), 50,000 allow 8. A study with one
comparison and a few feature families (10–20 families) is well served by
10,000 (24–49).
