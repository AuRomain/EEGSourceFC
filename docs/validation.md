# Validation and tests of v_static

## Tests

```matlab
results = runtests('tests');  table(results)
```

No private data are needed: the recordings are simulated on the
pipeline's forward model.

- `test_static_estimators`: every connectivity measure against its direct
  formula, on one piece and on several; the known couplings (lagged,
  zero-lag, independent); the split halves; the measures equal the
  whole-recording outputs of v_dynamic's windowed functions.
- `test_static_settings`: the static settings and their checks (the
  warning for settings of the windows; the scout function), the outputs
  folder, what changes the cache keys. `test_template`: this folder's template and the
  examples.
- `test_static_stage1`: stage 1 on simulated recordings in a BIDS layout:
  the outputs, the cache (skipped when nothing changed, redone when a
  setting changes), the same recording in EEGLAB, BrainVision, EDF and BDF
  formats, the individual bands of a recording without alpha, flat
  channels and channel subsets, the data rank (recordings with 20 and 50
  of their 64 dimensions removed; added to the reused outputs of earlier
  code), a boundary event, seven measures
  through stage 1 (the stored matrices equal `static_connectivity` of the
  stored region time series), and the five scout functions (mean, PCA,
  max, power, RMS) against their definitions and through stage 1.
- `test_static_statistics`: the features against their definitions, the
  graph metrics against the Brain Connectivity Toolbox (when installed),
  the reliability, and the enrichment, NBS and its figures on synthetic
  matrices with a planted network; and what the two new datasets showed
  (temporal channels named T3/T4, groups from the IDs, the sample-size
  report, the channel sets of a cap without 10-20 names).
- `test_static_integration`: both steps and the PSD analysis on a BIDS
  dataset (three file formats, `participants.tsv`), run by
  `run_full_pipeline` (the order of the stages, its `'from'` and `'to'`,
  its errors, the PSD analysis switched off, its outputs). Its test of v_dynamic's
  stage 2 is skipped here by design.
- `test_static_ground_truth`: 16 + 16 recordings with stronger DMN
  coupling in one group, and a null dataset without it; what both steps
  recover (below).
- `test_psd`: the PSD analysis (`docs/psd.md`): Welch's scaling and cuts
  (one piece equals MATLAB's `pwelch`), specparam and the band measures on
  synthetic spectra, the cluster test on a planted cluster, on noise and
  against FieldTrip's own `clusterstat`, the neighbours of channels and
  regions, the region signals against stage 1's saved scouts, the
  settings' checks, and the whole step on 8 + 8 simulated recordings with
  a planted occipital alpha increase that a score follows.
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
  AEC against ROInets (equal to 7 × 10⁻¹⁴). dwPLI against FieldTrip's
  `ft_connectivity_wpli` with `debias` (Vinck et al. 2011): given the same
  blocks of half the decorrelation time as trials, equal (5 × 10⁻¹³) when
  the pairs of close blocks, which this pipeline also leaves out, add
  nothing; with them, FieldTrip's estimator keeps a bias in 2 s windows
  (0.134 on pairs without lagged coupling) that leaving them out removes
  (0.016); against FieldTrip's usual dwPLI of 2 s epochs
  (`ft_freqanalysis`, `ft_connectivityanalysis` `wpli_debiased`), r 0.99
  over the edges (6 November 2026).
- `test_static_estimators`, `test_static_settings`, `test_static_stage1`,
  `test_static_integration`: the window mean (`S.static_matrices`, 6
  November 2026) equals that of a run with windows; window sizes from
  cycles or a fixed length, the overlaps, the AEC windows; no window fits:
  NaN and the reason; the settings and the keys; end to end, saved with its
  split halves, the windows not saved, the whole recording unchanged, both
  matrices analysed by the statistics. Windows of one measure
  (`S.measure_windows`): equal to the sliding-window functions at those
  sizes, the other measures unchanged, through stage 1 too; the settings
  checked, the warnings, the keys.
- `test_reference_statistics`: the group tests against `ttest2`, `fitlm`,
  PERMUTOOLS and PALM (statistic, degrees of freedom and permutation p;
  the z of each permutation on its own degrees of freedom is PALM's), the
  correlations against `partialcorr` and PERMUTOOLS, Hedges' g against
  PERMUTOOLS, `t_to_z` against PALM's `palm_gtoz` (its tests of v_dynamic's
  stage 2 are skipped here by design). The reference toolboxes are looked
  for in `Toolboxes/GitHub resources`; a test whose reference is missing
  is skipped, not failed.
- `test_same_code_as_v_dynamic`: every shared file is an unchanged copy of v_dynamic's,
  no other code file is here, and no function called here exists only in
  v_dynamic.

On 6 November 2026, late (the final code of version 1.0.0: the window
mean as the default, the windows of one measure, the demo over the whole
recording): 97 passed and 3 skipped by design of 101 (the windowed stage
2, and the duration and adjusted Rand index checks of
`test_reference_statistics`, which need the states) on R2026a (11.5
minutes, macOS); the remaining test, `test_same_code_as_v_dynamic`,
compares the shared files with v_dynamic's main branch, which receives
the same files (branch `static-windows`, 162 of 162 equal) once its
running analyses end. `run_static_demo` on this code: the numbers of
README, "Try it" (4.6 minutes). Earlier the same day, with the window
mean and before the demo was fixed: 96 passed and 3 skipped of 100
(16.5 minutes). The
code before the window mean became the default passed 86 and skipped 3
of 89 on R2026a (25.7 minutes) and on R2024b (54.8 minutes). (On 5 November, the commit of
the PSD analysis: 67 passed and 1 skipped of 68, 7.6 and 20.5 minutes.) Run
one release at a time on one machine: `test_static_integration` writes
to `tempdir/static_test_results`, which the run that finishes first
deletes.
`.github/workflows/tests.yml` runs them on GitHub
(Linux, Windows, macOS; R2024b and the latest release), fetching the
toolboxes with `tools/fetch_toolboxes.sh --static`.

## Validation

- **The code of v_dynamic.** Stage 1 and the statistics are v_dynamic's files
  (`tests/test_same_code_as_v_dynamic.m`), whose stage 1 was validated against
  independent implementations and on the 138 OpenNeuro recordings. The
  whole-recording matrices equal those of a v_dynamic run with windows on the
  same recording (`test_static_estimators`, 1e-10; dwPLI excepted, whose
  windowed version fits its blocks into the windows).
- **Same as Brainstorm.** `tools/compare_with_brainstorm.m` gives
  Brainstorm's connectivity process (`process_plv1n`, `bst_connectivity`,
  Hilbert transform) and `static_connectivity` the same analytic signal:
  on three OpenNeuro recordings (120 s, 24 regions, alpha 8–13 Hz) PLV,
  ciPLV, wPLI and PLI were equal (largest difference 1e-14). dwPLI has no
  usable Brainstorm counterpart on a continuous recording: Brainstorm's
  Hilbert dwPLI (not offered in its interface) gave 1.2 to 2.0, above the
  largest possible value, 1, because its sums are averaged before the
  debiasing. Against Brainstorm's Fourier dwPLI with 1 s epochs taken as
  trials (Vinck et al. 2011), this pipeline's dwPLI, whose observations
  are blocks of half the band's decorrelation time, agreed in level (mean
  0.021–0.031 against 0.021–0.028) and pattern (r 0.76–0.90 over the
  edges), as two estimators of the same quantity on 120 s do.
- **Simulations** (`test_static_ground_truth`, `simulation/run_static_simulation.m`;
  its run of 5 November 2026 on the current code: z statistics on the
  Welch df, max-statistic family test).
  16 + 16 recordings of 60 s, DMN coupling 0.30 against 0.45: the
  measured within-DMN wPLI followed each recording's planted coupling
  (r = 0.77); the group difference was found (g = 1.78, family-wise
  p = 0.0025, `p_bb` = 0.004), but not confined to the DMN: the
  salience–DMN block differed more (g = 2.29), 8 other blocks were
  significant, and the NBS component (37 edges, in the planted direction,
  none in the other) had 31% of its regions in the DMN (source leakage).
  The planted score correlated r = 0.67 with within-DMN wPLI (0.60
  planted) but was not significant after the correction over the family
  (family-wise p = 0.14). The null dataset (equal coupling, no
  correlation): no significant test in the run (0 of 214), no NBS
  component. (With the raw statistics and the family rule before 5
  November: 17 other blocks, an NBS component of 64 edges, the score's
  q 0.07.)
- **Usability** (3 November 2026). Twelve simulated recordings made as
  unlike each other as possible, in one folder without a BIDS layout:
  19 channels with the old 10-20 names (T3–T6), 32 with labels such as
  `EEG Fp1-REF`, 64 and 128 (10-05 positions), upper-case labels;
  128, 200, 250, 256, 500, 512, 1000 and 1024 Hz; EEGLAB, EDF, BDF,
  BrainVision and FIF files; ECG and VEOG channels, a flat channel, a
  boundary event, a 5 s stretch of NaN, an EEGLAB file of 2 s epochs; a
  semicolon-separated subject table with decimal commas and IDs 001–012.
  Both steps with the six bands, eLORETA and LCMV, wPLI and orthogonalised
  AEC, every family and NBS (stage 1 13 minutes, statistics 3 minutes, on
  R2024b): 11 recordings analysed, the epoched one stopped with the
  message explaining the 8-s minimum. With the individual bands and wMNE,
  the recordings in which restingIAF found no alpha peak (10 of 12: the
  simulation has no posterior alpha rhythm) kept the fixed bands; with
  the Destrieux atlas (148 regions) and no network table, the network
  features were skipped. The problems found this way (BDF files on Apple
  silicon, the individual bands of a recording without alpha, the message
  for short epochs) and by an independent review were fixed (CHANGES.md).
- **Not tested here**: Windows and Linux (the GitHub workflow runs the
  tests there), computing a leadfield for a new cap on Windows.

## Runs on real data

From an empty results folder, as a new user would run them (4 November
2026, R2026a, 8 workers, the code of v_dynamic edde89a):

- **OpenNeuro, with the settings of v_dynamic's run** (138 recordings, two caps,
  first-episode psychosis and controls; `examples/settings_openneuro.m`
  with v_dynamic's bands delta 1–4, theta 4–8, alpha 7–13, beta 12–30 Hz, PLV,
  wPLI and AEC, the `'window'` filter rule and v_dynamic's seed): stage 1 29
  minutes, statistics 15 minutes. Against v_dynamic's outputs of 22 and 29
  October: every whole-recording matrix identical (1656 of 1656: 138
  recordings × 4 bands × 3 measures; largest difference 0); the graph
  metrics of wPLI, their group comparisons and their correlations with the
  40 scores identical (t, g and r within 2e-15; p, max-T p and q equal);
  the NBS of alpha wPLI the same components with the same p. Only the
  adjustment over families differs, as it must: v_dynamic's runs adjusted over
  54 or 36 families (one band, with the state families), this one over
  30 (four bands); none was selected in either. After v_dynamic's next commit
  was copied here (the text of a message in a stage 1 file), `run_all` on
  these outputs reused all 138 recordings (48 s, no file rewritten), as
  `1- EEG source connectivity/equivalent_code.txt` intends.
- **OpenNeuro, the example as distributed** (six bands, wPLI and
  orthogonalised AEC, Brainstorm's filter width, every family, both
  comparisons, the correlation sets): stage 1 1.6
  hours (while the other run used the machine too), statistics 36
  minutes. All 138 recordings analysed (277–298 s each); 258 families
  (52,542 tests), none selected after the adjustment over families; one
  NBS component (beta2 wPLI, psychosis > controls, 16 edges, `p_both`
  0.03, not corrected over the 24 NBS analyses of bands × measures ×
  comparisons). Brainstorm's filter width against v_dynamic's `'window'` rule,
  same recordings and band edges: alpha wPLI r = 0.999 over the edges
  (median of the recordings), delta wPLI r = 0.93 (lowest 0.63), the band
  in which the two filters differ most.
- **OpenNeuro, the example on the code of 5–6 November** (R2026a, 8
  workers, beside three other analyses; stage 1 and the rest on two
  versions of 5 November with the same stage 1 code; the PSD analysis, the
  z statistics, the max-statistic family test, `S.n_perm` 50,000), with
  the connectivity over the whole recording (`S.static_matrices =
  {'full'}`, then the default; the example now computes the window mean,
  the default since 6 November 2026, which has not been run on these data
  with this folder; the PSD analysis does not depend on it). Stage 1 4.5 hours, the PSD analysis 3.9 hours (16 minutes for
  the spectra, the rest for the statistics of its 324 families), the
  statistics 1.9 hours (258 families). All 138 recordings analysed.
  - Statistics: 2 families selected, 2 significant tests, one result
    found with and without the chlorpromazine adjustment: in the 78
    patients, the strength of the left isthmus cingulate in theta
    (orthogonalised AEC) correlated with DISORGANIZATION (r = 0.57,
    `p_bb` 0.018 and 0.021). No group comparison significant. No NBS
    component: the smallest `p_both` was 0.071 (beta2 wPLI, psychosis >
    controls, 7 edges; on 4 November, with the NBS on the raw G, 0.03 and
    16 edges).
  - PSD analysis: 8 families selected, 10 significant tests, all
    correlations in the patients, each found with and without the
    chlorpromazine adjustment: absolute delta power at C2 and the
    aperiodic offset at C2 and CP2 with SANS_ANA (r −0.52 to −0.56,
    `p_bb` 0.003–0.021); relative and periodic alpha power of the right
    inferior temporal region with NEG_SYMPT (r −0.52, `p_bb` 0.03–0.04).
    No group comparison significant; 2,814 clusters, none significant.
  - These are exploratory correlations among many tests, each significant
    after the adjustment over families. They show that the final code runs
    the distributed example end to end; they do not validate the pipeline.
- **A new study: the lab's eight development recordings** (four patients
  and four controls, eyes-open rest; not shared). Set up
  from the template with only what this study differs in: EEGLAB files at
  512 Hz, the first 180 s, the lab's 61-channel cap, a semicolon-separated
  subject table with its own column names, patients − controls adjusted
  for age, correlations with eight cognitive (MCCB) domains in the
  patients. Both steps ran without an error, and without a warning from
  the pipeline: with the six fixed bands (stage 1 10.5 minutes,
  statistics 7.1 minutes; 86 families) and with the individual bands
  (12.2 and 6.8 minutes; 88 families, the IAF family included). The
  individual alpha frequency was found in all eight (9.8–11.2 Hz); two
  recordings were flagged because restingIAF's peak and specparam's
  differed by more than 0.5 Hz, as intended. With 4 + 4 recordings no
  family was selected (smallest family p 0.06): a check that the pipeline
  adapts to another study, not a result. Again on 5–6 November 2026 with
  the current code (the PSD analysis, the z statistics, the max-statistic
  family test): from an empty folder, then the PSD analysis and the
  statistics once more on the final code; no error. With age as covariate
  the MCCB correlations in the 4 patients are not testable, which the log
  says; 43 families with tests in the statistics and 54 in the PSD
  analysis, none selected; 45 PSD clusters, none significant.
- **Regularisation and scout function** (4 November 2026). The 138
  OpenNeuro recordings with the settings of v_dynamic's run (above), one change
  at a time against the reference (eLORETA λ 0.05, the mean): λ 0.01,
  0.1 and 0.5, and `S.scout_function = 'pca'`; both steps each time
  (stage 1 26–27 minutes, statistics 13–16 minutes, with other runs on
  the machine).

  | Against the reference | λ 0.01 | λ 0.1 | λ 0.5 | PCA |
  |---|---|---|---|---|
  | Each recording's pattern of edges, r (median of the recordings) | 0.78–0.90 | 0.94–0.98 | 0.57–0.78 | 0.39–0.63 |
  | Mean pattern over the recordings, r | 0.85–0.96 | 0.96–0.99 | 0.67–0.86 | 0.55–0.84 |
  | Each recording's mean connectivity, r across the recordings | 0.97–0.99 | 0.99–1.00 | 0.91–0.98 | 0.86–0.96 |
  | Level of PLV and AEC | ×0.91–0.98 | ×1.03–1.06 | ×1.17–1.28 | ×1.11–1.20 |
  | Level of wPLI | ×0.98–1.05 | ×0.98–1.01 | ×0.94–1.04 | ×1.00–1.03 |
  | Split-half r_sb of wPLI strength, alpha (reference 0.95) | 0.95 | 0.95 | 0.96 | 0.95 |
  | Same, delta (reference 0.80) | 0.81 | 0.79 | 0.75 | 0.70 |
  | Recordings identified from their other half, delta (reference 0.75) | 0.76 | 0.67 | 0.54 | 0.62 |
  | wPLI strength, psychosis − controls, alpha g (reference −0.13) | −0.13 | −0.13 | −0.13 | −0.16 |
  | Same, beta g (reference 0.23) | 0.22 | 0.22 | 0.15 | 0.15 |
  | Alpha correlations with p < 0.05, of 120 (reference 20) | 23 | 18 | 18 | 24 |

  (Four bands; PLV, wPLI and AEC for the matrices; wPLI for the rest.)
  A change of the inverse moves each recording's connectivity between
  neighbouring regions (leakage) more than it changes the group
  comparisons: no group comparison of strength was significant in any variant, and no
  NBS component. λ 0.05 lies in a flat range (0.01 to 0.1 give the same
  results); 0.5 regularises too much (lower reliability in delta and
  theta, higher PLV and AEC from leakage). The mean stays the default
  scout function and the PCA a sensitivity analysis: its matrices differ
  recording by recording, its group results agree. The nonlinear
  functions (`'max'`, `'power'`, `'rms'`) were compared on two recordings
  in v_dynamic (its `docs/validation.md`): power and RMS have no phase, and their
  wPLI is unrelated to the mean's.
- **Two datasets never seen before** (4 November 2026). Both pipelines
  from an empty results folder, each dataset on its own, set up from the
  template with only what the dataset differs in; statistics with age as
  covariate. No error in any stage; the problems found were fixed and the
  steps they affected run again on the fixed code (CHANGES.md: "Checks of
  a new cap", "Statistics on small or unusual datasets").
  - **NIRP** (2 patients with schizophrenia, 2 controls; eyes closed; EEGLAB
    files at 512 Hz; 27 channels with the old names T3–T6, no channel file:
    positions from the ASA 10-05 cap by label; IDs such as `18PZ006B`).
    Stage 1: all four recordings, 297 s each. Alpha power of the regions
    higher at the back than at the front in every recording (log10 +0.50 to
    +0.69, posterior minus frontal regions); alpha wPLI of the two halves of
    each recording r 0.91–0.92 (AEC 0.96, PLV 0.98), every recording
    identified from its other half; v_dynamic and v_static (different filters)
    alpha wPLI r 0.999–1.000. IAF 9.5–11.1 Hz, one recording flagged
    (restingIAF and specparam 0.7 Hz apart). v_dynamic's states (its stage 2) were
    not stable with 4 recordings, as expected. Statistics: with 2 + 2
    recordings and age as covariate nothing can be tested, which the log now
    says before the families.
  - **Geneva** (3 patients, 3 controls; eyes closed; ICA-cleaned EEGLAB
    files at 500 Hz, one subfolder per recording; a 129-channel ANT
    equidistant cap whose labels are Z1 … RR14, without a single 10-20 name;
    an .xlsx subject table without a group column: the group is the first
    letter of the ID). The cap: the lab's `.ced`, already on Brainstorm's
    ICBM152 head, in metres and projected on the scalp (median 0.3 mm) as a
    CSV. The leadfield was built by stage 1 itself (OpenMEEG, 40 s for the
    gain): electrodes 0.8 mm from the scalp (median), the sign and
    localisation checks passed at the electrodes nearest to Cz, Fpz, Oz, T7
    and T8 (6–16 mm away); built three times (alone, by v_dynamic, by v_static),
    the same leadfield (difference 0). Stage 1: all six recordings, 247–248
    s each. The recordings themselves have a weak and inconsistent alpha
    topography: at the scalp, outside the pipeline and with the files' own
    electrode positions, the posterior peak was at 6 Hz in one and 14 Hz in
    another, and the alpha power at the back exceeded the front's in two of
    six (by 0.05 and 0.13 log10). The regions' alpha was higher at the back
    in four (+0.21 to +0.44) and not in two (−0.11, −0.04); alpha wPLI of
    the two halves r 0.75 (median; P120, with no alpha topography at the
    scalp and the lowest alpha wPLI, 0.04 against 0.09–0.15, r 0.006: there
    is no alpha coupling to reproduce), 5 of 6 identified (AEC 0.95, all
    identified); v_dynamic and v_static alpha wPLI r 0.98–0.999. With the sets
    matched by position, the spectral profile found the IAF in four
    (9.5–11.9 Hz; the other two have no alpha peak on 3 posterior channels,
    as their topography suggests) and the muscle and eye indices of all six.
    Statistics: patients − controls with age (3 + 3) and the correlations
    with age (6) tested, the PANSS factors in the 3 patients with age not
    testable (said); 86 families with tests, none selected.
  - **Again on 5–6 November 2026**, with the current code (the PSD
    analysis, the z statistics, the max-statistic family test): each
    dataset from an empty folder, then the PSD analysis and the statistics
    once more on the final code; no error in any stage. NIRP: nothing
    testable with 2 + 2 recordings and age as covariate, said before the
    families of both analyses. Geneva: 86 families with tests in the
    statistics and 108 in the PSD analysis, none selected; 59 PSD
    clusters, none significant.
  - **Fixed** (found on these datasets): the leadfield's checks needed 10-20
    labels and measured the distance to the nearest scalp vertex; the groups
    from the IDs could not be named in comparisons (every family stopped); a
    comparison or correlation with too few recordings gave NaN or no rows
    without a word; the temporal channels ignored T3 and T4; the sample
    sizes counted no recording when stage 1 had been reused; on a cap
    without 10-20 names every posterior, frontal, frontopolar and temporal
    set was empty (no IAF; the alpha, muscle and eye measures NaN). Harmless
    messages: FieldTrip's "rounding sample index" (fractional event
    latencies in an EEGLAB file) and its adding its own copy of a few
    Brainstorm functions for specparam (none of them on the pipeline's path
    otherwise).
  - **To know when reading the outputs**: the band power of a region is not
    comparable between regions: with eLORETA and the mean of the region's
    sources, small medial temporal regions (parahippocampal, entorhinal,
    temporal pole) have the largest power in every recording (NIRP, Geneva
    and OpenNeuro alike) and large frontal ones the smallest; compare a
    region between recordings or groups.
- **Data rank** (4 November 2026; `stage1_summary.csv`, `data_rank`): the
  number of independent signals in the channels used, from the
  eigenvalues of their covariance (above 1e-12 of the largest). In real
  ICA-cleaned recordings the gap is wide (kept eigenvalues at least about
  1e-9 of the largest, removed ones about 1e-14), so the limit does not
  decide the rank. The 138 OpenNeuro recordings (60 channels): rank 9–49,
  median 33 (quartiles 29 and 37), 2 below `S.min_channels` (9 and 13,
  the warning's limit; 57 were below 32). NIRP: 14, 21, 21 and 23 of 27
  (one warning). Geneva: 117–128 of 129. The change of stage 1 computes
  nothing differently (checked in v_dynamic: an OpenNeuro recording recomputed,
  identical matrices; the NIRP and Geneva outputs reused, the rank added,
  nothing recomputed).
- **Other head models** (4 November 2026, `docs/anatomy.md`). Three heads
  made by `tools/make_anatomy.m` with Brainstorm 3.250107 (with its fix
  952475b) and compared with the shipped head through their leadfields
  (the 343 electrodes of the ASA 10-05 cap, the 68 Desikan-Killiany
  regions; for each region the topography of its sources' summed gain,
  whose overall sign is arbitrary: |r| between heads, and the ratio of
  its norms).
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
    which `make_anatomy` now reports: the four electrodes on them moved by
    12 to 34 mm (FCC3h, FCC4h, FCC2, FCC1); the cap otherwise by 2.0 mm
    (median). The cap goes through the linear part of the MNI
    normalization: through CAT12's non-linear deformation (fitted to the
    brain) the electrodes moved by 3.7 mm (median), nine by more than 10
    mm.

  | Heads | Regions: \|r\| median (min) | Regions with \|r\| < 0.9 | Norm ratio median (range) |
  |---|---|---|---|
  | ICBM152 (Brainstorm) vs shipped | 0.999 (0.926) | 0 of 68 | 0.91 |
  | Colin27_2016 vs shipped | 0.924 | 33 | 0.95 |
  | CAT12 on Colin27 vs shipped | 0.910 | 32 | 0.83 |
  | CAT12 on Colin27 vs Colin27_2016 | 0.969 (0.411) | 22 | 0.87 (0.68–1.14) |
  | CAT12 on Colin27 vs ICBM152 (Brainstorm) | 0.909 | 32 | 0.92 |

  The 60 electrodes of the OpenNeuro recordings give the same (0.999,
  0.920, 0.914, 0.967, 0.908). Another template changes the topographies
  of about half of the regions (|r| < 0.9) more than another segmentation
  of the same template does; CAT12 on the same MRI comes closer to that
  template than to the ICBM152, as it should.

  On the 138 OpenNeuro recordings, both steps with each head
  (`S.anatomy`; otherwise the settings of the reference run of the
  regularisation test above; stage 1 18 minutes and statistics 12 with
  each head), against the shipped head:

  | Against the shipped head | ICBM152 (Brainstorm) | Colin27_2016 | CAT12 on Colin27 |
  |---|---|---|---|
  | Each recording's pattern of edges, r (median of the recordings) | 0.97–0.98 | 0.45–0.63 | 0.40–0.60 |
  | Each edge across the recordings, r (median of the edges) | 0.99 | 0.44–0.75 | 0.35–0.73 |
  | Mean pattern over the recordings, r | 0.98–0.99 | 0.51–0.76 | 0.38–0.73 |
  | Each recording's mean connectivity, r across the recordings | ≥ 0.998 | 0.966–0.995 | 0.941–0.991 |
  | Level | ×0.98–1.00 | ×0.97–1.03 | ×0.93–1.00 |
  | Split-half r_sb of wPLI strength, alpha / delta (reference 0.95 / 0.80) | 0.95 / 0.80 | 0.95 / 0.83 | 0.96 / 0.83 |
  | wPLI strength, psychosis − controls, alpha / beta g (reference −0.13 / 0.23) | −0.12 / 0.24 | −0.15 / 0.23 | −0.13 / 0.22 |
  | Largest change of g over the 24 graph-metric tests | 0.02 | 0.07 | 0.12 |
  | Graph-metric tests with p < 0.05 (reference 0) | 0 | 0 | 0 |
  | Alpha correlations with p < 0.05, of 120 (reference 20) | 20 | 18 | 20 |
  | Smallest NBS p (reference 0.10) | 0.39 | 0.25 | 0.10 |

  (Four bands; PLV, wPLI and AEC for the matrices; wPLI for the rest;
  the run's tables `*_anatomy.csv` and `group_mean_pattern.csv`, kept
  by the lab.) Brainstorm's ICBM152, the same anatomy
  with other skulls, changes nothing that matters. Another head changes
  which pairs of regions are coupled, in each recording and in the group
  mean, as much as the PCA scout function does; each recording's overall
  connectivity, the reliability and the group results do not (no graph
  metric, correlation set or NBS component changes its conclusion). A
  head made from the same MRI as a template (CAT12 on Colin27) is as far
  from the shipped head as that template is. Edge-level results are
  therefore specific to the head model, and a study should report which
  one it used.
