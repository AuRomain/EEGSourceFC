# EEGSourceFC — static EEG source connectivity and power spectra

MATLAB pipeline for resting-state EEG: it estimates the activity of
cortical regions from the scalp recordings, measures the connectivity
between every pair of regions in several frequency bands, as the mean
of short windows (the default) or over the whole recording, and compares groups (and correlates with clinical
scores) on the edges and on the networks' graph metrics.

| Stage | What it does | Main function |
|---|---|---|
| 1. Source connectivity | Each recording goes through eLORETA (or LCMV, wMNE) on the ICBM152 template head, giving 68 Desikan-Killiany regions (or another atlas). Each continuous piece is band-passed in six bands (delta 1–4, theta 4–7, alpha 7–13, beta1 13–20, beta2 20–30, gamma 30–45 Hz; or bands anchored on each recording's individual alpha frequency). wPLI and orthogonalised AEC of every pair of regions (also PLV, PLI, ciPLV, AEC, debiased wPLI) and the band power of every region are computed as the mean of windows (the default), over the whole recording, or both (`S.static_matrices`). | `run_all` |
| 1b. Power spectra | Welch spectra (4 s) of every channel and region give absolute, relative and periodic (above the 1/f fit) power in eight bands, the 1/f exponent and offset, and the alpha frequency. These are compared between groups and correlated per channel and region, with cluster tests, and drawn as scalp and cortex maps. | `run_psd_analysis` |
| 3. Statistics | The network-based statistic (NBS) on the edges of every comparison, and where the edge differences concentrate (networks against random and spin-rotated labels). Group comparisons (permutation ANCOVA) and correlations of the graph metrics (strength, clustering coefficient, characteristic path length), network blocks, regions and band power. Each family is tested by its largest statistic; the families are adjusted over the whole run, bands included. | `run_statistics` |

`run_full_pipeline` runs the three stages in order. Stage 2 (brain network
states, connectivity over time) belongs to the companion pipeline
**EEGSourceStates** (folder `v_dynamic`), which shares this code
(section 9); it is skipped here. In the lab's project and in these
documents, this pipeline's folder is `v_static`.

**Contents:** 1 What you need · 2 Installation · 3 Try it · 4 Your study,
step by step · 5 Settings · 6 Outputs · 7 Limits · 8 Tests and validation ·
9 Reproducibility and shared code · More detail · How to cite · Licence

## 1. What you need

**Data**

- **EEG recordings, cleaned** (artefacts removed, bad channels removed or
  interpolated; the pipeline does not clean data), in any format
  FieldTrip reads: EEGLAB `.set`, EDF, BDF, BrainVision, FIF, ...
  Continuous recordings, with rejected stretches marked by boundary
  events, are best. Epoched files work only with epochs of at least 8 s
  (each epoch is filtered separately).
- **Channel positions: usually nothing to provide.** Standard 10-20,
  10-10 and 10-05 labels (also T3–T6, `EEG Fp1-REF`, any case) are placed
  at their template positions. Only a cap with other labels (e.g. EGI's
  E1…E128) needs a channel file (positions on the template head) and, once,
  OpenMEEG to compute its leadfield (**Step 1**).
- **A subject table** (.csv, .tsv or .xlsx): one row per recording, with
  the ID (as in the file names), the group, the covariates (age, sex, ...)
  and any scores to correlate. Without it only stage 1 runs.

**Software**

- MATLAB R2024b or later with the Signal Processing and the Statistics and
  Machine Learning toolboxes; Parallel Computing is optional (recordings
  then run in parallel, with the same results). Tested on macOS (Apple
  silicon) with R2024b and R2026a; `.github/workflows/tests.yml` runs the
  tests on Linux, Windows and macOS once the folder is on GitHub.
- Three free toolboxes, not redistributed here (versions tested):

| Toolbox | Used for |
|---|---|
| FieldTrip (commit 0beabe6, April 2023) | reading the files, eLORETA, LCMV |
| Brainstorm 3 (3.250107, commit 69bbf53) | band-pass filters, region values (only `toolbox/math` and `toolbox/core` are used; the whole folder for another head model) |
| restingIAF (commit b25db38) | individual alpha frequency (spectral profile, individual bands) |
| *included:* NBS 1.2 (`3- Analysis/fun/NBS`, unchanged) | network-based statistic |

  The graph metrics follow the Brain Connectivity Toolbox's definitions
  but do not need it.
- Only for a cap with non-standard labels: OpenMEEG, as command-line
  programs or as the Python package (`pip install openmeeg scipy numpy`,
  then `S.openmeeg_python` = that Python).
- Only for another head model (`tools/make_anatomy.m`): OpenMEEG (above),
  Brainstorm's whole folder with its fix 952475b (`tools/fetch_toolboxes.sh`
  applies it) and, for T1 MRIs, SPM12 and CAT12, which Brainstorm installs
  as its plugins the first time (`docs/anatomy.md`).

## 2. Installation

1. **Folders.** Put this folder and the toolboxes side by side:
   ```
   <project>/v_static/       this folder
   <project>/Toolboxes/      fieldtrip-toolbox, brainstorm3, restingIAF
   <project>/Datasets/       your recordings and subject table (or anywhere: S.data, S.subs_info)
   ```
2. **Toolboxes.** From `v_static/`, in a terminal (bash and git; on
   Windows, Git Bash):
   ```bash
   bash tools/fetch_toolboxes.sh --static
   ```
   It downloads the three toolboxes at the tested versions into
   `../Toolboxes` (FieldTrip is about 400 MB: up to half an hour). Copies
   you already have can be used instead: on the MATLAB path, or in
   `S.toolbox_paths` (fields `fieldtrip`, `brainstorm`, `restingiaf`: each
   toolbox's root folder, e.g. `S.toolbox_paths.brainstorm =
   '/apps/brainstorm3'`).
3. **Check.** In MATLAB, from `v_static/`:
   ```matlab
   addpath(pwd);  S = setup_toolboxes(pipeline_settings());
   ```
   It puts the code and the toolboxes on the path, stops with the list of
   what is missing, and warns when a toolbox differs from the tested
   version (`S.toolbox_check`: `'warn'`, `'error'` or `'off'`).

## 3. Try it

```matlab
run_static_demo                                % 16 + 16 simulated recordings, both steps, 5–10 minutes
results = runtests('tests');  table(results)   % the tests: 10–25 minutes, no private data needed
```

`run_static_demo` simulates recordings in which the default mode network
(DMN) regions share a lagged alpha rhythm, more strongly in one group,
with a score correlated with that coupling; runs both steps on the files
over the whole recording (`S.static_matrices = {'full'}`: the mean of
0.6 s alpha windows, the default, does not find this weak and stable
coupling, `docs/methods.md`, "Whole recording or windows"); and writes
what was recovered to
`simulated/static_demo/results/recovery_report.md`. In the run here
(4.7 minutes): the measured within-DMN wPLI followed each recording's
planted coupling (r = 0.90), the groups differed in it (g = 2.17,
`p_bb` = 0.0005 after the adjustment over families), and NBS found one
component in the planted direction (109 edges); 14 other network blocks
differed too (source leakage, **Limits**). The planted score correlation
(r = 0.60) was not found in the 16 recordings that have the score
(measured r = 0.36).

## 4. Your study, step by step

**Step 1 — Recordings.** Put the cleaned files in one folder (subfolders
are fine, e.g. a BIDS `sub-*/eeg/`). What is accepted:

| | |
|---|---|
| Files | EEGLAB, EDF, BDF, BrainVision, FIF (all tested), EGI `.mff` folders; continuous or epoched; converted to µV. One recording per subject ID (BIDS sessions: an ID pattern with the session). Events of bad data (`S.bad_events`, MNE's `BAD_` annotations) are removed; boundary events (EEGLAB `boundary`, BrainVision `New Segment`) and stretches of NaN cut a recording into pieces, filtered separately; a piece shorter than 6 s (7.6 s in delta and theta) after removing its edges is dropped |
| Channels | matched by label (case, `EEG ` prefixes, `-REF` suffixes and T3/T4/T5/T6 = T7/T8/P7/P8 handled); EOG, ECG and other channels left out; recordings may miss different channels; flat channels left out; 16 at least (fewer than 32: a warning). On a cap with other names (an ANT equidistant or EGI cap) the channel sets given by 10-20 names (posterior for the IAF, frontal, frontopolar, temporal) take the electrodes nearest to their positions (`cap_channels.m`) |
| Caps | `'ASA_10-05_343'` (any 10-20/10-10/10-05 labels; the default), `'EEGLAB_61'` (the lab's cap), or a file: a Brainstorm `channel_*.mat` or a CSV `label,x,y,z` in metres on the ICBM152 head, whose leadfield is computed once (OpenMEEG; the positions and the sign of the gain are checked, also for caps without 10-20 labels) |
| Head model | the ICBM152 template (default), another template (Colin27, ...) or one head per subject from T1 MRIs (CAT12): folders made by `tools/make_anatomy.m` with Brainstorm's processes, `S.anatomy` (`docs/anatomy.md`) |
| Atlases | `'Desikan_Killiany_68'` (default), a Brainstorm scout file on the ICBM152 `tess_cortex_pial_low` surface, or an atlas of the cortex by name (`'Destrieux'`, `'DKT'`, `'Schaefer_100_17net'`, ...: 22 on the ICBM152). The network features need a region table like `3- Analysis/DK68_networks.csv` (`S.network_table`); without one they are skipped |
| Rates, lengths | any; set `S.core_duration` so that every recording contributes the same number of seconds |

**Step 2 — Subject table.** One row per recording, any column names:

```
participant_id,group,age,sex,PANSS_pos
sub-01,control,24,f,n/a
sub-02,patient,27,m,18
```

Any delimiter, decimal point or comma; `n/a` or empty = missing; IDs are
read as text (`001` stays `001`, except in an .xlsx cell that holds the
number 1 formatted as 001: type such IDs as text). An ID that matches no
recording exactly is matched ignoring case, then leading zeros, when that
is unambiguous (the matches are printed). Numeric codes of a categorical
covariate (site 1/2/3) must be listed in `S.categorical_covariates`, or
they enter as one linear covariate.

**Step 3 — Settings.** `pipeline_settings.m` is the template: every field
has a working default; the fields marked EDIT describe a study. Write a
small function that starts from it and changes what differs, e.g.
`settings_mystudy.m` in this folder:

```matlab
function S = settings_mystudy()
S = pipeline_settings();
S.data          = '/data/mystudy/derivatives/cleaned';   % the recordings
S.file_pattern  = {'sub-*/eeg/*_eeg.set'};                % which files (may name subfolders)
S.id_regexp     = '^(sub-[^_]+)';                         % file name -> subject ID
S.core_duration = 180;                                    % seconds analysed per recording
S.montages      = {'ASA_10-05_343'};                      % standard labels
S.subs_info     = '/data/mystudy/participants.tsv';
S.id_column     = 'participant_id';
S.group_column  = 'group';
S.groups        = {'control', 'patient'};
S.comparisons   = struct('name', 'patient - control', ...
                         'A', {{'group', '==', 'patient'}}, 'B', {{'group', '==', 'control'}});
S.covariates    = {'age', 'sex'};                         % text columns become indicators
S.correlations  = struct('name', 'Symptoms', 'columns', '^PANSS_', ...
    'subjects', {{'group', '==', 'patient'}}, 'covariates', {{'age'}}, 'zero_as_missing', false);
end
```

A misspelled field stops the run and names the closest known one; a
step called without settings uses `pipeline_settings()` (a variable `S`
left in the workspace is never used). Two complete examples: `examples/settings_bids.m` (a BIDS dataset with a
standard cap) and `examples/settings_openneuro.m` (138 public recordings
in two folders, two caps, two comparisons, correlations with symptoms
and cognition).

**Step 4 — Run.**

```matlab
addpath(pwd);  run_full_pipeline(settings_mystudy());
```

1. It checks the settings and toolboxes.
2. It reads the subject table and lists the recordings without a row, the
   rows without a recording, and the sample size of each comparison
   (`check_subject_table`).
3. It runs stage 1, the power spectra (1b) and stage 3 in order.

Results that are already up to date are reused, so after an interruption
or a change of settings the same command redoes only what is missing or
changed. An error stops the run and says how to resume:
- `run_full_pipeline(S, 'from', 3)` restarts at the statistics;
- `'to', 1` stops after stage 1;
- `'from', '1b', 'to', '1b'` runs only the power spectra
  (`S.psd.run = false` leaves them out).

The statistics do not start without a recording with stage 1 outputs. At
the end it prints the time of each stage, the recordings that failed in
stage 1 (left out of the statistics), and the files to check (step 5).

Time on 8 workers, for the 138 OpenNeuro recordings of 300 s: stage 1
1.6 to 4.5 hours with the six bands over the whole recording (`'full'`;
the longer beside three other analyses; 29 minutes with four bands and a
shorter filter; the window mean adds the windows); with 50,000
permutations, the PSD analysis 3.9 hours and the statistics 1.9 hours,
also beside other analyses (the statistics took 36 minutes with 10,000).

Each stage can also be run on its own, so that its outputs can be checked
before the next one: `S = setup_toolboxes(settings_mystudy());` then
`run_all(S)`, `run_psd_analysis(S)` and `run_statistics(S)`.

**Step 5 — Check, then read the results.**

*Stage 1.* Open `results/Connectivity/stage1_summary.csv`. It has one row
per recording, with:
- the status (a recording that failed is listed with the reason and left
  out of the statistics);
- the seconds analysed and the pieces;
- the channels used, missing and flat;
- the noisiest channel's SD ratio (above 10: check the cleaning);
- `data_rank`: the number of independent signals left after cleaning
  (each removed ICA component or interpolated channel lowers it; below
  `S.min_channels` a warning; `run_all` prints its range);
- the head model (`anatomy`).

*Power spectra (1b).* In `results/PSD/<S.run_name>/`:
1. `recordings.csv`: every recording should be `ok`, with a specparam fit
   R² above about 0.9; the log names channels missing in some recordings
   (left out).
2. `Statistical analysis/summary.csv` and `clusters.csv`.
3. `figures/`: scalp maps of the channels; cortex maps of the regions in
   six views, and as 3-D `.fig` files to turn with the mouse in MATLAB.

The spectra are kept, so a change of band or test recomputes nothing.
Details: `docs/psd.md`.

*Stage 3.* In `results/Statistics/<S.run_name>/` (`main` by default; give
each variant of the analysis its own `S.run_name`):
1. The start of the log (`Statistical analysis/analysis_log.txt`): the
   recordings in each comparison and correlation set (those of stage 1,
   made or reused, in the subject table); a warning names any that cannot
   be tested (a comparison needs at least 2 per group and the covariates
   + 4 in all, a correlation the covariates + 4).
2. `Statistical analysis/summary.csv`: tests and significant results per
   band, measure and family. `families.csv`: the families selected by the
   adjustment over families (`family_selected`).
3. The tables in `Statistical analysis/<inverse>_<band>_<measure>/`:
   `static_GT_window_mean_group_comparison.csv` (strength, clustering,
   path length), `static_networks_window_mean_*`,
   `static_regions_window_mean_*`, `static_power_*` (with
   `S.static_matrices = {'full'}`: `full` instead of `window_mean`). A family is one band, measure and feature family with
   one comparison or correlation set. In each table:
   - `p_maxT_family` is the test's family-wise p over its family (its
     statistic against the family's largest in each permutation);
   - `family_p` tests the family as a whole; the families of the run are
     selected on it (Benjamini–Hochberg: `family_selected`);
   - `p_bb` = `p_maxT_family` × families / families selected
     (Benjamini–Bogomolov), and `significant` = 1 when `p_bb` ≤ `S.alpha`;
   - `q_fdr`, the FDR q within the family, describes; it decides nothing.
   Read Hedges' g and its CI with every p.
4. In `NBS/<inverse>_<band>_<measure>/`: `NBS_window_mean.csv` (each
   component, with `p_both`, the p to compare with alpha, since both
   directions are tested), `enrichment_window_mean.csv`, the figures.
5. `static_reliability_window_mean.csv` (`static_reliability.csv` for
   `'full'`): features with a split-half `r_sb` below about
   0.5 are measured too noisily for a group difference to be trusted.

## 5. Settings (main fields of `pipeline_settings.m`)

| Field | Default | Meaning |
|---|---|---|
| `S.data`, `S.file_pattern`, `S.id_regexp` | `../Datasets`, `*.set`, `*.edf`, `*.bdf`, `*.vhdr` | the recordings and how the subject ID is read from a file name |
| `S.srate`, `S.resample` | `[]`, `[]` | expected rate (checked), resampling |
| `S.core_start`, `S.core_duration` | 0, `[]` | the analysed stretch (s); `[]` = whole recording |
| `S.montages`, `S.atlases`, `S.inv` | `ASA_10-05_343`, `Desikan_Killiany_68`, `eloreta` | cap, atlas, inverse method(s) (`eloreta`, `lcmv`, `wmne`); the statistics analyse the first cap and atlas |
| `S.eloreta_lambda`, `S.lcmv_lambda`, `S.wmne` | 0.05, `'5%'`, SNR 3 | regularisation (eLORETA's is relative; a lower wMNE SNR regularises more); give each value its own `S.results`. On the OpenNeuro recordings eLORETA's 0.01 to 0.1 gave the same results and 0.5 too much regularisation (`docs/validation.md`) |
| `S.scout_function` | `'mean'` | each region's signal from its sources, Brainstorm's scout functions: `'mean'` (sign-flipped), `'pca'` (first principal component); also `'max'`, `'power'`, `'rms'` (not linear; power and RMS have no phase: `docs/methods.md`). The PCA changes each recording's matrices, not the group results (`docs/validation.md`) |
| `S.anatomy`, `S.anatomy_individual`, `S.anatomy_fallback` | `''` (the shipped ICBM152), `false`, `false` | another head for every recording (a folder), or one folder per subject ID; subjects without one stop the run unless the fallback lets them use the ICBM152 (`docs/anatomy.md`) |
| `S.bands` | delta 1–4, theta 4–7, alpha 7–13, beta1 13–20, beta2 20–30, gamma 30–45 Hz | any number of bands |
| `S.individual_bands` | `false` | `true`: bands anchored on each recording's individual alpha frequency (IAF ± 2 Hz for alpha; `S.bands_individual`, `S.iaf`; `docs/individual_alpha_frequency.md`) |
| `S.conn` | `wpli`, `aec_orth` | measures: `plv`, `wpli`, `pli`, `ciplv`, `dwpli`, `aec`, `aec_orth`, `aec_symorth`, `plv_symorth` |
| `S.static_matrices` | `{'window_mean'}` | how each measure is estimated: `'window_mean'`, in windows then averaged (the original pipeline's static matrix); `'full'`, over the whole recording; or both, each analysed on its own. Short windows bias phase measures upwards, the more so where the rhythm is stronger (`docs/methods.md`, "Whole recording or windows") |
| `S.window_cycles` or `S.window_length`, `S.window_overlap`; `S.aec_window`, `S.aec_step` | 6 cycles of the band centre (alpha 0.6 s), `[]`, 0.5; 6 s, 0.5 s | the windows of the window mean: phase measures in windows of a number of cycles of each band's centre or of a fixed length in seconds for every band, overlapping by the given fraction; AEC measures in windows of 6 s every 0.5 s (O'Neill et al. 2017) |
| `S.measure_windows` | none | windows of one measure instead of the common ones, e.g. `struct('wpli', struct('window_length', 2), 'aec_orth', struct('aec_window', 10, 'aec_step', 1))` (phase measures: `window_cycles`, `window_length`, `window_overlap`; AEC: `aec_window`, `aec_step`). Several window settings for one measure: separate runs, each with its own `S.results` |
| `S.filter_transition`, `S.edge_trim`, `S.min_segment` | `'brainstorm'`, 1, 6 | band-pass width (`docs/methods.md`), seconds dropped at each edge of a piece, shortest piece kept (6 s: the 6 s AEC windows fit in every piece; 4 is enough with `'full'` only) |
| `S.subs_info`, `S.id_column`, `S.group_column`, `S.groups` | `participants.tsv`, `participant_id`, `group`, control/patient | the subject table; `S.group_column = ''`: the group is the start of the ID (a column `group` is added for the comparisons) |
| `S.comparisons`, `S.covariates`, `S.categorical_covariates` | patient − control, age, none | comparisons (conditions on the table's columns), covariates (text columns become indicators; numeric codes too when listed as categorical) |
| `S.correlations`, `S.covariate_sets` | none | correlation sets: columns (regular expression), subjects, covariates |
| `S.static_measures`, `S.static_families` | every measure; global, networks, regions, power | what the statistics analyse |
| `S.n_perm`, `S.alpha`, `S.seed` | 10000, 0.05, 20261015 | permutations (a family can be selected only when 1/(`S.n_perm` + 1) ≤ `S.alpha` / families: the log says how many a run needs), level, random seed |
| `S.nbs`, `S.enrichment` | threshold 3.1, 5000 permutations, extent; 10000 | network-based statistic, network enrichment |
| `S.reliability_block`, `S.spectral_profile`, `S.figures`, `S.run_name` | 10 s, true, true, `main` | split halves, spectral profile, figures, output folder |
| `S.psd` | run, sensor and source levels, 8 bands (delta 1–4 … beta2 20–30 Hz), relative to 1–30 Hz, five measures, 4 s Welch segments, specparam 2–40 Hz, cluster test | the PSD analysis (`docs/psd.md`); its statistics use `S.comparisons`, `S.covariates`, `S.correlations`, `S.alpha`, and `S.psd.n_perm` permutations (`[]`, the default: `S.n_perm`) |

Settings of the network states (and of the windows, when the window mean
is not asked for) are not used here; one set to a non-default value gives
a warning (`check_settings:not_used`).

## 6. Outputs

`results/` (not shared; delete it to start afresh):

- `Connectivity/cmats/<recording>/static_cmats_<inverse>_<cap>_<atlas>_<band>.mat`:
  `static_full_cmats` (measures × regions × regions, in the order of
  `S.conn`), `static_halves` (odd and even blocks), `static_power`,
  `static_power_halves`, `info` (band, regions, channels, pieces, seconds,
  filter, settings); with the window mean (the default) also
  `static_cmats` (the window mean of each measure), `static_halves_window`
  (its odd and even blocks) and the windows of each measure in `info`
  (sizes and positions; the windowed matrices themselves are not saved). `Connectivity/stage1_summary.csv`; `iaf.csv` with
  individual bands.
- `Statistics/<run_name>/Statistical analysis/`: `summary.csv`,
  `families.csv`, `analysis_log.txt`, `spectral_profile/`, and per
  `<inverse>_<band>_<measure>/` the families' `_group_comparison.csv` and
  `_correlations.csv` (also .xlsx: n, means, g, t or G, p, `p_maxT_family`
  and `family_p` over the family, `q_fdr` within it, `p_bb` after the
  adjustment over families, `significant`),
  `features_window_mean.csv` (every feature of every recording),
  `static_reliability_window_mean.csv`, figures (with `'full'`:
  `features_full.csv`, `static_reliability.csv`, and the figures ending in
  `full`). Power does not depend on the measure:
  it is tested once per band and inverse method, in the first measure's
  folder.
- `Statistics/<run_name>/NBS/<inverse>_<band>_<measure>/`: `NBS_<matrix>.csv`
  and `.mat`, `enrichment_<matrix>.csv`, `edges_<comparison>_<matrix>.csv`
  (`<matrix>`: `window_mean` or `full`; every
  edge's group means, difference and G), figures of the edge statistics
  and of the components with `p_both` < alpha.
- `PSD/<run_name>/` (the PSD analysis; the spectra kept in `PSD/spectra/`):
  `recordings.csv`, `values_<level>_<measure>.csv`, `psd_values.mat`,
  `Statistical analysis/` (`summary.csv`, `families.csv`, `clusters.csv`,
  the families' tables), `figures/` (spectra, box plots, scalp and cortex
  maps), `psd_log.txt` (`docs/psd.md`).
- `provenance_*.mat`: settings, toolbox versions and summary of each step.

## 7. Limits

**For** resting-state scalp EEG, already cleaned, compared between two
groups (or any two sets of subjects defined by the subject table) with
covariates, and correlated with scores. **Validated on** adults'
eyes-open rest with 60–61-channel caps and on simulations (**Validation
and tests**).

- **Cleaned data only**: no artefact removal or interpolation (average
  reference, flat channels left out; a channel SD more than 10 times the
  median gets a warning).
- **Template head** by default: ICBM152 head, cortex and electrode
  positions for everyone; not for children. Another template, or each
  subject's head from a T1 MRI, can be used (`S.anatomy`,
  `tools/make_anatomy.m`); the caps are then template positions carried to
  each head through the MNI normalization, not digitised positions.
- **Source leakage**: a region's estimate mixes in its neighbours. wPLI,
  PLI, ciPLV and the orthogonalised measures ignore the zero-lag part, yet
  in simulation stronger coupling planted in the DMN also made the
  neighbouring network blocks differ (salience–DMN more than the DMN
  itself) and the NBS component spread beyond it. Read a network or edge
  result as "around these regions". Cleaning lowers the number of
  independent signals further: the OpenNeuro recordings (60 channels,
  ICA-cleaned) have a data rank of 9–49, median 33, so their 68 regions
  are mixtures of about 33 signals (`stage1_summary.csv`, `data_rank`).
- **Same length for everyone**: every estimate depends on the amount of
  data; set `S.core_duration` (the statistics warn when the groups differ
  in seconds analysed).
- **Many bands and measures** make many families: all are adjusted for
  together, so keep those of your question. NBS p values are corrected
  over edges, not over bands and measures.
- **Designs**: two-group contrasts with linear covariates, and
  correlations; no more than two groups at once, interactions or repeated
  measures. Not for MEG or intracranial data.

## 8. Tests and validation

After installing, run the tests from this folder to check your setup
(10–25 minutes; no data of yours needed):

```matlab
results = runtests('tests');  table(results)
```

Three tests are skipped by design (they test the windows and the states,
which this pipeline does not have), and a test whose reference toolbox
(PALM, PERMUTOOLS, ROInets) is not installed is skipped, not failed.

The methods were checked against reference toolboxes (Brainstorm,
FieldTrip, PALM, ROInets, PERMUTOOLS), on simulations with planted
effects, on 138 public recordings and on datasets the pipeline had never
seen: see `docs/validation.md`.

## 9. Reproducibility and shared code

Every random step is seeded (`S.seed`); parallel and serial runs agree up
to rounding; each step saves its settings, toolbox versions and summary;
re-running recomputes only what changed (a recording, a setting, the
code) and the statistics stop if stage 1's outputs were made with other
settings or code.

The code files are those of v_dynamic (the dynamic pipeline), unchanged:
`tests/test_same_code_as_v_dynamic.m` checks it, and a change to them is made
there and copied here (`docs/shared_code.md`).

## More detail

- `docs/methods.md`: filters (with each band's filter length and
  transient), measures, graph metrics, statistics and the adjustment over
  families, reliability.
- `docs/validation.md`: every test, the simulations, the usability test,
  the runs on real data.
- `docs/psd.md`: the PSD analysis: methods, outputs, checks, validation,
  and what changed from the lab's original scripts.
- `docs/individual_alpha_frequency.md`: the individual bands.
- `docs/anatomy.md`: another head model (a template, or each subject's
  head from a T1 MRI), made with Brainstorm by `tools/make_anatomy.m`.
- `docs/shared_code.md`: how the code is kept identical to v_dynamic's.
- `docs/references.md`: the papers to cite, with their DOIs.
- `CHANGES.md`: what changed and why.

## How to cite

Cite EEGSourceFC as given in `CITATION.cff` (on GitHub: "Cite this
repository"):

> Aubonnet R (2026). *EEGSourceFC: static EEG source connectivity and power
> spectra.* Version 1.0.0, MATLAB software.

Please also cite:
- the lab code it was built from: var-EEG-FC (Allouch et al. 2023) for
  stage 1, and dynCogPD (Duprez et al. 2022);
- the toolboxes it calls: FieldTrip, Brainstorm, OpenMEEG, NBS,
  restingIAF, ROInets;
- the methods you report: the inverse solution, the connectivity measure,
  the statistics.

`docs/references.md` gives each of them with its paper and DOI.

## Licence

GPL-3.0-or-later: the full text in `LICENSE`; `LICENSE.md` adds the
notices of the third-party code and data included here.
