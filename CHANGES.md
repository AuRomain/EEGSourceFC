# Origin of every file

Every code file of this folder is an unchanged copy of the same file of the
lab's brain-network-state pipeline (`../v_dynamic`, "v_dynamic"), whose
static-only mode (`S.dynamic = false`) this folder runs.
`tools/v_dynamic_files.txt` names the v_dynamic commit and lists every shared file with
its SHA-256; `tests/test_same_code_as_v_dynamic.m` checks them (and compares them
with v_dynamic when it is next to this folder); `tools/sync_from_v_dynamic.sh` copies
them again.

## Shared with v_dynamic (unchanged copies)

| Folder | Files |
|---|---|
| `1- EEG source connectivity/` | every file v_dynamic tracks there but its original 61-channel and ANT caps (read only by v_dynamic's `tools/make_corrected_inputs.m`; the original 61-channel positions are rotated): stage 1 (`run_all.m`, `run_pipeline.m`, reading, montage, head model, inverse solutions, scouts, band-pass, IAF, `get_all_connectivity_mats.m`, `static_connectivity.m`, `split_halves.m`, `+connectivity/`, `+wmne/`, `+ROInets/`) and `inputs/` (template head and cortex, the Desikan-Killiany atlas, the caps and their leadfields) |
| `1b- PSD analysis/` | every file v_dynamic tracks there: the PSD analysis (`run_psd_analysis.m`; spectra, measures, statistics with cluster tests, scalp and cortex maps; `docs/psd.md`) |
| `3- Analysis/` | `run_statistics.m`, `static_analysis.m`, `compute_NBS_analysis.m`, `spectral_profile_features.m`, `DK68_networks.csv`; `fun/`: the statistics (families, Welch t and G, ANCOVA and partial correlations with permutations, FDR, Benjamini–Bogomolov), the sample-size check (`check_sample_sizes.m`), the static features, reliability, network enrichment, spin permutations, the cluster test of the PSD analysis (`cluster_perm_test.m`), figures, NBS 1.2 (`NBSstats.m`, `get_components.m`) |
| top folder | `cap_channels.m`, `check_settings.m`, `settings_key.m`, `stage_settings.m`, `stage1_outputs.m`, `stage2_settings.m`, `analysis_dir.m`, `setup_toolboxes.m`, `toolbox_*`, `code_fingerprint.m`, `gt_window_metrics.m`, `read_subject_table.m`, `get_subject_groups.m`, `subject_ids.m`, `band_edges.m`, `bns_run_dir.m`, `file_tag.m`, `file_fingerprint.m`, `hash_text.m`, `run_full_pipeline.m` (the stages one after the other: 1, the PSD analysis 1b, 3; stage 2 is skipped here), `run_static_demo.m` |
| `simulation/` | `run_static_simulation.m`; `+sim/`: the simulator of recordings on the pipeline's forward model and the static dataset (`static_params.m`, `simulate_static_dataset.m`) |
| `tests/` | `test_static_estimators`, `test_static_settings`, `test_static_stage1`, `test_static_statistics`, `test_static_integration`, `test_static_ground_truth`, `test_anatomy`, `test_psd`, `test_reference_connectivity`, `test_reference_statistics`, their settings and recordings (`static_test_settings.m`, `static_test_recordings.m`); `tests/python/test_build_leadfield.py` |
| `tools/` | `fetch_toolboxes.sh` (`--static`: FieldTrip, Brainstorm with its fix 952475b, restingIAF), `brainstorm_base64_952475b.patch`, `build_leadfield_openmeeg.py`, `compare_with_brainstorm.m`, `make_anatomy.m` (other head models) |
| `docs/` | `individual_alpha_frequency.md`, `anatomy.md`, `psd.md` |

v_dynamic's stage 2 (`2- Dynamic connectivity/`: windows → states), its state
analyses and their tests are not here: `check_settings` sets
`S.dynamic = false` when that folder is absent and refuses `true`.

## This folder's own files

| File | What |
|---|---|
| `pipeline_settings.m` | The settings template for static analyses: `S.dynamic = false`; six bands (delta 1–4, theta 4–7, alpha 7–13, beta1 13–20, beta2 20–30, gamma 30–45 Hz) and the individual-band option (`S.individual_bands`, bands anchored on each recording's alpha frequency); wPLI and orthogonalised AEC, as the mean of their windows (`S.static_matrices`, `S.measure_windows`; `S.min_segment` 6 s so that the AEC windows fit); Brainstorm's filter width; the four feature families; the split-half blocks; outputs in `results/Statistics/main/`; `S.scout_function` (the mean) and the regularisation explained; the PSD analysis (`S.psd`, v_dynamic's defaults). Every other field is v_dynamic's default |
| `examples/settings_bids.m`, `examples/settings_openneuro.m` | A BIDS dataset with a standard cap; the OpenNeuro first-episode psychosis data |
| `tests/test_same_code_as_v_dynamic.m` | Every shared file equals its v_dynamic copy; no other code file |
| `tests/test_template.m` | The template and the examples pass `check_settings`, are static, and give the documented bands (fixed and individual) and PSD settings |
| `tools/sync_from_v_dynamic.sh` | Copies the shared files from v_dynamic and writes `tools/v_dynamic_files.txt` |
| `docs/methods.md` | Stage 1 step by step, the filters, the measures, the graph metrics, the statistics and the adjustment over families, reliability |
| `docs/validation.md` | Every test, the validation, the usability test, the runs on real data |
| `docs/shared_code.md` | How the code files are kept identical to v_dynamic's (`tools/sync_from_v_dynamic.sh`, `tests/test_same_code_as_v_dynamic.m`) |
| `README.md`, `CHANGES.md`, `LICENSE`, `LICENSE.md`, `.github/workflows/tests.yml`, `.gitattributes`, `.gitignore` | |

## History

The development history of this folder (33 commits, 3–6 November 2026) is kept
in the lab's archive; version 1.0.0 is published as a single commit. The
commit numbers of v_dynamic quoted below are those of its own history.

| Date | |
|---|---|
| 3 November 2026 | First version: v_dynamic's stage 1 and statistical engine adapted to whole-recording connectivity (its own `static_settings.m`, `static_connectivity.m`, `2- Statistics/` with its own features, edge analysis and reliability) |
| 3 November 2026 | PLV, ciPLV, wPLI, PLI and dwPLI compared with Brainstorm's connectivity process; seven measures through stage 1 |
| 3 November 2026 | Rebuilt from v_dynamic: what the first version added was moved into v_dynamic (its static-only mode, `S.dynamic = false`), with the same estimators and analyses, and this folder now holds v_dynamic's files unchanged. Changes for users of the first version: the statistics folder is `3- Analysis` (was `2- Statistics`); the template is `pipeline_settings.m` (was `static_settings.m`) and has six bands (was theta, alpha, beta); the demo is `run_static_demo` (was `run_demo`); the stage 1 files are `results/Connectivity/cmats/<recording>/static_cmats_*.mat` (was `results/Connectivity/static/`), with `static_full_cmats`, `static_halves`, `static_power`, `static_power_halves`; the statistics write to `results/Statistics/<run>/Statistical analysis/<inverse>_<band>_<measure>/` and `results/Statistics/<run>/NBS/<inverse>_<band>_<measure>/`; the global metrics are `Strength` (the mean node strength: (regions − 1) × the former mean connectivity), `ClustCoeff`, `CPL` (v_dynamic's `gt_window_metrics`, BCT's definitions as before); the reliability file is `static_reliability.csv`. Stage 1 outputs of the first version are recomputed (other keys) |
| 4 November 2026 | The open points of the audit (below), copied from v_dynamic edde89a |
| 4 November 2026 | The README rewritten as a step-by-step guide; `docs/methods.md` and `docs/validation.md`; copied from v_dynamic 71d7a86 (references to README sections); runs from scratch on real data (below) |
| 4 November 2026 | The scout function (`S.scout_function`: mean, PCA, max, power, RMS), copied from v_dynamic 6ebf912 |
| 4 November 2026 | Checks of a new cap's leadfield also for caps without 10-20 labels, copied from v_dynamic e16554a (found on a test dataset with an ANT equidistant cap) |
| 4 November 2026 | Statistics on small or unusual datasets (sample sizes reported, groups from the IDs usable in comparisons, T3/T4 temporal), copied from v_dynamic ae95090 |
| 4 November 2026 | The sensitivity of the results to eLORETA's regularisation and to the scout function, on the 138 OpenNeuro recordings (`docs/validation.md`); no code changed |
| 4 November 2026 | The sample-size check counts recordings whose stage 1 was reused, copied from v_dynamic f297a24 |
| 4 November 2026 | Channel sets (posterior for the IAF, frontal, frontopolar, temporal) matched by position on caps without 10-20 names, copied from v_dynamic be17811 |
| 4 November 2026 | The runs on the NIRP and Geneva datasets (`docs/validation.md`, "Two datasets never seen before"); 52 tests on R2026a and R2024b; no code changed |
| 4 November 2026 | The data rank in `stage1_summary.csv` (`data_rank`), copied from v_dynamic 3a7413e; 53 tests |
| 4 November 2026 | Other head models (`S.anatomy`, `S.anatomy_individual`, `S.anatomy_fallback`; `tools/make_anatomy.m`, `docs/anatomy.md`; Brainstorm's fix 952475b in `fetch_toolboxes.sh`), copied from v_dynamic a6c9da4; `test_anatomy`; 59 tests |
| 4 November 2026 | `run_full_pipeline.m`: stage 1 and the statistics with one command (README, "Steps 4 and 5 at once"), copied from v_dynamic 5ba01ca; `test_static_integration/test_full_pipeline`; 60 tests |
| 4 November 2026 | The README stands on its own: v_dynamic is named once (connectivity over time, section 9); how the code is kept identical to v_dynamic's moved to `docs/shared_code.md`; no code changed |
| 4 November 2026 | The other head models on the 138 OpenNeuro recordings (`docs/validation.md`, "Other head models"; README, Validation); the regularisation table's row "mean pattern over the recordings" held each recording's mean connectivity across the recordings: both rows now, with the mean pattern's own values (λ 0.01 0.85–0.96, PCA 0.55–0.84); no code changed |
| 4 November 2026 | New names: the dynamic pipeline v3 is v_dynamic (folder `v3_original_improved` → `v_dynamic`) and v2 is v_crosscheck, in every text; `tools/sync_from_v3.sh` → `tools/sync_from_v_dynamic.sh`, `tools/v3_files.txt` → `tools/v_dynamic_files.txt`, `tests/test_same_code_as_v3.m` → `tests/test_same_code_as_v_dynamic.m`; shared files copied from v_dynamic a18e943 (stage 1 changed only in a comment: `equivalent_code.txt` keeps the outputs of earlier runs up to date, checked on both OpenNeuro runs); nothing computed differently; 60 tests |
| 5 November 2026 | The full text of the licence in `LICENSE` (GNU GPL version 3 of 29 June 2007, as the Free Software Foundation published it); `LICENSE.md` and the README point to it; `docs/shared_code.md` and `tools/sync_from_v_dynamic.sh` list it among this folder's own files; no code changed |
| 5 November 2026 | The PSD analysis (stage 1b, `1b- PSD analysis/`, from the lab's scripts in `v_original/1b- PSD analysis`): the spectra of every channel and region, absolute, relative and periodic power in eight bands, the 1/f exponent and offset, the alpha frequency, the comparisons and correlations per channel and region with a cluster test (`3- Analysis/fun/cluster_perm_test.m`), scalp and cortex maps (`docs/psd.md`); `run_full_pipeline` runs it after stage 1 (`S.psd.run`); `S.psd` in `pipeline_settings.m` (section "PSD analysis", v_dynamic's defaults); copied from v_dynamic 8cb3d67 (its CHANGES.md, "PSD analysis: stage 1b"); `test_psd` (7), `test_static_integration/test_psd_step`, `test_template` checks `S.psd`; the statistics and stage 1 compute the same as before; 68 tests (67 passed, 1 skipped by design) on R2026a and R2024b |
| 5 November 2026 | The fixes of the newcomer review (below), copied from v_dynamic 3a8ad47 |
| 5 November 2026 | The max-statistic family test (below), copied from v_dynamic bd0bfcc |
| 5 November 2026 | `docs/psd.md` names the lab's original PSD scripts as not distributed, copied from v_dynamic 5195755; no code changed |
| 5 November 2026 | Methods checked against verified code (PALM, FieldTrip, MATLAB's own functions): the permutation statistics on z with the Welch df, as PALM; reference tests; copied from v_dynamic 6242a27 (below) |
| 5 November 2026 | Stage 1 comments no longer name the removed development tools, with the row of `equivalent_code.txt` that keeps the stage 1 outputs of 5 November accepted; copied from v_dynamic 20130e7; no code changed (below) |
| 5 November 2026 | `tests/test_same_code_as_v_dynamic.m` (this folder's own) allows the call of `duration_rank_summary` in the shared `tests/test_reference_statistics.m`, whose duration check skips itself here (no states); the guard test had failed since the copy of the methods on verified code (5 November) |
| 6 November 2026 | README, "Try it": the demo's numbers from its run of 5 November (4.7 minutes; `p_bb` 0.0005 instead of `q` 0.001, the rule before the max-statistic family test; 14 other blocks instead of 20; the score correlation, not found in 16 recordings); the tests take 10–25 minutes (24.5 for 89 tests beside three other analyses; 8.7 for 72 on 5 November), not 3–5; no code changed |
| 6 November 2026 | `1b- PSD analysis/psd_figures.m`: `savefig` without `'compact'`, which has no effect on R2026a and warned 12 times in the OpenNeuro PSD analysis; copied from v_dynamic f17e9d2 |
| 6 November 2026 | `docs/validation.md` and README on the final code: the simulations of `test_static_ground_truth` (run of 5 November: 8 other blocks, NBS 37 edges, the score's correlation not significant after the correction; null 0 of 214), the OpenNeuro example end to end (stage 1 4.5 h, PSD 3.9 h, statistics 1.9 h with 50,000 permutations; 2 and 10 significant correlations in the patients, no group difference, no NBS component), NIRP, Geneva and the lab's eight recordings again (no error), the times with 50,000 permutations, and the tests: 86 passed, 3 skipped by design of 89 on R2026a and R2024b; no code changed |
| 6 November 2026 | Public name **EEGSourceFC** (the folder stays `v_static` in the lab's project; v_dynamic is **EEGSourceStates**): README title and section "How to cite"; `CITATION.cff` (author Romain Aubonnet, GPL-3.0-or-later; the code it was built from, var-EEG-FC and dynCogPD, EEGSourceStates and the toolboxes as references; valid against the Citation File Format 1.2.0 schema; ORCID and affiliation still to add, version and DOI with the first release); `docs/references.md` (the papers of the base code, toolboxes, head model, atlas and methods, each checked against PubMed, Crossref or arXiv); README section 8, "Tests and validation", shortened to how to run the tests, pointing to `docs/validation.md`, which holds every detail it listed; no code changed |
| 6 November 2026 | Version 1.0.0, released 6 November 2026: `CITATION.cff` with the author's ORCID (0000-0002-5395-775X) and affiliation (Laboratory of Psychophysiology and Cognitive Neuroscience, Department of Systems Medicine, University of Rome Tor Vergata), the version and the date (and EEGSourceStates 1.0.0 among its references); README, "How to cite": the year and the version; still to add once hosted: the repository's URL and, optionally, a DOI |
| 6 November 2026 | Published as a single commit, without the development history (kept in the lab's archive): the commit numbers of this folder removed from this table and from `docs/validation.md`; no code changed |
| 6 November 2026, evening | Version 1.0.0 published again as a single commit, now with the two changes below (whole recording or windows; the window mean by default, windows of one measure); `docs/validation.md` and README: the tests on this code, and the OpenNeuro example's results labelled as those of the whole recording (`'full'`, the default when they were run); the abstract of `CITATION.cff` and the comment of `S.conn` in `pipeline_settings.m`: the window mean or the whole recording |
| 6 November 2026 | The commit numbers of v_dynamic quoted here updated: the co-author line was removed from v_dynamic's commit messages (and from this folder's), which gave every v_dynamic commit a new number; `equivalent_code.txt` copied again (only its notes column changed, which quotes them); no code changed |
| 6 November 2026, late | Version 1.0.0 published again as a single commit (the previous one kept in the lab's archive), now with the demo fixed. `run_static_demo` stopped at its report: it took the default window mean, but the report (`simulation/run_static_simulation.m`) reads the results of the whole recording (run of the demo as published: 5 minutes, nothing found, then `readtable` failed on `static_networks_full_group_comparison.csv`; the tests use `'full'` and did not see it). The demo now asks for `S.static_matrices = {'full'}` (the planted coupling is weak and stable: within-DMN g 0.42, family-wise p 0.99 as the mean of 0.6 s alpha windows in that run); `run_static_simulation` stops before simulating when `'full'` is not asked for (`run_static_simulation:full`); new test `test_simulation_needs_the_whole_recording` in `tests/test_static_settings.m`; copied from v_dynamic e43b2d3 (branch `static-windows`). README, "Try it": the demo over the whole recording; `docs/validation.md`: the tests on this code (97 passed and 3 skipped by design of 101) |

### Audit and usability test of the rebuild (3 November 2026)

An independent review of this folder and a run on unusual inputs (below
and README, **Validation**) found the following. The shared files were
fixed in v_dynamic (its CHANGES.md, last section) and copied here:

- band power (family `power`) was tested once per measure, so the same
  family entered the adjustment over families twice per band: now once
  per band and inverse method, in the folder of the first measure;
- BDF files failed on Apple silicon (the validated FieldTrip has no
  `read_24bit` mex file for it): `1- EEG source connectivity/read_24bit.m`
  reads them, and the format test includes BDF;
- an epoched file of short epochs stopped with an unclear message: the
  message now gives the longest piece and the reason, and the README the
  minimum epoch length;
- a recording without a row in the subject table stopped the spectral
  profile's statistics: it is left out with a warning;
- NBS printed a line per permutation into `analysis_log.txt`: quiet now;
- the demo needed restingIAF: FieldTrip and Brainstorm now;
- with the individual bands, a recording without an alpha peak failed
  entirely (the upper edge of `delta_iaf`, IAF − 6, was NaN and only the
  lower edge was checked): its anchored bands are now left out and its
  fixed bands computed, as documented
  (`test_static_stage1/test_individual_bands_without_alpha`);
- two tests needed v_dynamic's stage 2: skipped here;
- `docs/individual_alpha_frequency.md` listed the IAF of each private
  development recording: ranges only.

And in this folder's own files: `.gitattributes` (LF line endings on
every system, NBS kept byte for byte; the identity test compares bytes),
`core.autocrlf false` in the GitHub workflow, `tools/sync_from_v_dynamic.sh`
checks every file before removing or copying and replaces
`tools/v_dynamic_files.txt` only at the end, the original 61-channel and ANT caps
no longer copied, `test_same_code_as_v_dynamic` leaves a study's
`settings_<study>.m` in the top folder to the user, and the README's
filter table has rows for beta1 and beta2 and the delta filter's
attenuation at 0 Hz (recomputed at 256, 512 and 1000 Hz).

### Open points of the audit (4 November 2026)

Fixed in v_dynamic (commit edde89a, its CHANGES.md) and copied here:

- the settings of the windows and states, set to other values than their
  defaults in a static run, were ignored without a message: a warning
  (`check_settings:not_used`) now names them
  (`test_static_settings/test_settings_of_the_windows_are_reported`);
- the leadfield of a new cap on Windows: `where` instead of `which` to
  find the OpenMEEG programs, the Python command line quoted for
  `cmd.exe`, `python` by default (not run on Windows here).

In this folder's own tests: `test_same_code_as_v_dynamic/test_no_v_dynamic_file_missing`
checks that no function called here exists only in v_dynamic (with v_dynamic's code at
the end of the path, every direct dependency of every code file here must
resolve to this folder; the 15 files that `run_statistics.m` reaches only
with `S.dynamic = true`, and the stage 2 of a skipped test, are listed as
expected); it fails when a shared file is missing (checked by removing
`welch_G.m`). `test_template` checks that the template and the examples
give no warning.

### README, documentation and runs on real data (4 November 2026)

The README was rewritten as a step-by-step guide for a new reader (what
you need, installation, a first run, a study in six steps, the main
settings, the outputs, what it is for and its limits, validation, staying
in step with v_dynamic); the details moved to two files of this folder,
`docs/methods.md` and `docs/validation.md`. Copied from v_dynamic (commits
55d31b9 and 71d7a86, README restructured there too): the references to
README sections in messages and comments (`setup_toolboxes.m`,
`run_static_demo.m`, `get_ftHeadmodel_and_ftLeadfield.m`,
`tests/python/test_build_leadfield.py`, `docs/individual_alpha_frequency.md`).
The change of a message in a stage 1 file is listed in
`1- EEG source connectivity/equivalent_code.txt`, so stage 1 outputs of
the previous code are not recomputed. The runs from scratch on the
OpenNeuro data and on the eight development recordings are in
`docs/validation.md`, **Runs on real data**.

### Scout functions (4 November 2026)

Copied from v_dynamic (commit 6ebf912, its CHANGES.md): `S.scout_function`
chooses how each region's signal is formed from its sources, with
Brainstorm's scout functions: `'mean'` (the default, as before: the
sign-flipped mean), `'pca'` (the first principal component, Brainstorm's
current PCA), and `'max'`, `'power'`, `'rms'` (not linear: computed from
the source time series, region by region; power and RMS have no sign, so
phase measures on them give a warning). Outputs made with the mean keep
their key (`equivalent_code.txt`). `pipeline_settings.m` (this folder's)
has the new field and explains the regularisation settings. Three new
tests (`test_static_settings`, `test_static_stage1`).

Checked afterwards on the 138 OpenNeuro recordings (`docs/validation.md`,
"Regularisation and scout function"): eLORETA's λ 0.01, 0.1 and 0.5 and
the PCA, one at a time against λ 0.05 with the mean. λ 0.05 lies in a
flat range (0.01 to 0.1 give the same results), 0.5 regularises too much;
the PCA changes each recording's matrices but not the group results. The
mean stays the default.

### Checks of a new cap (4 November 2026)

Copied from v_dynamic (commit e16554a, its CHANGES.md). Found when testing on a
dataset recorded with ANT's 129-channel equidistant cap (labels Z1, LL1,
RA1, ...): `tools/build_leadfield_openmeeg.py` checked the sign and
localisation of a new cap's leadfield only at electrodes labelled Cz,
Fpz, Oz, T7, T8, so this cap got no check. It now uses the electrodes
nearest to these positions when the labels are absent, measures the
distance to the scalp surface (not to its nearest vertex), and says when
no label can show a rotated cap. Two new Python tests
(`tests/python/test_build_leadfield.py`, 4 in all).

### Statistics on small or unusual datasets (4 November 2026)

Copied from v_dynamic (commit ae95090, its CHANGES.md). Found by running both
pipelines from scratch on two test datasets never seen before (NIRP: 4
recordings, 27 channels labelled with T3/T4; Geneva: 6 recordings, a
subject table without a group column):

- `check_sample_sizes.m` (new, shared; added to `tools/sync_from_v_dynamic.sh`):
  before the families, the log lists the recordings of every comparison
  and correlation set, and a warning names those that cannot be tested.
  With 2 + 2 recordings and age as covariate, every comparison gave rows
  of NaN and the correlations no rows, without a word, and the summary
  counted the rows; it now counts the tests with a p-value.
- `read_subject_table.m`: with `S.group_column = ''` (groups from the
  start of the IDs) the table gets a column `group`, so the comparisons can
  name the groups; before, every family stopped on "Column group is not in
  the subject table".
- `spectral_profile_summary.m`: T3 and T4 count as temporal channels (the
  muscle index `hf_temporal` was NaN for caps with the old names).
- Three tests in `test_static_statistics`.
- v_dynamic f297a24: `check_sample_sizes.m` counted only the
  recordings of status `ok` in `stage1_summary.csv`; when `run_all` reused
  stage 1 (status `skipped (up to date)`), statistics run again on
  existing outputs reported 0 recordings and every set as not testable
  (found re-running NIRP's statistics). Made or reused now both count. The
  families themselves were not affected.
- v_dynamic be17811: `cap_channels.m` (new, shared; added to
  `tools/sync_from_v_dynamic.sh`), `setup_toolboxes.m`,
  `spectral_profile_summary.m`, `spectral_profile_features.m`. The
  posterior channels of the IAF (`S.iaf.channels`) and the frontal,
  frontopolar and temporal channels of the spectral profile were 10-20
  names: on the Geneva ANT cap (Z1 ... RR14) every set was empty, the IAF
  could not be found and those measures were NaN. A cap with none of a
  set's names now takes, for each name, its electrode nearest to the
  name's position on the ASA 10-05 cap (within 20 mm); `setup_toolboxes`
  so matches `S.iaf.channels` for cap files, with a message. Caps with
  10-20 names keep their channels and `S.iaf`; the spectral profile is
  computed again once (its code is part of its key). One test
  (`test_static_statistics`).

### Data rank (4 November 2026)

- `run_pipeline.m`, `run_all.m` (copied from v_dynamic 3a7413e): `data_rank` in the
  summary row of each recording, the number of independent signals in
  the channels used (eigenvalues of their covariance, after flat channels
  are left out and before the average reference, above 1e-12 of the
  largest); a warning below `S.min_channels`; `run_all` prints the range.
  Outputs of earlier code are reused (`equivalent_code.txt`): the
  recording is read once more and the rank added to the stored row,
  nothing recomputed. A limit of 32 was planned for the warning, but 57 of
  the 138 OpenNeuro recordings (rank 9–49, median 33) are below it; 16
  flags 2 of them. One test (`test_static_stage1/test_data_rank`).

### Newcomer review (5 November 2026)

A review of v_static and v_dynamic as a first-time user and reviewer
would read them (one reviewer per stage, and a usability reading), on the
OpenNeuro data. The fixes were made in v_dynamic (commits ebdf493,
96e1e13, 5bcfa6c; every row in its CHANGES.md) and copied here
(`tools/sync_from_v_dynamic.sh`, now with the new shared files
`check_subject_table.m`, `default_group_files.m`,
`3- Analysis/fun/family_figure_name.m`, `prepare_spectral_covariates.m`,
`prune_family_figures.m`, `resolution_note.m`,
`spectral_profile_measures.m`, `+connectivity/not_computed.m`,
`+connectivity/symmetric_orthogonalisation.m`). What changes here:

- **Stage 1.** Units converted to µV (FIF data are in V: band power was
  10^12 smaller); NaN and flat stretches detected on the cap's channels
  only; `S.bad_events` (MNE's `BAD_` annotations removed; EEGLAB's -99
  boundary is a cut); the symmetric orthogonalisation (`aec_symorth`,
  `plv_symorth`) checks the data rank, not `rank()` (float32 data passed
  it, and the measures were rounding noise), and a measure that cannot be
  computed is NaN for that recording only (`measures_nan`), instead of
  losing the recording; an empty or too-short split half is NaN (it was
  the whole recording); the IAF edge check on the frequency bins (it never
  fired at 250, 500 or 1000 Hz); adaptive memory for PLI, ciPLV and
  `aec_orth`; a ciPLV zero guard; EGI `.mff` folders; duplicate subject
  IDs (BIDS sessions) stop stage 1 before anything is computed.
- **Statistics.** A feature is never adjusted for itself (a perfect fit
  gave finite p-values; IAF and `IAF_cov` count as one);
  `S.categorical_covariates` (a site coded 1/2/3 entered linearly gave
  45% false positives); a covariate constant in one group is dropped;
  the sample-size check applies the covariates actually used; the
  smallest family q the permutations can reach is printed; spectral
  covariates work in a fresh PSD analysis; console and figure titles say
  "within-family q", and figures of tests not significant after the
  adjustment over families are removed; Spearman-Brown only for r > 0;
  maxT over |t| when the columns' N differ.
- **Subject table and settings.** The table is read before stage 1
  (`check_subject_table`: recordings without a row, rows without a
  recording, sample sizes); IDs matched ignoring case, then leading zeros
  when unambiguous; numeric `S.groups` gets a clear error; omitted
  settings get the template's values (seed 20261015; `save_scouts` false;
  the four families); `S.max_workers` at most the computer's workers;
  `S.toolbox_paths` may give a toolbox's root; fetched toolboxes keep LF
  line endings (Windows).
- `pipeline_settings.m` (this folder's own): `S.bad_events`,
  `S.categorical_covariates`, the ID pattern for BIDS sessions.

The stage 1 code fingerprint changes: existing stage 1 outputs are
recomputed.

### Max-statistic family test (5 November 2026)

Copied from v_dynamic bd0bfcc (its CHANGES.md, "Max-statistic family
test"). A family (band × measure × feature family × comparison or
correlation set) is tested as a whole by its largest statistic against
the largest of each permutation (`family_p`); each test gets its
family-wise p over the family (`p_maxT_family`); families are selected by
Benjamini–Hochberg over `family_p`, and a test is significant when `p_bb`
= `p_maxT_family` × F / R ≤ `S.alpha` (before: Simes family p and `q_bb` =
`q_fdr` × F / R, whose floor m / (n_perm + 1) left the 258 families of the
OpenNeuro run unselectable). `q_fdr` stays as a description.
`S.psd.n_perm` (default `[]` = `S.n_perm`) for the PSD statistics.

| Change | Before | Now |
|---|---|---|
| shared files (statistics, PSD statistics, settings check, demo, static simulation, tests) | | from v_dynamic bd0bfcc |
| `pipeline_settings.m` | | `S.psd.n_perm`; notes on the permutations a run needs |
| `examples/settings_openneuro.m` | `S.n_perm` 10000; the subject table described as the OpenNeuro participants and phenotype files | `S.n_perm` 50000 (258 families in the statistics, 324 in the PSD); the recordings and the subject table are the lab's (not distributed; the OpenNeuro phenotype files are not used) |
| `README.md`, `docs/methods.md` | `q_bb` | the new columns and rule; how many permutations a run needs |
| `README.md` | three numberings (steps 1–5 of the table, stages, steps 1–6 with "steps 4 and 5 at once", 4b and an item 0) | one: the stages (1, 1b, 3) in the table; steps 1–5 for a study (step 4 runs everything, step 5 checks and reads each stage), as in v_dynamic |
| `docs/validation.md` | a folder of the lab's computer for the anatomy tables | "the run's tables, kept by the lab" |

### Documentation: the lab's folders (5 November 2026)

Copied from v_dynamic 5195755: `docs/psd.md` names the lab's original
PSD scripts as not distributed, without a computer's path.

### Methods on verified code (5 November 2026)

Copied from v_dynamic 6242a27 (its CHANGES.md, "PALM and the Microstate
EEGlab toolbox as references" and the sections before it). The two-group
statistics (Welch t, the Welch-type G of the ANCOVA) are compared with
PALM (Winkler et al. 2014): the statistic and its degrees of freedom are
PALM's v and its Welch df. Each statistic is now converted to z on its
Welch df (PALM's `palm_gtoz`) before the maximum over the tests is taken,
in the data and in every permutation; with unequal groups the maximum of
the raw statistics favoured the smaller df (null FWER 0.114 at 0.05 for 8
against 30 recordings; with z 0.072, PALM's own result on the same data).
The result tables gain a `df` column; the NBS edges table gives `G` and
`z`, and the NBS threshold is on z. New shared files: `3- Analysis/fun/t_to_z.m`,
`3- Analysis/fun/hedges_g.m`, `tests/test_reference_connectivity.m`,
`tests/test_reference_statistics.m` (the PALM checks are skipped when
PALM is not in `<toolboxes>/GitHub resources`).

| Change | Before | Now |
|---|---|---|
| shared files (statistics, NBS, PSD statistics, settings check, tests) | | from v_dynamic 6242a27 (9b7f75d with the dates corrected) |
| `README.md` | | "Verified implementations (5 November 2026)" |
| `docs/methods.md` | maximum over the raw t or G | maximum over z on the Welch df; NBS threshold on z; edges `G` and `z` |
| `docs/validation.md` | | the reference checks (PSD against `pwelch` and FieldTrip's cluster statistic, connectivity against MATLAB and FieldTrip, statistics against PALM) |
| `tools/sync_from_v_dynamic.sh` | | the four new shared files |

### Stage 1 comments (5 November 2026)

Copied from v_dynamic 20130e7: `bandpass_pieces.m` names
`docs/settings.md` instead of two development tools removed in
v_dynamic 3a8ad47; `equivalent_code.txt` gets the row by which the new
stage 1 code accepts the outputs of the code before (the stage 1 files of
both folders are the same, so are their fingerprints). The PSD spectra
are computed once again, with the same values (their key holds the stage
1 fingerprint without the list).

### Whole recording or windows; dwPLI against FieldTrip (6 November 2026)

Asked for: static connectivity either over the whole recording or as the
mean of windows, the windows either a number of cycles of the band or a
fixed length, each with its overlap. Shared files from v_dynamic 286b9a5
(its branch `static-windows`, merged into v_dynamic once its running
analyses end). `S.static_matrices` chooses `'full'` (the default, as
before), `'window_mean'` or both; stage 1 then computes the windows of
v_dynamic (phase measures: `S.window_cycles` cycles of the band centre or
`S.window_length` s, `S.window_overlap`; AEC: `S.aec_window` s every
`S.aec_step` s), saves their mean and its split halves (not the windows),
and the statistics analyse each matrix on its own, with its reliability.
Outputs made before are still valid for `'full'` (the same values;
`equivalent_code.txt`).

| Change | Before | Now |
|---|---|---|
| shared files (stage 1, settings check and key, static analyses, tests) | from v_dynamic 504f25b | from v_dynamic 286b9a5: the window mean in a static-only run; dwPLI tested against FieldTrip (`ft_connectivity_wpli` with `debias`, and FieldTrip's dwPLI of 2 s epochs) |
| `pipeline_settings.m` | the connectivity over the whole recording only; no window settings | `S.static_matrices = {'full'}` with what each choice means, and the window settings (`window_cycles` 6, `window_length` `[]`, `window_overlap` 0.5, `aec_window` 6, `aec_step` 0.5) |
| `README.md` | "over the whole recording" | the two estimates, their settings and outputs |
| `docs/methods.md` | | step 5 (the window mean) and "Whole recording or windows": what each estimates, the demo data (planted within-DMN coupling: wPLI over the whole recording g 2.17; window mean of 0.6 s, 2 s, 6 s windows g 0.42, 0.93, 1.58), the short-window bias of phase measures (OpenNeuro alpha wPLI 0.15 over the whole recording, 0.55 as the mean of 0.6 s windows; its correlation with band power 0.48 and 0.67), AEC alike in both, bursts, windows longer than the shortest piece |
| `docs/validation.md` | | the FieldTrip checks of dwPLI and the tests of the window mean |
| `docs/references.md` | | O'Neill et al. 2017 (AEC windows) |

### The window mean by default; windows of one measure (6 November 2026)

Asked for: the window mean as the default, the whole recording still
selectable, and window settings that may differ between measures.
Several window settings for one measure in one run were not asked for
(separate runs, each with its own `S.results`). Shared files from
v_dynamic ae3cc80 (branch `static-windows`). Outputs made before with
`'full'` stay valid for `'full'`; a default run now computes the windows.

| Change | Before | Now |
|---|---|---|
| shared files (stage 1, settings check and key, tests) | from v_dynamic 286b9a5 | from v_dynamic ae3cc80: `S.static_matrices = {}` is the window mean; `S.measure_windows` (each measure's own windows, checked, in the stage 1 key only when given); the `window_vs_segment` warning reads each measure's windows and names the measures |
| `pipeline_settings.m` | `S.static_matrices = {'full'}`; `S.min_segment = 4` | `S.static_matrices = {'window_mean'}`; `S.measure_windows = struct()` with an example; `S.min_segment = 6`, so that the 6 s AEC windows fit in every piece kept and a default run gives no warning (user decision; on the OpenNeuro recordings, single pieces of about 5 minutes, 4 and 6 keep the same data) |
| `tests/test_template.m` | the default `'full'` | the default `'window_mean'`, `S.min_segment` 6, no `S.measure_windows` |
| `README.md`, `docs/methods.md`, `docs/validation.md` | the whole recording as the default | the window mean as the default, its output names, `S.measure_windows`, `S.min_segment` 6 and the piece lengths it needs (8 s, 9.6 s in delta and theta) |
