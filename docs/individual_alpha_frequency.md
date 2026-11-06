# Individual (IAF-anchored) frequency bands

Optional; off by default. Implemented 20 October 2026.

## Switching it on

In `pipeline_settings.m`, section **"Individual bands (optional): ON/OFF HERE"**:

```matlab
S.individual_bands = false;   % false: the fixed bands of S.bands; true: S.bands_individual
```

- **`false` (default):** the pipeline uses the fixed bands of `S.bands`, exactly as
  before. No IAF is estimated and restingIAF is not needed.
- **`true`:** `S.bands_individual` replaces `S.bands`. Each recording's IAF is
  estimated in stage 1, each band's edges are computed from it, and the IAF is analysed
  in stage 3.

An edge is either a number (Hz) or text relative to the IAF: `'iaf'`, `'iaf+2'`,
`'iaf-6'`, `'iaf+bins'`, `'iaf-bins'` (`band_edges.m`). Fixed and anchored bands can be
mixed. The default set gives the classical bands when IAF = 10 Hz:

| Band | Lower edge | Upper edge | At IAF = 10 Hz | At IAF = 9 Hz |
|---|---|---|---|---|
| delta_iaf | 1 Hz | IAF − 6 Hz | 1–4 | 1–3 |
| theta_iaf | IAF − 6 Hz | IAF − bins·bin_hz | 4–8 | 3–7 |
| alpha_iaf | IAF − bins·bin_hz | IAF + bins·bin_hz | 8–12 | 7–11 |
| beta1_iaf | IAF + bins·bin_hz | 20 Hz | 12–20 | 11–20 |
| beta2 | 20 Hz | 30 Hz | 20–30 | 20–30 |
| gamma | 30 Hz | 45 Hz | 30–45 | 30–45 |

- **`S.iaf.bins`** (default 2) and **`S.iaf.bin_hz`** (default 1 Hz) set the alpha band:
  IAF ± 2 bins of 1 Hz, so an IAF of 9 Hz gives 7–11 Hz. `bins = 3` gives a band as wide
  as 7–13 Hz; the default stays 2 after a comparison of the two on the OpenNeuro data
  (`docs/settings.md` of the full pipeline, v_dynamic: "Width of the IAF band",
  30 October 2026).
- **The upper edges of beta and gamma stay fixed:** beta peaks do not follow alpha
  harmonically (Haegens et al., 2014).
- **Narrow bands:** a band narrower than 1 Hz for a recording's IAF is left out for that
  recording (delta, 1 Hz to IAF − 6, for an IAF below 8 Hz); narrower than 1.5 Hz gives a
  warning. A band left out (no IAF, or too narrow) counts as done when the recording's IAF
  is stored with the same file and `S.iaf` (`iaf_<cap>_<atlas>.mat`): a rerun skips the
  recording (a band too narrow made it recompute at every run, until 5 November 2026).

## What happens

**Stage 1** (`run_pipeline.m`, `individual_alpha_frequency.m`), once per recording:

- **Channels:** the posterior channels of the cap (`S.iaf.channels`; on a cap file without
  these names, its electrodes nearest to their 10-05 positions, `cap_channels.m`), average-referenced,
  in the continuous pieces of at least 4 s.
- **IAF:** restingIAF (Corcoran et al., 2018), with its published defaults (search
  7–13 Hz, Savitzky–Golay frame 11 and order 5, at least 3 channels).
  - It gives the peak alpha frequency (PAF) and the centre of gravity (CoG).
  - The IAF used is `S.iaf.measure` (`'paf'` by default).
- **Checks** (`S.iaf.checks`), reported:
  - the IRASA oscillatory peak (Wen & Liu, 2016);
  - the specparam/FOOOF peak (Donoghue et al., 2020).

  Both are computed with FieldTrip's `ft_freqanalysis`; for specparam, FieldTrip calls
  Brainstorm's MATLAB implementation.
- **No estimate** (fewer than 3 channels with a peak): the recording is left out of the
  anchored bands only (`S.iaf.no_peak = 'exclude'`), or its CoG is used (`'cog'`). Old
  outputs of those bands are deleted, and the fixed bands are computed as usual.
- **Edge check** (`S.iaf.edge_check`): a "peak" within 0.5 Hz of the search range's edge,
  where the IRASA spectrum has no maximum inside the range, is not taken as the IAF. It is
  the flank of the beta or theta rhythm next to the range (see Validation).
  - Decided channel by channel on the frequency bins of the IRASA spectrum: a channel has
    no maximum inside the range when its largest value in the range is on the range's
    first (last) bin. The IAF is rejected when it lies within 0.5 Hz of that edge and at
    least half of the posterior channels have their maximum on that edge's bin
    (`R.irasa_edge`: that share).
  - Until 5 November 2026 the check compared the channels' median IRASA peak with 7 and
    13 Hz within 0.01 Hz. The bins of the 4 s pieces (padded to a power of 2) are
    0.244 Hz apart at 250, 500 and 1000 Hz, so the in-range edge bins are 7.08 and
    12.94 Hz and the check never fired; at 256 or 512 Hz (bins 0.25 Hz apart) it did not
    fire when exactly half of an even number of channels were on the edge bin (the
    median then halfway to the next bin). A simulated recording without alpha, resampled
    to 250, 500, 512 and 1000 Hz, with 10 or 11 posterior channels: its "peak" at
    12.85–12.88 Hz was rejected only at 512 Hz before, at every rate and channel count
    now. Of the 138 OpenNeuro recordings (1000 Hz), three have a PAF within 0.5 Hz of an
    edge; one IAF changes: a PAF of 7.19 Hz whose IRASA spectrum has its maximum on the
    7.08 Hz bin in 6 of the 10 posterior channels is now rejected (9 recordings without an
    IAF instead of 8); the two near 13 Hz (on the edge bin in 4 and 0 of 10) are kept.
- **Flags:** an IAF more than 0.5 Hz from either check's peak is flagged
  (`S.iaf.max_disagreement`). It is kept, but should be looked at.
- **Each band:** its edges from the IAF, the filter from the window rule, and windows of
  6 cycles of the band's centre, so the same number of cycles for every recording.
- **Saved:**
  - `stage1_summary.csv`: columns `iaf`, `iaf_paf`, `iaf_cog`, `iaf_irasa`,
    `iaf_specparam`, `iaf_channels`, `iaf_source`, `iaf_flag`, `bands_skipped`;
  - `Connectivity/iaf.csv`: one row per recording, with the subject ID;
  - each output's `info.band_hz` and `info.iaf`.

**Stage 2:**
- Window steps differ between recordings by design.
- `results.band_hz` and `results.iaf` give each recording's edges and IAF.

**Stage 3:**
- The IAF is compared between groups and correlated with the clinical scores, as a
  feature without states (`Statistical analysis/IAF`).
- The IAF columns are joined to the subject table (`read_subject_table`), so the IAF can
  also be a covariate of the other families: `S.covariates = {'age', 'IAF'}`. It is not a
  covariate of itself: the IAF family (and `IAF` of the family `spectral`) is then tested
  adjusted for age only (IAF, `IAF_cov` and the `IAF_PAF` or `IAF_CoG` that the IAF is
  count as one measure). Before 5 November 2026 the IAF was adjusted for itself there: a
  perfect fit, whose rounding-level residuals gave a meaningless p-value.

## Toolboxes

| Role | Toolbox | Where |
|---|---|---|
| IAF (PAF, CoG) | restingIAF (Corcoran et al., 2018), GPL-3 | github.com/corcorana/restingIAF, in `../Toolboxes/restingIAF` at commit b25db38; needed only when the switch is on |
| IRASA check | FieldTrip `ft_freqanalysis`, `cfg.method = 'irasa'` | `../Toolboxes/fieldtrip-toolbox` |
| specparam check | FieldTrip `ft_freqanalysis`, `cfg.output = 'fooof_peaks'` (Brainstorm's `process_fooof`) | `../Toolboxes/fieldtrip-toolbox` |

## Validation

**Simulation** (`simulation/`, `p.alpha_rms`): an alpha rhythm on the visual network with
a known IAF per subject (group means 10.4 and 9.7 Hz, as in Murphy & Öngür, 2019):

| Simulated alpha | restingIAF | IRASA | specparam |
|---|---|---|---|
| Clear (amplitude 2) | within 0.25 Hz of the truth in 4 of 4 | 4 of 4 | 4 of 4 |
| Weak (amplitude 1), next to strong beta at 13–20 Hz | 0 of 4 (the beta flank, 12.3–13.2 Hz) | 0 of 4 | 3 of 4 |
| None | 0 of 4 after the edge check (before it, a "peak" at 12.6–13.1 Hz) | — | — |

No method is safe in every case:
- restingIAF and IRASA look only within the search range, so a strong neighbouring rhythm
  can mislead them.
- specparam models separate peaks, but in the 8 development recordings it differed from
  the other two by up to 0.6 Hz.

Hence restingIAF (the reference method) sets the IAF, and disagreements are flagged.

**Development recordings** (eight, eyes open; a feasibility check run during
development): restingIAF peak 9.8–11.2 Hz (centre of gravity 9.2–10.9 Hz), from 4 to 11 of
the 11 posterior channels; specparam 9.5–11.8 Hz; IRASA 9.8–11.3 Hz. (Ranges only: the
values of each recording are not shared.)

**Tests:**
- `test_stage1/test_band_edges_and_individual_switch`: edge parsing, and the switch.
- `test_stage1/test_iaf_edge_check_at_any_rate`: the simulated recording without alpha is
  rejected at 250, 500, 512 and 1000 Hz with 10 and 11 posterior channels, and one with
  alpha is kept.
- `test_stage1/test_band_too_narrow_counts_as_done`: a band left out as too narrow is not
  recomputed at the next run.
- `test_simulation/test_individual_bands_in_stage1`:
  - the IAF is recovered within 0.3 Hz;
  - the anchored edges and windows are as specified;
  - a recording without alpha is left out of the anchored band only, and its old output
    is removed;
  - a second run is up to date.

**Still to do:**
- a split-half reliability of the IAF on the full sample (cf. Grandy et al., 2013);
- a simulation of whether anchored bands recover planted alpha connectivity better than
  fixed bands when the groups differ in IAF.

## Why individual bands

- **The alpha peak frequency differs between people.** Between subjects the SD was 2.8 Hz
  in an MEG source study; within a subject it varied by 0.9 Hz across conditions and
  regions (Haegens et al., 2014).
- **It is a stable trait** (Grandy et al., 2013).
- **It is lower in first-episode psychosis and in schizophrenia** (Murphy & Öngür, 2019;
  Catalano et al., 2024; Sponheim et al., 2024).
- **So fixed bands do not cover the same rhythm in patients and in controls.** Bands anchored on
  each subject's IAF avoid this (Doppelmayr et al., 1998; Klimesch, 1999).

## References (from PubMed)

- Klimesch W (1999). EEG alpha and theta oscillations reflect cognitive and memory performance: a review and analysis. *Brain Res Brain Res Rev* 29:169–195. [doi:10.1016/s0165-0173(98)00056-3](https://doi.org/10.1016/s0165-0173(98)00056-3)
- Doppelmayr M, Klimesch W, Pachinger T, Ripper B (1998). Individual differences in brain dynamics: important implications for the calculation of event-related band power. *Biol Cybern* 79:49–57. [doi:10.1007/s004220050457](https://doi.org/10.1007/s004220050457)
- Corcoran AW, Alday PM, Schlesewsky M, Bornkessel-Schlesewsky I (2018). Toward a reliable, automated method of individual alpha frequency (IAF) quantification. *Psychophysiology* 55:e13064. [doi:10.1111/psyp.13064](https://doi.org/10.1111/psyp.13064)
- Donoghue T, Haller M, Peterson EJ, et al. (2020). Parameterizing neural power spectra into periodic and aperiodic components. *Nat Neurosci* 23:1655–1665. [doi:10.1038/s41593-020-00744-x](https://doi.org/10.1038/s41593-020-00744-x)
- Wen H, Liu Z (2016). Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. *Brain Topogr* 29:13–26. [doi:10.1007/s10548-015-0448-0](https://doi.org/10.1007/s10548-015-0448-0)
- Haegens S, Cousijn H, Wallis G, Harrison PJ, Nobre AC (2014). Inter- and intra-individual variability in alpha peak frequency. *NeuroImage* 92:46–55. [doi:10.1016/j.neuroimage.2014.01.049](https://doi.org/10.1016/j.neuroimage.2014.01.049)
- Grandy TH, Werkle-Bergner M, Chicherio C, et al. (2013). Peak individual alpha frequency qualifies as a stable neurophysiological trait marker in healthy younger and older adults. *Psychophysiology* 50:570–582. [doi:10.1111/psyp.12043](https://doi.org/10.1111/psyp.12043)
- Murphy M, Öngür D (2019). Decreased peak alpha frequency and impaired visual evoked potentials in first episode psychosis. *NeuroImage Clin* 22:101693. [doi:10.1016/j.nicl.2019.101693](https://doi.org/10.1016/j.nicl.2019.101693)
- Catalano LT, Reavis EA, Wynn JK, Green MF (2024). Peak alpha frequency in schizophrenia, bipolar disorder, and healthy volunteers: associations with visual information processing and cognition. *Biol Psychiatry Cogn Neurosci Neuroimaging* 9:1132–1140. [doi:10.1016/j.bpsc.2024.06.004](https://doi.org/10.1016/j.bpsc.2024.06.004)
- Sponheim SR, Ramsay IS, Lynn PA, Vinogradov S (2024). Generalized slowing of resting-state neural oscillations in people with schizophrenia. *Biol Psychiatry Cogn Neurosci Neuroimaging* 10:693–698. [doi:10.1016/j.bpsc.2024.08.007](https://doi.org/10.1016/j.bpsc.2024.08.007)
