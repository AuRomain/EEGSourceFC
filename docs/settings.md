# Settings

Every setting of EEGSourceFC is a field of `S`, listed with its default
and its meaning in `pipeline_settings.m`; the main ones are in the README,
section 5. The code and its messages point here for these topics:

- **Caps** (a cap given as a file, its leadfield, the reference site):
  README, section 4, step 1 (the rows Channels and Caps).
- **The table of the filters** (each band's transition, filter length and
  transient) and the filter's transition rules (`S.filter_transition`,
  `S.filter_span`): `docs/methods.md`, "Filters".
- **Windows** (`S.window_cycles`, `S.window_length`, `S.window_overlap`,
  `S.measure_windows`) and the matrices (`S.static_matrices`):
  `docs/methods.md`, "Whole recording or windows".
- **The PSD analysis** (`S.psd`): `docs/psd.md`.
- **Individual bands** (`S.individual_bands`, `S.bands_individual`,
  `S.iaf`): `docs/individual_alpha_frequency.md`.
