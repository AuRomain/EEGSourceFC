# Another head model: templates and individual heads

By default stage 1 uses one head for every recording: Brainstorm's ICBM152
(scalp, outer and inner skull, the 15,002-vertex `tess_cortex_pial_low`
cortex) and the caps placed on it, all in `1- EEG source connectivity/inputs/`.
Three settings choose another head:

| Setting | Default | Meaning |
|---|---|---|
| `S.anatomy` | `''` | `''`: the shipped ICBM152. A folder: another head for every recording (a template such as Colin27). With `S.anatomy_individual`: the folder of the subjects' folders |
| `S.anatomy_individual` | `false` | `true`: each recording uses `<S.anatomy>/<ID>/`, ID = its subject ID (`S.id_regexp`, as in the subject table) |
| `S.anatomy_fallback` | `false` | with individual heads: `true` lets a subject without a folder use the shipped ICBM152 (with a warning naming them); `false` stops the run before anything is computed, naming every recording without a head |

The folders are made by `tools/make_anatomy.m`, a short script that calls
Brainstorm's own processes (below). Nothing else changes: the same
inverse methods, scout functions, measures and outputs; the leadfield of
each head is computed once, in its folder, with the same BEM as the
shipped one (OpenMEEG, conductivities brain/skull/scalp 1 / 0.0125 / 1).

## 1. Make the folders

```matlab
S = pipeline_settings();       % or your settings: the caps of S.montages are placed on each head

% another template for every recording (Brainstorm's list: 'ICBM152' is part
% of Brainstorm; 'Colin27_2016', 'ICBM152_2023b', 'FsAverage_2020', ... are
% downloaded once, about 20-80 MB)
make_anatomy(S, 'template', 'Colin27_2016', '/data/anatomy/Colin27_2016');

% individual heads from T1 MRIs (CAT12, about 30 minutes a subject, and 5 for the leadfield)
subjects = {'sub-01', '/data/mri/sub-01_T1w.nii.gz'
            'sub-02', '/data/mri/sub-02_T1w.nii.gz'};
make_anatomy(S, 't1', subjects, '/data/anatomy/individual');

% CAT12's server (neuro.uni-jena.de) cannot be reached: a CAT12 folder of your
% own, e.g. a release from github.com/ChristianGaser/cat12 (the folder of
% cat_version.m), anywhere but in Brainstorm's plugins folder
make_anatomy(S, 't1', subjects, '/data/anatomy/individual', 'cat12', '/data/toolboxes/cat12');

% individual heads from segmentations already made (CAT12 or FreeSurfer folders)
make_anatomy(S, 'segmentation', {'sub-03', '/data/fs/sub-03'}, '/data/anatomy/individual');

% a cap added to S.montages later: placed on existing folders
make_anatomy(S, 'caps', {'/data/anatomy/individual/sub-01', '/data/anatomy/individual/sub-02'});
```

`tools/` is put on the path by `addpath(fullfile(S.root, 'tools'))`. A
folder already complete is skipped (delete it to make it again). The IDs
must be those stage 1 reads from the file names (`S.id_regexp`): a
recording `sub-01_task-rest_eeg.set` with `S.id_regexp = '^([^_]+)'` looks
for `<S.anatomy>/sub-01/`.

What Brainstorm does (its processes, as in its tutorials' scripts):

| | Brainstorm |
|---|---|
| Template | `db_set_template` (its "Use template" menu; downloaded the first time); `process_generate_bem` only if the template has no skull surfaces |
| T1 MRI | `process_import_mri`; `process_segment_cat12` (CAT12: cortex of 15,000 vertices, the atlases Desikan-Killiany, Destrieux, DKT, ... with their sphere registration, non-linear MNI normalization; NAS/LPA/RPA from the MNI positions when none are given, Brainstorm's default in scripts); `process_generate_bem` (scalp, outer and inner skull of 1922 vertices, skull 4 mm) |
| Segmentation | `process_import_anatomy` (CAT12 or FreeSurfer folder); `process_mni_normalize` if it gave no MNI normalization; `process_generate_bem` |
| Caps | the pipeline's positions on the ICBM152 (`get_ftElec`) taken to the head through the two MNI normalizations (`cs_convert`, as "Add EEG positions" does with MNI coordinates in `channel_add_loc`), then projected on its scalp (`channel_project_scalp`, "Project electrodes on scalp"). With CAT12's non-linear normalization, its linear part: the deformation is fitted to the brain, and at the scalp it moved electrodes near the vertex by up to 3 cm (Colin27's T1) |

SPM12 and CAT12 are installed as Brainstorm plugins the first time (`'t1'`);
`'cat12', folder` gives Brainstorm a CAT12 folder of your own instead (its
custom plugin path, linked into SPM12's toolbox folder as Brainstorm links
the CAT12 it installs; not inside Brainstorm's plugins folder, where
Brainstorm takes any folder whose path starts with `plugins/cat12` for its
own). CAT12's surfaces are GIfTI files, which Brainstorm 3.250107 reads
with a Java class that Java 9 and later no longer have (MATLAB's Java on
Apple silicon is 11): the pinned Brainstorm needs its own fix of 13 October
2026 (commit 952475b), which `tools/fetch_toolboxes.sh` applies
(`tools/brainstorm_base64_952475b.patch`); without it `make_anatomy` stops
before CAT12 runs, saying so. The errors and warnings of Brainstorm's processes,
which it keeps in its report when called from a script, are printed, and an
error stops `make_anatomy` with Brainstorm's message.
Brainstorm runs headless in the same MATLAB (quit its interface first); its
user folder (plugins, downloaded templates) is
`<S.toolboxes>/brainstorm_user/.brainstorm`, kept for the next runs, never
`~/.brainstorm`; its database is a protocol deleted at the end; its folders
leave the MATLAB path when it stops. Then the leadfield of each cap is
computed in the folder (`make_anatomy(..., 'leadfield', false)` leaves it
to stage 1), with the builder's checks: closed and nested surfaces,
electrodes on the scalp, cap orientation, and the sign and localisation of
the gain. A folder already made is not made again (delete it to redo it),
but its missing leadfields are computed.

OpenMEEG's integrals are not finite for some exactly degenerate meshes:
the skulls of Brainstorm's current ICBM152 (1922 vertices each) have
two faces each with an exactly coplanar neighbour, and gave a head matrix
with 24 non-finite values. The
builder then computes the gain again with the vertices as Brainstorm writes
them for OpenMEEG (millimetres, 6 significant digits: `out_tess_tri`),
which removes them; the leadfield file's `bem_vertices` says which was used.
The pipeline's shipped head (642-vertex skulls) is not affected.

The vertex normals are those of Brainstorm's interface (MATLAB patch
normals, `tess_normals`): headless, Brainstorm sums the faces' normals
instead, which on the ICBM152 cortex differs from the interface's by a
median 5 degrees (up to 90, 23 vertices reversed). The sources are oriented
along these normals, so `make_anatomy` computes them as the interface does
(identical to the shipped cortex's).

## 2. What a folder holds

| File | |
|---|---|
| `tess_head.mat`, `tess_outerskull.mat`, `tess_innerskull.mat` | the BEM surfaces |
| `tess_cortex.mat` | the sources (its vertices and normals) and its atlases (`Atlas`), with the sphere registration (`Reg`) |
| `channel_<cap>.mat` | each cap of `S.montages` on this head (`<cap>` = `ASA_10-05_343`, `EEGLAB_61`, or the name of a cap file) |
| `subjectimage_T1.mat` | the MRI and its MNI normalization (for `make_anatomy(S, 'caps', ...)`) |
| `anatomy_info.txt` | what was done (source, Brainstorm version, MNI normalization, surfaces, atlases), and what to look at: scalp vertices standing out from their neighbours by more than 8 mm (a spike or fold of the surface; also a warning), how far the projection on the scalp moved the electrodes, naming those moved by more than 10 mm |
| `ftLeadfield_<cap>_normal.mat` | the leadfield, computed once (with `anatomy_key`, the fingerprint of the surfaces and cap it was made from: it is made again when they change) |

All in Brainstorm's subject coordinates (SCS, metres). Any folder with
these files works, e.g. a subject exported from a Brainstorm protocol.
Electrode positions digitised on a subject's head are not placed by
`make_anatomy`: import and align them in Brainstorm, then save the channel
file in the folder as `channel_<cap>.mat`, with the cap's labels in its
order (stage 1 stops if they differ).

## 3. Atlases on other heads

`S.atlases` can name an atlas of the cortex: `'Desikan_Killiany_68'` (the
cortex's `Desikan-Killiany`, which Brainstorm makes from FreeSurfer's or
CAT12's `aparc` / `aparc_DK40`), `'Destrieux'`, `'DKT'`,
`'Schaefer_100_17net'`, `'Brainnetome'`, `'HCP_MMP1'`, ... (the shipped
cortex has 22). On every head the regions are put in the order of the
shipped ICBM152's atlas of the same name, so the matrices of all
recordings have the same rows: a region missing on a head stops the run,
regions the shipped atlas does not have (`unknown`, ...) are left out.
Scout files (`S.atlases` = a file) are defined on one cortex: they work on
the shipped head or on the template whose cortex they were made on, not
with individual heads (`check_settings` stops). For the spin test and the
brain figures of stage 3, a template folder gives its cortex; individual
heads use the shipped ICBM152's (the same regions).

## 4. In the outputs

`stage1_summary.csv` has a column `anatomy`: `ICBM152 (shipped)`, the
folder's name (the template, or the subject ID), or `ICBM152 (shipped,
fallback)`. Each output's `info.anatomy` holds the name, folder and key of
the head. The head is part of each output's input key: outputs made with
another head are recomputed, and those made with the shipped head keep
their keys (a fallback subject reuses them).

## 5. Limits

- Individual heads need T1 MRIs of the subjects, and their segmentation
  should be checked in Brainstorm (the cortex inside the inner skull, the
  skull surfaces); `make_anatomy` does not show them.
- The caps are template positions warped through the MNI normalization,
  not digitised positions: the electrodes follow the head's shape, not
  where the cap was on that subject.
- The leadfield of each head takes a few minutes (OpenMEEG); stage 1
  computes them one at a time before the recordings.
- CAT12's surfaces can have defects. On Colin27's T1 its scalp had two
  bumps on top of the head (up to 39 mm), and the four electrodes on them
  were moved 12 to 34 mm by the projection: `anatomy_info.txt` names both.
  Check such a head in Brainstorm and remake or correct the surface (or
  give a FreeSurfer segmentation with `'segmentation'`).
- The cortex of a CAT12 head is its central surface (between the white
  matter and the pial surface), the one Brainstorm makes the default; the
  shipped ICBM152's is its pial surface (`tess_cortex_pial_low`).
- A CAT12 folder takes about 105 MB: the MRI with CAT12's deformation
  fields (89 MB), the leadfield of the 343-electrode cap (20 MB); a
  template folder about 25 MB.
- Timings on the development Mac (arm64, 28 cores): a Brainstorm template
  1.4-9 minutes with its leadfield; CAT12 on one T1 29 minutes, the caps
  and leadfield 5 more.
- How much the head changes the results: `docs/validation.md`, "Other head
  models".
