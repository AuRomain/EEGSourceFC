#!/usr/bin/env bash
# Copy into this folder (v_static) the files it shares with the full
# pipeline, v_dynamic, unchanged, and record their SHA-256 in
# tools/v_dynamic_files.txt (tests/test_same_code_as_v_dynamic.m checks them).
#
# Usage (from v_static/; bash, git):
#   bash tools/sync_from_v_dynamic.sh [path of v_dynamic]   (default ../v_dynamic)
#
# Shared: the stage 1 folder as v_dynamic's git tracks it (code, inputs/; not the
# original 61-channel and ANT caps, which only v_dynamic's tools/make_corrected_inputs.m
# reads: the original 61-channel positions are rotated, error E1), the PSD
# step ("1b- PSD analysis", as v_dynamic's git tracks it), the
# static analyses of "3- Analysis" (run_statistics runs only them when
# S.dynamic = false), the top-level helpers, the static simulation, the
# static tests, the anatomy test, four tools and the Brainstorm fix. Not shared, and written for this folder:
# pipeline_settings.m (the template: S.dynamic = false, six bands, wPLI and
# orthogonalised AEC, every feature family), examples/, README.md,
# CHANGES.md, LICENSE, LICENSE.md, .github/, tests/test_same_code_as_v_dynamic.m,
# tests/test_template.m and this script. A file of the previous list that
# is no longer shared is removed. Change the shared files in v_dynamic, test them
# there, then run this script and the tests here.
set -euo pipefail
V_DYNAMIC="${1:-../v_dynamic}"
if [ ! -f "$V_DYNAMIC/run_static_demo.m" ] || [ ! -d "$V_DYNAMIC/2- Dynamic connectivity" ]; then
  echo "$V_DYNAMIC is not v_dynamic" >&2; exit 1
fi
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

LIST=tools/v_dynamic_files.txt
NEW="$(mktemp)"
{
  git -C "$V_DYNAMIC" ls-files "1- EEG source connectivity" | grep -v -e '/inputs/channel_EEGLab_61_bs\.mat$' -e '/inputs/channel_ANT_Waveguard_65\.mat$'
  git -C "$V_DYNAMIC" ls-files "1b- PSD analysis"
  cat <<'FILES'
3- Analysis/run_statistics.m
3- Analysis/static_analysis.m
3- Analysis/compute_NBS_analysis.m
3- Analysis/spectral_profile_features.m
3- Analysis/DK68_networks.csv
3- Analysis/fun/atlas_regions.m
3- Analysis/fun/check_sample_sizes.m
3- Analysis/fun/circular_graph.m
3- Analysis/fun/cluster_perm_test.m
3- Analysis/fun/covariate_matrix.m
3- Analysis/fun/family_adjust.m
3- Analysis/fun/family_figure_name.m
3- Analysis/fun/fdr_bh.m
3- Analysis/fun/glass_brain.m
3- Analysis/fun/hedges_g.m
3- Analysis/fun/group_matrices.m
3- Analysis/fun/load_network_table.m
3- Analysis/fun/load_static.m
3- Analysis/fun/network_assignment.m
3- Analysis/fun/network_colours.m
3- Analysis/fun/network_enrichment.m
3- Analysis/fun/perm_ancova_maxT.m
3- Analysis/fun/perm_partial_corr_maxT.m
3- Analysis/fun/permu_ttest_permu_corr_BNS.m
3- Analysis/fun/prepare_spectral_covariates.m
3- Analysis/fun/prune_family_figures.m
3- Analysis/fun/resolution_note.m
3- Analysis/fun/select_subjects.m
3- Analysis/fun/spectral_profile_summary.m
3- Analysis/fun/spectral_profile_measures.m
3- Analysis/fun/spin_permutations.m
3- Analysis/fun/static_features.m
3- Analysis/fun/static_reliability.m
3- Analysis/fun/t_to_z.m
3- Analysis/fun/usable_covariates.m
3- Analysis/fun/welch_G.m
3- Analysis/fun/NBS/NBSstats.m
3- Analysis/fun/NBS/get_components.m
3- Analysis/fun/NBS/README.md
analysis_dir.m
band_edges.m
bns_run_dir.m
cap_channels.m
check_settings.m
check_subject_table.m
code_fingerprint.m
default_group_files.m
file_fingerprint.m
file_tag.m
get_subject_groups.m
gt_window_metrics.m
hash_text.m
read_subject_table.m
run_full_pipeline.m
run_static_demo.m
settings_key.m
setup_toolboxes.m
stage1_outputs.m
stage2_settings.m
stage_settings.m
subject_ids.m
toolbox_check.m
toolbox_manifest.m
toolbox_reference.json
docs/anatomy.md
docs/individual_alpha_frequency.md
docs/psd.md
simulation/run_static_simulation.m
simulation/+sim/band_noise.m
simulation/+sim/forward_model.m
simulation/+sim/peak_noise.m
simulation/+sim/pink_noise.m
simulation/+sim/roi_patterns.m
simulation/+sim/simulate_recording.m
simulation/+sim/simulate_static_dataset.m
simulation/+sim/static_params.m
simulation/+sim/write_eeg.m
tests/static_test_recordings.m
tests/test_anatomy.m
tests/test_psd.m
tests/test_reference_connectivity.m
tests/test_reference_statistics.m
tests/static_test_settings.m
tests/test_static_estimators.m
tests/test_static_ground_truth.m
tests/test_static_integration.m
tests/test_static_settings.m
tests/test_static_stage1.m
tests/test_static_statistics.m
tests/python/test_build_leadfield.py
tools/brainstorm_base64_952475b.patch
tools/build_leadfield_openmeeg.py
tools/compare_with_brainstorm.m
tools/fetch_toolboxes.sh
tools/make_anatomy.m
FILES
} > "$NEW"

while IFS= read -r f; do               # every file there before anything is removed or copied
  [ -f "$V_DYNAMIC/$f" ] || { echo "missing in $V_DYNAMIC: $f" >&2; rm "$NEW"; exit 1; }
done < "$NEW"

# remove the files of the previous list that are no longer shared
if [ -f "$LIST" ]; then
  grep -v '^#' "$LIST" | cut -d' ' -f3- | while IFS= read -r f; do
    if ! grep -qxF "$f" "$NEW" && [ -f "$f" ]; then echo "removed (no longer shared): $f"; rm "$f"; fi
  done
fi

REV="$(git -C "$V_DYNAMIC" rev-parse --short HEAD)"
DIRTY=""
if [ -n "$(git -C "$V_DYNAMIC" status --porcelain)" ]; then DIRTY=" (v_dynamic had uncommitted changes)"; fi
TMP="$(mktemp)"
{
  echo "# Files shared with v_dynamic, copied by tools/sync_from_v_dynamic.sh from"
  echo "# its commit $REV$DIRTY on $(date +%Y-%m-%d). SHA-256 and path; do not edit."
} > "$TMP"
n=0
while IFS= read -r f; do
  mkdir -p "$(dirname "$f")"
  cp -p "$V_DYNAMIC/$f" "$f"
  echo "$(sha "$f")  $f" >> "$TMP"
  n=$((n + 1))
done < "$NEW"
mv "$TMP" "$LIST"                       # the list is replaced only when every file was copied
rm "$NEW"
echo "$n files copied from $V_DYNAMIC ($REV$DIRTY); list and SHA-256 in $LIST"
