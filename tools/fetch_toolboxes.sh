#!/usr/bin/env bash
# Fetch the toolboxes of the pipeline at the versions validated with it
# (toolbox_manifest.m), into the folder setup_toolboxes looks in.
#
# Usage (from the pipeline folder, v_dynamic/ or v_static/; bash,
# git, curl, unzip and patch: on Windows, Git Bash):
#   bash tools/fetch_toolboxes.sh [destination] [--static] [--dry-run]
#   destination: default ../Toolboxes (S.toolboxes); --static: only the
#   toolboxes of a static-only run (S.dynamic = false: FieldTrip, Brainstorm,
#   restingIAF); --dry-run prints the commands only.
#
# Pinned: FieldTrip commit 0beabe6 (3 April 2023: the validated copy, found
# by content; about 400 MB to download, which can take half an hour),
# Brainstorm 3.250107 (commit 69bbf53; doc/version.txt gets the line that
# Brainstorm's zip adds, so the copy matches the validated one file by
# file; tools/brainstorm_base64_952475b.patch adds Brainstorm's own fix of
# 13 October 2026, commit 952475b, without which it cannot read the GIfTI
# surfaces of CAT12 with Java 9 or later: tools/make_anatomy.m), restingIAF (commit b25db38; only needed for individual bands and
# the spectral profile), and, with states (not with --static), BCT release
# 2019-03-03 (the zip of its website, checked by content), modkmeans.m and
# calc_fitmeas.m of the Microstate EEGlab toolbox (folder MST1.0 of commit
# a46ddd9) and CVIK (commit 22382ec). The last two are the published
# versions of copies that the lab had edited: toolbox_reference.json lists
# them as equivalent (same results; cvi_index.m bridges CVIK's two
# signatures). NBS 1.2 is not
# fetched: the two functions used are included in "3- Analysis/fun/NBS".
# Each git toolbox is fetched as one commit (no history). Afterwards,
# setup_toolboxes compares every toolbox with the validated files
# (toolbox_reference.json) and names any that differ: run in MATLAB
#   S = pipeline_settings(); S.toolboxes = '<destination>';
#   S = setup_toolboxes(S); disp(S.versions.toolboxes(:, {'toolbox', 'status'}))
set -euo pipefail

DEST="../Toolboxes"
DRY=0
STATIC=0
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --static) STATIC=1 ;;
    *) DEST="$a" ;;
  esac
done
FIELDTRIP_REV="${FIELDTRIP_REV:-0beabe62dd34bc64d3ccbdb0ebbb0d86d2ed7aff}"
BRAINSTORM_REV=69bbf5335fd967d89de90428fb22a1855526ce14
MICROSTATE_REV=a46ddd9459d2e2f8a2a7064841dcacb45766b329
CVIK_REV=22382ec1dd8e12025868a0b00d9e9a3a99fe6b4c
BCT_URL="https://drive.google.com/uc?export=download&id=1DmMvRnferBfGe057O-sZwB5jL4j8w1Hu"   # BCT.zip, 2019-03-03

run() {
  if [ "$DRY" = 1 ]; then printf '  %s\n' "$*"; else "$@"; fi
}

fetch_at() {   # url folder commit: that commit only, no history (status 1: the folder exists)
  local url="$1" dir="$2" rev="$3"
  if [ -d "$DEST/$dir" ]; then echo "$DEST/$dir exists: left as it is"; return 1; fi
  echo "== $dir ($rev)"
  run git init --quiet "$DEST/$dir"
  run git -C "$DEST/$dir" config core.autocrlf false   # the files' own bytes (Windows: no CRLF, same fingerprints)
  run git -C "$DEST/$dir" remote add origin "$url"
  run git -C "$DEST/$dir" fetch --quiet --depth 1 origin "$rev"
  run git -C "$DEST/$dir" -c advice.detachedHead=false checkout --quiet FETCH_HEAD
}

[ "$DRY" = 1 ] || mkdir -p "$DEST"
echo "Toolboxes into $DEST"

fetch_at https://github.com/fieldtrip/fieldtrip.git fieldtrip-toolbox "$FIELDTRIP_REV" || true
if fetch_at https://github.com/brainstorm-tools/brainstorm3.git brainstorm3 "$BRAINSTORM_REV"; then
  # the zip of this version ends doc/version.txt with its commit (no final newline)
  if [ "$DRY" = 1 ]; then
    printf '  %s\n' "append '% github $BRAINSTORM_REV' to $DEST/brainstorm3/doc/version.txt"
  else
    printf '\n%% github %s' "$BRAINSTORM_REV" >> "$DEST/brainstorm3/doc/version.txt"
  fi
  # Brainstorm's fix of the GIfTI reader for Java 9+ (commit 952475b), on this version
  PATCH="$(cd "$(dirname "$0")" && pwd)/brainstorm_base64_952475b.patch"
  if [ "$DRY" = 1 ]; then
    printf '  %s\n' "patch -d $DEST/brainstorm3 -p1 --binary --forward < $PATCH"
  else
    patch -d "$DEST/brainstorm3" -p1 --binary --forward --quiet < "$PATCH"
  fi
fi
fetch_at https://github.com/corcorana/restingIAF.git restingIAF b25db386c3fbbc9a40394b56cf94aaa7ad5824cc || true

if [ "$STATIC" = 1 ]; then          # static only: no clustering, no state graph metrics
  echo
  echo "Then check the versions in MATLAB (setup_toolboxes, see the header of this script)."
  exit 0
fi

if [ ! -f "$DEST/modkmeans/modkmeans.m" ]; then
  echo "== modkmeans (modkmeans.m, calc_fitmeas.m of the Microstate EEGlab toolbox, folder MST1.0)"
  TMP="$(mktemp -d)"
  DEST_SAVED="$DEST"; DEST="$TMP"
  fetch_at https://github.com/atpoulsen/Microstate-EEGlab-toolbox.git ms "$MICROSTATE_REV"
  DEST="$DEST_SAVED"
  run mkdir -p "$DEST/modkmeans"
  run cp "$TMP/ms/MST1.0/modkmeans.m" "$TMP/ms/MST1.0/calc_fitmeas.m" "$DEST/modkmeans/"
  rm -rf "$TMP"
fi

run mkdir -p "$DEST/GitHub resources"
fetch_at https://github.com/adanjoga/cvik-toolbox.git "GitHub resources/cvik-toolbox-master" "$CVIK_REV" || true

if [ ! -d "$DEST/BCT/2019_03_03_BCT" ]; then
  echo "== BCT (release 2019-03-03, the zip linked from https://sites.google.com/site/bctnet/)"
  TMP="$(mktemp -d)"
  if run curl -fsSL -o "$TMP/BCT.zip" "$BCT_URL" && run unzip -q "$TMP/BCT.zip" -d "$TMP/bct"; then
    run mkdir -p "$DEST/BCT"
    run mv "$TMP/bct/2019_03_03_BCT" "$DEST/BCT/"
  else
    echo "  BCT could not be downloaded: get BCT.zip (2019-03-03) from https://sites.google.com/site/bctnet/ and unzip it into $DEST/BCT" >&2
  fi
  rm -rf "$TMP"
fi

echo
echo "Then check the versions in MATLAB (setup_toolboxes, see the header of this script)."
