#!/usr/bin/env bash
# Render a locus of a bundle to a PNG with IGV desktop, headless. Optional helper; the bundle does not need IGV.
#
#   IGV=/path/to/IGV_Linux_2.19.8  scripts/igv_snapshot.sh SESSION.xml LOCUS OUT.png [WIDTHxHEIGHT]
#
# SESSION.xml comes from export_browser_tracks.py; LOCUS is e.g. scaffold_1:900-5400 or a feature name.
# Needs a display or xvfb-run (Linux). IGV writes preferences to $HOME/igv, so this script points HOME at a
# temporary directory (IGV_HOME_DIR, default: a fresh mktemp -d) to keep your real home untouched.
# IGV_EXPAND: comma-separated track names to show expanded (default: the evidence tracks). Expanding a dense
# real annotation track can make IGV allocate a huge blank area, so "annotation" is not in the default; add it for
# small annotations, e.g. IGV_EXPAND="annotation,long_read_models.bed.gz".
# Tested with IGV 2.19.8 (the "WithJava" Linux build).
set -euo pipefail
[ $# -ge 3 ] || { sed -n '2,9p' "$0"; exit 2; }
SESSION=$(realpath "$1"); LOCUS=$2; OUT=$(realpath -m "$3"); SIZE=${4:-1400x800}
: "${IGV:?set IGV to the IGV installation directory (containing igv.sh)}"
[ -x "$IGV/igv.sh" ] || { echo "error: $IGV/igv.sh not found" >&2; exit 2; }
export HOME=${IGV_HOME_DIR:-$(mktemp -d)}
BATCH=$(mktemp)
mkdir -p "$(dirname "$OUT")"
EXPAND=${IGV_EXPAND:-long_read_models.bed.gz,protein alignments (BED12 view),predictions.gff3.gz}
{
  echo new
  echo "load $SESSION"
  echo "goto $LOCUS"
  IFS=',' read -ra TRACKS <<< "$EXPAND"
  for tr in "${TRACKS[@]}"; do echo "expand $tr"; done
  echo "setSleepInterval 1500"
  echo "snapshotDirectory $(dirname "$OUT")"
  echo "snapshot $(basename "$OUT")"
  echo exit
} > "$BATCH"
RUN=("$IGV/igv.sh" -b "$BATCH")
if [ -z "${DISPLAY:-}" ]; then RUN=(xvfb-run -a -s "-screen 0 ${SIZE}x24" "${RUN[@]}"); fi
"${RUN[@]}" > "$BATCH.log" 2>&1 || { echo "IGV failed; see $BATCH.log" >&2; exit 1; }
[ -s "$OUT" ] || { echo "no snapshot written; see $BATCH.log" >&2; exit 1; }
echo "wrote $OUT"
