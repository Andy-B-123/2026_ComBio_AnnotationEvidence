#!/usr/bin/env bash
# Build examples/minimal/bundle from examples/minimal/inputs with the repository's own scripts.
# This is the whole "workflow": a few commands, no workflow manager.
#
#   bash examples/minimal/build.sh [OUT_DIR]
#
# Needs: python3 with PyYAML, and htslib (bgzip, tabix).
# Coverage: uses UCSC bedGraphToBigWig if present, otherwise the pre-built inputs/coverage.prebuilt.bw.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/../../scripts"
IN="$HERE/inputs"
OUT="${1:-$HERE/bundle}"
PY="${PYTHON:-python3}"
for tool in bgzip tabix; do
  command -v "$tool" >/dev/null || { echo "error: $tool (htslib) not found on PATH" >&2; exit 1; }
done

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
rm -rf "$OUT"; mkdir -p "$OUT/evidence"

# sort a GFF3 by seqid then start (## header lines stay on top), bgzip, tabix
gz_gff3() { { grep '^#' "$1"; grep -v '^#' "$1" | LC_ALL=C sort -t$'\t' -k1,1 -k4,4n -s; } | bgzip -c > "$2"; tabix -f -p gff "$2"; }
# same for a BED file (0-based start in column 2)
gz_bed()  { LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -s "$1" | bgzip -c > "$2"; tabix -f -p bed "$2"; }

echo "== identity: samples and sequence table"
cp "$IN/samples.tsv" "$OUT/samples.tsv"
"$PY" "$S/make_manifest.py" sequences --fasta "$IN/genome.fa" --alias-prefix 'lcl|' -o "$OUT/assembly_sequences.tsv"

echo "== annotation"
gz_gff3 "$IN/annotation.gff3" "$OUT/annotation.gff3.gz"

echo "== short-read RNA junctions (STAR SJ.out.tab -> rna_junctions)"
"$PY" "$S/aggregate_star_junctions.py" --samples "$IN/samples.tsv" --star-settings "not applicable (synthetic)" \
      --summary-json "$TMP/rna_short.json" -o "$OUT/evidence/rna_junctions.tsv.gz" "$IN"/star/*.SJ.out.tab

echo "== long-read models (per-read BED12 -> long_read_models)"
"$PY" "$S/collapse_long_reads.py" --samples "$IN/samples.tsv" --summary-json "$TMP/rna_long.json" \
      -o "$OUT/evidence/long_read_models.bed.gz" "$IN"/long_reads/*.bed12

echo "== protein alignments (miniprot --gff -> protein_alignments)"
"$PY" "$S/miniprot_to_evidence.py" "$IN/miniprot.gff" --aligner-version 0.18 --db "synthetic proteins" \
      --summary-json "$TMP/proteins.json" -o "$OUT/evidence/protein_alignments.gff3.gz"

echo "== predictions, repeats"
gz_gff3 "$IN/predictions.gff3" "$OUT/evidence/predictions.gff3.gz"
gz_bed  "$IN/repeats.bed"      "$OUT/evidence/repeats.bed.gz"

echo "== coverage (bedGraph -> bigWig)"
grep -v '^#' "$OUT/assembly_sequences.tsv" | cut -f1,2 > "$TMP/chrom.sizes"
if command -v bedGraphToBigWig >/dev/null; then
  LC_ALL=C sort -k1,1 -k2,2n "$IN/coverage.bedgraph" > "$TMP/coverage.sorted.bedgraph"
  bedGraphToBigWig "$TMP/coverage.sorted.bedgraph" "$TMP/chrom.sizes" "$OUT/evidence/coverage.bw"
  [ "${UPDATE_PREBUILT:-0}" = 1 ] && cp "$OUT/evidence/coverage.bw" "$IN/coverage.prebuilt.bw"
else
  echo "   bedGraphToBigWig not found: using inputs/coverage.prebuilt.bw"
  cp "$IN/coverage.prebuilt.bw" "$OUT/evidence/coverage.bw"
fi

echo "== manifest"
cp "$HERE/bundle_README.md" "$OUT/README.md"
"$PY" "$S/make_manifest.py" build "$OUT" --template "$HERE/manifest.template.yaml" --fasta "$IN/genome.fa" \
      --summary rna_short="$TMP/rna_short.json" --summary rna_long="$TMP/rna_long.json" \
      --summary proteins="$TMP/proteins.json" --created 2026-09-25

echo "== validate"
"$PY" "$S/validate_bundle.py" "$OUT"
