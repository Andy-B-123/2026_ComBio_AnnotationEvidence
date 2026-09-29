#!/usr/bin/env bash
# Build an annotation evidence bundle for YOUR genome with this repository's scripts.
# Copy this file, set the variables in the "INPUTS" block (or export them), and run it. Every step is optional except
# the sequence table, the annotation and the manifest; a step runs only if its input variable is set.
#
#   OUT=my_bundle GENOME=genome.fa ANNOTATION=annotation.gff3 SAMPLES=samples.tsv PROVENANCE=provenance.yaml \
#   SHORT_BAMS="S1=/data/S1.bam S2=/data/S2.bam" PROTEINS_FASTA=proteins.faa bash pipeline_template.sh
#
# Tools on PATH: python3 (+PyYAML), bgzip, tabix (htslib). Per step: samtools + bedtools (long reads, coverage),
# miniprot (proteins), bedGraphToBigWig from UCSC (bigWig coverage; otherwise a .bedgraph.gz is written).
# Nothing here is a workflow manager: read it top to bottom, it is the whole recipe. See docs/USER_GUIDE.md.
set -euo pipefail

# ------------------------------- INPUTS ---------------------------------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${SCRIPTS:-$HERE/../scripts}"            # where the repository scripts are
PY="${PYTHON:-python3}"
OUT="${OUT:?set OUT (the bundle directory to create)}"
GENOME="${GENOME:?set GENOME (assembly FASTA)}"
ANNOTATION="${ANNOTATION:?set ANNOTATION (the GFF3 being explained; .gff3 or .gff3.gz)}"
PROVENANCE="${PROVENANCE:?set PROVENANCE (your manifest template: copy schema/manifest.template.yaml and fill it in)}"
SAMPLES="${SAMPLES:?set SAMPLES (samples.tsv with a row per library; ids must match those below)}"

SHORT_BAMS="${SHORT_BAMS:-}"        # "id=path.bam ..." coordinate-sorted, spliced short-read RNA-seq (STAR BAMs incl. EGAPx)
SHORT_SJ="${SHORT_SJ:-}"            # alternative: "id=path.SJ.out.tab ..." STAR junction tables
LONG_BAMS="${LONG_BAMS:-}"          # "id=path.bam ..." long-read RNA-seq alignments (e.g. minimap2 -ax splice)
LONG_BED12="${LONG_BED12:-}"        # alternative: "id=path.bed12 ..." per-read BED12 you made yourself
PROTEINS_GFF="${PROTEINS_GFF:-}"    # miniprot --gff output, or
PROTEINS_FASTA="${PROTEINS_FASTA:-}"  # a protein FASTA to align with miniprot here
PREDICTIONS_GFF3="${PREDICTIONS_GFF3:-}"  # ab initio predictions (GFF3; column 2 = predictor tag listed in PROVENANCE)
REPEATS_BED="${REPEATS_BED:-}"      # optional repeat intervals (BED6)
COVERAGE="${COVERAGE:-no}"          # yes: merged, quantised coverage from SHORT_BAMS
STRIP_PREFIX="${STRIP_PREFIX:-}"    # e.g. 'lcl|' for EGAPx BAMs: removed from sequence names everywhere
PROTEINS_DB="${PROTEINS_DB:-protein set}"; PROTEINS_TAXID="${PROTEINS_TAXID:-}"; MINIPROT_VERSION="${MINIPROT_VERSION:-.}"
THREADS="${THREADS:-4}"
STAR_SETTINGS="${STAR_SETTINGS:-}"  # free text recorded in the junction table: the STAR --outSJfilter* values used upstream (default: see below)
# filters (recorded in the manifest; choose deliberately, see docs/SPECIFICATION.md section 8)
MIN_TOTAL_UNIQUE="${MIN_TOTAL_UNIQUE:-1}"   # keep junctions with >= N unique reads in total
MIN_SAMPLES="${MIN_SAMPLES:-1}"             # keep junctions seen in >= N samples
LONG_MAPQ="${LONG_MAPQ:-20}"                # long-read alignment MAPQ; primary alignments only
LONG_MIN_READS="${LONG_MIN_READS:-1}"       # keep long-read models with >= N reads
LONG_REGION="${LONG_REGION:-}"              # optional: restrict long-read BAMs to one sequence (needs an index)
PROTEIN_MIN_COVERAGE="${PROTEIN_MIN_COVERAGE:-0}"   # keep protein alignments covering >= N % of the protein
CREATED="${CREATED:-$(date +%F)}"
# ---------------------------------------------------------------------------------------------------------------------

for tool in bgzip tabix; do command -v "$tool" >/dev/null || { echo "error: $tool (htslib) not on PATH" >&2; exit 1; }; done
[ -e "$OUT/manifest.yaml" ] && { echo "error: $OUT/manifest.yaml exists; choose a new OUT" >&2; exit 1; }
TMP="$(mktemp -d)"; trap 'rm -rf "${TMP:?}"' EXIT
mkdir -p "$OUT/evidence"
SUMMARIES=()                                  # --summary arguments for make_manifest.py build
PFX=(); [ -n "$STRIP_PREFIX" ] && PFX=(--strip-prefix "$STRIP_PREFIX")

gz_gff3() { { grep '^#' "$1" || true; grep -v '^#' "$1" | LC_ALL=C sort -t$'\t' -k1,1 -k4,4n -s; } | bgzip -c > "$2"; tabix -f -p gff "$2"; }
gz_bed()  { LC_ALL=C sort -t$'\t' -k1,1 -k2,2n -s "$1" | bgzip -c > "$2"; tabix -f -p bed "$2"; }
opencat() { case "$1" in *.gz) gzip -dc "$1";; *) cat "$1";; esac; }

echo "== 1. identity: samples, sequence names and lengths"
cp "$SAMPLES" "$OUT/samples.tsv"
"$PY" "$SCRIPTS/make_manifest.py" sequences --fasta "$GENOME" ${STRIP_PREFIX:+--alias-prefix "$STRIP_PREFIX"} -o "$OUT/assembly_sequences.tsv"

echo "== 2. annotation (sorted, bgzip, tabix)"
opencat "$ANNOTATION" > "$TMP/annotation.gff3"
gz_gff3 "$TMP/annotation.gff3" "$OUT/annotation.gff3.gz"

if [ -n "$SHORT_BAMS" ] || [ -n "$SHORT_SJ" ]; then
  echo "== 3. short-read RNA junctions"
  SJ_FILES=()
  for pair in $SHORT_SJ; do SJ_FILES+=("$pair"); done
  for pair in $SHORT_BAMS; do
    id="${pair%%=*}"; bam="${pair#*=}"
    "$PY" "$SCRIPTS/bam_to_junctions.py" "$bam" --fasta "$GENOME" "${PFX[@]}" -o "$TMP/$id.SJ.out.tab" --summary-json "$TMP/b2j_$id.json" 2>/dev/null
    SJ_FILES+=("$id=$TMP/$id.SJ.out.tab"); SUMMARIES+=(--summary "rna_short=$TMP/b2j_$id.json")
  done
  if [ -z "$STAR_SETTINGS" ]; then
    if [ -n "$SHORT_BAMS" ] && [ -z "$SHORT_SJ" ]; then STAR_SETTINGS="not applicable: junctions extracted from BAMs by bam_to_junctions.py (unfiltered)"
    else STAR_SETTINGS="STAR defaults assumed: set STAR_SETTINGS to the --outSJfilter* values you used"; fi
  fi
  "$PY" "$SCRIPTS/aggregate_star_junctions.py" --samples "$SAMPLES" --min-total-unique "$MIN_TOTAL_UNIQUE" --min-samples "$MIN_SAMPLES" \
      --star-settings "$STAR_SETTINGS" "${PFX[@]}" --summary-json "$TMP/agg.json" -o "$OUT/evidence/rna_junctions.tsv.gz" "${SJ_FILES[@]}" 2>&1 | tail -1
  SUMMARIES+=(--summary "rna_short=$TMP/agg.json")
fi

if [ -n "$LONG_BAMS" ] || [ -n "$LONG_BED12" ]; then
  echo "== 4. long-read models (primary alignments, MAPQ >= $LONG_MAPQ) -> BED12 -> intron-chain models"
  BEDS=()
  for pair in $LONG_BED12; do BEDS+=("$pair"); done
  for pair in $LONG_BAMS; do
    id="${pair%%=*}"; bam="${pair#*=}"
    samtools view -b -F 0x900 -q "$LONG_MAPQ" "$bam" $LONG_REGION | bedtools bamtobed -bed12 -i - > "$TMP/$id.bed12"
    BEDS+=("$id=$TMP/$id.bed12")
  done
  "$PY" "$SCRIPTS/collapse_long_reads.py" --samples "$SAMPLES" --min-reads "$LONG_MIN_READS" "${PFX[@]}" \
      --summary-json "$TMP/lr.json" -o "$OUT/evidence/long_read_models.bed.gz" "${BEDS[@]}" 2>&1 | tail -1
  SUMMARIES+=(--summary "rna_long=$TMP/lr.json")
fi

if [ -n "$PROTEINS_GFF" ] || [ -n "$PROTEINS_FASTA" ]; then
  echo "== 5. protein-to-genome alignments (miniprot -> protein_match / match_part)"
  if [ -z "$PROTEINS_GFF" ]; then
    miniprot -It"$THREADS" --gff "$GENOME" "$PROTEINS_FASTA" > "$TMP/miniprot.gff" 2> "$TMP/miniprot.log"
    PROTEINS_GFF="$TMP/miniprot.gff"
  fi
  "$PY" "$SCRIPTS/miniprot_to_evidence.py" "$PROTEINS_GFF" --db "$PROTEINS_DB" ${PROTEINS_TAXID:+--taxid "$PROTEINS_TAXID"} \
      --aligner-version "$MINIPROT_VERSION" --min-pct-coverage "$PROTEIN_MIN_COVERAGE" "${PFX[@]}" \
      --summary-json "$TMP/prot.json" -o "$OUT/evidence/protein_alignments.gff3.gz" 2>&1 | tail -1
  SUMMARIES+=(--summary "proteins=$TMP/prot.json")
fi

if [ -n "$PREDICTIONS_GFF3" ]; then
  echo "== 6. ab initio predictions (kept as GFF3)"
  opencat "$PREDICTIONS_GFF3" > "$TMP/pred.gff3"; gz_gff3 "$TMP/pred.gff3" "$OUT/evidence/predictions.gff3.gz"
fi
if [ -n "$REPEATS_BED" ]; then
  echo "== 6b. repeats"; opencat "$REPEATS_BED" > "$TMP/rep.bed"; gz_bed "$TMP/rep.bed" "$OUT/evidence/repeats.bed.gz"
fi

if [ "$COVERAGE" = yes ] && [ -n "$SHORT_BAMS" ]; then
  echo "== 7. coverage: merged, spliced, quantised"
  BAMS=(); for pair in $SHORT_BAMS; do BAMS+=("${pair#*=}"); done
  EXT=bw; command -v bedGraphToBigWig >/dev/null || EXT=bedgraph.gz
  "$PY" "$SCRIPTS/make_coverage.py" --sequences "$OUT/assembly_sequences.tsv" "${PFX[@]}" --summary-json "$TMP/cov.json" \
      -o "$OUT/evidence/coverage.$EXT" "${BAMS[@]}" 2>&1 | tail -1
  SUMMARIES+=(--summary "coverage=$TMP/cov.json")
fi

echo "== 8. manifest (provenance from $PROVENANCE + computed facts)"
"$PY" "$SCRIPTS/make_manifest.py" build "$OUT" --template "$PROVENANCE" --fasta "$GENOME" --created "$CREATED" "${SUMMARIES[@]}"

echo "== 9. validate"
"$PY" "$SCRIPTS/validate_bundle.py" "$OUT"

echo "== done: $OUT   (try: python $SCRIPTS/summarise_gene.py $OUT <gene-id>)"
