# Design rationale

Status: **draft v0.1 — a proposal, not a standard.**
Companion documents: [SPECIFICATION.md](SPECIFICATION.md) (what the files are),
[RESEARCH_NOTES.md](RESEARCH_NOTES.md) (what was checked, with sources and measured numbers).

This page explains *why* each format was chosen. It is deliberately short.

## The problem in one paragraph

A genome annotation is a few megabytes of GFF3. The evidence that shaped it — hundreds of
gigabytes of RNA-seq alignments, protein-to-genome alignments, ab initio predictions — is
usually unavailable a few years later. A future reader (or an LLM-assisted tool) cannot
answer "why does this gene have this exon structure?" from the GFF3 alone. The proposal is a
*compact evidence shadow*: keep the **observations** that explain the structure, not the raw
data and not the pipeline.

## Design principles

1. **Observations, not pipelines.** Store what was seen (an intron crossed by reads in 14
   libraries; a protein aligned across two loci). Do not try to store how to re-derive it —
   that is what SRA, the raw reads and the workflow are for.
2. **Mature formats first.** BED, BED12, GFF3, bigWig, bgzip + tabix, TSV. The only new
   things are one column layout (the cross-sample junction table) and the manifest.
3. **One record per line, readable without a tool.** Both people and language models read
   text well and read `blockStarts`/CIGAR strings badly (see [LLM_GUIDE.md](LLM_GUIDE.md)).
4. **Coordinates are stated, not assumed.** Two conventions exist in the bundle, each tied
   to a file type, and never mixed inside a file (details below).
5. **Filters and definitions travel with the data.** Any threshold that removed records is
   written in the manifest *and* in the file header, because absence of a record can only be
   interpreted if you know what was filtered out.
6. **Size is a requirement, and it is measured.** See the table at the end.

## Format decisions

| Evidence class | File | Format | Coordinates |
|---|---|---|---|
| Annotation being explained | `annotation.gff3.gz` | GFF3 (sorted, bgzip, tabix) | 1-based, closed |
| Short-read RNA splice evidence | `evidence/rna_junctions.tsv.gz` | TSV with `##` header, tabix `-s1 -b2 -e3` | 1-based, closed (intron bases) |
| Long-read RNA structures | `evidence/long_read_models.bed.gz` | BED12 + 4 columns, tabix `-p bed` | 0-based, half-open (BED) |
| Protein-to-genome | `evidence/protein_alignments.gff3.gz` | GFF3 `protein_match` / `match_part` with `Target` | 1-based, closed |
| Ab initio / HMM predictions | `evidence/predictions.gff3.gz` | GFF3 (tool's own gene/mRNA/exon/CDS) | 1-based, closed |
| RNA-seq coverage (optional) | `evidence/coverage.bw` | bigWig | 0-based, half-open (bigWig) |
| Repeats (optional) | `evidence/repeats.bed.gz` | BED6 | 0-based, half-open |
| Datasets/samples | `samples.tsv` | TSV | – |
| Coordinate system identity | `assembly_sequences.tsv` | TSV (seqid, length, digest, aliases) | – |
| Provenance | `manifest.yaml` | YAML | – |

### RNA junctions: a plain TSV, not BED12

STAR's `SJ.out.tab` is the de facto per-sample junction format, but it is per-sample and its
columns cannot say "supported in 14 of 20 libraries". A BED12 junction track (TopHat/regtools
style) draws nicely, but its blocks are anchors, and BED has one score field — too few slots.
So the aggregate is a small TSV that **keeps STAR's coordinate convention** (first and last
*intron* base, 1-based), which is also exactly how a GFF3 intron is defined. A junction row
therefore compares to a GFF3 gene model by plain equality, with no off-by-one conversion. Columns:
`sample_count`, `group_count`, `total_unique_reads`, `max_unique_reads`,
`median_unique_reads`, `total_multimapping_reads`, `max_overhang`, `motif`.
`sample_count` versus `max_unique_reads` is what distinguishes "100 reads in one library" from
"supported independently in 20 libraries". Real data (six public libraries) show why this matters:
39 % of distinct junctions occurred in exactly one library and 35 % in all six.

### Long reads: BED12+, with the intron chain spelled out

Long reads observe *complete* structures, and BED12 is the standard, browser-native way to
carry a multi-block feature. Extra columns (BED12+4) hold `read_count`, `sample_count`,
`group_count` and `intron_chain`. The chain is redundant with the blocks, but it is written in
1-based intron coordinates so it can be compared to junction rows and GFF3 by eye; the
validator checks the redundancy. What BED12 cannot carry — per-read identity, poly-A, individual
TSS/TES positions — is deliberately dropped: models are collapsed by intron chain and report a
consensus start/end. An autoSql file is provided so the same records convert to bigBed.

### Protein alignments: GFF3, not PAF or BED12

* GFF3 is what the tools emit (miniprot `--gff`) and what NCBI publishes for alignment
  evidence (`protein_match`, `Target`, `Gap`, `pct_identity_gap` in its release directories).
  The `Target` attribute is part of the GFF3 specification.
* **BED12 loses the protein-side coordinates** of each block. Those coordinates are needed for the
  question "does the first half of this protein align to gene A and the second half to gene B?"
* PAF is compact, but is 0-based, protein CIGAR operators (`M I D F G N U V`) are understood by
  few tools, and EGAPx itself converts miniprot PAF to NCBI's ASN.1 immediately.
* Cost of GFF3: a parent line plus one child per aligned block. Region queries return the parent
  (which spans all its children), and `query_region.py` re-fetches the children.

### Coverage: optional, merged, quantised — and bigWig

Coverage is what makes *informative absence* possible: a missing junction means something very
different where the region is deeply covered than where it is not covered at all. But it is by
far the largest component. Measured on six real libraries (one 15 Mb chromosome):

| Representation | Size | vs 6 BAMs (408 MB) |
|---|---|---|
| Per-sample bigWig (mean of 6) | 1.20 MB each | 36–103× smaller each |
| Merged (sum) bigWig, full resolution | 3.61 MB | 113× |
| Merged, quantised to a geometric ladder (about 19 % steps) | 0.95 MB | 427× |
| Merged, 50-bp windows, quantised | 0.48 MB | 844× |

At full resolution bigWig and bgzip-compressed bedGraph are the **same size** (3.61 vs 3.60 MB);
after quantisation bgzipped bedGraph is 22–41 % *smaller* (0.74 vs 0.95 MB; 0.29 vs 0.48 MB),
because bigWig also stores an index and zoom levels. bigWig is chosen for random access and
browser support, not for size; `coverage.bedgraph.gz` (tabix-indexed) is an allowed alternative. Per-sample full-resolution bigWigs are what NCBI
publishes (145 files, about 10 GB, for one moth annotation), which shows the cost of not merging.
Recommendation: one merged track, quantised, optionally per strand; per-tissue tracks only when
tissue specificity is the point.

### Junctions from BAMs, single-exon models, browser views

* **BAM route.** EGAPx never writes `SJ.out.tab`, but every pipeline has BAMs, so `bam_to_junctions.py` derives STAR-format tables from them. The rules were
  found by testing candidates against STAR's own output and are exact on real data (unique and multimapper counts are *distinct read names*, so mates that both cross a junction
  count once; the overhang is the longest, over reads and mates, of the smaller of the two aligned runs beside the junction, with runs split at insertions, deletions and other
  introns). STAR's default output filters could not be reproduced exactly (a simple emulation of the documented defaults recovered 99.2–99.7 % of STAR's junctions in two libraries but also kept 9–12 % that STAR dropped), so the
  BAM-derived table is documented as unfiltered instead of presented as STAR-equivalent.
* **Single-exon long-read models** are clustered by start and end within a tolerance and always reported apart from spliced models, because they are a different kind of
  evidence (single-exon genes, but also fragments and unspliced RNA).
* **Browser views are derived, not canonical.** Testing in IGV showed that a GFF3 `protein_match` with `match_part` children is drawn as one unbroken bar, and that a
  junction table is not a browser format. The bundle keeps the GFF3 (protein-side coordinates are needed) and `export_browser_tracks.py` writes a BED12 view and IGV's native
  splice-junction BED.

### Predictions: keep the predictor's own GFF3

No attempt is made to standardise scores. Column 2 (`source`) carries a short tag that the
manifest maps to program, version and model. That is enough to answer "was this region called
mainly from experimental evidence or ab initio?".

### Provenance: a manifest, and two small tables

`manifest.yaml` is the entry point and sits at the bundle root. The initial sketch had a
`metadata/` folder and a `software.tsv`. Both were changed:

* **`software.tsv` was folded into the manifest.** Each evidence class lists its ordered
  `programs` (name, version, role, parameters). A second table would duplicate that and can
  disagree with it.
* **`samples.tsv` stays a table** — one row per library, so hundreds are fine — with a `group`
  column (tissue, stage, treatment) that makes `group_count` computable.
* **`assembly_sequences.tsv` was added.** The most common failure years later is not a wrong
  coordinate but a wrong *name*: in an EGAPx 1.0 run I inspected, the Gnomon evidence report,
  the STAR BAM headers and the miniprot PAF all use `lcl|NC_051855.1`, while the final GFF3
  uses `NC_051855.1`. This table (seqid, length, optional GA4GH refget digest, aliases) pins
  the coordinate system and lets a validator catch the mismatch.

## What is intentionally not done

* No pipeline, no workflow manager, no server, no custom binary format.
* No hard-coded biological filters. Defaults are minimal (`>=1` unique read overall) and are
  always recorded.
* No claim that a gene model is wrong. The tools report *support*, *absence of retained
  evidence* and *competing evidence* separately, and never turn the second into the third.
* No replacement for SRA, raw reads, BAM/CRAM or full workflow provenance.

## Sizes to expect

For one recent RefSeq annotation of *Helicoverpa armigera* (accessed 2026-09-25), the release
directory holds a 6.3 MB `genomic.gff.gz`, 145 per-run coverage bigWigs (about 10.3 GB) and a
50 MB transcript-alignment GFF3. The 145 RNA-seq runs behind it total 8.2 billion reads and
1.17 Tbp. The target for this proposal is a bundle within **tens of times** the size of the
annotation — for the junction, long-read, protein and prediction tables that is realistic; the
coverage track is the part to budget for. See [RESEARCH_NOTES.md](RESEARCH_NOTES.md).
