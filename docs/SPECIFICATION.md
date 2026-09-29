# Annotation evidence bundle — specification

**Version 0.1 (draft). A proposal and prototype, not an established standard.**
The key words MUST, SHOULD and MAY are used in their usual sense. Rationale for each choice is in
[DESIGN_RATIONALE.md](DESIGN_RATIONALE.md); sources and measurements are in [RESEARCH_NOTES.md](RESEARCH_NOTES.md).

An **annotation evidence bundle** is a directory that holds a genome annotation together with a compact,
region-queryable record of the *observations* that shaped it (RNA-seq splice junctions, long-read transcript
structures, protein-to-genome alignments, ab initio predictions, optionally coverage), plus the provenance
needed to interpret them years later. It is not a substitute for raw reads, BAM/CRAM or the workflow.

## 1. Layout

```text
<bundle>/
├── manifest.yaml                       REQUIRED  provenance; entry point
├── annotation.gff3.gz (+ .tbi)         REQUIRED  the annotation being explained (sorted, bgzip, tabix)
├── assembly_sequences.tsv              REQUIRED  names and lengths of the sequences all coordinates refer to
├── samples.tsv                         REQUIRED if any sample-based class is present
├── README.md                           SHOULD    a few lines for a human reader
└── evidence/
    ├── rna_junctions.tsv.gz (+ .tbi)         short-read RNA-seq junctions, collapsed across samples
    ├── long_read_models.bed.gz (+ .tbi)      long-read transcript models (BED12+4)
    ├── protein_alignments.gff3.gz (+ .tbi)   protein-to-genome alignments
    ├── predictions.gff3.gz (+ .tbi)          ab initio / HMM gene predictions
    ├── coverage.bw                           OPTIONAL merged RNA-seq coverage (bigWig)
    └── repeats.bed.gz (+ .tbi)               OPTIONAL repeat intervals (BED6)
```

File names above are conventional; the **manifest lists the actual files** and their `format`. Every evidence class
is optional, but each of `rna_short`, `rna_long`, `proteins`, `predictions` MUST appear in the manifest with
`present: true` or `present: false` and a `reason`. Compressed, indexed files use **bgzip** and a **tabix** index
(`.tbi`); readers SHOULD work without the index (by scanning) and MUST NOT require it.
The assembly FASTA is **not** part of the bundle.

## 2. Conventions

### 2.1 Coordinates

Two conventions exist, each tied to a file type. They are never mixed inside one file.

| Convention | Used by | Interval `[a, b]` covers |
|---|---|---|
| **1-based, closed** | GFF3 files; `rna_junctions` (columns `intron_start`, `intron_end`); the `intron_chain` column | bases `a … b` inclusive; the first base of a sequence is 1 |
| **0-based, half-open** | BED/BED12 (`chromStart`, `chromEnd`, `blockStarts`), bigWig, bedGraph | bases `a+1 … b` in 1-based terms; the first base is 0 and `chromEnd` is excluded |

Conversion: `bed_start = start − 1`, `bed_end = end` (and back: `start = bed_start + 1`, `end = bed_end`).
Length is `end − start + 1` (closed) or `chromEnd − chromStart` (BED).

**Introns** are always given as the **first and last intronic base** (1-based, closed). For consecutive exons
`[s1, e1]` and `[s2, e2]` the intron is `[e1 + 1, s2 − 1]`. This is how GFF3 defines an intron, and it is what
STAR reports in `SJ.out.tab`. (Many other tools instead report the flanking *exonic* bases; STAR's manual points
this out. Convert on input.) Checked on real data (six public libraries, one chromosome): 6 510 of 7 803 RefSeq introns (83.4 %)
match a junction row exactly with no conversion, and **0.0 %** match after shifting the junction by one base in any
direction (start ±1, end ±1, both −1, or the exonic-flank convention).

**BED12 blocks.** `blockStarts` are relative to `chromStart`; the first MUST be 0 and the last block MUST end at
`chromEnd`. Block `i` covers 1-based `chromStart + blockStarts[i] + 1 … chromStart + blockStarts[i] + blockSizes[i]`.
Do not do this arithmetic by hand: use `evidence_common.exons_from_bed12`, or the `intron_chain` column.

**Region strings** (`seqid:start-end`) in the tools and in `tabix` are 1-based, closed, whatever the file type.

### 2.2 Sequence names

`seqid` values MUST be identical, byte for byte, across the annotation, every evidence file, the bigWig chromosome
names and `assembly_sequences.tsv`. Do not carry pipeline-internal prefixes (`lcl|`) into a bundle; record them as
`aliases` in `assembly_sequences.tsv` if they matter. The validator reports a name that matches only an alias.

### 2.3 Strand, missing values, text

Strand is `+`, `-` or `.` (unknown / not applicable). Missing values are `.`. Files are UTF-8, tab-delimited,
LF line endings. Numbers use `.` as decimal separator and no thousands separators.

### 2.4 Self-describing headers

TSV and GFF3 files SHOULD begin with pragma lines `##ae.<key>=<value>` (application-specific directives are permitted
by GFF3; tabix and bedtools treat every `#` line as a comment). Keys: `profile`, `version`, `coordinates`,
`generated_by`, and the filters listed for each format. They repeat what the manifest says so that a file stays
interpretable if it is separated from its bundle. BED files carry no header (BED has no comment convention that all
browsers accept); their provenance is in the manifest only.

### 2.5 Sorting and indexing

Files MUST be sorted by `seqid` (each sequence contiguous) and then by start, so that `tabix` can index them:
TSV `tabix -s1 -b2 -e3`, GFF3 `tabix -p gff`, BED `tabix -p bed`.

## 3. `manifest.yaml`

Schema: [`schema/manifest.schema.yaml`](../schema/manifest.schema.yaml) (JSON Schema, checked by
`scripts/validate_bundle.py`). Commented skeleton: [`schema/manifest.template.yaml`](../schema/manifest.template.yaml).
Top-level keys:

| Key | Content |
|---|---|
| `spec_version` | `"0.1"` |
| `bundle` | `id`, `created` (YYYY-MM-DD), `title`, `description`, `license`, `created_by`, `notes` |
| `assembly` | `accession` (or `none`), `name`, `filename`, `checksum` (sha256 of the FASTA), `sequences_file`, `n_sequences`, `total_length`, `taxon` |
| `annotation` | `filename`, `checksum`, `pipeline`, `version`, `date`, `source` |
| `evidence.<class>` | see below |

Each evidence class (`rna_short`, `rna_long`, `proteins`, `predictions`, and optionally `coverage`, `repeats` or your own):

| Field | Meaning |
|---|---|
| `present` | boolean; if false, `reason` is required |
| `description`, `data_source` | what the evidence is and where it came from (`kind`: SRA, ENA, GenBank, RefSeq, UniProt, local, synthetic, other; `accessions`; `url`) |
| `n_datasets` | number of libraries / protein sets / predictor runs |
| `programs` | **ordered** steps, primary first: `name`, `version` (both required), `role`, `parameters`, and for predictors `source_tag` (the GFF3 column-2 value) and `model` |
| `databases` | reference databases with `name`, `version`, `taxon` (required for proteins) |
| `guided_by_annotation` | true if mapping used a junction database or GTF derived from the annotation being explained, or two-pass mode (junction support can then be circular) |
| `summarisation` | how raw records were reduced to the retained ones |
| `filters` | **every threshold that removed records**, see §8 |
| `stats` | counts reported by the producer scripts, keyed by script name (records in/out, records dropped by each filter) |
| `outputs` | list of `{file, format, coordinates, records, bytes, sha256, index}`; `format` is one of `rna_junctions`, `long_read_models`, `protein_alignments`, `predictions`, `coverage_bigwig`, `coverage_bedgraph`, `repeats_bed6`, `other` |

Provenance the pipeline cannot know (sources, versions, why something is absent) is written by hand into a template;
`scripts/make_manifest.py build` adds file sizes, checksums, record counts and the filters and statistics reported by
the producer scripts, so the numbers in the manifest are never typed.

## 4. `samples.tsv`

One row per library/dataset; a header line starting `#`. Columns (all present; `.` if unknown):
`sample_id` (unique; used by every evidence file), `evidence_class` (`rna_short` or `rna_long`), `group`
(tissue, stage or treatment — the unit for `group_count`), `description`, `run_accession`, `biosample`, `platform`,
`layout`, `strandedness`, `n_reads_input`, `n_reads_uniquely_mapped`, `notes`.

## 5. `assembly_sequences.tsv`

Header `#seqid  length  sha512t24u  aliases`. `sha512t24u` is the GA4GH refget sequence digest (sequence upper-cased,
SHA-512, first 24 bytes, base64url; the refget identifier is `SQ.` + digest). `aliases` is a comma-separated list of
other names for the sequence, or `.`. `scripts/make_manifest.py sequences` writes it from the FASTA.

## 6. `annotation.gff3.gz`

The annotation the evidence explains: GFF3, sorted, bgzip, tabix. Gene → transcript → exon/CDS with `ID`/`Parent`.
Features MAY span several lines with one `ID` if type, sequence and strand agree (NCBI does this for CDS). Partial
models SHOULD carry `partial=true`, which the tools report.

## 7. Evidence formats

### 7.1 `rna_junctions` — short-read splice junctions

Full column definition, the STAR mapping and worked examples: [`schema/rna_junctions.md`](../schema/rna_junctions.md).

TSV, header `##ae.*` lines then `#seqid …`. **One row per intron** (seqid, strand, first and last intron base):

| # | Column | Meaning |
|---|---|---|
| 1 | `seqid` | |
| 2–3 | `intron_start`, `intron_end` | first and last intron base, 1-based closed |
| 4 | `strand` | `+ - .` (`.` for non-canonical junctions, whose strand STAR does not infer) |
| 5 | `motif` | `GT/AG`, `GC/AG`, `AT/AC` or `non-canonical`, in **transcript orientation** (STAR codes 1/2, 3/4, 5/6 → these) |
| 6 | `sample_count` | supporting samples; a sample supports a junction if it has ≥ `sample_min_unique_reads` uniquely mapping reads on it |
| 7 | `group_count` | distinct `group` values among supporting samples; `.` if groups unknown |
| 8 | `total_unique_reads` | uniquely mapping reads, summed over supporting samples |
| 9 | `max_unique_reads` | the largest count in any one supporting sample |
| 10 | `median_unique_reads` | median over supporting samples |
| 11 | `total_multimapping_reads` | multimapper reads over supporting samples |
| 12 | `max_overhang` | the largest STAR spliced-alignment overhang in any supporting sample |

`sample_count` next to `max_unique_reads` separates *100 reads in one library* from *supported independently in 20*.
Header keys: `n_samples`, `n_groups`, `sample_min_unique_reads`, `min_samples`, `min_total_unique_reads`,
`upstream_star_filters` (free text: the `--outSJfilter*` values; STAR's `SJ.out.tab` is **already filtered**, and
annotated junctions bypass those filters), `seqid_prefix_stripped`. If mapping was guided by the annotation (a GTF at index
generation, or two-pass mode) say so in `guided_by_annotation`. Junction tables can also be **derived from BAMs**
(`scripts/bam_to_junctions.py`, needed for EGAPx, which never writes `SJ.out.tab`): the counts are identical to STAR's
(distinct read names; overhang = longest min of the two aligned runs beside the junction), but the result is unfiltered, so
it has many more low-support junctions than STAR's output; choose `min_total_unique_reads` accordingly.

### 7.2 `long_read_models` — long-read transcript models (BED12+4)

No header line. Reads with the same seqid, strand and **intron chain** are one model.

| # | Column | Meaning |
|---|---|---|
| 1–3 | `chrom`, `chromStart`, `chromEnd` | BED (0-based, half-open); start/end are the **median observed** start and end of the reads |
| 4 | `name` | unique model id (`LR000001`) |
| 5 | `score` | `min(read_count, 1000)` (BED shading only) |
| 6 | `strand` | |
| 7–8 | `thickStart`, `thickEnd` | equal to `chromStart`, `chromEnd` (no CDS claim) |
| 9 | `itemRgb` | `0,0,0` |
| 10–12 | `blockCount`, `blockSizes`, `blockStarts` | exon blocks, BED12 rules |
| 13 | `read_count` | reads collapsed into the model |
| 14 | `sample_count` | samples with ≥ 1 read |
| 15 | `group_count` | distinct groups; `.` if unknown |
| 16 | `intron_chain` | `1301-1700;1901-2500`: introns, 1-based closed, or `.` for a single-exon model; redundant with the blocks and **checked by the validator** |

**Single-exon models.** A long read with no intron has no intron chain, so it cannot be collapsed that way. Such reads are clustered instead: on
the same sequence and strand, reads whose start and end each lie within `--single-exon-tolerance` bases (default 50) of the cluster's first read form one
model, with one block and `intron_chain` `.`. They are **kept and counted** (single-exon genes and unspliced transcripts matter), but they are often fragments or
unspliced RNA, so consumers must report them separately from spliced models and never as intron or chain support. Record the tolerance in `filters`.
Gaps shorter than `min_intron_length` are not introns. An autoSql definition for `bedToBigBed` is in
[`schema/long_read_models.as`](../schema/long_read_models.as) (tested with UCSC `bedToBigBed` v2.8: `bedToBigBed -type=bed12+4 -as=long_read_models.as models.bed chrom.sizes models.bb`; note that
`bedToBigBed` prints "not in legal format" but still exits 0 when it rejects an autoSql file, so check that the output exists). Real long-read data contain many models supported by a single read; `min_reads` is therefore a filter worth choosing deliberately and recording.

### 7.3 `protein_alignments` — protein-to-genome (GFF3)

Per alignment one parent and one child per aligned block (one exon of the alignment):

* `protein_match` (SO:0000349): `ID` (unique), `Target=<protein> <start> <end>` (protein residues, 1-based, closed),
  strand in column 7, score in column 6 (aligner score), and attributes `rank`, `pct_identity`, `pct_positive`,
  `pct_coverage` (percent, 0–100; coverage = aligned protein span ÷ protein length), `query_length`, `frameshifts`,
  `stop_codons`, and optionally `db`, `taxid`. Column 2 is the aligner (`miniprot`).
* `match_part` (SO:0000039): `Parent`, `Target=<protein> <start> <end>` for that block, `pct_identity`.

The parent spans all its blocks (excluding any stop codon). Blocks of one alignment are contiguous in protein
coordinates: a codon split by an intron belongs to the downstream block (miniprot's rule). Attribute names are
lower-case, which GFF3 leaves to applications. Percentages, not fractions, so the unit cannot be confused.
`scripts/miniprot_to_evidence.py` converts `miniprot --gff`; by default only rank-1 alignments are kept (as NCBI does).

### 7.4 `predictions` — ab initio / HMM gene predictions (GFF3)

The predictor's own `gene`/`mRNA`/`exon`/`CDS` records. Column 2 is a short tag that the manifest maps, via `source_tag`
under `programs`, to program, version and model. Column 6 carries the predictor's confidence where it has one; no
attempt is made to standardise scores. IDs are unique within the file.

### 7.5 `coverage` — merged RNA-seq depth (optional)

A bigWig (or `coverage_bedgraph`: bgzip + tabix). Sequence names as §2.2. Values SHOULD be on the **depth scale** (reads
covering a base) so that a value means what it says, summed over libraries, spliced (introns not counted), and
SHOULD be quantised to a geometric ladder (e.g. `2^(k/4) − 1`, about 19 % steps) or binned to reduce size.
Record in `filters`: `quantisation`, `strand` (`both`, or one file per strand with `strand: plus|minus` on the output),
and whether multimappers were counted. `scripts/make_coverage.py` implements this basic recipe (spliced depth from BAMs, summed, snapped, bigWig or bedGraph). Coverage is what allows *informative absence* (§8); it is the largest component
and the reason this class is optional. Measured sizes are in the design rationale.

### 7.6 `repeats` — repeat intervals (optional)

BED6 (`name` = repeat class/family). Not interpreted by v0.1 tools beyond `query_region.py`.

### 7.7 Other classes

Any additional class MAY be added under `evidence:` with the same fields and `format: other`; document its format in `notes`.

## 8. Filters, and what absence means

Any threshold that removed records MUST be recorded in `filters` (and in the file header where the format has one).
Recommended keys:

| Class | Keys |
|---|---|
| `rna_short` | `min_sample_unique_reads`, `min_samples`, `min_total_unique_reads`, `upstream_star_filters`; for junctions taken from BAMs also `min_unique`, `min_total`, `min_overhang`, `star_outSJfilter` |
| `rna_long` | `min_intron_length`, `single_exon_alignments` (clustered with a tolerance, or excluded), `min_reads`, `min_samples`, `alignment_filters` |
| `proteins` | `max_rank`, `min_pct_identity`, `min_pct_coverage` |
| `coverage` | `quantisation`, `window_bp`, `strand`, `unique_only`, `multimappers` |

Consumers MUST keep three states apart:

1. **Support** — a record matches the structure (for an intron: identical first and last base).
2. **Absence of retained evidence** — the class exists but has no matching record. This is *not* evidence against the model:
   the record may have been filtered, the tissue may not have been sampled, expression may be too low, or the class may
   never have been able to observe the structure (junction tables cannot test transcript ends).
3. **Competing evidence** — a different record shares exactly one splice site with the structure (alternative
   acceptor/donor, exon skipping), or evidence spans two annotated genes. Report it with its support; do not turn it
   into a verdict.

Absence becomes *informative* only with supporting context, for example deep coverage of the flanking exons.

## 9. Conformance, versioning, size

* A bundle **conforms** if `scripts/validate_bundle.py` exits 0. It checks the schema, file presence, sizes and
  checksums, sequence names and bounds, per-format structure and internal consistency (junction counts, BED12 blocks
  against `intron_chain`, protein blocks against their parents), and sortedness.
* `spec_version` is `major.minor`. While the major version is 0, minor changes may break compatibility.
* **Size.** The retained text tables are small (on real data: 3 105 junctions from six public RNA-seq libraries are 38 KB, against 124 KB of
  gzipped per-sample tables, and the gap grows with the number of libraries).
  Coverage dominates; merge and quantise it. Aim for a bundle within tens of times the size of the annotation.

## 10. Making a bundle

The whole workflow is a handful of commands. [`examples/pipeline_template.sh`](../examples/pipeline_template.sh) is a parametrised, tested version of it
(read it top to bottom; every step is optional), [`examples/minimal/build.sh`](../examples/minimal/build.sh) builds the synthetic example, and
[USER_GUIDE.md](USER_GUIDE.md) explains inputs, outputs and how to look at the result:

```bash
make_manifest.py sequences --fasta genome.fa -o bundle/assembly_sequences.tsv
bam_to_junctions.py lib.bam --fasta genome.fa -o star/lib.SJ.out.tab                # per library; skip if you have STAR's SJ.out.tab
aggregate_star_junctions.py --samples samples.tsv -o bundle/evidence/rna_junctions.tsv.gz --summary-json rna.json star/*.SJ.out.tab
samtools view -b -F 0x900 -q 20 lib.bam | bedtools bamtobed -bed12 -i - > lib.bed12   # per long-read library
collapse_long_reads.py --samples samples.tsv -o bundle/evidence/long_read_models.bed.gz --summary-json lr.json *.bed12
miniprot -It8 --gff genome.fa proteins.faa | miniprot_to_evidence.py --db UniRef90 --aligner-version 0.18 \
    -o bundle/evidence/protein_alignments.gff3.gz --summary-json prot.json
make_coverage.py --sequences bundle/assembly_sequences.tsv -o bundle/evidence/coverage.bw --summary-json cov.json rna_*.bam    # optional
make_manifest.py build bundle/ --template provenance.yaml --fasta genome.fa \
    --summary rna_short=rna.json --summary rna_long=lr.json --summary proteins=prot.json --summary coverage=cov.json
validate_bundle.py bundle/
export_browser_tracks.py bundle/ bundle/browser --genome genome.fa                    # optional: IGV session and browser-friendly tracks
```

EGAPx note: EGAPx runs STAR with `--outSJtype None`, so it never writes `SJ.out.tab`; use `bam_to_junctions.py` on its STAR BAMs (in the
work directory). Its sequence names carry an `lcl|` prefix (BAMs, PAF, Gnomon report) that the bundle must not: pass `--strip-prefix 'lcl|'`.
