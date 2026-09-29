# `rna_junctions` — aggregated splice-junction table (spec 0.1)

The one format in this project that is not an existing standard: a plain TSV that collapses per-sample junction
tables (STAR `SJ.out.tab`) across samples, one row per intron. Everything else about it follows existing conventions.

## Coordinates

**1-based, closed. `intron_start` is the first intron base and `intron_end` the last.** For exons `[1001,1300]` and
`[1701,1900]` the intron is `1301 … 1700`. This equals STAR's `SJ.out.tab` columns 2–3 and a GFF3-derived intron, so
a junction row is compared to a gene model with plain equality. To turn a row into BED, use `chromStart = intron_start − 1`,
`chromEnd = intron_end`. Tools that report the flanking *exonic* bases (TopHat/regtools BED junctions) must be converted
on input (`intron_start = exon_end + 1`, `intron_end = next_exon_start − 1`).

## Columns

`#seqid  intron_start  intron_end  strand  motif  sample_count  group_count  total_unique_reads  max_unique_reads  median_unique_reads  total_multimapping_reads  max_overhang`

| # | Column | Type | Definition |
|---|---|---|---|
| 1 | `seqid` | text | as in `assembly_sequences.tsv` |
| 2 | `intron_start` | int ≥ 1 | first intron base |
| 3 | `intron_end` | int ≥ start | last intron base |
| 4 | `strand` | `+ - .` | transcript strand; `.` when STAR gives 0 (non-canonical junction) |
| 5 | `motif` | text | `GT/AG`, `GC/AG`, `AT/AC`, `non-canonical`; **transcript orientation** (a `-` strand junction that reads CT/AC on the forward strand is `GT/AG`) |
| 6 | `sample_count` | int ≥ 1 | number of supporting samples |
| 7 | `group_count` | int or `.` | distinct `group` values (`samples.tsv`) among supporting samples; `.` if no group table was given |
| 8 | `total_unique_reads` | int | uniquely mapping reads, summed over supporting samples |
| 9 | `max_unique_reads` | int | largest single-sample count among supporting samples |
| 10 | `median_unique_reads` | number | median over supporting samples (`2.5` is possible) |
| 11 | `total_multimapping_reads` | int | multi-mapping reads, summed over supporting samples |
| 12 | `max_overhang` | int | largest STAR "maximum spliced alignment overhang" among supporting samples |

**Supporting sample.** A sample supports a junction if it has at least `sample_min_unique_reads` (header; default 1)
uniquely mapping reads on it. A sample that does not support a junction contributes nothing to that row, including its
multimapping reads. With `--min-sample-unique 0` multimapper-only junctions are kept.

Consistency rules (checked by the validator): `sample_count ≥ 1`; `max ≤ total`; `0 ≤ median ≤ max`;
`total ≥ sample_count × sample_min_unique_reads`; `group_count ≤ sample_count`; `sample_count ≤ n_samples`; every row satisfies the
`min_samples` and `min_total_unique_reads` in the header; rows sorted by seqid then start; no duplicate
`(seqid, start, end, strand)`.

## STAR `SJ.out.tab` mapping

| STAR column | Meaning | Becomes |
|---|---|---|
| 1 | chromosome | `seqid` (optionally with a prefix removed) |
| 2, 3 | first, last intron base (1-based) | `intron_start`, `intron_end` (unchanged) |
| 4 | strand 0 undefined, 1 `+`, 2 `-` | `strand` (`.`, `+`, `-`) |
| 5 | motif 0 non-canonical; 1 GT/AG, 2 CT/AC; 3 GC/AG, 4 CT/GC; 5 AT/AC, 6 GT/AT | `motif`: 1,2 → `GT/AG`; 3,4 → `GC/AG`; 5,6 → `AT/AC`; 0 → `non-canonical` |
| 6 | annotated (0/1) | **not carried**: it depends on how the index was built, and in two-pass mode first-pass junctions are reported as annotated. Record `guided_by_annotation` in the manifest instead |
| 7 | unique reads | summed / max / median over supporting samples |
| 8 | multimapping reads | summed over supporting samples |
| 9 | maximum overhang | maximum |

## What STAR has already filtered

`SJ.out.tab` is **not raw**. STAR applies `--outSJfilterOverhangMin` (default `30 12 12 12`), `--outSJfilterCountUniqueMin`
(`3 1 1 1`), `--outSJfilterCountTotalMin` (`3 1 1 1`), `--outSJfilterDistToOtherSJmin` (`10 0 5 10`) and
`--outSJfilterIntronMaxVsReadN` (`50000 100000 200000`), the four values referring to non-canonical, GT/AG, GC/AG and AT/AC motifs;
annotated junctions bypass them. Record the values that were used in `upstream_star_filters`. A junction absent from this
table may have been seen and dropped by STAR before aggregation ever started.

## Choosing thresholds

The aggregator's defaults are deliberately minimal (`min_sample_unique_reads = 1`, `min_samples = 1`,
`min_total_unique_reads = 1`) and always recorded. Real data show why the columns matter: in six public libraries, 39 % of distinct
junctions occurred in exactly one library and 35 % in all six, and 12 % of one sample's `SJ.out.tab` rows had only multimapping reads.
A stricter `--min-samples 2` shrinks the table and removes single-library artefacts, but also removes real tissue-specific
junctions; that is a decision for the person making the bundle, and the manifest records it.

## Worked example

Rows from the synthetic example bundle (`examples/minimal`), eight libraries in three groups:

```text
#seqid      intron_start intron_end strand motif sample_count group_count total_unique_reads max_unique_reads median_unique_reads total_multimapping_reads max_overhang
scaffold_1  3601         4000       +      GT/AG 8            3           301                51               37                  9                        68
scaffold_1  3601         4800       +      GT/AG 1            1           100                100              100                 0                        65
scaffold_1  2701         3400       +      GT/AG 4            3           10                 4                2.5                 0                        46
```

Row 1 is an annotated intron seen in every library. Row 2 (an exon-skipping junction) has 100 reads, all from **one** library. Row 3 is
a low-count junction that joins the last exon of one annotated gene to the first exon of the next; it is seen in four libraries
from all three groups, so it is independently supported although each library contributes only 1–4 reads.
