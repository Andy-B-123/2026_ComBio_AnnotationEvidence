# Guide for LLMs and coding agents

Read this before interpreting an annotation evidence bundle. It is written for a language-model agent, but a
person new to the format will find it the fastest way in. The bundle is useful without an LLM; this guide only
says how to read it without fooling yourself. Format details are in [SPECIFICATION.md](SPECIFICATION.md).

## The one rule

**Keep three things apart, in every sentence you write about a gene model:**

| State | Means | Example wording |
|---|---|---|
| **Support** | a retained record matches the structure | "intron 2 is crossed by 154 uniquely mapping reads in 8 of 8 libraries (3 tissues)" |
| **Absence of retained evidence** | the evidence class exists but has no matching record | "no retained junction matches intron 2; this class cannot tell us whether the intron is real" |
| **Competing evidence** | a *different* record shares one splice site, or evidence spans two genes | "a junction sharing the donor but using an acceptor 10 bp upstream has 62 unique reads in 6 libraries; the annotated acceptor has none" |

**Absence of evidence is not evidence of absence, and it is never a finding that a gene model is wrong.**
"Not observed" and "contradicted" are different claims. Only competing evidence, weighed against the support for the
annotated structure, can suggest a problem, and even then it *suggests*; you have not shown it.

## Step 0 — orient yourself (always)

```bash
python scripts/validate_bundle.py BUNDLE        # is it a valid bundle? read any WARN lines
```

Then read `BUNDLE/manifest.yaml` and note, for each evidence class:

* `present`. A class with `present: false` says **nothing** about the genome. Do not say "no RNA-seq supports X" when there is no RNA-seq class.
* `filters`. Everything a record had to pass to be kept. A junction that is not in the table may have been seen and dropped.
* `n_datasets` and `data_source`. Eight libraries from one tissue cannot support or refute a gene expressed elsewhere.
* `guided_by_annotation`. If `true`, junction support may be partly circular: mapping was told about the annotation.
* `programs[].version` and `parameters` (for example the STAR `--outSJfilter*` values, protein aligner and database).

## Coordinates: the mistake to avoid

* GFF3 files and `rna_junctions.tsv.gz`: **1-based, closed**. BED (`long_read_models.bed.gz`, `repeats.bed.gz`) and bigWig: **0-based, half-open**.
* An intron is its **first and last intron base**. For exons 1001–1300 and 1701–1900 the intron is **1301–1700**.
  A junction row `1301 1700` therefore equals that GFF3 intron by plain equality. Do not add or subtract 1.
* **Never do BED12 arithmetic by hand.** `blockStarts` are relative to `chromStart`. Use the tools (below), or read the
  `intron_chain` column, which already lists introns in 1-based closed coordinates.
* Region strings (`scaffold_1:2500-3500`) are 1-based, closed, for every file type.
* Sequence names must match exactly. If a name looks off (`lcl|`, `chr1` vs `NC_…`), check `assembly_sequences.tsv` aliases.

## Looking at it

`scripts/export_browser_tracks.py` turns a bundle into an IGV session (junction arcs, protein blocks, coverage, long-read models). Pictures make support and
competition easy to see, but they are not evidence: every claim you make must still come from the tool output or the files, with numbers.

## The tools (use them instead of parsing files yourself)

```bash
python scripts/query_region.py BUNDLE scaffold_1:2500-3500            # everything overlapping a region
python scripts/query_region.py BUNDLE --gene geneA --flank 500 --format json
python scripts/summarise_gene.py BUNDLE geneA                          # per-intron evidence, competing junctions, cross-gene evidence
python scripts/summarise_gene.py BUNDLE geneA --format json            # same, machine-readable
```

Both print 1-based closed coordinates and put each class's filters next to its records. They use `tabix` if it is
installed and indexed, and scan otherwise (same results). Plain shell also works if you must:
`tabix evidence/rna_junctions.tsv.gz scaffold_1:2701-3400`, `zcat evidence/long_read_models.bed.gz | awk ...`.

## What each file means

| File | One record is | Read it as |
|---|---|---|
| `annotation.gff3.gz` | gene / transcript / exon / CDS | the structure being tested. `partial=true` marks models annotated without a start or stop |
| `evidence/rna_junctions.tsv.gz` | one **intron** seen in short-read RNA-seq, pooled over libraries | independent splicing observations; see "RNA splice support" |
| `evidence/long_read_models.bed.gz` | one **transcript structure** (intron chain, or a single-exon cluster) seen in long reads | a whole exon–intron structure on individual molecules; see below |
| `evidence/protein_alignments.gff3.gz` | one **protein aligned to the genome**: a `protein_match` parent and `match_part` blocks | homology-based support for exon–intron structure and coding potential |
| `evidence/predictions.gff3.gz` | a gene model from an ab initio / HMM predictor | *not* independent observation: the predictor's opinion, often trained on the same species |
| `evidence/coverage.bw` | merged RNA-seq depth (depth scale, quantised) | where transcription is seen, whatever the splicing |
| `evidence/repeats.bed.gz` | a repeat interval | context (repeat-derived models are a common false positive) |

## RNA splice support (`rna_junctions`)

Columns: `sample_count`, `group_count`, `total_unique_reads`, `max_unique_reads`, `median_unique_reads`,
`total_multimapping_reads`, `max_overhang`, `motif`.

* **Independence is `sample_count` and `group_count`, not read counts.** 100 reads in one library
  (`sample_count=1`, `max_unique_reads=100`) is far weaker evidence for a real, general junction than 20 reads spread over 8 libraries
  in 3 tissues. A single library can carry an artefact, a contaminant or a tissue-specific isoform.
* Compare `max_unique_reads` with `total_unique_reads`: if one library holds most of the reads, the support is not broad.
* `motif` `GT/AG`, `GC/AG`, `AT/AC` are canonical. `non-canonical` junctions are more often misalignments; treat them cautiously and look at overhang.
* A small `max_overhang` means every read crossing the junction had a short anchor on one side, which makes misalignment more likely; STAR itself imposes minimum overhangs for this reason.
* Multimapping reads are not counted as unique support; a junction with many multimappers and few unique reads may sit in a repeat.
* Depth matters: 3 reads is weak in a library of 200 million reads and strong in one of 5 million. `samples.tsv` has read counts when known.
* **Rules of thumb, not thresholds** (adapt to the data): a junction seen in several libraries and groups with unique reads in more than
  one is usually a real splice event; a junction seen once in one library is uncertain either way. The bundle deliberately does not
  hard-code such cut-offs; `filters` shows what the producer used.
* A junction table says **nothing about transcript ends**, and nothing about which junctions occurred on the same molecule.

## Long-read models (`long_read_models`)

* A model is a *collapsed* set of reads with an **identical intron chain**. `read_count`, `sample_count` and `group_count`
  say how many reads, libraries and groups saw it. Start and end are median observed ends, so they are approximate
  (long reads are often 5′-truncated).
* An **exact chain match** with an annotated transcript is the strongest single observation of that structure: introns seen together
  on one molecule. A model whose chain *contains* all annotated introns plus more is compatible but longer (read-through, a fusion or an
  alternative isoform). Say which.
* **Single-exon models** (`intron_chain` `.`, one block) are clusters of long reads with no intron, grouped by similar start and end. They are kept because single-exon
  genes and unspliced transcripts matter, but they are also what fragments, unspliced pre-mRNA and genomic contamination look like. Report them **separately** from
  spliced models: never count them as intron support or as support for an intron chain. They *can* support a single-exon annotated transcript (compare the overlap fractions
  `summarise_gene.py` prints) and can show an unannotated transcribed locus; a spliced model over a locus annotated without introns is competing evidence.
* `min_reads` and `min_samples` may have removed rare models: check `filters` before saying a structure is unseen.

## Protein support (`protein_alignments`)

* Percentages are **percent** (0–100): `pct_identity`, `pct_positive`, `pct_coverage`; `rank` 1 is the best placement of that protein.
* **Exact intron match** (identical first and last intron base) between an alignment and a gene model is good evidence for that splice
  site. Aligners can shift a boundary by a few bases to preserve frame, so also look at *near* matches (`summarise_gene.py` reports both).
* **Protein-side coordinates** (`Target=protein start end` per block) show which part of the protein each exon carries. That is what lets
  you ask whether the N-terminal half aligns to one gene and the C-terminal half to the next.
* Proteins are usually from *other species*, so identity is lower than for same-species transcripts; read `pct_identity` together with
  `pct_coverage` (high identity with low coverage may be a shared domain, not the whole gene) and check `databases` and `taxid`.
* A protein that aligns across two adjacent genes is evidence about the *protein family* having one gene there; it is not proof the
  annotated genes are one. Other proteins may align to each gene separately. Count them.

## Comparing evidence with a gene model

For each intron of the annotated transcript, `summarise_gene.py` reports exact matches per class, the class(es) that were checked and did
not match, and any junction sharing exactly one splice site. Interpret in this order:

1. Which classes are present at all (manifest)? Only those can be cited.
2. For each intron: which classes match exactly? (`supported_by`) Which exist but do not? (`not_observed_in`)
3. Is there **competing** junction evidence? Compare its `sample_count` and `total_unique_reads` with the annotated junction's.
   `competing_evidence: true` means a competitor is at least as well supported as the annotated intron in the same data (thresholds are
   parameters and are printed).
4. Terminal boundaries (start of the first exon, end of the last) can only be tested by long-read model ends, protein alignment ends and
   predictions, not by junctions.
5. Coverage: an unsupported intron next to deeply covered exons is *more* informative than one where the whole locus is uncovered.

## Two genes that might be one (or one that might be two)

Run `summarise_gene.py` on each gene and read "Evidence spanning this gene and a neighbouring annotated gene". Weigh:

| Suggests one gene | Does not |
|---|---|
| a long-read model whose chain covers both genes, with novel connecting intron(s) | most long-read reads cover only one gene (`single_gene` counts) |
| RNA junctions from the last exon of A to the first exon of B, in several libraries and groups | such junctions are rare compared with each gene's own junctions |
| a protein aligned across both, with protein coordinates continuing from A into B | other proteins align to A and to B separately |
| both models are `partial=true` (A lacks a stop, B lacks a start), same strand, adjacent | both are complete, with stop and start codons |
| an ab initio predictor calls one gene | other predictors call two |

Alternative explanations to state, not to dismiss: transcriptional read-through, overlapping UTRs, chimeric alignments,
a bicistronic or polycistronic locus, an assembly error. What the bundle cannot show: whether the fused reading frame is intact
(needs the genome sequence), whether one polypeptide is made, or expression in tissues that were not sampled.
For a **split** hypothesis (one annotated gene that may be two) the same tools apply: look for models and alignments that cover only one part,
for a gap in coverage, and for proteins that align to disjoint parts; but most of this is *absence* of spanning evidence,
so treat it with extra caution.

## What you must not conclude

* "No junction was found, so the intron is wrong / the gene is not real." *No.* See the filters, the tissues, the depth.
* "This isoform is not expressed." Only that it was not retained in these libraries after these filters.
* "The annotation is correct because it matches the evidence." The annotation may have been *built from* the same reads. If
  `guided_by_annotation` is true the agreement is partly circular.
* "The gene should be merged / split." The bundle can make it more or less likely; the decision needs the sequence, more tissues, or curation.
* Percentages or counts you have not read from a file or a tool output. Quote what the tool printed.
* Anything about a class with `present: false`.

## Answer template

> **Gene X** (scaffold:start–end, strand). Evidence classes present: … (n samples/datasets; filters …).
> **Supported:** intron 1 (exact match: 8 samples, 3 groups, 225 unique reads; long-read model LR000001, 84 reads; protein PROT_A_like 91 % identity).
> **Not observed:** intron 3 has no retained junction (class present, 8 samples); long reads and protein match it exactly. This is not evidence against it.
> **Competing:** a junction using an acceptor 10 bp upstream has 62 unique reads in 6 libraries; the annotated acceptor has none in this class.
> **Cannot say:** transcript ends (no long-read models within tolerance), expression outside the sampled tissues.
> **Caveats:** … Nothing here shows the model is right or wrong; it indicates where curation should look.

## A worked example

The synthetic bundle `examples/minimal/bundle` contains four genes built to show different situations
(`examples/minimal/README.md` walks through them):

* `geneA` and `geneB`: separate annotations, both `partial=true`; one long-read model (9 reads), one protein alignment, a junction seen in four libraries and one predictor span both; two other
  proteins and most long reads (75 of 84) cover only geneA, and another predictor calls two genes. → *suggestive, not conclusive*.
* `geneC` (minus strand): intron 2 has no RNA junction but long-read and protein support. → *absence in one class, support in others*.
* `geneD`: annotation follows ab initio predictions; RNA-seq favours an acceptor 10 bp upstream. → *competing evidence*.
* `geneE`: a single-exon gene: no junction can exist, and it is supported by a single-exon long-read model (20 reads, 3 libraries) and one predictor. → *evidence of a different kind, kept apart from spliced evidence*.
