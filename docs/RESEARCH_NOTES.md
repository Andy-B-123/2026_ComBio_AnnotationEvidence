# Research notes

Answers to the ten questions posed at the start of the project. All web sources were fetched and
read on **2026-09-25**; numbered references are listed at the end. Where something could not be
verified it says so. Measured numbers come from public data only and the commands used are given
so they can be repeated on other data.

---

## 1. Are there existing standards or projects that already address this?

**No standard for a compact, multi-class "evidence bundle" was found.** Several projects cover part of it:

| Precedent | What it does | Why it is not the same thing |
|---|---|---|
| AUGUSTUS / BRAKER **hints** [10][11] | GFF-style lines (`intron`, `exonpart`, `CDS`, …) with a `source=` code, made from RNA-seq and protein alignments | Designed to *drive prediction*; no sample provenance, no filter record, no aggregation semantics |
| **EVidenceModeler** [12] | Inputs are `gene_predictions.gff3`, `protein_alignments.gff3`, `transcript_alignments.gff3` plus a weights file | Same: an *input* to gene-model combination. Its file naming matches this project's classes, which is reassuring |
| **GeMoMa** [13] | Takes RNA-seq evidence as intron GFF files with supporting split-read counts (parameter `GeMoMa.r`) | Prediction input; homology-based annotation transfer |
| **NCBI RefSeq / EGAP** [1][2][3] | Publishes evidence alignments, per-run coverage and per-model support tables | Closest precedent; see Q2 and Q3 — but NCBI-specific, partly binary, and not designed as a reusable standard |
| **Snaptron / recount** [9] | Cross-study junction summaries (146 M junctions, >70 000 human samples) with tissue-specificity scoring | Shows cross-sample junction summarisation works at scale; human-only and not tied to an annotation |
| **GA4GH Sequence Annotation** [16] | A study group on interoperability of annotation models | About the *gene model*, not its evidence; not a ratified specification (see Q10) |

Not reviewed in depth: Ensembl/GENCODE transcript support levels, MAKER's AED score, RO-Crate and
BioCompute Objects as possible packaging wrappers (only that the RO-Crate project page exists [17]).

## 2. Which files does NCBI retain or distribute alongside RefSeq annotations?

From the FTP README [3] and a live directory listing for one insect release
(*Helicoverpa armigera*, `GCF_030705265.1-RS_2024_03`) [4]:

| Item | Where | Size in this release |
|---|---|---|
| Annotation (GFF3, GTF, GBFF, feature table) | release root | `genomic.gff.gz` 6.3 MB |
| Gnomon-only models (GFF3 + FASTA) | `Gnomon_models/` | 27.8 MB (three files) |
| **Same-species transcript alignments** (GFF3, `cDNA_match`/`EST_match`/`match` with `Target`, `Gap`, `pct_*`) | `Evidence_alignments/` | 50 MB |
| Alignments of the RefSeq transcripts (BAM) | `RefSeq_transcripts_alignments/` | 18 MB |
| **Per-run RNA-seq coverage** (bigWig) | `RNASeq_coverage_graphs/` | 145 files, **10.3 GB** |
| Run metadata (17 columns: platform, read/base count, alignment rate, spliced-alignment rate, BioSample, BioProject) | `*_rnaseq_runs.txt` | 28 KB |
| Per-gene expression counts | `*_gene_expression_counts.txt.gz` | 3.7 MB |
| RepeatMasker output | README lists `*_rm.out.gz`; not present in this release's listing | – |

The 145 runs are 8.18 billion reads and 1.174 Tbp (121 Illumina, 24 BGISEQ; 140 BioSamples, 16 BioProjects).

What I did **not** find in the release directory: protein-to-genome alignments, or a collapsed
splice-junction table. The README lists `protein_match` as a GFF3 alignment type, but the
evidence directory contains only transcript alignments. The Genome Data Viewer may display more
than the FTP listing shows; I did not check that.

## 3. How does EGAPx summarise or propagate RNA and protein evidence internally?

Source: the EGAPx repository at commit `f9a7392` (release v1.0) [1] and one real EGAPx 1.0 output directory.

* **Short reads.** STAR is run with `--alignSJoverhangMin 8 --outFilterMultimapNmax 200 … --outSJtype None`
  (`rnaseq_short/star_wnode/main.nf`), so **no `SJ.out.tab` is ever written**. SAM alignments are converted to NCBI
  ASN.1 (`convert_from_bam: sam2asn -filter 'pct_identity_gap >= 95' -collapse-identical`), then `rnaseq_collapse` merges alignments
  while recording sample information. EGAPx also forces local sequence ids (`convert_genome: force-local-ids`), which is why
  `lcl|` prefixes appear in its internals.
  NCBI's process page describes this as collapsing alignments "with the same splice structure and the same or
  similar start and end points" while recording samples and read counts, and filtering "very rare introns" [2].
* **Long reads.** minimap2 (`-x splice --eqx --MD --cs …`), then a filter (`matches`, `splices`, `exon_identity >= 0.85`,
  `rank == 1`, …) from `default_task_params.yaml`.
* **Proteins.** `miniprot` by default (`-p 0.4 --outs=0.4`), ProSplign optionally. miniprot writes **PAF**, which `paf2asn`
  converts to ASN.1 (`target_proteins/miniprot`, `paf2asn`).
* **Gene models.** Gnomon chains splice-compatible alignments into models and adds ab initio parts where needed.
  `gnomon_evidence_summary` then scores every model against every line of evidence.
* **What the user receives.** The README documents `GNOMON/new.gnomon_report.txt` (one row per model *per line of
  evidence*: identity, base/CDS coverage, precise/approximate splice-site support, "core support") and
  `GNOMON/new.gnomon_quality_report.txt` (one row per model: minimal supporting alignments, ab initio percentage,
  number of supporting short-read BioSamples, non-consensus introns).

In an EGAPx 1.0 test run with short and long reads, the published output directory held the annotation as GFF3, GTF and ASN.1
(`annotated_genome.asn`), the two Gnomon reports, and a large binary ASN.1 file, `stats/rnaseq_long/align.asn`, whose strings show it holds **miniprot
protein alignments** (`miniprot_score_ms`, `paf2asn:` placement ids, `XP_` query ids) despite the directory name. The `filtered_protein_alignments/`
directory of that run held only long-read alignment statistics.
The evidence is therefore *retained*, but as NCBI-specific binary ASN.1 that ordinary genomics tools cannot query.

**Model-centric versus observation-centric.** The Gnomon reports answer "what supported model X?". They cannot answer
"what other evidence sits at this locus that no model used?" — for example a long read that spans two annotated
genes. That is the gap an observation-centric, region-queryable bundle fills. The two are complementary, and
the Gnomon reports can be listed in a bundle as an optional extra.

Seqids differ between EGAPx internals and its final GFF3: the Gnomon report, BAM headers and PAF use `lcl|NC_051855.1`;
the GFF3 uses `NC_051855.1`.

## 4. Are there existing standard representations of aggregated splice-junction evidence?

There is a de facto per-sample format and several ad hoc aggregates, but **no accepted cross-sample standard**.

* **STAR `SJ.out.tab`** [5]: nine columns — seqid; first and last *intron* base (1-based); strand (0/1/2); motif (0–6);
  annotated flag; unique-read count; multimapper count; maximum overhang. Three properties matter for
  aggregation: (a) STAR "defines the junction start/end as intronic bases, while many other software define them as
  exonic bases"; (b) the file is **already filtered** by `--outSJfilter*` (defaults: non-canonical junctions need ≥3
  unique reads and ≥30 bp overhang; canonical ones ≥1 read and ≥12 bp), and annotated junctions bypass those filters;
  (c) in 2-pass mode, junctions found in the first pass are reported as "annotated". All three must be recorded in the
  manifest, so the aggregator does.
* **regtools** `junctions extract` [14]: BED12 output from any spliced BAM (tested with HISAT2, TopHat2, STAR, minimap2).
* **Snaptron / recount** [9]: cross-sample junction summaries; I could not retrieve the column definitions, so no
  claims are made about them.
* **GeMoMa, AUGUSTUS hints** [10][13]: junction-like evidence as GFF with counts or multiplicities.

**Coordinate check on real data.** STAR junctions from six public libraries were compared with the 7 803 annotated introns
(all transcripts) of the RefSeq annotation on the same chromosome: 6 510 (83.4 %) are identical (same first and last intron base,
no conversion, no strand disagreement in any of them), and **0.0 %** are identical after shifting by one base in any direction or using
the exonic-flank convention. So STAR's `SJ.out.tab` intron and a GFF3-derived intron are the same object.

Real STAR data (six public libraries, one chromosome; Q5): 12 377 rows collapse to 3 711 distinct junctions;
1 440 (39 %) occur in one library and 1 290 (35 %) in all six; 12 % of one sample's rows have zero unique reads
(multimappers only); overhang 12–63 bp; no strand conflicts between samples. These are why the aggregate keeps
sample count, group count, maximum, median and multimapper columns, and why the minimum-unique-read filter defaults
to a recorded, minimal value instead of being hidden.

## 5. How compact is bigWig compared with bedGraph and BAM?

Method (one 15.4 Mb chromosome of *Hermetia illucens*; six public paired-end libraries ERR1801985–ERR1801990,
mapped with STAR 2.7.11b; coverage from `bedtools genomecov -bg -split`, i.e. spliced coverage, all alignments
including multimappers; `bedGraphToBigWig` from UCSC tools; `bgzip` from htslib 1.18):

| Representation | Size | vs 6 BAMs (407.8 MB) |
|---|---|---|
| 6 STAR BAMs (44–127 MB each) | 407.8 MB | 1× |
| Per-sample bedGraph, plain (mean) | 6.8 MB | – |
| Per-sample bedGraph, bgzip (mean) | 1.35 MB | – |
| Per-sample bigWig (mean) | 1.20 MB | 36–103× per sample |
| Merged (sum) bedGraph, plain | 18.0 MB | 23× |
| Merged bedGraph, bgzip | 3.60 MB | 113× |
| Merged bigWig | **3.61 MB** | 113× |
| Merged, quantised to the ladder 2^(k/4)−1 (about 19 % steps): bedGraph.gz / bigWig | 0.74 / 0.95 MB | 550× / 427× |
| Merged, 50-bp windows + quantised: bedGraph.gz / bigWig | 0.29 / 0.48 MB | 1 419× / 844× |

Reading: **at full resolution bigWig is the same size as bgzipped bedGraph** (3.61 vs 3.60 MB), and after quantisation
bgzipped bedGraph is 22–41 % smaller, because bigWig also stores an index and zoom levels. (Quantised values are stored on the depth scale, snapped to the nearest 2^(k/4)−1, so a value still means "about this many reads".) The bigWig advantage is random access,
zoom summaries and native browser support, not size. **Merging and quantising are what save space**, by two to three orders of
magnitude. Caveats: one chromosome, six libraries, all alignments counted; the deepest position had a depth of 453 433, so
extrapolating to a genome is unsafe. For comparison, NCBI publishes one full-resolution bigWig per run (145 files, 54–91 MB in the
ones listed, 10.3 GB total) [4].

## 6. Is BED12 sufficient for long-read structural support?

For **exon–intron structure, yes**: blocks carry it exactly, it is browser-native, tabix-indexable and convertible to
bigBed with an autoSql file [6]. BED requires 0-based, half-open coordinates, blockStarts relative to `chromStart` with the
first equal to 0, and a score in 0–1000 [6]. Limits:

* BED12 has no place for read counts or sample support → extra columns (BED12+4).
* Per-read information (poly-A, individual TSS/TES, read identity) is lost by design. Models are collapsed by
  **intron chain** and report a consensus start/end.
* Single-exon reads have no intron chain and cannot be collapsed that way; they are clustered by similar start and end, counted, and kept apart from spliced models.
* Relative `blockStarts` are error-prone for people and models, hence the explicit `intron_chain` column.

## 7. Is GFF3 the most interoperable representation for protein-to-genome evidence?

**As an interchange format, yes; as a compact query format, BED12+ is simpler — but loses protein coordinates.**

* The GFF3 specification defines `Target` ("target_id start end [strand]") and `Gap` for protein-to-nucleotide
  alignments, states that discontinuous alignments may be several rows or one row plus `Gap`, and reserves only
  attributes that start with a capital letter [7].
* NCBI publishes alignment evidence as GFF3 (`protein_match`, `cDNA_match`, `match`; `Target`, `Gap`, `pct_coverage`,
  `pct_identity_gap`) [3][15], and miniprot emits GFF3 natively (`mRNA`, `CDS`, `stop_codon`, `Rank`, `Identity`,
  `Positive`, `Target`, with `##PAF` lines carrying protein length) [8].
* Sequence Ontology terms: `protein_match` SO:0000349, `match_part` SO:0000039.
* PAF is compact but its protein CIGAR (`M I D F G N U V`) is unfamiliar to most tools, and EGAPx converts it to ASN.1.
* Limit found in testing: **IGV 2.19.8 draws a `protein_match` with `match_part` children as one unbroken bar**, so the block structure is invisible; the project therefore exports a
  BED12 view for browsers while keeping the GFF3 as the canonical record.
* Decisive reason for GFF3: it keeps per-block `Target` ranges, which are needed to see that one protein's N-terminal half
  aligns to one gene and its C-terminal half to the next.

## 8. How should evidence from many samples be collapsed without erasing useful information?

Keep, per junction: **how many** samples support it, **how many groups** (tissues/stages), the **total**, the **maximum** and
the **median** of unique reads over supporting samples, multimapper reads, and the largest overhang. Define
"supporting sample" explicitly (default ≥1 unique read; recorded). Keep `samples.tsv` so counts can be interpreted.

What is lost, and stated as such: *which* samples supported a junction (only counts survive), depth normalised to library size,
and which junctions co-occurred on one read (only long reads retain that). An optional sparse per-sample table would restore the
first at a cost of about one row per (junction, supporting sample) — 12 377 rows against 3 711 here. It is left out of v0.1.

## 9. What metadata are essential for interpretability 5–10 years later?

The manifest schema encodes these: the **assembly identity** (accession, checksum, and the exact sequence names and lengths);
the **annotation** version, date and pipeline; for each evidence class the **data source and accessions**, the ordered
**programs with versions and parameters**, any **reference database and version** (proteins), the **summarisation procedure**
and **every filter**; whether mapping was **guided by an annotation** or two-pass (a circularity risk); the **coordinate
convention and format version** of each file; the **spec version**; **checksums**; **licence and contact**. For samples:
group, layout, strandedness, platform and read counts. Everything else is deliberately left to SRA and the workflow.

## 10. Could this be expressed through GA4GH standards?

Partly — as *building blocks*, not as a ready-made standard [16]:

* **BED v1.0** is a listed GA4GH product ("stores information on the location of genomic features"); the long-read,
  repeat and (via bigBed) other tables build on it.
* **refget Sequences / Sequence Collections** give digests for sequences and for whole reference sets, "to replace or
  complement the human-readable identifiers currently used for reference genomes" [18]. This is the right way to pin
  assembly identity; `assembly_sequences.tsv` carries the per-sequence digest.
* **Experiments Metadata Checklist (Expmeta)** is a reference point for the minimum sample/experiment metadata.
* **htsget, DRS, Service Info** are candidate access-layer standards for a future service (see [FUTURE_MCP.md](FUTURE_MCP.md)).
* **RNAget** targets expression matrices, not junctions.
* **Sequence Annotation (SA)** is a *study group* "evaluating interoperability gaps" in annotation, with no ratified
  specification. The older `ga4gh-schemas` repository (retired 2018-01-24) had a `sequence_annotations.proto` with a
  GFF3/SO-based `Feature` and a bigWig-like `Continuous` message [19].

Conclusion: nothing to adopt wholesale today; nothing conflicting either. SA is the natural place to raise the evidence-bundle
idea.

---

## Limits of this review

Not verified: Ensembl/GENCODE support tags; how JBrowse2 and UCSC load each file (IGV 2.19.8 was tested, see below); Genome Data Viewer tracks; Snaptron column definitions; RO-Crate and BioCompute suitability.
The EGAPx observations are from one repository commit and one run and may not hold for other versions.

## References

1. EGAPx repository, commit `f9a7392`, release v1.0 — <https://github.com/ncbi/egapx> (`README.md`; `nf/subworkflows/ncbi/rnaseq_short/{star_wnode,rnaseq_collapse}`; `nf/subworkflows/ncbi/target_proteins/{miniprot,paf2asn}`; `nf/subworkflows/ncbi/gnomon/gnomon_evidence_summary`; `ui/assets/default_task_params.yaml`)
2. NCBI Eukaryotic Genome Annotation Pipeline, annotation process — <https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process/>
3. NCBI genomes FTP README — <https://ftp.ncbi.nlm.nih.gov/genomes/all/README.txt> ; NCBI GFF3 README — <https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt>
4. Release directory listing, *H. armigera* — <https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/29058/GCF_030705265.1-RS_2024_03/>
5. STAR manual §5.5 and §17.18 — <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf> ; Dobin et al., *Bioinformatics* 29:15–21 (2013), doi:10.1093/bioinformatics/bts635
6. UCSC BED format — <https://genome.ucsc.edu/FAQ/FAQformat.html> ; bigWig — <https://genome.ucsc.edu/goldenPath/help/bigWig.html> ; Kent et al., *Bioinformatics* 26:2204–2207 (2010), doi:10.1093/bioinformatics/btq351
7. GFF3 specification — <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md> ; Sequence Ontology `so-simple.obo`
8. miniprot README and manual page — <https://github.com/lh3/miniprot> ; Li, *Bioinformatics* 39:btad014 (2023), doi:10.1093/bioinformatics/btad014
9. Wilks et al., Snaptron, *Bioinformatics* (2018), doi:10.1093/bioinformatics/btx547 — <https://snaptron.cs.jhu.edu/>
10. AUGUSTUS, "Using hints" — <https://github.com/Gaius-Augustus/Augustus/blob/master/docs/RUNNING-AUGUSTUS.md>
11. BRAKER — <https://github.com/Gaius-Augustus/BRAKER>
12. EVidenceModeler wiki — <https://github.com/EVidenceModeler/EVidenceModeler/wiki>
13. GeMoMa documentation — <https://www.jstacs.de/index.php/GeMoMa-Docs>
14. regtools `junctions extract` — <https://regtools.readthedocs.io/en/latest/commands/junctions-extract/>
15. (see 3) NCBI GFF3 README, "Alignments" section
16. GA4GH products — <https://www.ga4gh.org/our-products/> ; Sequence Annotation — <https://www.ga4gh.org/product/sequence-annotation/>
17. RO-Crate — <https://www.researchobject.org/ro-crate/>
18. GA4GH refget — <https://ga4gh.github.io/refget/> (sequences: `/sequences/`, collections: `/seqcols/`)
19. Retired GA4GH schemas — <https://github.com/ga4gh/ga4gh-schemas> (archived; `sequence_annotations.proto`)
20. tabix — Li, *Bioinformatics* 27:718–719 (2011), doi:10.1093/bioinformatics/btq671 ; <https://samtools.github.io/hts-specs/tabix.pdf>
