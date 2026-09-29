# Minimal synthetic evidence bundle

This directory is an **annotation evidence bundle** (spec v0.1, draft): a compact, standard-format
record of the evidence behind an annotation. **All data here are synthetic.** Start with
`manifest.yaml`, which describes everything else.

```
manifest.yaml              provenance: sources, programs, versions, filters, checksums
samples.tsv                one row per library/dataset (group = tissue/stage)
assembly_sequences.tsv     sequence names, lengths, GA4GH refget digests, aliases
annotation.gff3.gz         the annotation the evidence explains (+ .tbi)
evidence/
  rna_junctions.tsv.gz         short-read splice junctions collapsed across samples (+ .tbi)
  long_read_models.bed.gz      long-read transcript models, BED12+4 (+ .tbi)
  protein_alignments.gff3.gz   protein-to-genome alignments (+ .tbi)
  predictions.gff3.gz          ab initio / HMM gene predictions (+ .tbi)
  coverage.bw                  merged RNA-seq coverage (bigWig)
  repeats.bed.gz               repeat intervals, BED6 (+ .tbi)
```

**Coordinates.** GFF3 files and the junction table are **1-based, closed**. BED files and bigWig are
**0-based, half-open**. Junction rows give the first and last *intron* base, which is also how a GFF3
intron is defined. Each file's format is defined in `docs/SPECIFICATION.md` of the repository.

**What absence means.** A missing record is *not* evidence that a gene model is wrong: the evidence
was filtered, sampled and summarised (see `filters` in the manifest), and RNA-seq only sees what was
expressed in the sampled tissues. See `docs/LLM_GUIDE.md`.

Try it (from the repository root):

```bash
python scripts/validate_bundle.py examples/minimal/bundle
python scripts/query_region.py examples/minimal/bundle scaffold_1:2500-3500
python scripts/summarise_gene.py examples/minimal/bundle geneA
```
