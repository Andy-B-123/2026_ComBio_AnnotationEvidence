---
name: annotation-evidence
description: Explain the evidence behind a genome annotation from an annotation evidence bundle (splice junctions, long-read models, protein and ab initio alignments), compare competing gene models against it, or build a compact bundle from BAMs, protein alignments and a GFF3. Use when asked why a gene has its exon structure, whether introns or boundaries are supported, whether two genes may be one, or to package annotation evidence.
---

# Annotation evidence

Prototype instructions; not tied to any one LLM. Everything is defined in the repository:
`docs/SPECIFICATION.md` (formats), `docs/LLM_GUIDE.md` (how to interpret; **read it first**), `schema/` (manifest, junction table),
`scripts/` (tools). A bundle is a directory with `manifest.yaml`. Keep this file short and follow those documents rather than restating them.

## The rule

Always keep **support** (an exact matching record), **absence of retained evidence** (class present, no record; never a finding
that a model is wrong) and **competing evidence** (a different record shares one splice site, or evidence spans two genes) apart.
Report numbers exactly as a tool printed them. Say which evidence classes exist and what was filtered (manifest `filters`).
Coordinates: GFF3 and junction tables are 1-based closed; BED and bigWig are 0-based half-open; introns are first..last intron base.

## 1. "Examine this bundle and explain the evidence for gene ABC"

1. `python scripts/validate_bundle.py BUNDLE` (note any WARN), then read `BUNDLE/manifest.yaml`: classes present, `n_datasets`, `filters`, `guided_by_annotation`.
2. `python scripts/summarise_gene.py BUNDLE ABC` (add `--format json` for structured use).
3. Answer with the template at the end of `docs/LLM_GUIDE.md`: supported introns (which classes, how many samples/groups/reads), introns with no
   retained evidence, competing junctions, terminal boundaries, and what the bundle cannot tell you.
4. For a region rather than a gene: `python scripts/query_region.py BUNDLE seqid:start-end`.

## 2. "Produce a compact bundle from these BAMs, protein alignments and GFF3"

Follow `docs/USER_GUIDE.md`; the recipe is `examples/pipeline_template.sh` (set its variables, or run its steps by hand).
Short-read RNA: STAR `SJ.out.tab`, or coordinate-sorted STAR BAMs → `bam_to_junctions.py` (EGAPx never writes `SJ.out.tab`; use its STAR BAMs and
`--strip-prefix 'lcl|'`) → `aggregate_star_junctions.py`. Long reads: filtered BAM → `bedtools bamtobed -bed12` → `collapse_long_reads.py` (single-exon
reads are kept as separate, counted models). Proteins: `miniprot --gff` → `miniprot_to_evidence.py`. Optional coverage: `make_coverage.py`. Then write a
provenance template from `schema/manifest.template.yaml`, run `make_manifest.py sequences` and `build`, and finish with `validate_bundle.py` (must exit 0).
To look at the result: `export_browser_tracks.py`, then open the IGV session (`scripts/igv_snapshot.sh` renders a PNG headless).
Ask the user for what a script cannot know (data sources, tool versions, parameters, why a class is absent); never invent provenance, and
record every filter you applied. Do not put private sequence names (`lcl|…`) in the bundle.

## 3. "Compare two competing gene models against the retained evidence"

1. Put the competing model in its own GFF3 (gene, mRNA, exon, CDS), e.g. `examples/minimal/competing_model_AB.gff3`.
2. Run `summarise_gene.py BUNDLE MODEL_A` and `summarise_gene.py BUNDLE MODEL_B --annotation other_models.gff3` (or one run per model file).
3. Compare, per model: which introns have exact support in which class; competing junctions; the whole-chain support
   (`full_chain`: how many long-read models have exactly this chain, which proteins contain all its introns); cross-gene evidence
   and the counts of evidence that touches only one gene.
4. State where the models differ and how strong the evidence is for each difference, in support / absence / competing terms.
   **Do not declare a winner** unless the difference is a matter of counted evidence you can quote; say what further data would decide it
   (sequence, more tissues, curation).

## Never

* infer that a model is wrong from a missing record; * describe a class with `present: false`; * do BED12 block arithmetic by hand;
* quote a number you did not read from a file or tool output.
