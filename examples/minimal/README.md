# Minimal synthetic example

A tiny, deterministic evidence bundle: two scaffolds (30 kb and 12 kb), five genes, eight short-read libraries in three tissues, three long-read
libraries, four proteins, two ab initio predictors. **Everything is synthetic** — no reads were sequenced. The sequence is random apart from
what the gene models require (splice motifs at every intron end; reading frames that run through the exons), and the protein alignments are
genuine `miniprot 0.18` output on that synthetic genome. Total bundle size: about 45 KB, of which 26 KB is the fixed overhead of one tiny bigWig.

```text
inputs/            raw synthetic inputs (genome, STAR-format junction tables, per-read BED12, proteins, miniprot GFF3, ...)
make_synthetic.py  writes inputs/ deterministically (--run-miniprot re-runs miniprot)
build.sh           inputs/ -> bundle/ using the repository scripts (the whole "workflow")
manifest.template.yaml   hand-written provenance; entries for programs that were NOT run say so
bundle/            the built bundle (validated; a test checks that build.sh reproduces it)
competing_model_AB.gff3  a merged geneA+geneB model, to score against the same evidence
```

## The five genes

| Gene | Situation | What the bundle shows |
|---|---|---|
| **geneA** (+) and **geneB** (+) | annotated separately; A has no stop codon and B no start codon (both `partial=true`) | one long-read model (9 reads, 2 libraries), one protein alignment (466 aa; aa 1–200 on A, 201–466 on B), a junction 2701–3400 seen in 4 libraries of 3 tissues (10 unique reads), and the AUGUSTUS prediction all span both. **But** 75 of the 84 long reads over geneA and 113 of the 122 over geneB touch only one gene, two other proteins align to each gene separately, and Helixer predicts two genes. → *suggestive, not conclusive* |
| geneB, intron 1 | a junction 3601–4800 skips exon 2 | **100 unique reads, all in one library** (gut), against 301 reads in 8 libraries for the annotated intron. Single-library support is not broad support |
| **geneC** (−) | minus strand | intron 2 has **no RNA junction** but exact long-read, protein and AUGUSTUS support: absence in one class, support in others, no verdict |
| **geneD** (+) | annotation follows two ab initio predictors | the annotated intron 7301–7800 has no RNA, long-read or protein support (`prediction_only`), while a junction with an acceptor 10 bp upstream (7301–7790) has **62 unique reads in 6 libraries**: competing evidence |
| **geneE** (+) | a single-exon gene | no junction can exist; supported by **one single-exon long-read model** (20 reads, 3 libraries), coverage and one predictor. Four other single-exon long-read loci (3 reads each) are unannotated |

## Commands

```bash
python scripts/validate_bundle.py examples/minimal/bundle
python scripts/query_region.py    examples/minimal/bundle scaffold_1:2500-3500
python scripts/summarise_gene.py  examples/minimal/bundle geneA          # the split-gene candidate
python scripts/summarise_gene.py  examples/minimal/bundle geneD          # prediction_only + COMPETING EVIDENCE
python scripts/summarise_gene.py  examples/minimal/bundle geneC          # minus strand; RNA absent for intron 2
python scripts/summarise_gene.py  examples/minimal/bundle geneE          # single-exon gene: single-exon long-read evidence

# score a competing model that is not in the annotation against the same evidence
python scripts/summarise_gene.py  examples/minimal/bundle geneAB --annotation examples/minimal/competing_model_AB.gff3
```

The last command shows the merged model's five introns are all individually supported, but the whole chain is seen on only 9 long reads and one
protein (against 75 and 105 reads for the two separate models). That contrast, not any single number, is what a curator has to weigh.

## Looking at it in IGV

```bash
python scripts/export_browser_tracks.py examples/minimal/bundle browser --genome examples/minimal/inputs/genome.fa --locus scaffold_1:900-5400
IGV=/path/to/IGV_Linux_2.19.8 IGV_EXPAND="annotation,long_read_models.bed.gz,protein alignments (BED12 view),predictions.gff3.gz" \
    scripts/igv_snapshot.sh browser/igv_session.xml scaffold_1:900-5400 figure.png       # or open browser/igv_session.xml in IGV
```

The `.fai` index of `genome.fa` is needed by IGV (`samtools faidx`). Pictures of this example are in [`docs/img/`](../../docs/img) and in the user guide.

## Rebuilding

```bash
python examples/minimal/make_synthetic.py                 # add --run-miniprot to re-run the aligner (needs miniprot on PATH)
bash   examples/minimal/build.sh /tmp/rebuilt             # needs bgzip + tabix; uses bedGraphToBigWig if present, else inputs/coverage.prebuilt.bw
```

`make_synthetic.py` is the place to see exactly how each situation above was constructed. The predictions, junction tables, long reads and coverage
are hand-authored to look like tool output; the manifest says which programs were not run.
