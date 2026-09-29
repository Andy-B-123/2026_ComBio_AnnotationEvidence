# Test fixtures

* `real_star/*.SJ.out.tab` — excerpts (a few dozen rows each, chosen to include multimapper-only rows,
  minus-strand rows and non-GT/AG motifs) of STAR 2.7.11b `SJ.out.tab` files from three public
  *Hermetia illucens* paired-end RNA-seq runs (ENA ERR1801985, ERR1801986, ERR1801987), mapped to
  chromosome NC_051855.1 of assembly GCF_905115235.1 with default `--outSJfilter*` settings and no GTF.
* `miniprot_DPP3.gff` — output of `miniprot --gff` 0.18 on miniprot's own public test data
  (`test/DPP3-hs.gen.fa.gz` and `test/DPP3-mm.pep.fa.gz`, https://github.com/lh3/miniprot, MIT licence).
