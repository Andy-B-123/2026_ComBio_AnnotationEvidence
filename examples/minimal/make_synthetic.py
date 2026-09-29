#!/usr/bin/env python3
"""Generate the raw inputs of the minimal synthetic example (deterministic; standard library only).

Nothing here is real biology. A random genome is built so that the annotated gene structures are
internally consistent (splice motifs at intron ends, reading frames through exons), and evidence
files are written *as if* produced by common tools. Only ``miniprot.gff`` can be real tool output:
run with --run-miniprot to align the synthetic proteins to the synthetic genome.

Story (see README.md):
  geneA + geneB  annotated separately, both partial (A lacks a stop, B lacks a start), but one
                 long-read model, one protein alignment and some RNA-seq junctions span both.
                 => a POSSIBLE split gene; the evidence is suggestive, not conclusive.
  geneC          minus strand; RNA-seq supports intron 1 but not intron 2; protein, long reads
                 and AUGUSTUS support both.
  geneD          annotation follows the ab initio predictors; RNA-seq favours an alternative
                 acceptor 10 bp upstream of the annotated one.
  geneE          a single-exon gene: no junction can exist, so it is supported by single-exon
                 long-read models (20 reads, 3 libraries), coverage and one predictor.
  scaffold_1:20000-25400  four unannotated single-exon long-read loci (3 reads each).
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import evidence_common as ec  # noqa: E402

SEQ_LEN = {"scaffold_1": 30000, "scaffold_2": 12000}

# gene models, 1-based closed coordinates. CDS blocks listed in ascending genomic order.
GENES = [
    dict(id="geneA", seq="scaffold_1", strand="+", partial="3-prime",
         exons=[(1001, 1300), (1701, 1900), (2501, 2700)], cds=[(1101, 1300), (1701, 1900), (2501, 2700)]),
    dict(id="geneB", seq="scaffold_1", strand="+", partial="5-prime",
         exons=[(3401, 3600), (4001, 4300), (4801, 5200)], cds=[(3401, 3600), (4001, 4300), (4801, 5101)]),
    dict(id="geneC", seq="scaffold_2", strand="-", partial=None,
         exons=[(1001, 1400), (2001, 2300), (3201, 3500)], cds=[(1041, 1400), (2001, 2300), (3201, 3449)]),
    dict(id="geneD", seq="scaffold_2", strand="+", partial=None,
         exons=[(7001, 7300), (7801, 8200)], cds=[(7001, 7300), (7801, 8100)]),
    dict(id="geneE", seq="scaffold_1", strand="+", partial=None,          # single exon: no junction can exist
         exons=[(8001, 8900)], cds=[(8001, 8900)]),
]

# RNA-seq samples and groups (tissues); long-read samples
SHORT = [("S1", "head"), ("S2", "head"), ("S3", "head"), ("S4", "gut"), ("S5", "gut"),
         ("S6", "gonad"), ("S7", "gonad"), ("S8", "gonad")]
LONG = [("L1", "head"), ("L2", "gut"), ("L3", "gonad")]

# junctions: (seqid, intron_start, intron_end, strand) -> {sample: (unique, multi)}
def _all(counts, multi=0):
    return {s: (c, multi if i % 3 == 0 else 0) for i, ((s, _), c) in enumerate(zip(SHORT, counts))}

JUNCTIONS = {
    ("scaffold_1", 1301, 1700, "+"): _all([22, 31, 18, 40, 25, 33, 29, 27], 2),      # A intron 1
    ("scaffold_1", 1901, 2500, "+"): _all([15, 20, 12, 28, 19, 22, 17, 21], 1),      # A intron 2
    ("scaffold_1", 1901, 2490, "+"): {"S2": (1, 0), "S5": (2, 0), "S6": (1, 0)},      # weak alternative acceptor in A
    ("scaffold_1", 2701, 3400, "+"): {"S1": (2, 0), "S3": (3, 0), "S5": (1, 0), "S7": (4, 0)},   # A -> B cross-gene junction
    ("scaffold_1", 3601, 4000, "+"): _all([30, 42, 25, 51, 36, 44, 38, 35], 3),      # B intron 1
    ("scaffold_1", 4301, 4800, "+"): _all([28, 39, 24, 47, 33, 41, 35, 32], 0),      # B intron 2
    ("scaffold_1", 3601, 4800, "+"): {"S4": (100, 0)},                                # B exon-skipping, one library only
    ("scaffold_1", 15000, 15500, "."): {"S6": (3, 1)},                                # isolated non-canonical junction
    ("scaffold_2", 1401, 2000, "-"): _all([18, 22, 15, 26, 20, 24, 19, 21], 0),      # C intron 1 (minus strand)
    ("scaffold_2", 7301, 7790, "+"): {"S1": (9, 0), "S2": (12, 0), "S3": (8, 1), "S4": (15, 0), "S5": (11, 0), "S6": (7, 0)},  # D competing acceptor
}
# long-read models: (seqid, strand, exons) -> {sample: reads}
LR_MODELS = [
    ("scaffold_1", "+", [(1001, 1300), (1701, 1900), (2501, 2700)], {"L1": 30, "L2": 20, "L3": 25}),          # A only
    ("scaffold_1", "+", [(3401, 3600), (4001, 4300), (4801, 5200)], {"L1": 40, "L2": 35, "L3": 30}),          # B only
    ("scaffold_1", "+", [(1001, 1300), (1701, 1900), (2501, 2700), (3401, 3600), (4001, 4300), (4801, 5200)],
     {"L1": 6, "L2": 3}),                                                                                     # A+B in one read
    ("scaffold_1", "+", [(3401, 3600), (4801, 5200)], {"L2": 8}),                                             # B skipping exon 2
    ("scaffold_2", "-", [(1001, 1400), (2001, 2300), (3201, 3500)], {"L1": 12, "L2": 10, "L3": 9}),          # C
    ("scaffold_1", "+", [(8001, 8900)], {"L1": 8, "L2": 7, "L3": 5}),                                         # E: single-exon gene
]
# coverage (sum over the 8 short-read libraries): (seqid, start, end, depth) 1-based closed
COVERAGE = [
    ("scaffold_1", 1001, 1300, 800), ("scaffold_1", 1301, 1700, 8), ("scaffold_1", 1701, 1900, 700), ("scaffold_1", 1901, 2500, 7),
    ("scaffold_1", 2501, 2700, 750), ("scaffold_1", 2701, 3400, 14), ("scaffold_1", 3401, 3600, 900), ("scaffold_1", 3601, 4000, 9),
    ("scaffold_1", 4001, 4300, 850), ("scaffold_1", 4301, 4800, 8), ("scaffold_1", 4801, 5200, 900),
    ("scaffold_2", 1001, 1400, 600), ("scaffold_2", 1401, 2000, 5), ("scaffold_2", 2001, 2300, 550), ("scaffold_2", 2301, 3200, 4),
    ("scaffold_2", 3201, 3500, 650), ("scaffold_2", 7001, 7300, 40), ("scaffold_2", 7301, 7790, 1), ("scaffold_2", 7791, 8200, 45),   # coverage follows the RNA-supported acceptor (7790), not the annotated one (7800)
    ("scaffold_1", 8001, 8900, 120),
]

CODE_B1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
CODE_B2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
CODE_B3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
CODE_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: aa for a, b, c, aa in zip(CODE_B1, CODE_B2, CODE_B3, CODE_AA)}
SENSE_CODONS = sorted(c for c, aa in CODON.items() if aa != "*")
AMINO = "ACDEFGHIKLMNPQRSTVWY"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(dna: str) -> str:
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna) - len(dna) % 3, 3))


def place(genome, seqid, strand, blocks_asc, sense):
    """Write the spliced sense sequence into the genome across blocks (reverse-complemented on '-')."""
    pos = 0
    ordered = blocks_asc if strand == "+" else list(reversed(blocks_asc))
    for s, e in ordered:
        n = e - s + 1
        seg = sense[pos:pos + n]
        genome[seqid][s - 1:e] = list(seg if strand == "+" else revcomp(seg))
        pos += n
    assert pos == len(sense)


def set_site(genome, seqid, start, end, strand):
    """Impose splice dinucleotides on the intron [start, end] (GT..AG, or CT..AC seen on the + strand of a - gene)."""
    donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
    genome[seqid][start - 1:start + 1] = list(donor)
    genome[seqid][end - 2:end] = list(acceptor)


def make_cds(rng, nt_len, start_codon, stop_codon):
    codons = [rng.choice(SENSE_CODONS) for _ in range(nt_len // 3)]
    if start_codon:
        codons[0] = "ATG"
    if stop_codon:
        codons[-1] = "TAA"
    return "".join(codons)


def build_genome(rng):
    genome = {s: [rng.choice("ACGT") for _ in range(n)] for s, n in SEQ_LEN.items()}
    rng_e = random.Random(20260926)
    cds_seq = {}
    for g in GENES:
        n = sum(e - s + 1 for s, e in g["cds"])
        assert n % 3 == 0, (g["id"], n)
        source = rng_e if g["id"] == "geneE" else rng                  # geneE has its own stream: adding it must not change the others
        cds_seq[g["id"]] = make_cds(source, n, start_codon=g["partial"] != "5-prime", stop_codon=g["partial"] != "3-prime")
        place(genome, g["seq"], g["strand"], g["cds"], cds_seq[g["id"]])
        for s, e in ec.introns_from_exons(g["exons"]):
            set_site(genome, g["seq"], s, e, g["strand"])
    # extra splice sites for the alternative / cross-gene junctions
    for seqid, s, e, strand in JUNCTIONS:
        if strand in "+-":
            set_site(genome, seqid, s, e, strand)
    return {k: "".join(v) for k, v in genome.items()}, cds_seq


def mutate(rng, prot, rate):
    p = list(prot)
    for i in range(len(p)):
        if rng.random() < rate:
            p[i] = rng.choice(AMINO.replace(p[i], ""))
    return "".join(p)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def wrap(seq, n=60):
    return "\n".join(seq[i:i + n] for i in range(0, len(seq), n))


def gff_lines_for_gene(g, source, gid, tid, score=".", extra="", cds_phase=True):
    strand = g["strand"]
    ex, cds = g["exons"], g["cds"]
    gs, ge = ex[0][0], ex[-1][1]
    partial = ";partial=true" if g.get("partial") else ""
    lines = [f"{g['seq']}\t{source}\tgene\t{gs}\t{ge}\t{score}\t{strand}\t.\tID={gid}{extra}{partial}",
             f"{g['seq']}\t{source}\tmRNA\t{gs}\t{ge}\t{score}\t{strand}\t.\tID={tid};Parent={gid}{extra}{partial}"]
    for s, e in ex:
        lines.append(f"{g['seq']}\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={tid}")
    order = cds if strand == "+" else list(reversed(cds))
    done = 0
    phases = {}
    for s, e in order:
        phases[(s, e)] = (3 - done % 3) % 3
        done += e - s + 1
    for s, e in cds:
        lines.append(f"{g['seq']}\t{source}\tCDS\t{s}\t{e}\t.\t{strand}\t{phases[(s, e)]}\tParent={tid}")
    return lines


def annotation_gff3():
    out = ["##gff-version 3"] + [f"##sequence-region {s} 1 {n}" for s, n in SEQ_LEN.items()]
    rows = []
    for g in GENES:
        rows += gff_lines_for_gene(g, "synthetic", g["id"], g["id"] + ".t1", extra=f";Name={g['id']}")
    return "\n".join(out + rows) + "\n"


def predictions_gff3():
    by = {g["id"]: g for g in GENES}
    merged = dict(id="AB", seq="scaffold_1", strand="+", partial=None,
                  exons=by["geneA"]["exons"] + by["geneB"]["exons"], cds=by["geneA"]["cds"] + by["geneB"]["cds"])
    rows = []
    rows += gff_lines_for_gene(merged, "AUGUSTUS", "augustus_g1", "augustus_g1.t1", score="0.62")
    rows += gff_lines_for_gene(by["geneC"], "AUGUSTUS", "augustus_g2", "augustus_g2.t1", score="0.91")
    rows += gff_lines_for_gene(by["geneD"], "AUGUSTUS", "augustus_g3", "augustus_g3.t1", score="0.55")
    rows += gff_lines_for_gene(by["geneA"], "Helixer", "helixer_g1", "helixer_g1.t1", score="0.88")
    rows += gff_lines_for_gene(by["geneB"], "Helixer", "helixer_g2", "helixer_g2.t1", score="0.84")
    rows += gff_lines_for_gene(by["geneD"], "Helixer", "helixer_g3", "helixer_g3.t1", score="0.71")
    rows += gff_lines_for_gene(by["geneE"], "Helixer", "helixer_g4", "helixer_g4.t1", score="0.79")
    keyed = sorted(rows, key=lambda r: (r.split("\t")[0], int(r.split("\t")[3]), 0 if r.split("\t")[2] == "gene" else 1, r.split("\t")[1]))
    header = ["##gff-version 3"] + ec.pragma_lines(profile="predictions", version=ec.SPEC_VERSION, coordinates="1-based-closed",
                                                  note="hand-authored for the synthetic example; no predictor was run")
    return "\n".join(header + keyed) + "\n"


def star_files(rng):
    files = {}
    for sid, _ in SHORT:
        lines = []
        for (seqid, s, e, strand), per in sorted(JUNCTIONS.items()):
            if sid not in per:
                continue
            uniq, multi = per[sid]
            code = {"+": (1, 1), "-": (2, 2), ".": (0, 0)}[strand]
            lines.append((seqid, s, e, code[0], code[1], 0, uniq, multi, rng.randint(24, 70)))
        lines.sort(key=lambda r: (r[0], r[1], r[2]))
        files[sid] = "".join("\t".join(map(str, r)) + "\n" for r in lines)
    return files


def long_read_files(rng):
    """Per-read BED12 as written by `bedtools bamtobed -bed12` (one line per alignment)."""
    out = {sid: [] for sid, _ in LONG}
    n = 0
    for seqid, strand, exons, per in LR_MODELS:
        for sid, reads in per.items():
            for _ in range(reads):
                ex = [list(x) for x in exons]
                r = rng.random()
                if r > 0.7:                           # 5' (or 3') truncated / extended reads
                    j = rng.choice([10, 25, 40] if len(exons) == 1 else [30, 60, 100])
                    if strand == "+" or rng.random() < 0.5:
                        ex[0][0] = min(ex[0][1] - 20, ex[0][0] + j)
                    else:
                        ex[-1][1] = max(ex[-1][0] + 20, ex[-1][1] - j)
                cs, ce, sizes, starts = ec.bed12_from_exons([tuple(x) for x in ex])
                n += 1
                out[sid].append("\t".join(map(str, (seqid, cs, ce, f"read{n:05d}/ccs", 60, strand, cs, ce, "0,0,0", len(sizes),
                                                    ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ",")))
                                )
    # four unannotated single-exon loci, one read per library each (become single-exon models)
    for sid, _ in LONG:
        for k in range(4):
            s = 20000 + 1500 * k
            out[sid].append("\t".join(map(str, ("scaffold_1", s, s + 900, f"single{sid}{k}/ccs", 60, "+", s, s + 900, "0,0,0", 1, "900,", "0,"))))
    return {sid: "\n".join(v) + "\n" for sid, v in out.items()}


def coverage_bedgraph():
    import math
    rows = []
    for seqid, s, e, d in COVERAGE:
        v = 2 ** (round(4 * math.log2(1 + d)) / 4) - 1        # depth snapped to the ladder 2^(k/4) - 1
        rows.append(f"{seqid}\t{s - 1}\t{e}\t{v:.2f}")
    return "\n".join(rows) + "\n"


def samples_tsv():
    hdr = ["sample_id", "evidence_class", "group", "description", "run_accession", "biosample", "platform", "layout",
           "strandedness", "n_reads_input", "n_reads_uniquely_mapped", "notes"]
    rows = [hdr]
    for i, (sid, grp) in enumerate(SHORT, 1):
        rows.append([sid, "rna_short", grp, f"synthetic short-read library {i}", ".", ".", "ILLUMINA", "paired",
                     "unstranded", ".", ".", "synthetic"])
    for i, (sid, grp) in enumerate(LONG, 1):
        rows.append([sid, "rna_long", grp, f"synthetic long-read library {i}", ".", ".", "PACBIO", "single",
                     "stranded", ".", ".", "synthetic; Iso-Seq/Kinnex-like, refined full-length reads"])
    return "#" + "\t".join(rows[0]) + "\n" + "".join("\t".join(r) + "\n" for r in rows[1:])


def repeats_bed():
    return "scaffold_1\t6000\t6300\tLINE/L1\t0\t+\nscaffold_1\t22000\t22150\tSimple_repeat\t0\t.\nscaffold_2\t9500\t9800\tDNA/hAT\t0\t-\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", default=os.path.join(HERE, "inputs"))
    ap.add_argument("--run-miniprot", action="store_true", help="align proteins.faa to genome.fa with miniprot --gff (needs miniprot)")
    args = ap.parse_args()
    rng = random.Random(20260925)
    genome, cds = build_genome(rng)
    o = args.outdir
    write(f"{o}/genome.fa", "".join(f">{s}\n{wrap(q)}\n" for s, q in genome.items()))
    write(f"{o}/annotation.gff3", annotation_gff3())
    ab = translate(cds["geneA"] + cds["geneB"]).rstrip("*")
    prot = {
        "PROT_AB_ortholog": (mutate(rng, ab, 0.12), "fused ortholog spanning geneA and geneB (synthetic)"),
        "PROT_A_like": (mutate(rng, translate(cds["geneA"]), 0.10), "N-terminal domain protein, geneA only (synthetic)"),
        "PROT_B_like": (mutate(rng, translate(cds["geneB"]).rstrip("*"), 0.10), "C-terminal domain protein, geneB only (synthetic)"),
        "PROT_C_homolog": (mutate(rng, translate(cds["geneC"]).rstrip("*"), 0.08), "geneC homolog (synthetic)"),
    }
    write(f"{o}/proteins.faa", "".join(f">{k} {d}\n{wrap(p)}\n" for k, (p, d) in prot.items()))
    for sid, text in star_files(rng).items():
        write(f"{o}/star/{sid}.SJ.out.tab", text)
    for sid, text in long_read_files(rng).items():
        write(f"{o}/long_reads/{sid}.bed12", text)
    write(f"{o}/predictions.gff3", predictions_gff3())
    write(f"{o}/coverage.bedgraph", coverage_bedgraph())
    write(f"{o}/repeats.bed", repeats_bed())
    write(f"{o}/samples.tsv", samples_tsv())
    if args.run_miniprot:
        res = subprocess.run(["miniprot", "--gff", f"{o}/genome.fa", f"{o}/proteins.faa"], capture_output=True, text=True, check=True)
        write(f"{o}/miniprot.gff", res.stdout)
        print(res.stderr.strip().splitlines()[-1])
    print(f"wrote synthetic inputs to {o}")


if __name__ == "__main__":
    main()
