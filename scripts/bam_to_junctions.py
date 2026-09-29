#!/usr/bin/env python3
"""Extract splice junctions from a spliced short-read RNA-seq alignment (BAM/CRAM/SAM) in STAR's SJ.out.tab format.

Use this when the STAR run did not write SJ.out.tab (EGAPx runs STAR with --outSJtype None), or when the BAM came from
elsewhere. The output has STAR's nine columns, so ``aggregate_star_junctions.py`` reads it unchanged.

Counting rules (checked against STAR's own SJ.out.tab: identical counts and overhangs on all 1 755 and all 2 218
junctions of two real libraries):
  * unique reads      = distinct READ NAMES with NH:i:1 that cross the junction (mates that both cross count once);
  * multimapping      = distinct read names with NH > 1 that cross it (secondary alignments must be in the BAM);
  * max overhang      = the largest, over reads and mates, of min(left run, right run), where a "run" is the aligned bases
                        (M, =, X) next to the junction, stopping at an insertion, deletion or another intron;
  * strand and motif  = from the genome sequence, as STAR does (needs --fasta): 1/2 GT/AG, 3/4 GC/AG, 5/6 AT/AC on
                        the + / - strand; 0 (strand 0) for anything else;
  * column 6 (annotated) is always 0: the BAM does not know an annotation.

The output is UNFILTERED: STAR applies --outSJfilter* before writing SJ.out.tab, so a BAM-derived table has more low-support
junctions. Filter with --min-unique / --min-total / --min-overhang here, or downstream (aggregate_star_junctions.py
--min-total-unique, --min-samples). Whatever you use, record it in the manifest.

Records that are unmapped or QC-fail (flags 0x4, 0x200) are skipped; everything else counts, including secondary
alignments (that is how multimappers are seen). Sequence names are looked up in the FASTA after --strip-prefix.

Examples:
  bam_to_junctions.py sample.bam --fasta genome.fa -o star/sample.SJ.out.tab
  bam_to_junctions.py egapx.bam --fasta genome.fa --strip-prefix 'lcl|' -o star/sample.SJ.out.tab
  samtools view -h sample.bam chr1 | bam_to_junctions.py - --fasta genome.fa -o chr1.SJ.out.tab       # SAM text on stdin
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
# (donor, acceptor) on the forward strand -> (strand code, motif code)
MOTIFS = {("GT", "AG"): (1, 1), ("CT", "AC"): (2, 2), ("GC", "AG"): (1, 3), ("CT", "GC"): (2, 4),
          ("AT", "AC"): (1, 5), ("GT", "AT"): (2, 6)}


class Fasta:
    """Random access to an (optionally .fai-indexed) FASTA without external libraries."""

    def __init__(self, path: str):
        self.path = path
        self.index: dict[str, tuple[int, int, int, int]] = {}          # name -> (length, offset, bases/line, bytes/line)
        fai = path + ".fai"
        if os.path.exists(fai):
            for line in open(fai):
                f = line.split("\t")
                self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        else:
            self._scan()
        self.fh = open(path, "rb")

    def _scan(self):
        name = None
        length = offset = bases = width = 0
        pos = 0
        with open(self.path, "rb") as fh:
            for raw in fh:
                if raw.startswith(b">"):
                    if name:
                        self.index[name] = (length, offset, bases, width)
                    name, length, offset, bases, width = raw[1:].split()[0].decode(), 0, pos + len(raw), 0, 0
                elif name:
                    if not bases:
                        bases, width = len(raw.rstrip(b"\r\n")), len(raw)
                    length += len(raw.rstrip(b"\r\n"))
                pos += len(raw)
        if name:
            self.index[name] = (length, offset, bases, width)

    def fetch(self, name: str, start: int, end: int) -> str:
        """Bases start..end (1-based, closed), upper-case."""
        length, offset, bases, width = self.index[name]
        if start < 1 or end > length:
            raise ec.EvidenceError(f"{name}:{start}-{end} is outside the sequence (length {length})")
        s0, e0 = start - 1, end
        first = offset + (s0 // bases) * width + s0 % bases
        last = offset + ((e0 - 1) // bases) * width + (e0 - 1) % bases
        self.fh.seek(first)
        return self.fh.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def sam_lines(path: str):
    """Yield SAM text lines from a BAM/CRAM/SAM path, or from stdin for '-'."""
    if path == "-":
        yield from sys.stdin
        return
    if path.endswith(".sam"):
        with open(path) as fh:
            yield from fh
        return
    if not ec.have_tool("samtools"):
        raise ec.EvidenceError("samtools is needed to read BAM/CRAM; alternatively pipe SAM text on stdin ('-') or give a .sam file")
    proc = subprocess.Popen(["samtools", "view", "-h", path], stdout=subprocess.PIPE, text=True)
    yield from proc.stdout
    if proc.wait() != 0:
        raise ec.EvidenceError(f"samtools view failed on {path}")


def read_junctions(lines, strip_prefix: str | None, exclude_flags: int):
    """Yield (chrom, {(start, end): [unique_names, multi_names, max_overhang]}), one table per reference.

    For a coordinate-sorted input (``@HD SO:coordinate``) a reference is emitted as soon as the next one starts, which
    bounds memory; otherwise every reference is kept until the end of the input. Only records with an N in the CIGAR are read.
    """
    sorted_input = False
    tables: dict[str, dict] = {}
    cur = None
    done: set[str] = set()
    for line in lines:
        if line.startswith("@"):
            if line.startswith("@HD") and "SO:coordinate" in line:
                sorted_input = True
            continue
        f = line.split("\t", 12)
        cigar = f[5]
        if "N" not in cigar or int(f[1]) & exclude_flags:
            continue
        chrom = f[2]
        if strip_prefix and chrom.startswith(strip_prefix):
            chrom = chrom[len(strip_prefix):]
        if chrom != cur:
            if sorted_input:
                if chrom in done:
                    raise ec.EvidenceError("input declares SO:coordinate but a reference reappears: not coordinate-sorted")
                if cur is not None:
                    yield cur, tables.pop(cur)
                    done.add(cur)
            cur = chrom
        table = tables.setdefault(chrom, {})
        nh = 1
        for tag in f[11].split("\t") if len(f) > 11 else ():
            if tag.startswith("NH:i:"):
                nh = int(tag[5:])
                break
        name = f[0]
        ref, cur_run = int(f[3]), 0
        runs: list[int] = []
        juncs: list[tuple[int, int, int]] = []
        for n_s, op in CIGAR.findall(cigar):
            n = int(n_s)
            if op in "M=X":
                cur_run += n
                ref += n
            elif op == "I":
                runs.append(cur_run)
                cur_run = 0
            elif op == "D":
                runs.append(cur_run)
                cur_run = 0
                ref += n
            elif op == "N":
                runs.append(cur_run)
                juncs.append((ref, ref + n - 1, len(runs) - 1))
                cur_run = 0
                ref += n
        runs.append(cur_run)
        for s, e, ri in juncs:
            rec = table.setdefault((s, e), [set(), set(), 0])
            rec[0 if nh == 1 else 1].add(name)
            rec[2] = max(rec[2], min(runs[ri], runs[ri + 1]))
    for chrom, table in tables.items():
        yield chrom, table


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("alignments", help="BAM/CRAM/SAM file, or '-' for SAM text on stdin (with header)")
    ap.add_argument("--fasta", required=True, help="assembly FASTA (a .fai index is used if present); needed for motif and strand")
    ap.add_argument("-o", "--output", required=True, help="output SJ.out.tab-format file ('-' for stdout)")
    ap.add_argument("--strip-prefix", help="remove this prefix from alignment sequence names (e.g. 'lcl|')")
    ap.add_argument("--min-unique", type=int, default=0, help="keep junctions with >= N unique reads (default 0)")
    ap.add_argument("--min-total", type=int, default=1, help="keep junctions with >= N reads, unique + multimapping (default 1)")
    ap.add_argument("--min-overhang", type=int, default=0, help="keep junctions whose maximum overhang is >= N (default 0)")
    ap.add_argument("--exclude-flags", type=lambda x: int(x, 0), default=0x204, help="skip records with any of these SAM flag bits (default 0x204)")
    ap.add_argument("--summary-json", help="write provenance/statistics as JSON")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    fasta = Fasta(args.fasta)
    out_rows = []
    stats = {"junctions_seen": 0, "junctions_written": 0, "dropped_by_thresholds": 0, "noncanonical": 0, "sequences_not_in_fasta": 0}
    missing = set()
    for chrom, table in read_junctions(sam_lines(args.alignments), args.strip_prefix, args.exclude_flags):
        if chrom not in fasta.index:
            missing.add(chrom)
            continue
        rows = []
        for (s, e), (uniq, multi, ovh) in table.items():
            stats["junctions_seen"] += 1
            u, m = len(uniq), len(multi)
            if u < args.min_unique or u + m < args.min_total or ovh < args.min_overhang:
                stats["dropped_by_thresholds"] += 1
                continue
            donor, acceptor = fasta.fetch(chrom, s, s + 1), fasta.fetch(chrom, e - 1, e)
            strand, motif = MOTIFS.get((donor, acceptor), (0, 0))
            stats["noncanonical"] += motif == 0
            rows.append((chrom, s, e, strand, motif, 0, u, m, ovh))
        out_rows += sorted(rows, key=lambda r: (r[1], r[2]))
    stats["junctions_written"] = len(out_rows)
    stats["sequences_not_in_fasta"] = len(missing)
    if missing:
        ec.warn(f"{len(missing)} alignment sequence(s) not found in the FASTA and skipped: {', '.join(sorted(missing)[:5])}"
                " (use --strip-prefix for names such as 'lcl|NC_...')")
    text = "".join("\t".join(map(str, r)) + "\n" for r in out_rows)
    if args.output == "-":
        sys.stdout.write(text)
    else:
        with open(args.output, "w") as fh:
            fh.write(text)
    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    ec.write_summary_json(
        args.summary_json, "bam_to_junctions.py", __version__,
        parameters={"input": os.path.basename(args.alignments), "counting": "distinct read names; NH==1 unique, NH>1 multimapping",
                    "strip_prefix": args.strip_prefix, "exclude_flags": hex(args.exclude_flags)},
        filters={"min_unique": args.min_unique, "min_total": args.min_total, "min_overhang": args.min_overhang,
                 "star_outSJfilter": "not applied"}, stats=stats)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
