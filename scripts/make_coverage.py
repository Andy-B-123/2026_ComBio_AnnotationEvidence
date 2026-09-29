#!/usr/bin/env python3
"""Make the optional merged RNA-seq coverage track from one or more BAMs.

Basic, deliberately simple recipe (the "coverage" evidence class):
  1. per BAM: depth of SPLICED coverage (introns are not counted) with ``bedtools genomecov -bg -split``;
  2. sum the depths over all BAMs (``bedtools unionbedg``), so a value means "about this many reads cover the base";
  3. optionally average in fixed windows (--window), then snap the depth to a geometric ladder ``2^(k/L) - 1``
     (--ladder L, default 4, about 19 % steps; 0 = keep exact values). Values stay on the DEPTH scale;
  4. write a bigWig (needs UCSC ``bedGraphToBigWig``), or a bedGraph if the output ends in .bedgraph / .bedgraph.gz.

On six real libraries (one 15 Mb chromosome, 408 MB of BAM) this gave a 3.6 MB bigWig at full resolution, 0.95 MB with the
default ladder and 0.48 MB with the ladder plus 50-bp windows (docs/DESIGN_RATIONALE.md). Not strand-specific. By default every
alignment is counted, including multimappers; --unique-only counts NH:i:1 alignments only (needs samtools).

Examples:
  make_coverage.py --sequences bundle/assembly_sequences.tsv -o bundle/evidence/coverage.bw a.bam b.bam c.bam
  make_coverage.py --sizes chrom.sizes --window 50 --strip-prefix 'lcl|' -o cov.bedgraph.gz *.bam
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"


def snap(depth: float, ladder: int) -> float:
    """Snap a depth to the geometric ladder 2^(k/ladder) - 1 (k integer); 0 stays 0. ladder 0 = unchanged."""
    if ladder <= 0 or depth <= 0:
        return depth
    return 2 ** (round(ladder * math.log2(1 + depth)) / ladder) - 1


def summed_records(lines, strip_prefix: str | None):
    """Yield (chrom, start0, end0, summed_depth) from ``bedtools unionbedg`` output (chrom start end v1 v2 ...)."""
    for line in lines:
        f = line.rstrip("\n").split("\t")
        chrom = f[0][len(strip_prefix):] if strip_prefix and f[0].startswith(strip_prefix) else f[0]
        yield chrom, int(f[1]), int(f[2]), sum(float(x) for x in f[3:])


def windowed(records, size: int, lengths: dict[str, int]):
    """Turn (chrom, start0, end0, depth) records (sorted) into per-window mean depth records."""
    cur, acc = None, {}

    def flush():
        L = lengths[cur]
        for b in sorted(acc):
            lo, hi = b * size, min((b + 1) * size, L)
            yield cur, lo, hi, acc[b] / (hi - lo)

    for chrom, s, e, d in records:
        if chrom != cur:
            if cur is not None:
                yield from flush()
            cur, acc = chrom, {}
        b = s // size
        while b * size < e:
            lo, hi = max(s, b * size), min(e, (b + 1) * size)
            acc[b] = acc.get(b, 0.0) + (hi - lo) * d
            b += 1
    if cur is not None:
        yield from flush()


def quantised(records, ladder: int, lengths: dict[str, int]):
    """Snap values, clip to sequence length, merge adjacent equal intervals; yields bedGraph lines."""
    prev = None
    for chrom, s, e, d in records:
        if chrom not in lengths:
            continue
        e = min(e, lengths[chrom])
        if e <= s:
            continue
        v = f"{snap(d, ladder):.2f}"
        if v == "0.00":
            if prev:
                yield "\t".join(map(str, prev))
                prev = None
            continue
        if prev and prev[0] == chrom and prev[2] == s and prev[3] == v:
            prev[2] = e
        else:
            if prev:
                yield "\t".join(map(str, prev))
            prev = [chrom, s, e, v]
    if prev:
        yield "\t".join(map(str, prev))


def read_sizes(args) -> dict[str, int]:
    if args.sizes:
        return {l.split()[0]: int(l.split()[1]) for l in open(args.sizes) if l.strip() and not l.startswith("#")}
    header, rows = ec.read_tsv_table(args.sequences)
    return {r["seqid"]: int(r["length"]) for r in rows}


def require(tool: str, why: str):
    if not ec.have_tool(tool):
        raise ec.EvidenceError(f"{tool} not found on PATH ({why})")


def tool_version(cmd: list[str]) -> str:
    res = subprocess.run(cmd, capture_output=True, text=True)
    return ((res.stdout or "") + (res.stderr or "")).strip().splitlines()[0] if (res.stdout or res.stderr) else "unknown"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bams", nargs="+", help="coordinate-sorted BAMs (spliced short-read RNA-seq)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sizes", help="two-column chrom.sizes file")
    src.add_argument("--sequences", help="assembly_sequences.tsv of the bundle")
    ap.add_argument("-o", "--output", required=True, help=".bw (bigWig), .bedgraph or .bedgraph.gz")
    ap.add_argument("--ladder", type=int, default=4, help="steps per doubling of (1 + depth); 0 keeps exact depths (default 4)")
    ap.add_argument("--window", type=int, default=0, help="average depth in windows of N bases (default 0 = off)")
    ap.add_argument("--unique-only", action="store_true", help="count only uniquely mapping (NH:i:1) alignments; needs samtools")
    ap.add_argument("--strip-prefix", help="remove this prefix from BAM sequence names (e.g. 'lcl|')")
    ap.add_argument("--summary-json", help="write provenance as JSON (read by make_manifest.py)")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    lengths = read_sizes(args)
    require("bedtools", "needed for genomecov and unionbedg")
    if args.unique_only:
        require("samtools", "needed for --unique-only")
    out = args.output
    kind = "bigwig" if out.endswith((".bw", ".bigwig")) else ("bedgraph.gz" if out.endswith(".bedgraph.gz") else "bedgraph")
    if kind == "bigwig":
        require("bedGraphToBigWig", "UCSC tool; alternatively write .bedgraph.gz")
    outdir = os.path.dirname(os.path.abspath(out))
    with tempfile.TemporaryDirectory(dir=outdir, prefix=".cov_") as tmp:
        files = []
        for i, bam in enumerate(args.bams):
            bg = os.path.join(tmp, f"s{i}.bedgraph")
            with open(bg, "w") as fh:
                if args.unique_only:
                    p1 = subprocess.Popen(["samtools", "view", "-b", "-h", "-d", "NH:1", bam], stdout=subprocess.PIPE)
                    subprocess.run(["bedtools", "genomecov", "-ibam", "stdin", "-bg", "-split"], stdin=p1.stdout, stdout=fh, check=True)
                    p1.stdout.close()
                    if p1.wait() != 0:
                        raise ec.EvidenceError(f"samtools view failed on {bam}")
                else:
                    subprocess.run(["bedtools", "genomecov", "-ibam", bam, "-bg", "-split"], stdout=fh, check=True)
            files.append(bg)
        if len(files) == 1:
            union = subprocess.Popen(["cat", files[0]], stdout=subprocess.PIPE, text=True)
        else:
            union = subprocess.Popen(["bedtools", "unionbedg", "-i", *files], stdout=subprocess.PIPE, text=True)
        recs = summed_records(union.stdout, args.strip_prefix)
        if args.window > 0:
            recs = windowed(recs, args.window, lengths)
        lines = list(quantised(recs, args.ladder, lengths))
        if union.wait() != 0:
            raise ec.EvidenceError("bedtools unionbedg failed")
        bg_out = os.path.join(tmp, "merged.bedgraph")
        with open(bg_out, "w") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))
        if kind == "bigwig":
            sizes = os.path.join(tmp, "chrom.sizes")
            with open(sizes, "w") as fh:
                fh.write("".join(f"{k}\t{v}\n" for k, v in lengths.items()))
            subprocess.run(["bedGraphToBigWig", bg_out, sizes, out], check=True)
        elif kind == "bedgraph.gz":
            ec.bgzip_and_index(bg_out, out, "bed")
        else:
            shutil.copyfile(bg_out, out)
    stats = {"input_bams": len(args.bams), "intervals_written": len(lines), "output_bytes": os.path.getsize(out)}
    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    ec.write_summary_json(
        args.summary_json, "make_coverage.py", __version__,
        parameters={"bedtools": tool_version(["bedtools", "--version"]), "counting": "spliced (-split), summed over BAMs",
                    "n_bams": len(args.bams)},
        filters={"quantisation": (f"depth snapped to 2^(k/{args.ladder}) - 1, k integer; values are depth, not log depth" if args.ladder > 0 else "none (exact depth)"),
                 "window_bp": args.window or "none", "strand": "both", "unique_only": args.unique_only,
                 "multimappers": "excluded" if args.unique_only else "counted"}, stats=stats)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ec.EvidenceError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
