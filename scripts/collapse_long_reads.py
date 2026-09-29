#!/usr/bin/env python3
"""Collapse per-sample long-read alignments (BED12) into ``long_read_models`` (BED12+4).

Reads with the same seqid, strand and INTRON CHAIN (identical intron coordinates) become one
transcript model. The model records how many reads, samples and groups support it. Its start and
end are the median observed start and end (an observed value, ``median_low``), so a few truncated
or over-extended reads do not move the boundaries.

Input: one BED12 file per sample, one line per alignment, e.g. from a filtered BAM:
    samtools view -b -F 0x900 -q 20 sample.bam | bedtools bamtobed -bed12 -i - > sample.bed12
(filter secondary/supplementary alignments and low mapping quality BEFORE this step; that choice
belongs in the manifest.) Alignment gaps shorter than --min-intron are merged back into their
block (treated as deletions, not introns).

Single-exon alignments have no intron chain. They are clustered instead: reads on the same seqid and strand whose
start AND end each lie within --single-exon-tolerance bases (default 50) of the cluster's first (leftmost) read form
one single-exon model (one block, intron_chain "."). They are counted and kept because single-exon genes and
unspliced transcripts matter, but they are often fragments or unspliced RNA, so they are always reported separately
from spliced models. Use --no-single-exon to drop them (the count is still recorded).

Output columns: BED12 + read_count, sample_count, group_count, intron_chain.
  chromStart/chromEnd 0-based half-open; intron_chain uses 1-based closed intron coordinates,
  e.g. "1301-1700;1901-2500", identical in convention to rna_junctions rows and GFF3 introns.

Example:
  collapse_long_reads.py --samples samples.tsv -o evidence/long_read_models.bed.gz  BC01=bc01.bed12 BC02=bc02.bed12
"""
from __future__ import annotations

import argparse
import os
import statistics
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"


def merge_short_gaps(exons, min_intron: int):
    """Merge blocks separated by a gap of fewer than ``min_intron`` bases (1-based closed exons)."""
    out = [list(exons[0])]
    for s, e in exons[1:]:
        gap = s - out[-1][1] - 1
        if gap < min_intron:
            out[-1][1] = e
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def read_alignments(path: str, strip_prefix: str | None):
    with ec.open_text(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ec.EvidenceError(f"{path}:{n}: expected BED12 (>= 12 columns), got {len(f)}")
            seqid = f[0][len(strip_prefix):] if strip_prefix and f[0].startswith(strip_prefix) else f[0]
            cs = int(f[1])
            sizes = [int(x) for x in f[10].strip(",").split(",")]
            starts = [int(x) for x in f[11].strip(",").split(",")]
            yield seqid, f[5], ec.exons_from_bed12(cs, sizes, starts)


def cluster_single_exon(reads, tol: int):
    """Cluster single-exon reads (seqid, strand, start, end, sample). Deterministic: reads are processed in
    (seqid, strand, start, end, sample) order and join the first open cluster whose anchor (its first read)
    has start and end within ``tol`` bases; a cluster closes once a read starts more than ``tol`` after its anchor."""
    clusters = []
    open_by_key: dict[tuple, list] = {}
    for seqid, strand, start, end, sample in sorted(reads):
        opened = open_by_key.setdefault((seqid, strand), [])
        opened[:] = [c for c in opened if start - c["anchor"][0] <= tol]
        for c in opened:
            if abs(end - c["anchor"][1]) <= tol:
                break
        else:
            c = {"seqid": seqid, "strand": strand, "anchor": (start, end), "starts": [], "ends": [], "samples": defaultdict(int)}
            opened.append(c)
            clusters.append(c)
        c["starts"].append(start)
        c["ends"].append(end)
        c["samples"][sample] += 1
    return clusters


def collapse(inputs: dict[str, str], groups, min_intron: int, strip_prefix: str | None,
             single_exon: bool = True, single_exon_tolerance: int = 50):
    models: dict[tuple, dict] = {}
    single_reads = []
    stats = {"input_files": len(inputs), "input_alignments": 0, "single_exon_alignments": 0}
    for sample, path in inputs.items():
        for seqid, strand, exons in read_alignments(path, strip_prefix):
            stats["input_alignments"] += 1
            exons = merge_short_gaps(exons, min_intron)
            introns = tuple(ec.introns_from_exons(exons))
            if not introns:
                stats["single_exon_alignments"] += 1
                if single_exon:
                    single_reads.append((seqid, strand, exons[0][0], exons[-1][1], sample))
                continue
            m = models.setdefault((seqid, strand, introns), {"starts": [], "ends": [], "samples": defaultdict(int)})
            m["starts"].append(exons[0][0])
            m["ends"].append(exons[-1][1])
            m["samples"][sample] += 1
    rows = []
    for (seqid, strand, introns), m in models.items():
        start = statistics.median_low(m["starts"])
        end = statistics.median_low(m["ends"])
        exons = []
        pos = start
        for i_s, i_e in introns:
            exons.append((pos, i_s - 1))
            pos = i_e + 1
        exons.append((pos, end))
        known = {groups.get(s, ".") for s in m["samples"]} - {"."} if groups is not None else set()
        rows.append({"seqid": seqid, "strand": strand, "exons": exons, "introns": list(introns),
                     "read_count": sum(m["samples"].values()), "sample_count": len(m["samples"]),
                     "group_count": len(known) if known else "."})
    n_spliced = len(rows)
    for c in cluster_single_exon(single_reads, single_exon_tolerance):
        known = ({groups.get(s, ".") for s in c["samples"]} - {"."}) if groups is not None else set()
        rows.append({"seqid": c["seqid"], "strand": c["strand"],
                     "exons": [(statistics.median_low(c["starts"]), statistics.median_low(c["ends"]))], "introns": [],
                     "read_count": sum(c["samples"].values()), "sample_count": len(c["samples"]),
                     "group_count": len(known) if known else "."})
    rows.sort(key=lambda r: (r["seqid"], r["exons"][0][0], r["exons"][-1][1], r["strand"], r["introns"]))
    stats["spliced_models"] = n_spliced
    stats["single_exon_models"] = len(rows) - n_spliced
    stats["single_exon_reads_in_models"] = sum(r["read_count"] for r in rows if not r["introns"])
    stats["models"] = len(rows)
    return rows, stats


def bed_line(row, model_id: str) -> str:
    cs, ce, sizes, starts = ec.bed12_from_exons(row["exons"])
    return "\t".join(str(x) for x in (
        row["seqid"], cs, ce, model_id, min(row["read_count"], 1000), row["strand"], cs, ce, "0,0,0",
        len(sizes), ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ",",
        row["read_count"], row["sample_count"], row["group_count"], ec.format_chain(row["introns"])))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bed12", nargs="+", metavar="BED12", help="per-sample BED12 alignments; SAMPLE=path sets the sample id")
    ap.add_argument("-o", "--output", required=True, help="output path; '.gz' means bgzip + tabix (needs htslib)")
    ap.add_argument("--samples", help="samples.tsv (sample_id, group, evidence_class) to define groups")
    ap.add_argument("--min-intron", type=int, default=30, help="gaps shorter than this are not introns (default 30)")
    ap.add_argument("--no-single-exon", action="store_true", help="do not build single-exon models (they are still counted)")
    ap.add_argument("--single-exon-tolerance", type=int, default=50,
                    help="single-exon reads whose start and end are each within N bases of a cluster's first read are one model (default 50)")
    ap.add_argument("--min-reads", type=int, default=1, help="keep models with >= N reads (default 1)")
    ap.add_argument("--min-samples", type=int, default=1, help="keep models seen in >= N samples (default 1)")
    ap.add_argument("--id-prefix", default="LR", help="model id prefix (default LR)")
    ap.add_argument("--strip-prefix", help="remove this prefix from seqids (e.g. 'lcl|')")
    ap.add_argument("--summary-json", help="write provenance/statistics as JSON (read by make_manifest.py)")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    inputs: dict[str, str] = {}
    for item in args.bed12:
        sid, sep, path = item.partition("=")
        if not sep or os.path.exists(item):
            sid, path = os.path.basename(item).split(".")[0], item
        if sid in inputs:
            raise ec.EvidenceError(f"duplicate sample id '{sid}'")
        if not os.path.exists(path):
            raise ec.EvidenceError(f"file not found: {path}")
        inputs[sid] = path

    groups = None
    if args.samples:
        header, rows = ec.read_tsv_table(args.samples)
        groups = {r["sample_id"]: r.get("group", ".") for r in rows if r.get("evidence_class", "rna_long") in ("rna_long", ".")}
        missing = [s for s in inputs if s not in groups]
        if missing:
            ec.warn(f"sample(s) not listed as rna_long in {args.samples}: {', '.join(missing[:5])}")

    rows, stats = collapse(inputs, groups, args.min_intron, args.strip_prefix, not args.no_single_exon, args.single_exon_tolerance)
    kept = [r for r in rows if r["read_count"] >= args.min_reads and r["sample_count"] >= args.min_samples]
    stats["output_models"] = len(kept)
    lines = [bed_line(r, f"{args.id_prefix}{i:06d}") for i, r in enumerate(kept, 1)]

    if args.output.endswith(".gz"):
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=os.path.dirname(os.path.abspath(args.output))) as tmp:
            tmp.write("\n".join(lines) + ("\n" if lines else ""))
        try:
            ec.bgzip_and_index(tmp.name, args.output, "bed")
        finally:
            os.unlink(tmp.name)
    else:
        with open(args.output, "w") as fh:
            fh.write("\n".join(lines) + ("\n" if lines else ""))

    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    ec.write_summary_json(
        args.summary_json, "collapse_long_reads.py", __version__,
        parameters={"collapse_key": "seqid+strand+intron_chain (single-exon: start/end clusters)", "boundaries": "median observed start/end",
                    "n_samples": len(inputs), "strip_prefix": args.strip_prefix},
        filters={"min_intron_length": args.min_intron,
                 "single_exon_alignments": ("excluded (counted in stats)" if args.no_single_exon
                                            else f"clustered as single-exon models (start and end within {args.single_exon_tolerance} bp)"),
                 "min_reads": args.min_reads, "min_samples": args.min_samples}, stats=stats)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
