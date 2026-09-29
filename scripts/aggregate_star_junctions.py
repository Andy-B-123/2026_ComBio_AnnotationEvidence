#!/usr/bin/env python3
"""Collapse per-sample STAR ``SJ.out.tab`` files into one ``rna_junctions`` table.

Each output row is one intron (seqid, strand, first and last intron base, 1-based closed) with
how many samples and groups support it and summary statistics of the unique-read counts, so that
"100 reads in one library" can be told apart from "supported independently in 20 libraries".

Definitions (recorded in the output header and, via make_manifest.py, in manifest.yaml):
  * A sample SUPPORTS a junction if it has >= --min-sample-unique uniquely mapping reads on it.
    Samples that do not support a junction contribute nothing to that junction's row.
  * sample_count / group_count : number of supporting samples / distinct groups among them.
  * total, max, median unique reads: over supporting samples only.
  * total_multimapping_reads: multimapper reads summed over supporting samples.
  * max_overhang: largest STAR overhang seen in any supporting sample.
  * After aggregation a junction is kept only if sample_count >= --min-samples and
    total_unique_reads >= --min-total-unique.
NOTE: SJ.out.tab is already filtered by STAR (--outSJfilter*, see docs/SPECIFICATION.md);
record those settings in the manifest.

Example:
  aggregate_star_junctions.py --samples samples.tsv -o evidence/rna_junctions.tsv.gz  star/*.SJ.out.tab
"""
from __future__ import annotations

import argparse
import os
import statistics
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"

STAR_STRAND = {0: ".", 1: "+", 2: "-"}
# Motif is reported in transcript orientation: CT/AC on the minus strand is GT/AG, etc.
STAR_MOTIF = {0: "non-canonical", 1: "GT/AG", 2: "GT/AG", 3: "GC/AG", 4: "GC/AG", 5: "AT/AC", 6: "AT/AC"}
STAR_MOTIF_STRAND = {0: ".", 1: "+", 2: "-", 3: "+", 4: "-", 5: "+", 6: "-"}
SUFFIXES = (".SJ.out.tab", "_SJ.out.tab", "SJ.out.tab", ".tab")


def sample_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    for suf in SUFFIXES:
        if base.endswith(suf) and len(base) > len(suf):
            return base[: -len(suf)]
    return base


def parse_star_line(line: str, where: str):
    f = line.rstrip("\n").split("\t")
    if len(f) != 9:
        raise ec.EvidenceError(f"{where}: expected 9 tab-separated columns (STAR SJ.out.tab), got {len(f)}")
    try:
        start, end = int(f[1]), int(f[2])
        strand_code, motif_code = int(f[3]), int(f[4])
        unique, multi, overhang = int(f[6]), int(f[7]), int(f[8])
    except ValueError as exc:
        raise ec.EvidenceError(f"{where}: non-numeric field ({exc})") from exc
    if strand_code not in STAR_STRAND or motif_code not in STAR_MOTIF:
        raise ec.EvidenceError(f"{where}: unknown strand/motif code {strand_code}/{motif_code}")
    if end < start or start < 1:
        raise ec.EvidenceError(f"{where}: invalid intron coordinates {start}-{end}")
    return f[0], start, end, strand_code, motif_code, unique, multi, overhang


def aggregate(sj_files: dict[str, str], groups: dict[str, str] | None, min_sample_unique: int,
              strip_prefix: str | None):
    """Return (rows, stats). ``sj_files`` maps sample_id -> path; rows are dicts sorted for tabix."""
    acc: dict[tuple, dict] = {}
    stats = {"input_files": len(sj_files), "input_rows": 0, "rows_below_sample_threshold": 0,
             "motif_conflicts": 0, "strand_motif_inconsistent": 0}
    for sample, path in sj_files.items():
        with ec.open_text(path) as fh:
            for n, line in enumerate(fh, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                seqid, start, end, sc, mc, uniq, multi, ovh = parse_star_line(line, f"{path}:{n}")
                stats["input_rows"] += 1
                if strip_prefix and seqid.startswith(strip_prefix):
                    seqid = seqid[len(strip_prefix):]
                if uniq < min_sample_unique:
                    stats["rows_below_sample_threshold"] += 1
                    continue
                if STAR_MOTIF_STRAND[mc] != STAR_STRAND[sc]:
                    stats["strand_motif_inconsistent"] += 1
                key = (seqid, start, end, STAR_STRAND[sc])
                a = acc.get(key)
                if a is None:
                    a = acc[key] = {"motif": STAR_MOTIF[mc], "uniq": [], "multi": 0, "ovh": 0, "groups": set()}
                elif a["motif"] != STAR_MOTIF[mc]:
                    stats["motif_conflicts"] += 1
                a["uniq"].append(uniq)
                a["multi"] += multi
                a["ovh"] = max(a["ovh"], ovh)
                if groups is not None:
                    a["groups"].add(groups.get(sample, "."))
    rows = []
    for (seqid, start, end, strand), a in acc.items():
        u = a["uniq"]
        known = {g for g in a["groups"] if g != "."}
        gcount = len(known) if known else "."
        rows.append({"seqid": seqid, "intron_start": start, "intron_end": end, "strand": strand,
                     "motif": a["motif"], "sample_count": len(u), "group_count": gcount,
                     "total_unique_reads": sum(u), "max_unique_reads": max(u),
                     "median_unique_reads": f"{statistics.median(u):g}",
                     "total_multimapping_reads": a["multi"], "max_overhang": a["ovh"]})
    rows.sort(key=lambda r: (r["seqid"], r["intron_start"], r["intron_end"], r["strand"]))
    stats["distinct_junctions_before_filters"] = len(rows)
    return rows, stats


def write_table(rows, handle, header_pragmas):
    for line in header_pragmas:
        handle.write(line + "\n")
    handle.write("#" + "\t".join(ec.JUNCTION_COLUMNS) + "\n")
    for r in rows:
        handle.write("\t".join(str(r[c]) for c in ec.JUNCTION_COLUMNS) + "\n")


def load_groups(samples_path: str, sample_ids) -> dict[str, str]:
    header, rows = ec.read_tsv_table(samples_path)
    if "sample_id" not in header:
        raise ec.EvidenceError(f"{samples_path}: no 'sample_id' column")
    groups = {}
    for r in rows:
        if r.get("evidence_class", "rna_short") in ("rna_short", "."):
            groups[r["sample_id"]] = r.get("group", ".") or "."
    missing = [s for s in sample_ids if s not in groups]
    if missing:
        ec.warn(f"{len(missing)} sample(s) not found in {samples_path} as rna_short: {', '.join(missing[:5])}"
                + (" ..." if len(missing) > 5 else ""))
    return groups


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("sj", nargs="+", metavar="SJ", help="SJ.out.tab files; use SAMPLE=path to set the sample id "
                    "(default: file name without .SJ.out.tab)")
    ap.add_argument("-o", "--output", required=True,
                    help="output path; '.gz' means bgzip + tabix index (needs htslib), otherwise plain text")
    ap.add_argument("--samples", help="samples.tsv (columns sample_id, group, evidence_class) to define groups")
    ap.add_argument("--min-sample-unique", type=int, default=1,
                    help="unique reads a sample needs to count as supporting a junction (default 1)")
    ap.add_argument("--min-samples", type=int, default=1, help="keep junctions with >= N supporting samples (default 1)")
    ap.add_argument("--min-total-unique", type=int, default=1,
                    help="keep junctions with >= N unique reads in total (default 1)")
    ap.add_argument("--strip-prefix", help="remove this prefix from seqids (e.g. 'lcl|')")
    ap.add_argument("--star-settings", default=None,
                    help="free text recorded in the header: the STAR outSJfilter* settings used upstream")
    ap.add_argument("--summary-json", help="write run statistics as JSON")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    sj_files: dict[str, str] = {}
    for item in args.sj:
        sid, sep, path = item.partition("=")
        if not sep or os.path.exists(item):
            sid, path = sample_id_from_path(item), item
        if sid in sj_files:
            raise ec.EvidenceError(f"duplicate sample id '{sid}'")
        if not os.path.exists(path):
            raise ec.EvidenceError(f"file not found: {path}")
        sj_files[sid] = path

    groups = load_groups(args.samples, sj_files) if args.samples else None
    rows, stats = aggregate(sj_files, groups, args.min_sample_unique, args.strip_prefix)
    kept = [r for r in rows if r["sample_count"] >= args.min_samples and r["total_unique_reads"] >= args.min_total_unique]
    stats.update(output_rows=len(kept), dropped_by_min_samples_or_total=len(rows) - len(kept))
    if stats["strand_motif_inconsistent"]:
        ec.warn(f"{stats['strand_motif_inconsistent']} rows had a strand that disagrees with the intron motif")
    if stats["motif_conflicts"]:
        ec.warn(f"{stats['motif_conflicts']} junctions had different motif codes in different samples (first kept)")

    n_groups = len({g for s, g in (groups or {}).items() if s in sj_files and g != "."}) or None
    pragmas = ["##ae.profile=rna_junctions"] + ec.pragma_lines(
        version=ec.SPEC_VERSION, coordinates="1-based-closed;first-and-last-intron-base",
        generated_by=f"aggregate_star_junctions.py {__version__}", source_format="STAR SJ.out.tab",
        n_samples=len(sj_files), n_groups=n_groups if n_groups is not None else ".",
        sample_min_unique_reads=args.min_sample_unique, min_samples=args.min_samples,
        min_total_unique_reads=args.min_total_unique,
        seqid_prefix_stripped=args.strip_prefix or ".", upstream_star_filters=args.star_settings or ".",
        motif="transcript-orientation")

    out = args.output
    if out.endswith(".gz"):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False, dir=os.path.dirname(os.path.abspath(out))) as tmp:
            write_table(kept, tmp, pragmas)
        try:
            ec.bgzip_and_index(tmp.name, out, "tsv")
        finally:
            os.unlink(tmp.name)
    else:
        with open(out, "w") as fh:
            write_table(kept, fh, pragmas)

    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    ec.write_summary_json(
        args.summary_json, "aggregate_star_junctions.py", __version__,
        parameters={"source_format": "STAR SJ.out.tab", "n_samples": len(sj_files), "strip_prefix": args.strip_prefix,
                    "upstream_star_filters": args.star_settings},
        filters={"min_sample_unique_reads": args.min_sample_unique, "min_samples": args.min_samples,
                 "min_total_unique_reads": args.min_total_unique}, stats=stats)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
