#!/usr/bin/env python3
"""Show every piece of retained evidence that overlaps a genomic region (or a gene) of a bundle.

All coordinates in the output are 1-based, closed, whatever the storage convention of the file
(BED files are converted). Junction and intron coordinates are the first and last INTRON base.
Each evidence class is printed with the filters recorded in the manifest, because a missing
record can only be interpreted if you know what was filtered out.

Uses tabix when the file is indexed and tabix is installed, otherwise scans the file.

Examples:
  query_region.py bundle/ scaffold_1:2500-3500
  query_region.py bundle/ --gene geneA --flank 500 --format json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

CLASSES = ("annotation", "rna_short", "rna_long", "proteins", "predictions", "repeats", "coverage")


def spans(intervals) -> str:
    return ", ".join(f"{s}-{e}" for s, e in intervals) if intervals else "-"


def gene_dict(g: ec.Gene) -> dict:
    return {"id": g.id, "seqid": g.seqid, "start": g.start, "end": g.end, "strand": g.strand, "source": g.source,
            "names": sorted(g.names), "partial": g.attrs.get("partial") == "true",
            "transcripts": [{"id": t.id, "type": t.type, "start": t.start, "end": t.end, "exons": t.exons, "cds": t.cds,
                             "introns": t.introns, "partial": t.partial} for t in g.transcripts]}


def collect(b: ec.Bundle, seqid: str, start: int, end: int | None, classes, bins: int, use_tabix: bool) -> dict:
    out: dict = {"bundle": (b.manifest.get("bundle") or {}).get("id"), "spec_version": b.manifest.get("spec_version"),
                 "region": {"seqid": seqid, "start": start, "end": end, "coordinates": "1-based-closed"}, "classes": {}}
    lens = b.sequences()
    if lens and seqid not in lens:
        raise ec.EvidenceError(f"'{seqid}' is not a sequence of this assembly (see assembly_sequences.tsv)")
    if end is None:
        end = lens.get(seqid)
    out["region"]["end"] = end
    ev = b.classes()

    def cls_meta(name):
        e = ev.get(name) or {}
        return {"present": bool(e.get("present")), "filters": e.get("filters") or {}}

    if "annotation" in classes:
        out["classes"]["annotation"] = {"present": True, "genes": [gene_dict(g) for g in b.annotation_genes(seqid, start, end)]}
    if "rna_short" in classes:
        m = cls_meta("rna_short")
        m["n_datasets"] = (ev.get("rna_short") or {}).get("n_datasets")
        m["junctions"] = b.junctions(seqid, start, end, use_tabix) if m["present"] else []
        out["classes"]["rna_short"] = m
    if "rna_long" in classes:
        m = cls_meta("rna_long")
        m["models"] = b.long_read_models(seqid, start, end, use_tabix) if m["present"] else []
        out["classes"]["rna_long"] = m
    if "proteins" in classes:
        m = cls_meta("proteins")
        m["alignments"] = b.protein_alignments(seqid, start, end) if m["present"] else []
        out["classes"]["proteins"] = m
    if "predictions" in classes:
        m = cls_meta("predictions")
        m["genes"] = [gene_dict(g) for g in b.prediction_genes(seqid, start, end)] if m["present"] else []
        out["classes"]["predictions"] = m
    if "repeats" in classes:
        m = cls_meta("repeats")
        m["repeats"] = b.repeats(seqid, start, end) if m["present"] else []
        out["classes"]["repeats"] = m
    if "coverage" in classes:
        m = cls_meta("coverage")
        m["bins"] = []
        if m["present"] and bins > 0 and end:
            iv = b.coverage_intervals(seqid, start, end)
            if iv is None:
                m["note"] = "coverage is present but cannot be read: install pyBigWig or the UCSC bigWigToBedGraph tool"
            else:
                size = max(1, -(-(end - start + 1) // bins))
                for s in range(start, end + 1, size):
                    e = min(end, s + size - 1)
                    m["bins"].append({"start": s, "end": e, "mean_depth": round(ec.mean_depth(iv, s, e), 2)})
        out["classes"]["coverage"] = m
    return out


def fmt_filters(f: dict) -> str:
    return "; ".join(f"{k}={v}" for k, v in f.items()) if f else "none recorded"


def render_text(d: dict) -> str:
    r = d["region"]
    L = [f"# bundle: {d['bundle']} (spec {d['spec_version']})",
         f"# region: {r['seqid']}:{r['start']}-{r['end']}   ALL COORDINATES BELOW ARE 1-BASED, CLOSED",
         "# intron / junction coordinates are the first and last INTRON base; absence of a record is not evidence against a model", ""]
    c = d["classes"]
    if "annotation" in c:
        L.append(f"== annotation: {len(c['annotation']['genes'])} gene(s)")
        for g in c["annotation"]["genes"]:
            L.append(f"gene {g['id']}  {g['seqid']}:{g['start']}-{g['end']} ({g['strand']})" + ("  partial=true" if g["partial"] else ""))
            for t in g["transcripts"]:
                L.append(f"  {t['type']} {t['id']}  exons: {spans(t['exons'])}")
                L.append(f"    CDS: {spans(t['cds'])}   introns: {spans(t['introns'])}")
        L.append("")
    if "rna_short" in c:
        s = c["rna_short"]
        if not s["present"]:
            L.append("== rna_short: not present in this bundle\n")
        else:
            L.append(f"== rna_short: {len(s['junctions'])} junction(s) from {s['n_datasets']} sample(s)   [filters: {fmt_filters(s['filters'])}]")
            for j in s["junctions"]:
                gc = "?" if j["group_count"] is None else j["group_count"]
                L.append(f"{j['seqid']}:{j['start']}-{j['end']} ({j['strand']}) {j['motif']}  samples={j['sample_count']} groups={gc} "
                         f"unique total/max/median={j['total_unique_reads']}/{j['max_unique_reads']}/{j['median_unique_reads']:g} "
                         f"multimapping={j['total_multimapping_reads']} max_overhang={j['max_overhang']}")
            L.append("")
    if "rna_long" in c:
        s = c["rna_long"]
        if not s["present"]:
            L.append("== rna_long: not present in this bundle\n")
        else:
            L.append(f"== rna_long: {len(s['models'])} model(s)   [filters: {fmt_filters(s['filters'])}]")
            for m in s["models"]:
                gc = "?" if m.get("group_count") is None else m["group_count"]
                L.append(f"{m['id']}  {m['seqid']}:{m['start']}-{m['end']} ({m['strand']})  reads={m.get('read_count')} "
                         f"samples={m.get('sample_count')} groups={gc}")
                L.append(f"    exons: {spans(m['exons'])}   introns: {spans(m['introns'])}")
            L.append("")
    if "proteins" in c:
        s = c["proteins"]
        if not s["present"]:
            L.append("== proteins: not present in this bundle\n")
        else:
            L.append(f"== proteins: {len(s['alignments'])} alignment(s)   [filters: {fmt_filters(s['filters'])}]")
            for a in s["alignments"]:
                at = a["attrs"]
                L.append(f"{a['id']}  {a['protein']} aa {a['protein_start']}-{a['protein_end']} -> {a['seqid']}:{a['start']}-{a['end']} ({a['strand']})  "
                         f"identity={at.get('pct_identity')}% coverage={at.get('pct_coverage')}% rank={at.get('rank')} db={at.get('db', '.')}")
                L.append("    blocks: " + ", ".join(f"{b['start']}-{b['end']} (aa {b['target_start']}-{b['target_end']})" for b in a["blocks"]))
            L.append("")
    if "predictions" in c:
        s = c["predictions"]
        if not s["present"]:
            L.append("== predictions: not present in this bundle\n")
        else:
            L.append(f"== predictions: {len(s['genes'])} gene model(s)")
            for g in s["genes"]:
                for t in g["transcripts"]:
                    L.append(f"{g['source']} {t['id']}  {g['seqid']}:{t['start']}-{t['end']} ({g['strand']})  exons: {spans(t['exons'])}")
            L.append("")
    if "repeats" in c and c["repeats"]["present"]:
        L.append(f"== repeats: {len(c['repeats']['repeats'])} interval(s)")
        for x in c["repeats"]["repeats"]:
            L.append(f"{x['name']}  {x['seqid']}:{x['start']}-{x['end']} ({x['strand']})")
        L.append("")
    if "coverage" in c and c["coverage"]["present"]:
        s = c["coverage"]
        L.append(f"== coverage: mean depth per bin   [{fmt_filters(s['filters'])}]")
        if s.get("note"):
            L.append(s["note"])
        for b in s["bins"]:
            L.append(f"{b['start']}-{b['end']}  {b['mean_depth']}")
        L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bundle", help="bundle directory")
    ap.add_argument("region", nargs="?", help="seqid:start-end (1-based, closed; commas allowed), or just seqid")
    ap.add_argument("--gene", help="use the span of this gene (ID, Name or locus_tag) instead of a region")
    ap.add_argument("--flank", type=int, default=0, help="extend the region/gene span by N bases on each side")
    ap.add_argument("--classes", default=",".join(CLASSES), help=f"comma-separated subset of: {','.join(CLASSES)}")
    ap.add_argument("--format", choices=("text", "json"), default="text")
    ap.add_argument("--coverage-bins", type=int, default=10, help="number of coverage bins to report (0 = off; default 10)")
    ap.add_argument("--no-tabix", action="store_true", help="scan files instead of using tabix (for testing)")
    args = ap.parse_args(argv)

    b = ec.Bundle(args.bundle)
    if args.gene:
        g = ec.find_gene(b.annotation_path, args.gene)
        seqid, start, end = g.seqid, g.start, g.end
    elif args.region:
        seqid, start, end = ec.parse_region(args.region)
    else:
        ap.error("give a REGION or --gene")
    start = max(1, start - args.flank)
    if end is not None:
        end += args.flank
    classes = [c.strip() for c in args.classes.split(",") if c.strip()]
    bad = [c for c in classes if c not in CLASSES]
    if bad:
        ap.error(f"unknown class(es) {bad}; choose from {list(CLASSES)}")
    d = collect(b, seqid, start, end, classes, args.coverage_bins, not args.no_tabix)
    print(json.dumps(d, indent=1) if args.format == "json" else render_text(d))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
