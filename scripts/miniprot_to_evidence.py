#!/usr/bin/env python3
"""Convert ``miniprot --gff`` output to the ``protein_alignments`` GFF3 profile.

miniprot writes one ``mRNA`` line per alignment and one ``CDS`` line per aligned block (plus a
``stop_codon`` line), with the protein length in a preceding ``##PAF`` line. This script retypes
them so that they are not mistaken for annotated genes:

    mRNA  -> protein_match  (SO:0000349)  with Target, rank, pct_identity, pct_positive, pct_coverage, ...
    CDS   -> match_part     (SO:0000039)  one per aligned block, with its own Target range
    stop_codon lines are dropped (not an aligned residue); the parent span is the CDS span.

Input must come from ``miniprot --gff`` (with ``##PAF`` lines) to get query length and coverage.
By default only the best-ranked alignment of each protein (rank 1) is kept, as NCBI does; use
--max-rank 0 to keep all. Identity and positives are converted from fractions to percent.

Example:
  miniprot -Iut16 --gff genome.fa proteins.faa | miniprot_to_evidence.py --db UniRef90 --taxid 7147 \\
      --aligner-version 0.18 -o evidence/protein_alignments.gff3.gz
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from urllib.parse import quote

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"


def pct(x: str | None) -> str:
    return "." if x in (None, "") else f"{float(x) * 100:.2f}"


def convert(lines, db, taxid, max_rank, min_pct_identity, min_pct_coverage):
    """Yield (parent_feature_dict, [child dicts]) for every kept alignment."""
    stats = {"alignments_in": 0, "kept": 0, "dropped_rank": 0, "dropped_identity": 0, "dropped_coverage": 0,
             "no_query_length": 0}
    pending_qlen = None
    current = None                      # dict for the alignment being read
    out = []

    def flush():
        nonlocal current
        if current is None:
            return
        stats["alignments_in"] += 1
        parts = current["parts"]
        if not parts:
            current = None
            return
        ident = float(current["attrs"].get("Identity", "nan")) * 100
        qs, qe = ec.parse_target(current["attrs"].get("Target", ""))[1:]
        qlen = current["qlen"]
        cov = 100.0 * (qe - qs + 1) / qlen if (qlen and qs and qe) else None
        if max_rank and int(current["attrs"].get("Rank", 1)) > max_rank:
            stats["dropped_rank"] += 1
        elif ident == ident and ident < min_pct_identity:
            stats["dropped_identity"] += 1
        elif cov is not None and cov < min_pct_coverage:
            stats["dropped_coverage"] += 1
        else:
            if qlen is None:
                stats["no_query_length"] += 1
            stats["kept"] += 1
            out.append((current, cov))
        current = None

    for line in lines:
        if line.startswith("##PAF"):
            f = line.rstrip("\n").split("\t")
            pending_qlen = int(f[2])
            continue
        if line.startswith("#") or not line.strip():
            continue
        feat = ec.parse_gff3_line(line)
        if feat.type == "mRNA":
            flush()
            current = {"feat": feat, "attrs": feat.attrs, "qlen": pending_qlen, "parts": []}
            pending_qlen = None
        elif feat.type == "CDS" and current is not None:
            current["parts"].append(feat)
    flush()
    return out, stats


def render(alignments, source, db, taxid):
    """Return sorted GFF3 lines (parents before children at equal start)."""
    keyed = []
    for i, (cur, cov) in enumerate(sorted(alignments, key=lambda a: (a[0]["feat"].seqid, min(p.start for p in a[0]["parts"]),
                                                                       max(p.end for p in a[0]["parts"]), a[0]["feat"].attrs.get("ID", ""))), 1):
        m, at, parts = cur["feat"], cur["attrs"], cur["parts"]
        aid = f"PA{i:06d}"
        protein, qs, qe = ec.parse_target(at.get("Target", ""))
        span_s, span_e = min(p.start for p in parts), max(p.end for p in parts)
        attrs = [f"ID={aid}", f"Target={protein} {qs} {qe}", f"rank={at.get('Rank', '1')}",
                 f"pct_identity={pct(at.get('Identity'))}", f"pct_positive={pct(at.get('Positive'))}",
                 f"pct_coverage={'.' if cov is None else f'{cov:.2f}'}",
                 f"query_length={cur['qlen'] if cur['qlen'] else '.'}",
                 f"frameshifts={at.get('Frameshift', '0')}", f"stop_codons={at.get('StopCodon', '0')}"]
        if db:
            attrs.append(f"db={quote(db, safe='()/.-_:')}")
        if taxid:
            attrs.append(f"taxid={quote(str(taxid), safe='')}")
        score = "." if m.score is None else f"{m.score:g}"
        keyed.append(((m.seqid, span_s, 0, span_e),
                      "\t".join([m.seqid, source, "protein_match", str(span_s), str(span_e), score, m.strand, ".", ";".join(attrs)])))
        for p in parts:
            _, ts, te = ec.parse_target(p.attrs.get("Target", ""))
            cattrs = [f"Parent={aid}", f"Target={protein} {ts} {te}", f"pct_identity={pct(p.attrs.get('Identity'))}"]
            keyed.append(((m.seqid, p.start, 1, p.end),
                          "\t".join([m.seqid, source, "match_part", str(p.start), str(p.end), "." if p.score is None else f"{p.score:g}",
                                     p.strand, ".", ";".join(cattrs)])))
    keyed.sort(key=lambda kv: kv[0])
    return [ln for _, ln in keyed]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("gff", nargs="?", default="-", help="miniprot --gff output (default: stdin)")
    ap.add_argument("-o", "--output", required=True, help="output path; '.gz' means bgzip + tabix (needs htslib)")
    ap.add_argument("--db", help="name of the protein database/set, recorded as db=... on every alignment")
    ap.add_argument("--taxid", help="NCBI taxid of the protein set's source taxon, recorded as taxid=...")
    ap.add_argument("--aligner-version", default=".", help="miniprot version, recorded in the header")
    ap.add_argument("--source", default="miniprot", help="GFF3 column 2 (default miniprot)")
    ap.add_argument("--max-rank", type=int, default=1, help="keep alignments with rank <= N; 0 keeps all (default 1)")
    ap.add_argument("--min-pct-identity", type=float, default=0.0)
    ap.add_argument("--min-pct-coverage", type=float, default=0.0)
    ap.add_argument("--strip-prefix", help="remove this prefix from seqids (e.g. 'lcl|')")
    ap.add_argument("--summary-json", help="write provenance/statistics as JSON (read by make_manifest.py)")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    fh = sys.stdin if args.gff == "-" else ec.open_text(args.gff)
    lines = fh.readlines()
    if args.strip_prefix:
        lines = [ln[len(args.strip_prefix):] if ln.startswith(args.strip_prefix) else ln for ln in lines]
    alignments, stats = convert(lines, args.db, args.taxid, args.max_rank, args.min_pct_identity, args.min_pct_coverage)
    body = render(alignments, args.source, args.db, args.taxid)
    header = ["##gff-version 3"] + ec.pragma_lines(
        profile="protein_alignments", version=ec.SPEC_VERSION, coordinates="1-based-closed",
        generated_by=f"miniprot_to_evidence.py {__version__}", source_format="miniprot --gff",
        aligner=f"miniprot {args.aligner_version}", db=args.db or ".", taxid=args.taxid or ".",
        max_rank=args.max_rank, min_pct_identity=args.min_pct_identity, min_pct_coverage=args.min_pct_coverage)
    text = "\n".join(header + body) + "\n"

    if args.output.endswith(".gz"):
        with tempfile.NamedTemporaryFile("w", suffix=".gff3", delete=False, dir=os.path.dirname(os.path.abspath(args.output))) as tmp:
            tmp.write(text)
        try:
            ec.bgzip_and_index(tmp.name, args.output, "gff3")
        finally:
            os.unlink(tmp.name)
    else:
        with open(args.output, "w") as out:
            out.write(text)

    for k, v in stats.items():
        print(f"{k}: {v}", file=sys.stderr)
    if stats["no_query_length"]:
        ec.warn("some alignments had no ##PAF line: query_length and pct_coverage are '.'; run miniprot with --gff, not --gff-only")
    ec.write_summary_json(
        args.summary_json, "miniprot_to_evidence.py", __version__,
        parameters={"aligner": f"miniprot {args.aligner_version}", "source_format": "miniprot --gff", "db": args.db, "taxid": args.taxid},
        filters={"max_rank": args.max_rank, "min_pct_identity": args.min_pct_identity,
                 "min_pct_coverage": args.min_pct_coverage}, stats=stats)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
