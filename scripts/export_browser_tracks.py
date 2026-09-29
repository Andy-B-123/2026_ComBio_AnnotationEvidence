#!/usr/bin/env python3
"""Write browser-friendly views of a bundle for IGV (desktop) and other genome browsers.

Most bundle files load directly (annotation, long-read models, predictions, repeats, coverage). Two do not,
so this script derives them (the bundle itself is not changed):

  junctions.bed            RNA-seq junctions in the TopHat junction format that IGV draws as a native splice-junction
                           track: one BED12 feature per junction with two anchor blocks around the intron (--anchor bases
                           each) and a `track name=junctions` header. Arc height follows the score = total unique reads
                           (not capped at 1000, which is fine for IGV but not for the UCSC browser). The name is
                           s<samples>/g<groups>/u<unique reads>; hover to read it.
  protein_alignments.bed   BED12, one feature per alignment, exon blocks visible. IGV draws a GFF3
                           protein_match with match_part children as one unbroken bar, so blocks are lost; the
                           protein-side coordinates stay in the bundle's GFF3 and in the BED12 name.
  igv_session.xml          an IGV session that loads everything (paths relative to the session file).

Coordinates: BED is 0-based half-open, converted from the bundle's 1-based closed values (the intron's first and
last base are the bases just inside the two anchors).

Example:
  export_browser_tracks.py bundle/ browser/ --genome genome.fa --locus scaffold_1:900-5400
  then open browser/igv_session.xml in IGV (File > Open Session)
"""
from __future__ import annotations

import argparse
import os
import sys
from xml.sax.saxutils import quoteattr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"


JUNCTION_HEADER = 'track name=junctions description="TopHat junctions"'


def junction_bed12(b: ec.Bundle, anchor: int = 10) -> list[str]:
    """Junction rows in TopHat junction BED12 form (data lines only; the track header is written by the caller)."""
    path = b.output_of_format("rna_junctions")
    if not path:
        return []
    lens = b.sequences()
    rows = []
    for line in ec.open_text(path):
        if line.startswith("#") or not line.strip():
            continue
        j = ec.parse_junction_line(line)
        cs = max(0, j["start"] - 1 - anchor)                    # anchor block ends at the last exon base = intron_start - 1 (0-based half-open)
        ce = min(lens.get(j["seqid"], j["end"] + anchor), j["end"] + anchor)
        left, right = j["start"] - 1 - cs, ce - j["end"]
        gc = "?" if j["group_count"] is None else j["group_count"]
        rows.append((j["seqid"], cs, "\t".join(map(str, (
            j["seqid"], cs, ce, f"s{j['sample_count']}/g{gc}/u{j['total_unique_reads']}", j["total_unique_reads"],
            j["strand"] if j["strand"] in "+-" else ".", cs, ce, "255,0,0", 2, f"{left},{right}", f"0,{j['end'] - cs}")))))
    return [r[2] for r in sorted(rows, key=lambda r: (r[0], r[1]))]


def protein_bed12(b: ec.Bundle) -> list[str]:
    path = b.output_of_format("protein_alignments")
    if not path:
        return []
    rows = []
    for a in ec.build_protein_alignments(ec.iter_gff3(path)):
        cs, ce, sizes, starts = ec.bed12_from_exons(a["exons"])
        ident = a["attrs"].get("pct_identity", ".")
        try:
            score = min(1000, int(round(float(ident) * 10)))
        except ValueError:
            score = 0
        name = f"{a['protein']}:{a['protein_start']}-{a['protein_end']} id={ident}%"
        rows.append((a["seqid"], cs, "\t".join(map(str, (
            a["seqid"], cs, ce, name, score, a["strand"], cs, ce, "0,0,0", len(sizes),
            ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ",")))))
    return [r[2] for r in sorted(rows, key=lambda r: (r[0], r[1]))]


def write_lines(path: str, lines: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + ("\n" if lines else ""))


def session_xml(genome: str, locus: str, resources: list[tuple[str, str]]) -> str:
    res = "\n".join(f"        <Resource name={quoteattr(name)} path={quoteattr(p)}/>" for name, p in resources)
    return (f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
            f'<Session genome={quoteattr(genome)} hasGeneTrack="false" hasSequenceTrack="true" locus={quoteattr(locus)} version="8">\n'
            f'    <Resources>\n{res}\n    </Resources>\n</Session>\n')


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bundle")
    ap.add_argument("outdir", help="directory for the derived files and the IGV session (created if needed)")
    ap.add_argument("--genome", help="assembly FASTA (indexed with .fai) or an IGV genome id; written into the session")
    ap.add_argument("--anchor", type=int, default=10, help="bases of exon drawn on each side of a junction (default 10)")
    ap.add_argument("--locus", help="initial locus for the session, e.g. scaffold_1:900-5400 (default: first sequence)")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)

    b = ec.Bundle(args.bundle)
    os.makedirs(args.outdir, exist_ok=True)
    rel = lambda p: os.path.relpath(p, args.outdir)
    resources: list[tuple[str, str]] = []
    resources.append(("annotation", rel(b.annotation_path)))
    order = [("coverage_bigwig", None), ("rna_junctions", None), ("long_read_models", None),
             ("protein_alignments", None), ("predictions", None), ("repeats_bed6", None)]
    written = []
    for fmt, derived in order:
        if fmt == "rna_junctions":
            lines = junction_bed12(b, args.anchor)
            if lines:
                write_lines(os.path.join(args.outdir, "junctions.bed"), [JUNCTION_HEADER] + lines)
                resources.append(("junctions", "junctions.bed"))
                written.append(("junctions.bed", len(lines)))
        elif fmt == "protein_alignments":
            lines = protein_bed12(b)
            if lines:
                write_lines(os.path.join(args.outdir, "protein_alignments.bed"), lines)
                resources.append(("protein alignments (BED12 view)", "protein_alignments.bed"))
                written.append(("protein_alignments.bed", len(lines)))
        else:
            p = b.output_of_format(fmt)
            if p:
                resources.append((os.path.basename(p), rel(p)))
    lens = b.sequences()
    locus = args.locus or (f"{next(iter(lens))}:1-{min(5000, next(iter(lens.values())))}" if lens else "")
    genome = args.genome
    if genome and os.path.exists(genome):
        genome = rel(os.path.abspath(genome))
    with open(os.path.join(args.outdir, "igv_session.xml"), "w") as fh:
        fh.write(session_xml(genome or "", locus, resources))
    for name, n in written:
        print(f"wrote {name}: {n} records", file=sys.stderr)
    print(f"wrote igv_session.xml ({len(resources)} tracks)", file=sys.stderr)
    if not genome:
        ec.warn("no --genome given: set the genome in IGV (Genomes > Load Genome from File) before opening the session")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
