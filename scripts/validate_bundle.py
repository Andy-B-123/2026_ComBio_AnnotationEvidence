#!/usr/bin/env python3
"""Validate an annotation-evidence bundle against docs/SPECIFICATION.md.

Checks (ERROR fails the run, WARN does not):
  * manifest.yaml against schema/manifest.schema.yaml (needs `jsonschema`; otherwise basic key checks)
  * listed files exist, and their size and sha256 match the manifest; tabix indexes are present
  * every seqid is in assembly_sequences.tsv, and no coordinate exceeds the sequence length
  * each file follows its format: columns, value ranges, internal consistency
    (e.g. BED12 blocks vs the intron_chain column, junction counts, protein block Target ranges)
  * files are sorted the way tabix requires
Exit status 1 if there is any ERROR.

Example:  validate_bundle.py examples/minimal/bundle
"""
from __future__ import annotations

import argparse
import datetime
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

SCHEMA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "schema", "manifest.schema.yaml")
MOTIFS = {"GT/AG", "GC/AG", "AT/AC", "non-canonical"}


class Report:
    def __init__(self, quiet=False):
        self.errors, self.warnings, self.passed = [], [], []
        self.quiet = quiet

    def error(self, msg):
        self.errors.append(msg)
        print(f"ERROR  {msg}")

    def warn(self, msg):
        self.warnings.append(msg)
        print(f"WARN   {msg}")

    def ok(self, msg):
        self.passed.append(msg)
        if not self.quiet:
            print(f"ok     {msg}")


class SortCheck:
    """tabix needs each seqid contiguous and starts non-decreasing within a seqid."""

    def __init__(self, rep, label):
        self.rep, self.label, self.seen, self.cur, self.last, self.bad = rep, label, set(), None, 0, 0

    def add(self, seqid, start, where):
        if seqid != self.cur:
            if seqid in self.seen:
                self._fail(f"{self.label}: seqid '{seqid}' is not contiguous ({where})")
            self.seen.add(seqid)
            self.cur, self.last = seqid, 0
        if start < self.last:
            self._fail(f"{self.label}: not sorted by start at {where} ({start} after {self.last}); tabix needs sorted input")
        self.last = max(self.last, start)

    def _fail(self, msg):
        self.bad += 1
        if self.bad <= 3:
            self.rep.error(msg)


def _plain(obj):
    if isinstance(obj, dict):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_plain(v) for v in obj]
    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()
    return obj


def check_manifest(b: ec.Bundle, rep: Report) -> None:
    m = _plain(b.manifest)
    try:
        import jsonschema
    except ImportError:
        rep.warn("jsonschema not installed; only checking that required top-level keys exist")
        for k in ("spec_version", "bundle", "assembly", "annotation", "evidence"):
            if k not in m:
                rep.error(f"manifest.yaml: missing required key '{k}'")
        return
    schema = ec.load_yaml(SCHEMA)
    errs = sorted(jsonschema.Draft202012Validator(schema).iter_errors(m), key=lambda e: list(e.absolute_path))
    for e in errs[:15]:
        rep.error(f"manifest.yaml: {'/'.join(str(p) for p in e.absolute_path) or '(root)'}: {e.message[:160]}")
    if not errs:
        rep.ok("manifest.yaml conforms to the schema")
    if m.get("spec_version") != ec.SPEC_VERSION:
        rep.warn(f"manifest spec_version {m.get('spec_version')} differs from this validator's {ec.SPEC_VERSION}")


def check_files(b: ec.Bundle, rep: Report) -> None:
    ann = b.manifest.get("annotation", {})
    entries = [(ann.get("filename"), ann.get("checksum", {}).get("value"), ann.get("bytes"))]
    for cls in b.classes():
        for o in b.outputs(cls):
            entries.append((o.get("file"), o.get("sha256"), o.get("bytes")))
    for rel, sha, size in entries:
        if not rel:
            continue
        p = b.path(rel)
        if not os.path.exists(p):
            rep.error(f"listed file is missing: {rel}")
            continue
        if size is not None and os.path.getsize(p) != size:
            rep.error(f"{rel}: size {os.path.getsize(p)} differs from manifest {size}")
        if sha and ec.sha256_file(p) != sha:
            rep.error(f"{rel}: sha256 differs from manifest (file changed after the manifest was written?)")
        if rel.endswith(".gz") and not (os.path.exists(p + ".tbi") or os.path.exists(p + ".csi")):
            rep.warn(f"{rel}: no .tbi/.csi index; region queries will scan the whole file")
    rep.ok(f"checked existence, size and sha256 of {sum(1 for e in entries if e[0])} file(s)")


def load_sequences(b: ec.Bundle, rep: Report):
    lens = b.sequences()
    aliases = {}
    p = b.path((b.manifest.get("assembly") or {}).get("sequences_file", "assembly_sequences.tsv"))
    if not lens:
        rep.error("assembly_sequences.tsv missing or empty: cannot check seqids or coordinate bounds")
        return {}, {}
    for r in ec.read_tsv_table(p)[1]:
        for a in (r.get("aliases") or ".").split(","):
            if a and a != ".":
                aliases[a] = r["seqid"]
    rep.ok(f"{len(lens)} sequences in assembly_sequences.tsv")
    return lens, aliases


class Coord:
    """Bounds/seqid checker that reports each kind of problem once per file."""

    def __init__(self, rep, label, lens, aliases):
        self.rep, self.label, self.lens, self.aliases, self.reported = rep, label, lens, aliases, set()

    def check(self, seqid, start, end, where):
        if not self.lens:
            return
        if seqid not in self.lens:
            key = ("seqid", seqid)
            if key not in self.reported:
                self.reported.add(key)
                hint = f" (it is an alias of '{self.aliases[seqid]}': use the primary name)" if seqid in self.aliases else ""
                self.rep.error(f"{self.label}: seqid '{seqid}' is not in assembly_sequences.tsv{hint} ({where})")
            return
        if start < 1 or end < start or end > self.lens[seqid]:
            key = ("bounds", seqid)
            if key not in self.reported:
                self.reported.add(key)
                self.rep.error(f"{self.label}: interval {start}-{end} on {seqid} (length {self.lens[seqid]}) is invalid ({where})")


def data_lines(path):
    with ec.open_text(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.strip() and not line.startswith("#"):
                yield n, line.rstrip("\n")


def check_pragma(path, profile, rep, label):
    pr = ec.read_pragmas(path)
    if pr.get("profile") != profile:
        rep.warn(f"{label}: header ##ae.profile is '{pr.get('profile')}', expected '{profile}'")
    return pr


def check_rna_junctions(b, o, lens, aliases, rep):
    p, label = b.path(o["file"]), o["file"]
    pr = check_pragma(p, "rna_junctions", rep, label)
    cols_ok = False
    with ec.open_text(p) as fh:
        for line in fh:
            if line.startswith("#") and not line.startswith("##"):
                cols_ok = line[1:].rstrip("\n").split("\t") == list(ec.JUNCTION_COLUMNS)
                break
    if not cols_ok:
        rep.error(f"{label}: column header must be '#{chr(9).join(ec.JUNCTION_COLUMNS)}'")
    cc, sc = Coord(rep, label, lens, aliases), SortCheck(rep, label)
    n_samples = int(pr["n_samples"]) if pr.get("n_samples", ".").isdigit() else None
    min_u = int(pr.get("sample_min_unique_reads", 1))
    min_tot = int(pr.get("min_total_unique_reads", 0))
    min_s = int(pr.get("min_samples", 1))
    seen, n, bad = set(), 0, 0
    for ln, line in data_lines(p):
        n += 1
        try:
            r = ec.parse_junction_line(line)
        except (ec.EvidenceError, ValueError) as exc:
            rep.error(f"{label}:{ln}: {exc}")
            bad += 1
            if bad > 5:
                break
            continue
        w = f"line {ln}"
        cc.check(r["seqid"], r["start"], r["end"], w)
        sc.add(r["seqid"], r["start"], w)
        key = (r["seqid"], r["start"], r["end"], r["strand"])
        if key in seen:
            rep.error(f"{label}:{ln}: duplicate junction {key}")
        seen.add(key)
        problems = []
        if r["strand"] not in "+-.":
            problems.append("strand must be + - or .")
        if r["motif"] not in MOTIFS:
            problems.append(f"motif '{r['motif']}' not one of {sorted(MOTIFS)}")
        if r["sample_count"] < 1:
            problems.append("sample_count < 1")
        if r["max_unique_reads"] > r["total_unique_reads"]:
            problems.append("max_unique_reads > total_unique_reads")
        if not (0 <= r["median_unique_reads"] <= r["max_unique_reads"]):
            problems.append("median_unique_reads outside 0..max")
        if r["total_unique_reads"] < r["sample_count"] * min_u:
            problems.append("total_unique_reads < sample_count * sample_min_unique_reads")
        if r["sample_count"] < min_s or r["total_unique_reads"] < min_tot:
            problems.append("row violates the min_samples / min_total_unique_reads filter in the header")
        if r["group_count"] is not None and r["group_count"] > r["sample_count"]:
            problems.append("group_count > sample_count")
        if n_samples and r["sample_count"] > n_samples:
            problems.append("sample_count > n_samples")
        for pmsg in problems:
            rep.error(f"{label}:{ln}: {pmsg}")
    if not bad:
        rep.ok(f"{label}: {n} junction rows")
    manifest_filters = b.filters("rna_short")
    for key, hkey in (("min_samples", "min_samples"), ("min_total_unique_reads", "min_total_unique_reads"),
                      ("min_sample_unique_reads", "sample_min_unique_reads")):
        if key in manifest_filters and hkey in pr and str(manifest_filters[key]) != pr[hkey]:
            rep.warn(f"{label}: header {hkey}={pr[hkey]} differs from manifest filters.{key}={manifest_filters[key]}")


def check_long_reads(b, o, lens, aliases, rep):
    p, label = b.path(o["file"]), o["file"]
    cc, sc, ids, n = Coord(rep, label, lens, aliases), SortCheck(rep, label), set(), 0
    for ln, line in data_lines(p):
        n += 1
        w = f"line {ln}"
        try:
            f = line.split("\t")
            if len(f) != 16:
                raise ec.EvidenceError(f"expected 16 columns (BED12+4), found {len(f)}")
            r = ec.parse_long_read_line(line)
            cs, ce = int(f[1]), int(f[2])
            sizes = [int(x) for x in f[10].strip(",").split(",")]
            starts = [int(x) for x in f[11].strip(",").split(",")]
        except (ec.EvidenceError, ValueError) as exc:
            rep.error(f"{label}:{ln}: {exc}")
            continue
        cc.check(r["seqid"], r["start"], r["end"], w)
        sc.add(r["seqid"], r["start"], w)
        errs = []
        if r["id"] in ids:
            errs.append(f"duplicate model id {r['id']}")
        ids.add(r["id"])
        if int(f[9]) != len(sizes) or len(sizes) != len(starts):
            errs.append("blockCount does not match blockSizes/blockStarts")
        if starts[0] != 0 or starts[-1] + sizes[-1] != ce - cs:
            errs.append("BED12 blocks must start at chromStart (first blockStart 0) and end at chromEnd")
        if any(starts[i] + sizes[i] >= starts[i + 1] for i in range(len(starts) - 1)):
            errs.append("blocks overlap, touch or are out of order")
        if r["strand"] not in "+-.":
            errs.append("strand must be + - or .")
        if not (cs <= int(f[6]) <= int(f[7]) <= ce):
            errs.append("thickStart/thickEnd outside chromStart..chromEnd")
        if f[15] != ec.format_chain(r["introns"]):
            errs.append(f"intron_chain '{f[15]}' does not match blocks ('{ec.format_chain(r['introns'])}')")
        if not r["introns"] and int(f[9]) != 1:
            errs.append("a model without introns must have exactly one block")
        if r["read_count"] < r["sample_count"]:
            errs.append("read_count < sample_count")
        if r["group_count"] is not None and r["group_count"] > r["sample_count"]:
            errs.append("group_count > sample_count")
        for e in errs:
            rep.error(f"{label}:{ln}: {e}")
        if r["score"] != min(r["read_count"], 1000):
            rep.warn(f"{label}:{ln}: score {r['score']} is not min(read_count, 1000)")
    rep.ok(f"{label}: {n} long-read models")


def check_protein(b, o, lens, aliases, rep):
    p, label = b.path(o["file"]), o["file"]
    check_pragma(p, "protein_alignments", rep, label)
    feats = [(n, ec.parse_gff3_line(l)) for n, l in data_lines(p)]
    cc, sc = Coord(rep, label, lens, aliases), SortCheck(rep, label)
    parents = {}
    for n, f in feats:
        cc.check(f.seqid, f.start, f.end, f"line {n}")
        sc.add(f.seqid, f.start, f"line {n}")
        if f.type == "protein_match":
            if "ID" not in f.attrs or "Target" not in f.attrs:
                rep.error(f"{label}:{n}: protein_match needs ID and Target")
                continue
            if f.attrs["ID"] in parents:
                rep.error(f"{label}:{n}: duplicate ID {f.attrs['ID']}")
            parents[f.attrs["ID"]] = f
            for k in ("pct_identity", "pct_coverage", "pct_positive"):
                v = f.attrs.get(k)
                if v is None:
                    rep.warn(f"{label}:{n}: attribute {k} missing")
                elif v != "." and not 0 <= float(v) <= 100:
                    rep.error(f"{label}:{n}: {k}={v} is outside 0-100 (percent, not a fraction?)")
        elif f.type != "match_part":
            rep.warn(f"{label}:{n}: unexpected feature type '{f.type}' (expected protein_match / match_part)")
    kids = {}
    for n, f in feats:
        if f.type == "match_part":
            par = f.attrs.get("Parent")
            if par not in parents:
                rep.error(f"{label}:{n}: match_part Parent '{par}' has no protein_match")
                continue
            pf = parents[par]
            if not (pf.start <= f.start and f.end <= pf.end) or f.strand != pf.strand:
                rep.error(f"{label}:{n}: match_part lies outside its parent or has a different strand")
            kids.setdefault(par, []).append(f)
    for pid, pf in parents.items():
        blocks = sorted(kids.get(pid, []), key=lambda x: x.start, reverse=(pf.strand == "-"))
        if not blocks:
            rep.error(f"{label}: protein_match {pid} has no match_part")
            continue
        prev_end = 0
        for f in blocks:
            _, ts, te = ec.parse_target(f.attrs.get("Target", ""))
            if ts is None or te < ts:
                rep.error(f"{label}: {pid}: bad Target on a match_part")
            elif ts < prev_end:
                rep.warn(f"{label}: {pid}: match_part Target ranges overlap or run backwards ({ts} after {prev_end})")
            prev_end = te if te else prev_end
        gs = min(x.start for x in blocks), max(x.end for x in blocks)
        if (gs[0], gs[1]) != (pf.start, pf.end):
            rep.warn(f"{label}: {pid}: parent span {pf.start}-{pf.end} differs from its blocks {gs[0]}-{gs[1]}")
    rep.ok(f"{label}: {len(parents)} protein alignments, {sum(len(v) for v in kids.values())} blocks")


def check_gff_models(b, o, lens, aliases, rep, profile, source_tags=None):
    p, label = b.path(o["file"]), o["file"]
    if profile:
        check_pragma(p, profile, rep, label)
    feats = [(n, ec.parse_gff3_line(l)) for n, l in data_lines(p)]
    cc, sc, ids = Coord(rep, label, lens, aliases), SortCheck(rep, label), {}
    for n, f in feats:
        cc.check(f.seqid, f.start, f.end, f"line {n}")
        sc.add(f.seqid, f.start, f"line {n}")
        if "ID" in f.attrs:
            prev = ids.get(f.attrs["ID"])
            # GFF3 allows one feature to span several lines with the same ID (e.g. the CDS lines of a
            # transcript) provided type, sequence and strand agree; anything else is a duplicate.
            if prev is not None and (prev.type, prev.seqid, prev.strand) != (f.type, f.seqid, f.strand):
                rep.error(f"{label}:{n}: duplicate ID {f.attrs['ID']}")
            ids.setdefault(f.attrs["ID"], f)
        if f.strand not in "+-.":
            rep.error(f"{label}:{n}: strand must be + - or .")
    for n, f in feats:
        for par in f.attrs.get("Parent", "").split(","):
            if par and par not in ids:
                rep.error(f"{label}:{n}: Parent '{par}' not found in this file")
            elif par and not (ids[par].start <= f.start and f.end <= ids[par].end):
                rep.error(f"{label}:{n}: feature lies outside its parent {par}")
    if source_tags is not None:
        used = {f.source for _, f in feats}
        missing = used - set(source_tags)
        if missing:
            rep.warn(f"{label}: source tag(s) {sorted(missing)} are not listed as source_tag under programs in the manifest")
    rep.ok(f"{label}: {len(feats)} features")


def check_bed6(b, o, lens, aliases, rep):
    p, label = b.path(o["file"]), o["file"]
    cc, sc, n = Coord(rep, label, lens, aliases), SortCheck(rep, label), 0
    for ln, line in data_lines(p):
        f = line.split("\t")
        n += 1
        if len(f) < 6:
            rep.error(f"{label}:{ln}: expected >= 6 columns")
            continue
        s, e = int(f[1]) + 1, int(f[2])
        cc.check(f[0], s, e, f"line {ln}")
        sc.add(f[0], s, f"line {ln}")
    rep.ok(f"{label}: {n} BED records")


def check_samples(b, rep):
    sp = b.path("samples.tsv")
    if not os.path.exists(sp):
        rep.warn("samples.tsv not found")
        return
    header, rows = ec.read_tsv_table(sp)
    for req in ("sample_id", "evidence_class", "group"):
        if req not in header:
            rep.error(f"samples.tsv: missing required column '{req}'")
    ids = [r.get("sample_id") for r in rows]
    if len(ids) != len(set(ids)):
        rep.error("samples.tsv: sample_id values are not unique")
    for cls, ev in b.classes().items():
        if ev.get("present") and cls in ("rna_short", "rna_long"):
            n = sum(1 for r in rows if r.get("evidence_class") == cls)
            if ev.get("n_datasets") is not None and ev["n_datasets"] != n:
                rep.warn(f"samples.tsv has {n} '{cls}' rows but the manifest says n_datasets={ev['n_datasets']}")
    rep.ok(f"samples.tsv: {len(rows)} rows")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bundle", help="bundle directory (contains manifest.yaml)")
    ap.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")
    args = ap.parse_args(argv)
    rep = Report(args.quiet)
    b = ec.Bundle(args.bundle)
    check_manifest(b, rep)
    check_files(b, rep)
    lens, aliases = load_sequences(b, rep)
    check_samples(b, rep)
    ann = {"file": b.manifest["annotation"]["filename"]}
    check_gff_models(b, ann, lens, aliases, rep, None)
    for cls in b.classes():
        for o in b.outputs(cls):
            if not os.path.exists(b.path(o["file"])):
                continue
            fmt = o.get("format")
            if fmt == "rna_junctions":
                check_rna_junctions(b, o, lens, aliases, rep)
            elif fmt == "long_read_models":
                check_long_reads(b, o, lens, aliases, rep)
            elif fmt == "protein_alignments":
                check_protein(b, o, lens, aliases, rep)
            elif fmt == "predictions":
                tags = [p.get("source_tag") for p in b.classes()[cls].get("programs", []) if p.get("source_tag")]
                check_gff_models(b, o, lens, aliases, rep, "predictions", tags or None)
            elif fmt == "repeats_bed6":
                check_bed6(b, o, lens, aliases, rep)
    print(f"\n{len(rep.errors)} error(s), {len(rep.warnings)} warning(s), {len(rep.passed)} check group(s) passed")
    return 1 if rep.errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
