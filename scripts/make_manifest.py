#!/usr/bin/env python3
"""Build the provenance files of an evidence bundle.

Two subcommands:

  sequences   FASTA -> assembly_sequences.tsv  (seqid, length, GA4GH refget digest, aliases)
  build       hand-written provenance template + the files in the bundle -> manifest.yaml

``build`` never invents provenance. You write what only you know (data sources, programs and
versions, databases, why an evidence class is absent) in a template, see schema/manifest.template.yaml.
The script adds what can be computed: file sizes, sha256, record counts, sample counts, the
sequence table summary, and the filters/statistics reported by the producer scripts
(aggregate_star_junctions.py, collapse_long_reads.py, miniprot_to_evidence.py --summary-json).

Examples:
  make_manifest.py sequences --fasta genome.fa --alias-prefix 'lcl|' -o bundle/assembly_sequences.tsv
  make_manifest.py build bundle/ --template manifest.template.yaml --fasta genome.fa \\
      --summary rna_short=agg.json --summary rna_long=lr.json --summary proteins=prot.json
"""
from __future__ import annotations

import argparse
import base64
import copy
import datetime
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"

COORDINATES = {"rna_junctions": "1-based-closed", "long_read_models": "0-based-half-open",
               "protein_alignments": "1-based-closed", "predictions": "1-based-closed",
               "repeats_bed6": "0-based-half-open", "coverage_bigwig": "0-based-half-open",
               "coverage_bedgraph": "0-based-half-open", "annotation_gff3": "1-based-closed"}
SAMPLE_CLASS = {"rna_short": "rna_short", "rna_long": "rna_long"}


def sha512t24u(hasher) -> str:
    """GA4GH refget sequence digest (base64url of the first 24 bytes of SHA-512)."""
    return base64.urlsafe_b64encode(hasher.digest()[:24]).decode("ascii")


def fasta_records(path):
    """Yield (name, length, sha512t24u) streaming through a FASTA (plain or .gz)."""
    name, length, h = None, 0, None
    with ec.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, length, sha512t24u(h)
                name, length, h = line[1:].split()[0], 0, hashlib.sha512()
            elif name is not None:
                s = "".join(line.split()).upper()
                length += len(s)
                h.update(s.encode("ascii"))
    if name is not None:
        yield name, length, sha512t24u(h)


def cmd_sequences(args) -> int:
    aliases = {}
    if args.aliases:
        for line in open(args.aliases):
            if line.strip() and not line.startswith("#"):
                sid, _, rest = line.rstrip("\n").partition("\t")
                aliases[sid] = rest.strip()
    rows = []
    for name, length, digest in fasta_records(args.fasta):
        al = aliases.get(name) or (args.alias_prefix + name if args.alias_prefix else ".")
        rows.append((name, length, digest, al or "."))
    if not rows:
        raise ec.EvidenceError(f"no sequences found in {args.fasta}")
    with open(args.output, "w") as out:
        out.write("##ae.profile=assembly_sequences\n##ae.version=" + ec.SPEC_VERSION + "\n")
        out.write("#seqid\tlength\tsha512t24u\taliases\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} sequences, {sum(r[1] for r in rows)} bp -> {args.output}", file=sys.stderr)
    return 0


def count_records(path: str) -> int | None:
    if path.endswith((".bw", ".bigWig", ".bigwig")):
        return None
    n = 0
    with ec.open_text(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def file_facts(bundle: str, rel: str, fmt: str | None) -> dict:
    p = os.path.join(bundle, rel)
    if not os.path.exists(p):
        raise ec.EvidenceError(f"file listed in the template does not exist: {p}")
    facts = {"bytes": os.path.getsize(p), "sha256": ec.sha256_file(p)}
    rec = count_records(p)
    if rec is not None:
        facts["records"] = rec
    for ext in (".tbi", ".csi"):
        if os.path.exists(p + ext):
            facts["index"] = rel + ext
            break
    if fmt in COORDINATES:
        facts["coordinates"] = COORDINATES[fmt]
    return facts


def merge_summary(ev: dict, summary: dict) -> None:
    """Fold one producer summary (--summary-json of a script) into an evidence class.

    Filters are merged, the program is appended once, and statistics are stored per tool
    (``stats: {tool: {...}}``). Several summaries of the same tool, for instance one per sample, add up their counts.
    """
    ev.setdefault("filters", {}).update(summary.get("filters") or {})
    programs = ev.setdefault("programs", [])
    if not any(p.get("name") == summary["tool"] for p in programs):
        programs.append({"name": summary["tool"], "version": summary["version"], "role": "summarisation",
                         "parameters": json.dumps(summary.get("parameters") or {}, sort_keys=True)})
    cur = ev.setdefault("stats", {}).setdefault(summary["tool"], {})
    for k, v in (summary.get("stats") or {}).items():
        both_numbers = isinstance(v, (int, float)) and isinstance(cur.get(k), (int, float))
        cur[k] = cur[k] + v if both_numbers else v


def cmd_build(args) -> int:
    bundle = args.bundle
    m = copy.deepcopy(ec.load_yaml(args.template))
    m["spec_version"] = ec.SPEC_VERSION
    b = m.setdefault("bundle", {})
    b["created"] = args.created or datetime.date.today().isoformat()
    b["generator"] = f"make_manifest.py {__version__}"

    asm = m.setdefault("assembly", {})
    seqfile = os.path.join(bundle, "assembly_sequences.tsv")
    if os.path.exists(seqfile):
        header, rows = ec.read_tsv_table(seqfile)
        asm["sequences_file"] = "assembly_sequences.tsv"
        asm["n_sequences"] = len(rows)
        asm["total_length"] = sum(int(r["length"]) for r in rows)
    if args.fasta:
        asm.setdefault("filename", os.path.basename(args.fasta))
        asm["checksum"] = {"algorithm": "sha256", "value": ec.sha256_file(args.fasta)}

    ann = m.setdefault("annotation", {})
    ann.update(file_facts(bundle, ann["filename"], "annotation_gff3"))
    ann["checksum"] = {"algorithm": "sha256", "value": ann.pop("sha256")}

    samples = []
    sp = os.path.join(bundle, "samples.tsv")
    if os.path.exists(sp):
        samples = ec.read_tsv_table(sp)[1]

    summaries: dict[str, list[dict]] = {}
    for item in args.summary or []:
        cls, _, path = item.partition("=")
        with open(path) as fh:
            summaries.setdefault(cls, []).append(json.load(fh))

    for cls, ev in (m.get("evidence") or {}).items():
        if not isinstance(ev, dict) or not ev.get("present"):
            continue
        if "n_datasets" not in ev and cls in SAMPLE_CLASS:
            ev["n_datasets"] = sum(1 for r in samples if r.get("evidence_class") == SAMPLE_CLASS[cls])
        for out in ev.get("outputs") or []:
            out.update(file_facts(bundle, out["file"], out.get("format")))
        for s in summaries.get(cls, []):
            merge_summary(ev, s)

    header = ("# annotation-evidence manifest, spec %s. Provenance fields come from the template;\n"
              "# file facts (bytes, sha256, records) and producer filters were computed by make_manifest.py.\n"
              % ec.SPEC_VERSION)
    import yaml
    text = yaml.safe_dump(m, sort_keys=False, default_flow_style=False, width=110, allow_unicode=True)
    out_path = args.output or os.path.join(bundle, "manifest.yaml")
    with open(out_path, "w") as fh:
        fh.write(header + text)
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sequences", help="FASTA -> assembly_sequences.tsv")
    s.add_argument("--fasta", required=True)
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--alias-prefix", help="add PREFIX+seqid as an alias for every sequence (e.g. 'lcl|')")
    s.add_argument("--aliases", help="TSV: seqid <TAB> comma-separated aliases (e.g. from an NCBI assembly report)")
    s.set_defaults(func=cmd_sequences)
    b = sub.add_parser("build", help="template + bundle files -> manifest.yaml")
    b.add_argument("bundle", help="bundle directory")
    b.add_argument("--template", required=True, help="hand-written provenance YAML (see schema/manifest.template.yaml)")
    b.add_argument("--fasta", help="assembly FASTA, to record its sha256 (the FASTA itself is not part of the bundle)")
    b.add_argument("--summary", action="append", metavar="CLASS=JSON",
                   help="producer summary for an evidence class, e.g. rna_short=agg.json (repeatable, also several per class)")
    b.add_argument("--created", help="creation date YYYY-MM-DD (default: today)")
    b.add_argument("-o", "--output", help="default: BUNDLE/manifest.yaml")
    b.set_defaults(func=cmd_build)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
