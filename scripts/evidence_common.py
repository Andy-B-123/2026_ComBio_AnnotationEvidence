"""Shared helpers for the annotation-evidence scripts.

Standard library only, except PyYAML for reading ``manifest.yaml`` (imported lazily).
``tabix``/``bgzip`` (htslib) are used when present; everything also works without them.

COORDINATE CONVENTIONS (see docs/SPECIFICATION.md, section 2)
    * Text tables defined by this project and GFF3 : 1-based, closed   [start, end]
    * BED-family files (BED6/BED12, bigWig, bedGraph): 0-based, half-open [chromStart, chromEnd)
    * All functions in this module that return or accept ``start``/``end`` without saying
      otherwise use 1-based closed coordinates. BED values are always called ``chrom_start``/``chrom_end``.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict, namedtuple
from urllib.parse import unquote

SPEC_VERSION = "0.1"
PRAGMA_PREFIX = "##ae."

# tabix arguments per file kind
TABIX_PRESETS = {
    "tsv": ["-s1", "-b2", "-e3"],     # 1-based closed, columns seqid,start,end
    "gff3": ["-p", "gff"],            # 1-based closed, columns 1,4,5
    "bed": ["-p", "bed"],             # 0-based half-open, columns 1,2,3
}


class EvidenceError(Exception):
    """A problem with a bundle, a file or an argument (message is meant for the user)."""


# --------------------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------------------
def closed_to_bed(start: int, end: int) -> tuple[int, int]:
    """1-based closed [start, end] -> BED 0-based half-open (chrom_start, chrom_end)."""
    if start < 1 or end < start:
        raise ValueError(f"invalid 1-based closed interval {start}-{end}")
    return start - 1, end


def bed_to_closed(chrom_start: int, chrom_end: int) -> tuple[int, int]:
    """BED 0-based half-open [chrom_start, chrom_end) -> 1-based closed (start, end)."""
    if chrom_start < 0 or chrom_end <= chrom_start:
        raise ValueError(f"invalid BED interval {chrom_start}-{chrom_end}")
    return chrom_start + 1, chrom_end


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """True if two 1-based closed intervals share at least one base."""
    return a_start <= b_end and b_start <= a_end


def introns_from_exons(exons) -> list[tuple[int, int]]:
    """Introns (1-based closed, first..last intronic base) between consecutive exons.

    ``exons`` are 1-based closed intervals in any order. ``intron_start = exon_end + 1`` and
    ``intron_end = next_exon_start - 1``. Adjacent or overlapping exons produce no intron.
    """
    ex = sorted((int(s), int(e)) for s, e in exons)
    out = []
    for (_, e1), (s2, _) in zip(ex, ex[1:]):
        if s2 - e1 >= 2:
            out.append((e1 + 1, s2 - 1))
    return out


def exons_from_bed12(chrom_start: int, block_sizes, block_starts) -> list[tuple[int, int]]:
    """BED12 blocks -> exons as 1-based closed intervals.

    ``block_starts`` are relative to ``chrom_start`` (BED rule). Block i covers
    0-based half-open [chrom_start + start_i, chrom_start + start_i + size_i).
    """
    if len(block_sizes) != len(block_starts):
        raise ValueError("blockSizes and blockStarts differ in length")
    return [(chrom_start + st + 1, chrom_start + st + sz) for sz, st in zip(block_sizes, block_starts)]


def bed12_from_exons(exons) -> tuple[int, int, list[int], list[int]]:
    """Exons (1-based closed) -> (chrom_start, chrom_end, blockSizes, blockStarts) for BED12."""
    ex = sorted((int(s), int(e)) for s, e in exons)
    chrom_start = ex[0][0] - 1
    chrom_end = ex[-1][1]
    sizes = [e - s + 1 for s, e in ex]
    starts = [s - 1 - chrom_start for s, _ in ex]
    return chrom_start, chrom_end, sizes, starts


def format_chain(introns) -> str:
    """Intron chain as text: '1301-1700;1901-2500' (1-based closed); '.' if there are no introns."""
    return ";".join(f"{s}-{e}" for s, e in introns) if introns else "."


def parse_chain(text: str) -> list[tuple[int, int]]:
    if text in (".", ""):
        return []
    return [tuple(int(x) for x in part.split("-")) for part in text.split(";")]


_REGION_RE = re.compile(r"^(?P<seqid>.+?):(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?$")


def parse_region(text: str) -> tuple[str, int, int | None]:
    """Parse ``seqid``, ``seqid:start`` or ``seqid:start-end`` (1-based closed, commas allowed).

    Returns (seqid, start, end); ``end`` is None for a whole sequence or an open end.
    The split is on the LAST colon so that sequence names containing ':' still work.
    """
    text = text.strip()
    m = _REGION_RE.match(text)
    if not m:
        return text, 1, None
    start = int(m.group("start").replace(",", ""))
    end = m.group("end")
    end = int(end.replace(",", "")) if end else None
    if start < 1 or (end is not None and end < start):
        raise EvidenceError(f"invalid region '{text}': need 1 <= start <= end (1-based, closed)")
    return m.group("seqid"), start, end


# --------------------------------------------------------------------------------------
# File access
# --------------------------------------------------------------------------------------
def open_text(path: str):
    """Open a plain or gzip/bgzip text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_pragmas(path: str) -> dict[str, str]:
    """Return the ``##ae.key=value`` header lines of a TSV/GFF3 file as a dict (keys without prefix)."""
    out: dict[str, str] = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith(PRAGMA_PREFIX):
                key, _, value = line[len(PRAGMA_PREFIX):].rstrip("\n").partition("=")
                out[key] = value
    return out


def pragma_lines(**kv) -> list[str]:
    return [f"{PRAGMA_PREFIX}{k}={v}" for k, v in kv.items() if v is not None]


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def have_tool(name: str) -> bool:
    return shutil.which(name) is not None


def bgzip_and_index(plain_path: str, out_gz: str, kind: str) -> None:
    """Compress ``plain_path`` with bgzip to ``out_gz`` and tabix-index it. Needs htslib."""
    if not (have_tool("bgzip") and have_tool("tabix")):
        raise EvidenceError("bgzip/tabix (htslib) not found on PATH; install htslib, or write an uncompressed file")
    with open(out_gz, "wb") as fh:
        subprocess.run(["bgzip", "-c", plain_path], stdout=fh, check=True)
    subprocess.run(["tabix", "-f", *TABIX_PRESETS[kind], out_gz], check=True)


def _line_interval(line: str, kind: str) -> tuple[str, int, int]:
    f = line.split("\t", 5)
    if kind == "gff3":
        return f[0], int(f[3]), int(f[4])
    if kind == "bed":
        return f[0], int(f[1]) + 1, int(f[2])
    return f[0], int(f[1]), int(f[2])          # tsv


def fetch_lines(path: str, kind: str, seqid: str, start: int = 1, end: int | None = None,
                use_tabix: bool = True) -> list[str]:
    """Data lines of ``path`` overlapping the 1-based closed region (comment lines are skipped).

    Uses ``tabix`` when the file has a .tbi/.csi index and tabix is installed; otherwise scans
    the whole file. Both paths apply the same overlap rule, and tests check that they agree.
    """
    indexed = os.path.exists(path + ".tbi") or os.path.exists(path + ".csi")
    if use_tabix and indexed and have_tool("tabix"):
        region = seqid if (start == 1 and end is None) else f"{seqid}:{start}-{end if end is not None else 2**31 - 1}"
        if ":" in seqid:
            region = "{" + seqid + "}" + region[len(seqid):]
        res = subprocess.run(["tabix", path, region], capture_output=True, text=True)
        if res.returncode != 0:
            raise EvidenceError(f"tabix failed on {path}: {res.stderr.strip()}")
        return [ln for ln in res.stdout.split("\n") if ln]
    hi = end if end is not None else 2**62
    out = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            line = line.rstrip("\n")
            sid, s, e = _line_interval(line, kind)
            if sid == seqid and overlaps(s, e, start, hi):
                out.append(line)
    return out


def read_tsv_table(path: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a header-carrying TSV (``samples.tsv``, ``assembly_sequences.tsv``).

    Lines starting with '##' and blank lines are skipped. The first other line is the header; a
    single leading '#' on its first column name is removed. Returns (header, rows-as-dicts).
    """
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = [c[1:] if (i == 0 and c.startswith("#")) else c for i, c in enumerate(f)]
                continue
            rows.append(dict(zip(header, f + ["."] * (len(header) - len(f)))))
    return header or [], rows


# --------------------------------------------------------------------------------------
# GFF3
# --------------------------------------------------------------------------------------
Feature = namedtuple("Feature", "seqid source type start end score strand phase attrs")


def parse_attrs(col9: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if col9 in ("", "."):
        return out
    for part in col9.strip().split(";"):
        if part:
            k, _, v = part.partition("=")
            out[k.strip()] = unquote(v.strip())
    return out


def parse_gff3_line(line: str) -> Feature:
    f = line.rstrip("\n").split("\t")
    if len(f) != 9:
        raise EvidenceError(f"GFF3 line has {len(f)} columns, expected 9: {line[:80]!r}")
    score = None if f[5] == "." else float(f[5])
    return Feature(f[0], f[1], f[2], int(f[3]), int(f[4]), score, f[6], f[7], parse_attrs(f[8]))


def iter_gff3(path: str):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            yield parse_gff3_line(line)


def fetch_gff3(path: str, seqid: str, start: int = 1, end: int | None = None) -> list[Feature]:
    """Features overlapping a region, with parents' children completed.

    A region query only returns the child lines that overlap the region. Parent lines span all
    their children, so we widen the query to the union of returned extents until it is stable.
    """
    lo, hi = start, end
    for _ in range(4):
        feats = [parse_gff3_line(ln) for ln in fetch_lines(path, "gff3", seqid, lo, hi)]
        if not feats:
            return []
        new_lo = min(lo, min(f.start for f in feats))
        new_hi = None if hi is None else max(hi, max(f.end for f in feats))
        if new_lo == lo and new_hi == hi:
            return feats
        lo, hi = new_lo, new_hi
    return feats


# ---- gene-model builder (annotation and predictions) ------------------------------------
class Transcript:
    def __init__(self, feat: Feature, gene_id: str | None):
        self.id = feat.attrs.get("ID", f"{feat.seqid}:{feat.start}-{feat.end}")
        self.type = feat.type
        self.source = feat.source
        self.gene_id = gene_id
        self.seqid, self.start, self.end, self.strand = feat.seqid, feat.start, feat.end, feat.strand
        self.attrs = feat.attrs
        self.exons: list[tuple[int, int]] = []
        self.cds: list[tuple[int, int]] = []

    def finalise(self):
        self.exons = sorted(set(self.exons))
        self.cds = sorted(set(self.cds))
        if not self.exons and self.cds:
            self.exons = list(self.cds)          # CDS-only annotations (e.g. some predictors)
        if self.exons:
            self.start, self.end = self.exons[0][0], self.exons[-1][1]
        return self

    @property
    def introns(self):
        return introns_from_exons(self.exons)

    @property
    def partial(self) -> bool:
        return self.attrs.get("partial", "").lower() == "true"


class Gene:
    def __init__(self, feat: Feature):
        self.id = feat.attrs.get("ID", f"{feat.seqid}:{feat.start}-{feat.end}")
        self.attrs = feat.attrs
        self.source = feat.source
        self.seqid, self.start, self.end, self.strand = feat.seqid, feat.start, feat.end, feat.strand
        self.transcripts: list[Transcript] = []

    @property
    def names(self):
        keys = ("ID", "Name", "gene", "gene_id", "locus_tag", "Alias")
        return {v for k in keys for v in self.attrs.get(k, "").split(",") if v}


def build_gene_models(features) -> dict[str, Gene]:
    """Assemble gene -> transcript -> exon/CDS from GFF3 features (gene/mRNA/exon/CDS hierarchy).

    A feature is a gene if its type ends in 'gene'. Any other feature that has exon or CDS
    children is a transcript. Transcripts without a gene parent get a synthetic single-transcript gene.
    """
    features = list(features)                       # may be a generator; it is traversed several times
    by_id: dict[str, Feature] = {}
    children: dict[str, list[Feature]] = defaultdict(list)
    for f in features:
        if "ID" in f.attrs:
            by_id[f.attrs["ID"]] = f
        for p in f.attrs.get("Parent", "").split(","):
            if p:
                children[p].append(f)
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    for f in features:
        if f.type.endswith("gene") and "ID" in f.attrs:
            genes[f.attrs["ID"]] = Gene(f)
    for f in features:
        if f.type in ("exon", "CDS"):
            for p in f.attrs.get("Parent", "").split(","):
                if p in by_id and by_id[p].type not in ("exon", "CDS") and not by_id[p].type.endswith("gene"):
                    tf = by_id[p]
                    t = transcripts.get(p)
                    if t is None:
                        gene_id = next((g for g in tf.attrs.get("Parent", "").split(",") if g in genes), None)
                        t = transcripts[p] = Transcript(tf, gene_id)
                    (t.exons if f.type == "exon" else t.cds).append((f.start, f.end))
    for tid, t in transcripts.items():
        t.finalise()
        if t.gene_id is None:
            synth = Feature(t.seqid, t.source, "gene", t.start, t.end, None, t.strand, ".", {"ID": f"gene:{tid}", "Name": tid})
            genes.setdefault(synth.attrs["ID"], Gene(synth))
            t.gene_id = synth.attrs["ID"]
        genes[t.gene_id].transcripts.append(t)
    for g in genes.values():
        if g.transcripts:
            g.start = min(g.start, min(t.start for t in g.transcripts))
            g.end = max(g.end, max(t.end for t in g.transcripts))
    return genes


def find_gene(annotation_path: str, name: str) -> Gene:
    """Look a gene up by ID/Name/gene/locus_tag/Alias (linear scan of the annotation)."""
    hit = None
    for f in iter_gff3(annotation_path):
        if f.type.endswith("gene") and name in {v for k in ("ID", "Name", "gene", "gene_id", "locus_tag", "Alias")
                                                for v in f.attrs.get(k, "").split(",")}:
            hit = f
            break
    if hit is None:
        raise EvidenceError(f"gene '{name}' not found in {annotation_path}")
    feats = fetch_gff3(annotation_path, hit.seqid, hit.start, hit.end)
    genes = build_gene_models(feats)
    gid = hit.attrs.get("ID")
    if gid not in genes:
        raise EvidenceError(f"gene '{name}' has no transcript with exons in {annotation_path}")
    return genes[gid]


# --------------------------------------------------------------------------------------
# Evidence record parsers (all return dicts with 1-based closed start/end)
# --------------------------------------------------------------------------------------
JUNCTION_COLUMNS = ("seqid", "intron_start", "intron_end", "strand", "motif", "sample_count", "group_count",
                    "total_unique_reads", "max_unique_reads", "median_unique_reads",
                    "total_multimapping_reads", "max_overhang")
_JUNCTION_INT = {"sample_count", "total_unique_reads", "max_unique_reads", "total_multimapping_reads", "max_overhang"}


def _opt_int(x: str):
    return None if x == "." else int(x)


def parse_junction_line(line: str) -> dict:
    f = line.rstrip("\n").split("\t")
    if len(f) != len(JUNCTION_COLUMNS):
        raise EvidenceError(f"rna_junctions line has {len(f)} columns, expected {len(JUNCTION_COLUMNS)}: {line[:80]!r}")
    d = dict(zip(JUNCTION_COLUMNS, f))
    rec = {"kind": "junction", "seqid": d["seqid"], "start": int(d["intron_start"]), "end": int(d["intron_end"]),
           "strand": d["strand"], "motif": d["motif"], "group_count": _opt_int(d["group_count"]),
           "median_unique_reads": float(d["median_unique_reads"])}
    for k in _JUNCTION_INT:
        rec[k] = int(d[k])
    return rec


LONG_READ_EXTRA = ("read_count", "sample_count", "group_count", "intron_chain")


def parse_long_read_line(line: str) -> dict:
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        raise EvidenceError(f"BED12 line has {len(f)} columns, expected >= 12: {line[:80]!r}")
    chrom_start, chrom_end = int(f[1]), int(f[2])
    sizes = [int(x) for x in f[10].strip(",").split(",")]
    starts = [int(x) for x in f[11].strip(",").split(",")]
    exons = exons_from_bed12(chrom_start, sizes, starts)
    rec = {"kind": "long_read_model", "seqid": f[0], "start": chrom_start + 1, "end": chrom_end, "id": f[3],
           "score": int(float(f[4])), "strand": f[5], "exons": exons, "introns": introns_from_exons(exons)}
    if len(f) >= 16:
        rec.update(read_count=int(f[12]), sample_count=int(f[13]), group_count=_opt_int(f[14]), intron_chain=f[15])
    return rec


def build_protein_alignments(features) -> list[dict]:
    """Group ``protein_match`` parents and their ``match_part`` children into alignment dicts."""
    parents: dict[str, Feature] = {}
    parts: dict[str, list[Feature]] = defaultdict(list)
    for f in features:
        if f.type == "protein_match" and "ID" in f.attrs:
            parents[f.attrs["ID"]] = f
        elif f.type == "match_part":
            for p in f.attrs.get("Parent", "").split(","):
                parts[p].append(f)
    out = []
    for pid, pf in parents.items():
        protein, pstart, pend = parse_target(pf.attrs.get("Target", ""))
        blocks = []
        for c in sorted(parts.get(pid, []), key=lambda x: x.start):
            _, ts, te = parse_target(c.attrs.get("Target", ""))
            blocks.append({"start": c.start, "end": c.end, "target_start": ts, "target_end": te,
                           "pct_identity": _opt_float(c.attrs.get("pct_identity"))})
        exons = [(b["start"], b["end"]) for b in blocks]
        out.append({"kind": "protein_alignment", "id": pid, "seqid": pf.seqid, "start": pf.start, "end": pf.end,
                    "strand": pf.strand, "source": pf.source, "score": pf.score, "protein": protein,
                    "protein_start": pstart, "protein_end": pend, "blocks": blocks, "exons": exons,
                    "introns": introns_from_exons(exons), "attrs": pf.attrs})
    return sorted(out, key=lambda a: (a["start"], a["end"], a["id"]))


def _opt_float(x):
    return None if x in (None, "", ".") else float(x)


def parse_target(text: str):
    """GFF3 Target 'protein_id start end [strand]' -> (id, start, end). Missing -> (None, None, None)."""
    parts = text.split()
    if len(parts) < 3:
        return None, None, None
    return parts[0], int(parts[1]), int(parts[2])


# --------------------------------------------------------------------------------------
# bigWig (optional: pyBigWig, or UCSC bigWigToBedGraph)
# --------------------------------------------------------------------------------------
def bigwig_backend() -> str | None:
    try:
        import pyBigWig  # noqa: F401
        return "pyBigWig"
    except ImportError:
        return "bigWigToBedGraph" if have_tool("bigWigToBedGraph") else None


def bigwig_intervals(path: str, seqid: str, start: int, end: int) -> list[tuple[int, int, float]] | None:
    """Covered intervals of a bigWig within a 1-based closed region, as 1-based closed (start, end, value).

    Returns None if no reader is available (pyBigWig / UCSC bigWigToBedGraph). Uncovered bases are absent.
    """
    backend = bigwig_backend()
    if backend is None:
        return None
    s0, e0 = start - 1, end
    if backend == "pyBigWig":
        import pyBigWig
        bw = pyBigWig.open(path)
        try:
            if seqid not in bw.chroms():
                return []
            e0 = min(e0, bw.chroms()[seqid])
            return [(a + 1, b, float(v)) for a, b, v in (bw.intervals(seqid, s0, e0) or [])]
        finally:
            bw.close()
    res = subprocess.run(["bigWigToBedGraph", f"-chrom={seqid}", f"-start={s0}", f"-end={e0}", path, "/dev/stdout"],
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise EvidenceError(f"bigWigToBedGraph failed on {path}: {res.stderr.strip()}")
    out = []
    for line in res.stdout.splitlines():
        f = line.split("\t")
        out.append((int(f[1]) + 1, int(f[2]), float(f[3])))
    return out


def mean_depth(intervals, start: int, end: int) -> float:
    """Mean value over [start, end] (1-based closed) counting uncovered bases as 0."""
    total = 0.0
    for s, e, v in intervals:
        lo, hi = max(s, start), min(e, end)
        if hi >= lo:
            total += v * (hi - lo + 1)
    return total / (end - start + 1)


# --------------------------------------------------------------------------------------
# Bundle
# --------------------------------------------------------------------------------------
def load_yaml(path: str):
    try:
        import yaml
    except ImportError as exc:                      # pragma: no cover
        raise EvidenceError("PyYAML is required to read manifest.yaml (pip install pyyaml)") from exc
    with open(path) as fh:
        return yaml.safe_load(fh)


EVIDENCE_KIND = {"rna_junctions": "tsv", "long_read_models": "bed", "protein_alignments": "gff3",
                 "predictions": "gff3", "repeats_bed6": "bed", "annotation_gff3": "gff3"}


class Bundle:
    """A bundle directory: manifest.yaml plus the files it lists."""

    def __init__(self, root: str):
        self.root = os.path.abspath(root)
        mpath = os.path.join(self.root, "manifest.yaml")
        if not os.path.exists(mpath):
            raise EvidenceError(f"{root} is not an evidence bundle: manifest.yaml not found")
        self.manifest = load_yaml(mpath) or {}

    def path(self, rel: str) -> str:
        return os.path.join(self.root, rel)

    @property
    def annotation_path(self) -> str:
        return self.path(self.manifest["annotation"]["filename"])

    def outputs(self, cls: str) -> list[dict]:
        """Output file records of one evidence class (empty if the class is absent)."""
        ev = (self.manifest.get("evidence") or {}).get(cls) or {}
        if not ev.get("present"):
            return []
        return list(ev.get("outputs") or [])

    def classes(self) -> dict[str, dict]:
        return {k: v for k, v in (self.manifest.get("evidence") or {}).items() if isinstance(v, dict)}

    def output_of_format(self, fmt: str) -> str | None:
        for cls in self.classes():
            for o in self.outputs(cls):
                if o.get("format") == fmt:
                    return self.path(o["file"])
        return None

    def sequences(self) -> dict[str, int]:
        rel = (self.manifest.get("assembly") or {}).get("sequences_file", "assembly_sequences.tsv")
        p = self.path(rel)
        out = {}
        if os.path.exists(p):
            with open_text(p) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    f = line.rstrip("\n").split("\t")
                    out[f[0]] = int(f[1])
        return out

    # ---- typed region access (all coordinates 1-based closed) -------------------------------
    def _rows(self, fmt: str, kind: str, seqid: str, start: int, end: int | None, use_tabix: bool = True):
        path = self.output_of_format(fmt)
        return [] if path is None else fetch_lines(path, kind, seqid, start, end, use_tabix)

    def annotation_genes(self, seqid: str, start: int = 1, end: int | None = None) -> list[Gene]:
        hi = end if end is not None else 2**62
        genes = build_gene_models(fetch_gff3(self.annotation_path, seqid, start, end)).values()
        return sorted((g for g in genes if g.transcripts and overlaps(g.start, g.end, start, hi)), key=lambda g: (g.start, g.id))

    def junctions(self, seqid: str, start: int = 1, end: int | None = None, use_tabix: bool = True) -> list[dict]:
        return [parse_junction_line(ln) for ln in self._rows("rna_junctions", "tsv", seqid, start, end, use_tabix)]

    def long_read_models(self, seqid: str, start: int = 1, end: int | None = None, use_tabix: bool = True) -> list[dict]:
        return [parse_long_read_line(ln) for ln in self._rows("long_read_models", "bed", seqid, start, end, use_tabix)]

    def protein_alignments(self, seqid: str, start: int = 1, end: int | None = None) -> list[dict]:
        path = self.output_of_format("protein_alignments")
        if path is None:
            return []
        hi = end if end is not None else 2**62
        return [a for a in build_protein_alignments(fetch_gff3(path, seqid, start, end)) if overlaps(a["start"], a["end"], start, hi)]

    def prediction_genes(self, seqid: str, start: int = 1, end: int | None = None) -> list[Gene]:
        path = self.output_of_format("predictions")
        if path is None:
            return []
        hi = end if end is not None else 2**62
        genes = build_gene_models(fetch_gff3(path, seqid, start, end)).values()
        return sorted((g for g in genes if g.transcripts and overlaps(g.start, g.end, start, hi)), key=lambda g: (g.start, g.source, g.id))

    def repeats(self, seqid: str, start: int = 1, end: int | None = None) -> list[dict]:
        out = []
        for ln in self._rows("repeats_bed6", "bed", seqid, start, end):
            f = ln.split("\t")
            out.append({"kind": "repeat", "seqid": f[0], "start": int(f[1]) + 1, "end": int(f[2]), "name": f[3], "strand": f[5]})
        return out

    def coverage_intervals(self, seqid: str, start: int, end: int):
        path = self.output_of_format("coverage_bigwig")
        return None if path is None else bigwig_intervals(path, seqid, start, end)

    def filters(self, cls: str) -> dict:
        return (self.manifest.get("evidence", {}).get(cls) or {}).get("filters") or {}


def write_summary_json(path: str | None, tool: str, version: str, parameters: dict, filters: dict, stats: dict) -> None:
    """Provenance sidecar written by every producer script and read by make_manifest.py.

    ``parameters`` are the settings of this summarisation step; ``filters`` are the thresholds that
    removed records (these end up in manifest.yaml under ``evidence.<class>.filters``).
    """
    if not path:
        return
    import json
    with open(path, "w") as fh:
        json.dump({"tool": tool, "version": version, "parameters": parameters, "filters": filters, "stats": stats},
                  fh, indent=1, sort_keys=True)
        fh.write("\n")


def warn(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)
