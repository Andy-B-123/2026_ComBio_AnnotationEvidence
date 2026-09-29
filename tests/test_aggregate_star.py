"""aggregate_star_junctions.py: semantics, and agreement with an independent computation on REAL STAR output."""
import glob
import os
import statistics
import subprocess
import sys

import pytest

import aggregate_star_junctions as agg
import evidence_common as ec
from conftest import DATA, SCRIPTS, needs_htslib

REAL = sorted(glob.glob(os.path.join(DATA, "real_star", "*.SJ.out.tab")))


def write(path, rows):
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def run(files, groups=None, **kw):
    sj = {agg.sample_id_from_path(f): f for f in files}
    return agg.aggregate(sj, groups, kw.get("min_sample_unique", 1), kw.get("strip_prefix"))


# ---- semantics on small hand-made inputs ---------------------------------------------------
def test_counts_distinguish_one_library_from_many(tmp_path):
    """100 reads in one library vs 20 reads spread over 4 libraries."""
    j1, j2 = ("c", 100, 200, 1, 1, 0), ("c", 300, 400, 1, 1, 0)
    files = [write(tmp_path / f"s{i}.SJ.out.tab", ([j1 + (100 if i == 0 else 0, 0, 40)] if i == 0 else []) + [j2 + (5, 1, 30)])
             for i in range(4)]
    files[0] = write(tmp_path / "s0.SJ.out.tab", [j1 + (100, 0, 40), j2 + (5, 1, 30)])
    rows, _ = run(files)
    by = {(r["intron_start"], r["intron_end"]): r for r in rows}
    one, many = by[(100, 200)], by[(300, 400)]
    assert (one["sample_count"], one["total_unique_reads"], one["max_unique_reads"]) == (1, 100, 100)
    assert (many["sample_count"], many["total_unique_reads"], many["max_unique_reads"]) == (4, 20, 5)
    assert many["total_multimapping_reads"] == 4 and many["median_unique_reads"] == "5"


def test_star_codes_are_translated(tmp_path):
    rows_in = [("c", 10, 20, 1, 1, 0, 3, 0, 20), ("c", 30, 40, 2, 2, 0, 3, 0, 20), ("c", 50, 60, 1, 3, 0, 3, 0, 20),
               ("c", 70, 80, 2, 4, 0, 3, 0, 20), ("c", 90, 99, 1, 5, 0, 3, 0, 20), ("c", 110, 120, 2, 6, 0, 3, 0, 20),
               ("c", 130, 140, 0, 0, 0, 3, 0, 20)]
    rows, _ = run([write(tmp_path / "a.SJ.out.tab", rows_in)])
    got = [(r["strand"], r["motif"]) for r in rows]
    assert got == [("+", "GT/AG"), ("-", "GT/AG"), ("+", "GC/AG"), ("-", "GC/AG"), ("+", "AT/AC"), ("-", "AT/AC"), (".", "non-canonical")]


def test_coordinates_pass_through_unchanged(tmp_path):
    """STAR already uses first/last intron base, 1-based: no conversion must happen."""
    rows, _ = run([write(tmp_path / "a.SJ.out.tab", [("c", 1301, 1700, 1, 1, 0, 2, 0, 30)])])
    assert (rows[0]["intron_start"], rows[0]["intron_end"]) == (1301, 1700)


def test_min_sample_unique_and_multimapper_only_rows(tmp_path):
    f1 = write(tmp_path / "a.SJ.out.tab", [("c", 10, 20, 1, 1, 0, 0, 7, 20), ("c", 30, 40, 1, 1, 0, 1, 0, 20)])
    f2 = write(tmp_path / "b.SJ.out.tab", [("c", 30, 40, 1, 1, 0, 4, 2, 25)])
    rows, st = run([f1, f2])                      # default: a sample needs >= 1 unique read
    assert [(r["intron_start"], r["sample_count"], r["total_unique_reads"], r["total_multimapping_reads"]) for r in rows] == [(30, 2, 5, 2)]
    assert st["rows_below_sample_threshold"] == 1
    rows0, _ = run([f1, f2], min_sample_unique=0)  # keep multimapper-only rows
    assert {r["intron_start"] for r in rows0} == {10, 30}


def test_groups_are_counted(tmp_path):
    files = [write(tmp_path / f"{s}.SJ.out.tab", [("c", 10, 20, 1, 1, 0, 2, 0, 20)]) for s in "abc"]
    rows, _ = run(files, groups={"a": "head", "b": "head", "c": "gut"})
    assert rows[0]["sample_count"] == 3 and rows[0]["group_count"] == 2
    rows, _ = run(files, groups=None)
    assert rows[0]["group_count"] == "."


def test_strip_prefix(tmp_path):
    rows, _ = run([write(tmp_path / "a.SJ.out.tab", [("lcl|chr1", 10, 20, 1, 1, 0, 2, 0, 20)])], strip_prefix="lcl|")
    assert rows[0]["seqid"] == "chr1"


@pytest.mark.parametrize("line", ["c\t10\t20\t1\t1\t0\t2\t0", "c\t10\tx\t1\t1\t0\t2\t0\t20", "c\t10\t20\t3\t1\t0\t2\t0\t20", "c\t30\t20\t1\t1\t0\t2\t0\t20"])
def test_malformed_input_is_rejected(tmp_path, line):
    p = tmp_path / "bad.SJ.out.tab"
    p.write_text(line + "\n")
    with pytest.raises(ec.EvidenceError):
        run([str(p)])


def test_global_filters_via_cli(tmp_path):
    f1 = write(tmp_path / "a.SJ.out.tab", [("c", 10, 20, 1, 1, 0, 1, 0, 20), ("c", 30, 40, 1, 1, 0, 9, 0, 20)])
    f2 = write(tmp_path / "b.SJ.out.tab", [("c", 30, 40, 1, 1, 0, 1, 0, 20)])
    out = tmp_path / "o.tsv"
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "aggregate_star_junctions.py"), "--min-samples", "2", "-o", str(out), f1, f2],
                   check=True, capture_output=True)
    text = out.read_text()
    assert "##ae.min_samples=2" in text
    data = [l for l in text.splitlines() if not l.startswith("#")]
    assert len(data) == 1 and data[0].split("\t")[1] == "30"


# ---- real STAR output: compare with an independent re-implementation ---------------------------
def independent(files):
    acc = {}
    for f in files:
        for line in open(f):
            c = line.split("\t")
            if int(c[6]) < 1:
                continue
            key = (c[0], int(c[1]), int(c[2]), {"0": ".", "1": "+", "2": "-"}[c[3]])
            a = acc.setdefault(key, {"u": [], "m": 0, "o": 0})
            a["u"].append(int(c[6]))
            a["m"] += int(c[7])
            a["o"] = max(a["o"], int(c[8]))
    return {k: (len(a["u"]), sum(a["u"]), max(a["u"]), statistics.median(a["u"]), a["m"], a["o"]) for k, a in acc.items()}


def test_real_star_fixture_is_present():
    assert len(REAL) == 3


def test_real_star_matches_independent_computation():
    rows, stats = run(REAL)
    want = independent(REAL)
    got = {(r["seqid"], r["intron_start"], r["intron_end"], r["strand"]):
           (r["sample_count"], r["total_unique_reads"], r["max_unique_reads"], float(r["median_unique_reads"]),
            r["total_multimapping_reads"], r["max_overhang"]) for r in rows}
    assert got == want
    assert stats["strand_motif_inconsistent"] == 0 and stats["motif_conflicts"] == 0
    assert any(r["motif"] == "AT/AC" for r in rows) and any(r["strand"] == "-" for r in rows)


def test_real_star_multimapper_only_rows_are_dropped_by_default():
    n_multi_only = sum(1 for f in REAL for l in open(f) if int(l.split("\t")[6]) == 0)
    assert n_multi_only > 0                      # the fixture really contains them
    _, stats = run(REAL)
    assert stats["rows_below_sample_threshold"] == n_multi_only


@needs_htslib
def test_bgzip_output_is_tabix_indexed_and_queryable(tmp_path):
    out = tmp_path / "rna_junctions.tsv.gz"
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "aggregate_star_junctions.py"), "-o", str(out), *REAL], check=True, capture_output=True)
    assert os.path.exists(str(out) + ".tbi")
    first = [l for l in ec.open_text(str(out)) if not l.startswith("#")][0].split("\t")
    lines = ec.fetch_lines(str(out), "tsv", first[0], int(first[1]), int(first[2]))
    assert any(l.split("\t")[1] == first[1] for l in lines)
    assert ec.read_pragmas(str(out))["profile"] == "rna_junctions"
