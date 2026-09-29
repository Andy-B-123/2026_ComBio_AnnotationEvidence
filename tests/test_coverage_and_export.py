"""make_coverage.py and export_browser_tracks.py."""
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

import evidence_common as ec
import export_browser_tracks as ebt
import make_coverage as mc
from conftest import EXAMPLE, SCRIPTS


# ---- quantisation: the ladder 2^(k/L) - 1 ----------------------------------------------------------
def test_snap_examples_and_identity_cases():
    assert mc.snap(0, 4) == 0 and mc.snap(5.5, 0) == 5.5
    assert mc.snap(750, 4) == pytest.approx(723.08, abs=0.01)
    assert mc.snap(800, 4) == pytest.approx(860.08, abs=0.01)
    assert mc.snap(1, 4) == pytest.approx(1.0)                          # 2^1 - 1 is on the ladder


def test_snap_is_idempotent_monotone_and_within_half_a_step():
    xs = [0.5, 1, 2, 3, 7, 10, 33, 100, 750, 4321, 100000]
    for L in (2, 4, 8):
        snapped = [mc.snap(x, L) for x in xs]
        assert all(abs(mc.snap(v, L) - v) < 1e-6 * max(1, v) for v in snapped)                # ladder values map to themselves
        assert snapped == sorted(snapped)
        half_step = 2 ** (1 / (2 * L)) - 1
        assert all(abs((s + 1) / (x + 1) - 1) <= half_step + 1e-9 for s, x in zip(snapped, xs))


def test_summed_records_adds_columns_and_strips_prefix():
    lines = ["lcl|c\t0\t10\t1\t2\t0.5\n", "lcl|c\t10\t20\t0\t1\t0\n"]
    assert list(mc.summed_records(lines, "lcl|")) == [("c", 0, 10, 3.5), ("c", 10, 20, 1.0)]


def test_windowed_means_include_uncovered_bases_and_clip_the_last_window():
    recs = [("c", 0, 10, 4.0), ("c", 30, 35, 10.0)]
    out = list(mc.windowed(iter(recs), 20, {"c": 35}))
    assert out == [("c", 0, 20, 2.0), ("c", 20, 35, (5 * 10.0) / 15)]         # first window: 10 bases at 4 over 20; second: 5 bases at 10 over 15


def test_quantised_merges_adjacent_equal_values_drops_zeros_and_clips():
    recs = [("c", 0, 10, 1.0), ("c", 10, 20, 1.0), ("c", 20, 30, 0.0), ("c", 30, 60, 3.0), ("d", 0, 5, 1.0)]
    lines = list(mc.quantised(iter(recs), 4, {"c": 50}))
    assert lines == ["c\t0\t20\t1.00", "c\t30\t50\t3.00"]                    # merged; zero dropped; clipped at 50; unknown sequence dropped
    assert list(mc.quantised(iter([("c", 0, 5, 2.0)]), 0, {"c": 50})) == ["c\t0\t5\t2.00"]


# ---- pipeline on a tiny BAM -----------------------------------------------------------------------------
tools = pytest.mark.skipif(not (shutil.which("bedtools") and shutil.which("samtools")), reason="bedtools/samtools not on PATH")


def tiny_bam(tmp_path, reads):
    head = "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:c\tLN:300\n"
    body = "".join(f"{n}\t0\tc\t{pos}\t60\t{cig}\t*\t0\t0\t*\t*\t{tag}\n" for n, pos, cig, tag in reads)
    sam = tmp_path / "t.sam"
    sam.write_text(head + body)
    bam = tmp_path / "t.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    return str(bam)


def run_cov(*args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_coverage.py"), *args], capture_output=True, text=True)


@tools
def test_spliced_coverage_skips_introns_and_sums_over_bams(tmp_path):
    b1 = tiny_bam(tmp_path, [("a", 81, "20M100N30M", "NH:i:1"), ("b", 81, "20M100N30M", "NH:i:1")])
    (tmp_path / "2").mkdir()
    b2 = tiny_bam(tmp_path / "2", [("c", 91, "10M", "NH:i:1")])
    (tmp_path / "s").write_text("c\t300\n")
    out = tmp_path / "cov.bedgraph"
    r = run_cov("--sizes", str(tmp_path / "s"), "--ladder", "0", "-o", str(out), b1, b2)
    assert r.returncode == 0, r.stderr
    got = [l.split("\t") for l in out.read_text().splitlines()]
    # reads a,b cover 0-based [80,100) and [200,230); read c adds [90,100) from the second BAM; the intron [100,200) is NOT covered
    assert got == [["c", "80", "90", "2.00"], ["c", "90", "100", "3.00"], ["c", "200", "230", "2.00"]]


@tools
def test_unique_only_excludes_multimappers_and_records_filters(tmp_path):
    import json
    bam = tiny_bam(tmp_path, [("a", 81, "20M", "NH:i:1"), ("m", 81, "20M", "NH:i:3")])
    (tmp_path / "s").write_text("c\t300\n")
    for flag, depth in (([], "2.00"), (["--unique-only"], "1.00")):
        out = tmp_path / f"o{depth}.bedgraph"
        js = tmp_path / f"{depth}.json"
        assert run_cov("--sizes", str(tmp_path / "s"), "--ladder", "0", "-o", str(out), "--summary-json", str(js), bam, *flag).returncode == 0
        assert out.read_text().split("\t")[3].strip() == depth
        d = json.load(open(js))
        assert d["filters"]["multimappers"] == ("excluded" if flag else "counted") and d["filters"]["strand"] == "both"


@tools
def test_bigwig_and_bgzip_outputs_read_back_identically(tmp_path):
    bam = tiny_bam(tmp_path, [("a", 81, "20M100N30M", "NH:i:1")])
    (tmp_path / "s").write_text("c\t300\n")
    gz = tmp_path / "c.bedgraph.gz"
    assert run_cov("--sizes", str(tmp_path / "s"), "--ladder", "0", "-o", str(gz), bam).returncode == 0 or not (shutil.which("bgzip") and shutil.which("tabix"))
    if gz.exists():
        assert os.path.exists(str(gz) + ".tbi") and ec.open_text(str(gz)).read().splitlines()[0] == "c\t80\t100\t1.00"
    if shutil.which("bedGraphToBigWig") and ec.bigwig_backend():
        bw = tmp_path / "c.bw"
        assert run_cov("--sizes", str(tmp_path / "s"), "--ladder", "0", "-o", str(bw), bam).returncode == 0
        assert ec.bigwig_intervals(str(bw), "c", 1, 300) == [(81, 100, 1.0), (201, 230, 1.0)]     # 1-based closed on read-back


def test_missing_tools_give_a_clear_error(tmp_path):
    env = dict(os.environ, PATH="/nonexistent")
    (tmp_path / "s").write_text("c\t300\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_coverage.py"), "--sizes", str(tmp_path / "s"), "-o", str(tmp_path / "x.bw"), "x.bam"],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 2 and "bedtools not found" in r.stderr


# ---- browser export ------------------------------------------------------------------------------------
def test_junction_bed12_is_a_tophat_junction_with_the_intron_between_two_anchor_blocks():
    b = ec.Bundle(EXAMPLE)
    rows = ebt.junction_bed12(b, anchor=10)
    assert len(rows) == 10
    f = next(r for r in rows if r.split("\t")[3].startswith("s8/g3/u225")).split("\t")
    # intron 1301..1700 (1-based) = BED [1300, 1700); anchors are the 10 exon bases on each side
    assert f[1:3] == ["1290", "1710"] and f[4] == "225" and f[5] == "+" and f[9:12] == ["2", "10,10", "0,410"]
    starts = [int(x) for x in f[11].split(",")]
    sizes = [int(x) for x in f[10].split(",")]
    exons = ec.exons_from_bed12(int(f[1]), sizes, starts)
    assert ec.introns_from_exons(exons) == [(1301, 1700)]                # the gap between the anchors is exactly the intron
    assert rows == sorted(rows, key=lambda r: (r.split("\t")[0], int(r.split("\t")[1])))


def test_junction_anchors_are_clipped_to_the_sequence_ends(tmp_path):
    rows = ebt.junction_bed12(ec.Bundle(EXAMPLE), anchor=10**6)
    assert all(int(r.split("\t")[1]) >= 0 and int(r.split("\t")[2]) <= 30000 for r in rows if r.startswith("scaffold_1"))


def test_protein_bed12_keeps_the_blocks_and_the_protein_coordinates_in_the_name():
    b = ec.Bundle(EXAMPLE)
    rows = ebt.protein_bed12(b)
    assert len(rows) == 4
    ab = next(r for r in rows if r.startswith("scaffold_1") and "PROT_AB_ortholog" in r).split("\t")
    assert ab[3] == "PROT_AB_ortholog:1-466 id=87.77%" and ab[9] == "6"
    exons = ec.exons_from_bed12(int(ab[1]), [int(x) for x in ab[10].strip(",").split(",")], [int(x) for x in ab[11].strip(",").split(",")])
    gff = next(a for a in ec.build_protein_alignments(ec.iter_gff3(b.output_of_format("protein_alignments"))) if a["protein"] == "PROT_AB_ortholog")
    assert exons == gff["exons"]
    assert [r.split("\t")[:2] for r in rows] == sorted([r.split("\t")[:2] for r in rows], key=lambda x: (x[0], int(x[1])))


def test_export_writes_a_loadable_session_with_existing_relative_paths(tmp_path):
    out = tmp_path / "browser"
    genome = tmp_path / "genome.fa"
    genome.write_text(">scaffold_1\nACGT\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "export_browser_tracks.py"), EXAMPLE, str(out), "--genome", str(genome),
                        "--locus", "scaffold_1:900-5400"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    root = ET.parse(out / "igv_session.xml").getroot()
    assert root.tag == "Session" and root.get("locus") == "scaffold_1:900-5400"
    assert os.path.exists(os.path.normpath(os.path.join(out, root.get("genome"))))
    paths = [res.get("path") for res in root.iter("Resource")]
    assert len(paths) == 7 and all(os.path.exists(os.path.normpath(os.path.join(out, p))) for p in paths)
    assert "junctions.bed" in paths and "protein_alignments.bed" in paths
    assert (out / "junctions.bed").read_text().startswith('track name=junctions description="TopHat junctions"\n')
