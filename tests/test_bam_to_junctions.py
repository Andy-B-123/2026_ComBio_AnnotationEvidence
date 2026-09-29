"""bam_to_junctions.py: STAR-equivalent junction counts from alignments."""
import os
import random
import shutil
import subprocess
import sys

import pytest

import aggregate_star_junctions as agg
import bam_to_junctions as b2j
import evidence_common as ec
from conftest import DATA, SCRIPTS

REGION = os.path.join(DATA, "bam_region")


def run(*args, stdin=None):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, "bam_to_junctions.py"), *args], capture_output=True, text=True, input=stdin)


def load(path):
    out = {}
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        out[(c[0], int(c[1]), int(c[2]))] = tuple(c[3:5]) + tuple(c[6:9])      # strand, motif, unique, multi, overhang
    return out


# ---- real data: identical to STAR's own SJ.out.tab ---------------------------------------------------
def test_real_region_reproduces_every_star_junction_exactly(tmp_path):
    out = tmp_path / "o.tab"
    r = run(os.path.join(REGION, "region.sam"), "--fasta", os.path.join(REGION, "region.fa"), "-o", str(out))
    assert r.returncode == 0, r.stderr
    star, mine = load(os.path.join(REGION, "expected.SJ.out.tab")), load(str(out))
    assert len(star) == 51
    wrong = {k: (v, mine.get(k)) for k, v in star.items() if mine.get(k) != v}
    assert wrong == {}                                     # strand, motif, unique, multimapping AND overhang all identical
    assert len(mine) > len(star)                           # unfiltered: junctions that STAR's --outSJfilter* removed are present
    motifs = {v[1] for v in star.values()}
    assert {"1", "2", "3", "4", "6"} <= motifs and any(v[3] != "0" for v in star.values()) and any(v[2] == "0" for v in star.values())


def test_real_region_mates_that_both_cross_a_junction_count_once():
    lines = [l for l in open(os.path.join(REGION, "region.sam")) if not l.startswith("@")]
    names = [l.split("\t")[0] for l in lines]
    assert len(names) > len(set(names))                    # the fixture contains pairs with two spliced records
    fa = b2j.Fasta(os.path.join(REGION, "region.fa"))
    tables = dict(b2j.read_junctions(open(os.path.join(REGION, "region.sam")), None, 0x204))
    for (s, e), (uniq, multi, ovh) in tables["NC_051855.1"].items():
        assert isinstance(uniq, set) and not (uniq & multi)


def test_output_feeds_the_aggregator_unchanged(tmp_path):
    out = tmp_path / "S1.SJ.out.tab"
    assert run(os.path.join(REGION, "region.sam"), "--fasta", os.path.join(REGION, "region.fa"), "-o", str(out)).returncode == 0
    rows, st = agg.aggregate({"S1": str(out)}, None, 1, None)
    assert st["strand_motif_inconsistent"] == 0 and rows


# ---- hand-made cases --------------------------------------------------------------------------------
def make_fasta(tmp_path, seqs, width=10):
    p = tmp_path / "g.fa"
    with open(p, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name} test\n" + "\n".join(s[i:i + width] for i in range(0, len(s), width)) + "\n")
    return str(p)


def sam(records, so="coordinate", contigs=(("c", 1000),)):
    head = f"@HD\tVN:1.4\tSO:{so}\n" + "".join(f"@SQ\tSN:{n}\tLN:{l}\n" for n, l in contigs)
    body = ""
    for r in records:
        name, flag, chrom, pos, cigar = r[:5]
        tags = r[5] if len(r) > 5 else "NH:i:1"
        body += "\t".join([name, str(flag), chrom, str(pos), "60", cigar, "*", "0", "0", "*", "*", tags]) + "\n"
    return head + body


def genome_with(intron_sites, length=1000):
    """A 'c' sequence with chosen donor/acceptor dinucleotides: {(start, end): (donor, acceptor)} for 1-based intron bases."""
    g = list("A" * length)
    for (s, e), (d, a) in intron_sites.items():
        g[s - 1:s + 1] = list(d)
        g[e - 2:e] = list(a)
    return "".join(g)


def call(tmp_path, records, genome, extra=(), **kw):
    fa = make_fasta(tmp_path, genome if isinstance(genome, dict) else {"c": genome})
    inp = tmp_path / "in.sam"
    inp.write_text(sam(records, **kw))
    out = tmp_path / "out.tab"
    r = run(str(inp), "--fasta", fa, "-o", str(out), *extra)
    assert r.returncode == 0, r.stderr
    return [l.rstrip("\n").split("\t") for l in open(out)], r


def test_mates_crossing_the_same_junction_count_once_and_unique_vs_multimapper_are_distinct_names(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    recs = [("r1", 99, "c", 81, "20M100N30M"), ("r1", 147, "c", 90, "11M100N40M"),      # pair, both mates cross: ONE read
            ("r2", 0, "c", 81, "20M100N30M"),
            ("m1", 0, "c", 81, "20M100N30M", "NH:i:2"), ("m1", 256, "c", 81, "20M100N30M", "NH:i:2"),   # two alignments of one multimapper read: ONE read
            ("m2", 0, "c", 81, "20M100N30M", "NH:i:3")]
    rows, _ = call(tmp_path, recs, g)
    # unique: r1 (both mates) and r2 = 2; multimapping: m1 (two alignments) and m2 = 2; overhang: max(min(20,30), min(11,40)) = 20
    assert rows == [["c", "101", "200", "1", "1", "0", "2", "2", "20"]]


def test_overhang_is_the_longest_min_of_the_two_runs_next_to_the_junction(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    recs = [("a", 0, "c", 91, "10M100N40M"), ("b", 0, "c", 71, "30M100N25M"), ("c", 0, "c", 61, "40M100N5M")]
    rows, _ = call(tmp_path, recs, g)
    assert rows[0][8] == "25"                                           # min(10,40)=10, min(30,25)=25, min(40,5)=5


def test_runs_are_split_at_insertions_and_deletions(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    # 60M... a deletion 4 bases before the junction makes the left run 4, not 64; an insertion after makes the right run 6
    recs = [("a", 0, "c", 41, "56M4D4M100N6M2I34M")]
    rows, _ = call(tmp_path, recs, g)
    assert rows[0][8] == "4"                                            # min(4, 6): M-runs stop at D and I
    recs = [("a", 0, "c", 41, "60M100N40M")]
    rows, _ = call(tmp_path, recs, g)
    assert rows[0][8] == "40"
    # an insertion that decides the minimum, on the right and on the left of the junction
    rows, _ = call(tmp_path, [("a", 0, "c", 41, "60M100N6M2I34M")], g)
    assert rows[0][8] == "6"                                            # right run stops at the insertion: min(60, 6)
    g2 = genome_with({(131, 230): ("GT", "AG")})
    rows, _ = call(tmp_path, [("a", 0, "c", 71, "30M2I30M100N40M")], g2)
    assert rows[0][1:3] == ["131", "230"] and rows[0][8] == "30"        # left run: the 30M next to the junction, not 30+2+30


def test_soft_clips_and_deletions_shift_positions_but_do_not_count_as_anchor(tmp_path):
    g = genome_with({(153, 252): ("GT", "AG")})
    # 5S ignored; 30M covers 101..130; 2D skips 131-132; 20M covers 133..152; the intron is 153..252
    rows, _ = call(tmp_path, [("a", 0, "c", 101, "5S30M2D20M100N25M3S")], g)
    assert rows == [["c", "153", "252", "1", "1", "0", "1", "0", "20"]]           # left run is 20 (stops at the deletion), right run 25


@pytest.mark.parametrize("donor,acceptor,strand,motif", [("GT", "AG", "1", "1"), ("CT", "AC", "2", "2"), ("GC", "AG", "1", "3"),
                                                          ("CT", "GC", "2", "4"), ("AT", "AC", "1", "5"), ("GT", "AT", "2", "6"),
                                                          ("GA", "AG", "0", "0"), ("AG", "GT", "0", "0")])
def test_motif_and_strand_codes_follow_star(tmp_path, donor, acceptor, strand, motif):
    g = genome_with({(101, 200): (donor, acceptor)})
    rows, _ = call(tmp_path, [("a", 0, "c", 81, "20M100N30M")], g)
    assert rows[0][3:6] == [strand, motif, "0"]


def test_motif_is_read_from_the_forward_strand_of_the_fasta_case_insensitively(tmp_path):
    g = genome_with({(101, 200): ("gt", "ag")}).replace("A", "a")
    rows, _ = call(tmp_path, [("a", 0, "c", 81, "20M100N30M")], g)
    assert rows[0][3:5] == ["1", "1"]


def test_unmapped_and_qcfail_records_are_skipped_but_secondary_and_duplicates_count(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    recs = [("u", 4, "c", 81, "20M100N30M"), ("q", 512, "c", 81, "20M100N30M"), ("s", 256, "c", 81, "20M100N30M", "NH:i:2"), ("d", 1024, "c", 81, "20M100N30M")]
    rows, _ = call(tmp_path, recs, g)
    assert rows[0][6:8] == ["1", "1"]                                   # d unique; s multimapping
    rows, _ = call(tmp_path, recs, g, extra=("--exclude-flags", "0x604"))     # also skip duplicates
    assert rows[0][6:8] == ["0", "1"]


def test_thresholds(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG"), (301, 400): ("GT", "AG")})
    recs = [("a", 0, "c", 81, "20M100N30M"), ("b", 0, "c", 281, "5M100N30M", "NH:i:2")]
    assert len(call(tmp_path, recs, g)[0]) == 2
    assert len(call(tmp_path, recs, g, extra=("--min-unique", "1"))[0]) == 1
    assert len(call(tmp_path, recs, g, extra=("--min-overhang", "10"))[0]) == 1
    assert len(call(tmp_path, recs, g, extra=("--min-total", "2"))[0]) == 0


def test_strip_prefix_and_missing_sequences(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    recs = [("a", 0, "lcl|c", 81, "20M100N30M"), ("b", 0, "other", 81, "20M100N30M")]
    rows, r = call(tmp_path, recs, g, extra=("--strip-prefix", "lcl|"), contigs=(("lcl|c", 1000), ("other", 1000)))
    assert [x[0] for x in rows] == ["c"] and "not found in the FASTA" in r.stderr


def test_unsorted_input_is_fully_processed_and_a_false_sorted_claim_is_an_error(tmp_path):
    g = {"c": genome_with({(101, 200): ("GT", "AG")}), "d": genome_with({(101, 200): ("CT", "AC")})}
    recs = [("a", 0, "c", 81, "20M100N30M"), ("b", 0, "d", 81, "20M100N30M"), ("c", 0, "c", 81, "20M100N30M")]
    rows, _ = call(tmp_path, recs, g, so="unsorted", contigs=(("c", 1000), ("d", 1000)))
    assert {(r[0], r[6]) for r in rows} == {("c", "2"), ("d", "1")}
    fa = make_fasta(tmp_path, g)
    inp = tmp_path / "bad.sam"
    inp.write_text(sam(recs, so="coordinate", contigs=(("c", 1000), ("d", 1000))))
    r = run(str(inp), "--fasta", fa, "-o", str(tmp_path / "x.tab"))
    assert r.returncode == 2 and "not coordinate-sorted" in r.stderr


def test_sam_on_stdin(tmp_path):
    g = genome_with({(101, 200): ("GT", "AG")})
    fa = make_fasta(tmp_path, {"c": g})
    r = run("-", "--fasta", fa, "-o", "-", stdin=sam([("a", 0, "c", 81, "20M100N30M")]))
    assert r.returncode == 0 and r.stdout.strip().split("\t")[:5] == ["c", "101", "200", "1", "1"]


@pytest.mark.skipif(not shutil.which("samtools"), reason="samtools not on PATH")
def test_bam_input_gives_the_same_result_as_sam(tmp_path):
    fa = os.path.join(REGION, "region.fa")
    bam = tmp_path / "r.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), os.path.join(REGION, "region.sam")], check=True)
    o1, o2 = tmp_path / "1.tab", tmp_path / "2.tab"
    assert run(str(bam), "--fasta", fa, "-o", str(o1)).returncode == 0
    assert run(os.path.join(REGION, "region.sam"), "--fasta", fa, "-o", str(o2)).returncode == 0
    assert open(o1).read() == open(o2).read()


def test_summary_json_records_filters(tmp_path):
    import json
    g = genome_with({(101, 200): ("GT", "AG")})
    fa = make_fasta(tmp_path, {"c": g})
    inp = tmp_path / "i.sam"
    inp.write_text(sam([("a", 0, "c", 81, "20M100N30M")]))
    js = tmp_path / "s.json"
    assert run(str(inp), "--fasta", fa, "-o", str(tmp_path / "o.tab"), "--min-unique", "1", "--summary-json", str(js)).returncode == 0
    d = json.load(open(js))
    assert d["filters"]["min_unique"] == 1 and d["filters"]["star_outSJfilter"] == "not applied" and d["stats"]["junctions_written"] == 1


# ---- FASTA access ------------------------------------------------------------------------------------
@pytest.mark.parametrize("width", [7, 10, 60])
@pytest.mark.parametrize("indexed", [False, True])
def test_fasta_fetch_matches_python_slicing_across_line_breaks(tmp_path, width, indexed):
    rng = random.Random(3)
    seqs = {"one": "".join(rng.choice("ACGT") for _ in range(137)), "two": "".join(rng.choice("ACGT") for _ in range(41))}
    fa = make_fasta(tmp_path, seqs, width=width)
    if indexed and shutil.which("samtools"):
        subprocess.run(["samtools", "faidx", fa], check=True)
    f = b2j.Fasta(fa)
    for name, s in seqs.items():
        for a, b in [(1, 1), (1, len(s)), (len(s), len(s)), (2, 3), (width, width + 1), (width + 1, width * 2 + 2)] + [tuple(sorted((rng.randint(1, len(s)), rng.randint(1, len(s))))) for _ in range(40)]:
            if b <= len(s):
                assert f.fetch(name, a, b) == s[a - 1:b], (name, a, b, width)
    with pytest.raises(ec.EvidenceError):
        f.fetch("two", 40, 42)
