"""collapse_long_reads.py and miniprot_to_evidence.py."""
import os
import subprocess
import sys

import pytest

import collapse_long_reads as clr
import evidence_common as ec
import miniprot_to_evidence as mpe
from conftest import DATA, SCRIPTS


# ---- long reads ----------------------------------------------------------------------------
def bed(chrom, exons, strand="+", name="r"):
    cs, ce, sizes, starts = ec.bed12_from_exons(exons)
    return "\t".join(map(str, (chrom, cs, ce, name, 60, strand, cs, ce, "0,0,0", len(sizes),
                               ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ",")))


def collapse(tmp_path, per_sample, groups=None, min_intron=30):
    inputs = {}
    for sid, lines in per_sample.items():
        p = tmp_path / f"{sid}.bed12"
        p.write_text("\n".join(lines) + "\n")
        inputs[sid] = str(p)
    return clr.collapse(inputs, groups, min_intron, None)


def test_same_intron_chain_collapses_across_samples(tmp_path):
    ex = [(1001, 1300), (1701, 1900)]
    rows, st = collapse(tmp_path, {"L1": [bed("c", ex)] * 3, "L2": [bed("c", ex)] * 2}, groups={"L1": "head", "L2": "gut"})
    assert len(rows) == 1
    r = rows[0]
    assert (r["read_count"], r["sample_count"], r["group_count"]) == (5, 2, 2)
    assert r["introns"] == [(1301, 1700)]


def test_different_ends_same_chain_use_median_observed_ends(tmp_path):
    chain = [(1301, 1700)]
    starts, ends = [1001, 1101, 1101, 1101, 1201], [1850, 1900, 1900, 1900, 1950]
    reads = [bed("c", [(s, 1300), (1701, e)]) for s, e in zip(starts, ends)]
    rows, _ = collapse(tmp_path, {"L1": reads})
    r = rows[0]
    assert r["introns"] == chain
    # median of the starts is 1101 (neither the minimum 1001 nor the maximum 1201); median of the ends is 1900
    assert r["exons"] == [(1101, 1300), (1701, 1900)]


def test_different_chain_or_strand_is_a_different_model(tmp_path):
    a = [(1001, 1300), (1701, 1900)]
    b = [(1001, 1300), (1751, 1900)]
    rows, _ = collapse(tmp_path, {"L1": [bed("c", a), bed("c", b), bed("c", a, strand="-")]})
    assert len(rows) == 3


def single(tmp_path, per_sample, tol=50, keep=True, groups=None):
    inputs = {}
    for sid, lines in per_sample.items():
        p = tmp_path / f"{sid}.bed12"
        p.write_text("\n".join(lines) + "\n")
        inputs[sid] = str(p)
    return clr.collapse(inputs, groups, 30, None, keep, tol)


def test_single_exon_reads_are_clustered_counted_and_kept_apart_from_spliced_models(tmp_path):
    reads = {"L1": [bed("c", [(100, 900)]), bed("c", [(110, 905)]), bed("c", [(1001, 1300), (1701, 1900)])],
             "L2": [bed("c", [(130, 890)]), bed("c", [(100, 900)])]}
    rows, st = single(tmp_path, reads, groups={"L1": "head", "L2": "gut"})
    se = [r for r in rows if not r["introns"]]
    assert len(se) == 1 and len(rows) == 2
    r = se[0]
    assert (r["read_count"], r["sample_count"], r["group_count"]) == (4, 2, 2)
    assert r["exons"] == [(100, 900)]                                    # median_low (2nd of 4) of starts 100,100,110,130 and ends 890,900,900,905
    assert (st["single_exon_alignments"], st["single_exon_models"], st["single_exon_reads_in_models"], st["spliced_models"]) == (4, 1, 4, 1)


def test_single_exon_reads_beyond_tolerance_or_on_the_other_strand_are_separate_models(tmp_path):
    reads = {"L1": [bed("c", [(100, 900)]), bed("c", [(300, 900)]), bed("c", [(100, 1200)]), bed("c", [(100, 900)], strand="-")]}
    rows, _ = single(tmp_path, reads)
    assert len(rows) == 4


def test_single_exon_tolerance_is_measured_from_the_first_read_of_a_cluster(tmp_path):
    reads = {"L1": [bed("c", [(100, 900)]), bed("c", [(140, 900)]), bed("c", [(190, 900)])]}
    rows, _ = single(tmp_path, reads, tol=50)
    assert sorted(r["read_count"] for r in rows) == [1, 2]                # 190 is 90 from the anchor at 100: not chained through 140
    rows, _ = single(tmp_path, reads, tol=100)
    assert [r["read_count"] for r in rows] == [3]


def test_single_exon_tolerance_boundary_is_inclusive_for_start_and_end(tmp_path):
    # membership is checked through the cluster coordinates, not just the sizes: (100,150) together and 151 alone, not 100 alone and (150,151) together
    rows, _ = single(tmp_path, {"L1": [bed("c", [(100, 900)]), bed("c", [(150, 900)]), bed("c", [(151, 900)])]}, tol=50)
    assert {r["read_count"]: r["exons"][0] for r in rows} == {2: (100, 900), 1: (151, 900)}          # start 150 is exactly 50 from the anchor: joins; 151 does not
    rows, _ = single(tmp_path, {"L1": [bed("c", [(100, 900)]), bed("c", [(100, 950)]), bed("c", [(100, 951)])]}, tol=50)
    assert {r["read_count"]: r["exons"][0] for r in rows} == {2: (100, 900), 1: (100, 951)}          # end 950 is exactly 50 from the anchor's end: joins; 951 does not


def test_single_exon_clustering_does_not_depend_on_input_order(tmp_path):
    a = [bed("c", [(s, s + 800)]) for s in (100, 120, 500, 510, 900)]
    r1, _ = single(tmp_path, {"L1": a[:3], "L2": a[3:]})
    r2, _ = single(tmp_path, {"L2": a[3:], "L1": a[:3]})
    key = lambda rows: sorted((r["exons"][0], r["read_count"], r["sample_count"]) for r in rows)
    assert key(r1) == key(r2)


def test_no_single_exon_option_drops_the_models_but_keeps_the_count(tmp_path):
    rows, st = single(tmp_path, {"L1": [bed("c", [(100, 900)])] * 3}, keep=False)
    assert rows == [] and st["single_exon_alignments"] == 3 and st["single_exon_models"] == 0


def test_single_exon_bed_line_has_one_block_and_a_dot_chain():
    row = {"seqid": "c", "strand": "+", "exons": [(101, 900)], "introns": [], "read_count": 4, "sample_count": 2, "group_count": 2}
    f = clr.bed_line(row, "LR1").split("\t")
    assert f[1:3] == ["100", "900"] and f[9:12] == ["1", "800,", "0,"] and f[12:] == ["4", "2", "2", "."]
    rec = ec.parse_long_read_line("\t".join(f))
    assert rec["introns"] == [] and rec["exons"] == [(101, 900)]


def test_short_gaps_are_merged_not_called_introns(tmp_path):
    read = bed("c", [(1001, 1300), (1311, 1500), (1901, 2000)])        # 10-base gap = deletion, 400-base gap = intron
    rows, _ = collapse(tmp_path, {"L1": [read]}, min_intron=30)
    assert rows[0]["introns"] == [(1501, 1900)]


def test_bed_output_line_is_valid_bed12_plus_4():
    row = {"seqid": "c", "strand": "+", "exons": [(1001, 1300), (1701, 1900)], "introns": [(1301, 1700)],
           "read_count": 7, "sample_count": 2, "group_count": 2}
    f = clr.bed_line(row, "LR000001").split("\t")
    assert len(f) == 16 and f[1:3] == ["1000", "1900"] and f[9] == "2" and f[10] == "300,200," and f[11] == "0,700,"
    assert f[12:] == ["7", "2", "2", "1301-1700"]
    assert ec.parse_long_read_line("\t".join(f))["introns"] == [(1301, 1700)]


def test_score_is_capped_at_1000():
    row = {"seqid": "c", "strand": "+", "exons": [(1, 10), (50, 60)], "introns": [(11, 49)], "read_count": 5000, "sample_count": 1, "group_count": "."}
    assert clr.bed_line(row, "x").split("\t")[4] == "1000"


# ---- miniprot ---------------------------------------------------------------------------------
def real_lines():
    return open(os.path.join(DATA, "miniprot_DPP3.gff")).read().splitlines(True)


def test_real_miniprot_output_is_converted_block_for_block():
    alns, stats = mpe.convert(real_lines(), None, None, 1, 0, 0)
    assert stats["kept"] == 1 and stats["no_query_length"] == 0
    body = mpe.render(alns, "miniprot", "test db", "10090")
    feats = [ec.parse_gff3_line(l) for l in body]
    parents = [f for f in feats if f.type == "protein_match"]
    parts = [f for f in feats if f.type == "match_part"]
    assert len(parents) == 1 and len(parts) == 17                       # the real alignment has 17 CDS blocks
    orig = [ec.parse_gff3_line(l) for l in real_lines() if not l.startswith("#")]
    assert [(p.start, p.end) for p in parts] == [(f.start, f.end) for f in orig if f.type == "CDS"]
    assert (parents[0].start, parents[0].end) == (1, 27033)             # stop_codon excluded from the span
    a = parents[0].attrs
    assert a["Target"] == "ENSMUSP00000025851.4 1 738" and a["query_length"] == "738"
    assert a["pct_identity"] == "92.82" and a["pct_positive"] == "96.34" and a["pct_coverage"] == "100.00"
    assert a["db"] == "test db" and a["taxid"] == "10090"           # parsed (URL-decoded) value
    assert "db=test%20db" in body[0]                                # stored URL-escaped, as GFF3 requires for reserved characters


def test_protein_side_coordinates_are_contiguous():
    alns, _ = mpe.convert(real_lines(), None, None, 1, 0, 0)
    feats = [ec.parse_gff3_line(l) for l in mpe.render(alns, "miniprot", None, None)]
    al = ec.build_protein_alignments(feats)[0]
    t = [(b["target_start"], b["target_end"]) for b in al["blocks"]]
    assert t[0][0] == 1 and t[-1][1] == 738 and all(t[i][1] + 1 == t[i + 1][0] for i in range(len(t) - 1))
    assert al["introns"][0] == (271, 2972)                              # between CDS 1-270 and 2973-3062


def test_rank_and_identity_filters():
    lines = real_lines()
    second = [l.replace("Rank=1", "Rank=2").replace("MP000001", "MP000002") for l in lines if not l.startswith("##gff")]
    both = lines + second
    _, st = mpe.convert(both, None, None, 1, 0, 0)
    assert st["kept"] == 1 and st["dropped_rank"] == 1
    _, st = mpe.convert(both, None, None, 0, 0, 0)
    assert st["kept"] == 2
    _, st = mpe.convert(lines, None, None, 1, 99.0, 0)
    assert st["dropped_identity"] == 1


def test_missing_paf_lines_give_unknown_coverage_and_a_warning():
    lines = [l for l in real_lines() if not l.startswith("##PAF")]
    alns, st = mpe.convert(lines, None, None, 1, 0, 0)
    assert st["no_query_length"] == 1
    feats = [ec.parse_gff3_line(l) for l in mpe.render(alns, "miniprot", None, None)]
    assert feats[0].attrs["pct_coverage"] == "." and feats[0].attrs["query_length"] == "."


def test_cli_writes_gff3_with_header(tmp_path):
    out = tmp_path / "p.gff3"
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "miniprot_to_evidence.py"), os.path.join(DATA, "miniprot_DPP3.gff"),
                    "-o", str(out), "--aligner-version", "0.18"], check=True, capture_output=True)
    head = out.read_text().splitlines()[:4]
    assert head[0] == "##gff-version 3" and "##ae.profile=protein_alignments" in head
    assert ec.read_pragmas(str(out))["aligner"] == "miniprot 0.18"
