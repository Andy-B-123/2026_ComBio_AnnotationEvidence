"""Coordinate-convention tests.

Conventions under test (docs/SPECIFICATION.md section 2):
  GFF3 and the project's text tables : 1-based, closed
  BED / BED12 / bigWig               : 0-based, half-open
  STAR SJ.out.tab                    : 1-based, first and last INTRON base (same as a GFF3 intron)
"""
import random

import pytest

import evidence_common as ec
from conftest import needs_htslib


# ---- 1-based closed <-> BED 0-based half-open ------------------------------------------
@pytest.mark.parametrize("start,end,bed", [(1, 1, (0, 1)), (1, 10, (0, 10)), (1301, 1700, (1300, 1700)), (5, 5, (4, 5))])
def test_closed_to_bed(start, end, bed):
    assert ec.closed_to_bed(start, end) == bed
    assert ec.bed_to_closed(*bed) == (start, end)


def test_length_is_preserved_by_conversion():
    for s, e in [(1, 1), (10, 20), (1301, 1700)]:
        cs, ce = ec.closed_to_bed(s, e)
        assert ce - cs == e - s + 1


@pytest.mark.parametrize("bad", [(0, 5), (5, 4), (-1, 3)])
def test_closed_to_bed_rejects_invalid(bad):
    with pytest.raises(ValueError):
        ec.closed_to_bed(*bad)


@pytest.mark.parametrize("bad", [(-1, 5), (5, 5), (6, 5)])
def test_bed_to_closed_rejects_invalid(bad):
    with pytest.raises(ValueError):
        ec.bed_to_closed(*bad)


# ---- GFF3 exons -> introns (identical to STAR's intron convention) ---------------------------
def test_introns_from_exons_matches_star_convention():
    # exon 1001-1300, exon 1701-1900: intron is the bases 1301..1700 inclusive.
    assert ec.introns_from_exons([(1001, 1300), (1701, 1900)]) == [(1301, 1700)]
    # STAR SJ.out.tab reports the same intron as start=1301, end=1700 (first and last intron base)


def test_introns_order_independent_and_adjacent_exons():
    assert ec.introns_from_exons([(2501, 2700), (1001, 1300), (1701, 1900)]) == [(1301, 1700), (1901, 2500)]
    assert ec.introns_from_exons([(1, 10), (11, 20)]) == []          # touching exons: no intron
    assert ec.introns_from_exons([(1, 10), (12, 20)]) == [(11, 11)]  # smallest possible intron, one base


# ---- BED12 blocks -----------------------------------------------------------------------
def test_bed12_blocks_to_exons_and_back():
    exons = [(1001, 1300), (1701, 1900), (2501, 2700)]
    cs, ce, sizes, starts = ec.bed12_from_exons(exons)
    assert (cs, ce) == (1000, 2700)
    assert sizes == [300, 200, 200] and starts == [0, 700, 1500]
    assert ec.exons_from_bed12(cs, sizes, starts) == exons
    # BED rule: first blockStart is 0 and the last block ends exactly at chromEnd
    assert starts[0] == 0 and starts[-1] + sizes[-1] == ce - cs


def test_bed12_introns_equal_gff3_introns():
    exons = [(1001, 1300), (1701, 1900)]
    cs, ce, sizes, starts = ec.bed12_from_exons(exons)
    via_bed = ec.introns_from_exons(ec.exons_from_bed12(cs, sizes, starts))
    assert via_bed == ec.introns_from_exons(exons) == [(1301, 1700)]


def test_intron_chain_text_roundtrip():
    chain = [(1301, 1700), (1901, 2500)]
    assert ec.format_chain(chain) == "1301-1700;1901-2500"
    assert ec.parse_chain("1301-1700;1901-2500") == chain
    assert ec.format_chain([]) == "." and ec.parse_chain(".") == []


def test_bed12_parse_line_gives_one_based_coordinates():
    line = "s1\t1000\t2700\tLR1\t5\t+\t1000\t2700\t0,0,0\t3\t300,200,200,\t0,700,1500,\t5\t2\t1\t1301-1700;1901-2500"
    rec = ec.parse_long_read_line(line)
    assert (rec["start"], rec["end"]) == (1001, 2700)                 # BED chromStart+1 .. chromEnd
    assert rec["exons"] == [(1001, 1300), (1701, 1900), (2501, 2700)]
    assert rec["introns"] == [(1301, 1700), (1901, 2500)]
    assert ec.format_chain(rec["introns"]) == "1301-1700;1901-2500"


# ---- regions ------------------------------------------------------------------------------
def test_parse_region_forms():
    assert ec.parse_region("scaffold_1:1001-2000") == ("scaffold_1", 1001, 2000)
    assert ec.parse_region("scaffold_1:1,001-2,000") == ("scaffold_1", 1001, 2000)
    assert ec.parse_region("scaffold_1:1001") == ("scaffold_1", 1001, None)
    assert ec.parse_region("scaffold_1") == ("scaffold_1", 1, None)
    assert ec.parse_region("NC_051855.1:5-9") == ("NC_051855.1", 5, 9)


def test_parse_region_splits_on_last_colon():
    assert ec.parse_region("HLA-A*01:01:5-9") == ("HLA-A*01:01", 5, 9)


@pytest.mark.parametrize("bad", ["chr1:0-10", "chr1:10-5"])
def test_parse_region_rejects_invalid(bad):
    with pytest.raises(ec.EvidenceError):
        ec.parse_region(bad)


def test_overlap_boundaries_closed_intervals():
    assert ec.overlaps(1, 10, 10, 20)          # share base 10
    assert not ec.overlaps(1, 10, 11, 20)      # adjacent, no shared base
    assert ec.overlaps(5, 5, 5, 5)             # single base
    assert not ec.overlaps(5, 5, 6, 6)


# ---- tabix must agree with the pure-Python scan on every boundary ---------------------------
FILES = [("rna_junctions.tsv.gz", "tsv"), ("long_read_models.bed.gz", "bed"), ("predictions.gff3.gz", "gff3"),
         ("repeats.bed.gz", "bed"), ("protein_alignments.gff3.gz", "gff3")]


@needs_htslib
@pytest.mark.parametrize("name,kind", FILES)
def test_tabix_and_scan_agree_at_all_edges(example_bundle, name, kind):
    path = f"{example_bundle}/evidence/{name}"
    edges = set()
    for ln in ec.fetch_lines(path, kind, "scaffold_1", 1, None, use_tabix=False):
        s, e = ec._line_interval(ln, kind)[1:]
        edges.update({s - 1, s, s + 1, e - 1, e, e + 1})
    rng = random.Random(1)
    regions = [(a, b) for a in sorted(edges) for b in sorted(edges) if 1 <= a <= b][:400]
    regions += [(a, a) for a in sorted(edges) if a >= 1]
    regions += [(rng.randint(1, 30000), rng.randint(1, 30000)) for _ in range(50)]
    for a, b in regions:
        a, b = min(a, b), max(a, b)
        scan = ec.fetch_lines(path, kind, "scaffold_1", a, b, use_tabix=False)
        tbx = ec.fetch_lines(path, kind, "scaffold_1", a, b, use_tabix=True)
        assert sorted(scan) == sorted(tbx), (name, a, b)


@needs_htslib
def test_bed_query_boundaries_are_one_based(example_bundle):
    """LR000001 is BED 1000-2700, i.e. 1-based 1001..2700: it overlaps base 1001 but not base 1000, and 2700 but not 2701."""
    path = f"{example_bundle}/evidence/long_read_models.bed.gz"
    ids = lambda a, b, **k: {ln.split("\t")[3] for ln in ec.fetch_lines(path, "bed", "scaffold_1", a, b, **k)}
    for tabix in (True, False):
        assert "LR000001" in ids(1001, 1001, use_tabix=tabix)
        assert "LR000001" not in ids(1000, 1000, use_tabix=tabix)
        assert "LR000001" in ids(2700, 2700, use_tabix=tabix)
        assert "LR000001" not in ids(2701, 2701, use_tabix=tabix)


def test_junction_boundaries_intron_bases(example_bundle):
    """The cross-gene junction covers the intron bases 2701..3400: annotated exon bases 2700 and 3401 are outside it."""
    import evidence_common as m
    b = m.Bundle(example_bundle)
    covering = lambda a, z: {(j["start"], j["end"]) for j in b.junctions("scaffold_1", a, z, use_tabix=False)}
    assert (2701, 3400) in covering(2701, 2701)
    assert (2701, 3400) in covering(3400, 3400)
    assert (2701, 3400) not in covering(2700, 2700)
    assert (2701, 3400) not in covering(3401, 3401)
