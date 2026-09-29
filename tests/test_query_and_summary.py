"""query_region.py and summarise_gene.py on the synthetic example (whose four genes are built to show specific situations)."""
import json
import os
import subprocess
import sys

import pytest

import evidence_common as ec
import query_region as qr
import summarise_gene as sg
from conftest import EXAMPLE, SCRIPTS


def cli(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], capture_output=True, text=True)


@pytest.fixture(scope="module")
def bundle():
    return ec.Bundle(EXAMPLE)


def report(bundle, gene, **kw):
    return sg.analyse(bundle, gene, kw.get("tol", 5), kw.get("end_tol", 50), kw.get("ratio", 1.0), kw.get("min_comp", 2))


# ---- query_region ------------------------------------------------------------------------------
def test_query_region_text_mentions_every_class():
    r = cli("query_region.py", EXAMPLE, "scaffold_1:2500-3500")
    assert r.returncode == 0, r.stderr
    out = r.stdout
    assert "1-BASED, CLOSED" in out
    for token in ("geneA", "geneB", "scaffold_1:2701-3400", "LR000002", "PA000002", "AUGUSTUS augustus_g1", "== coverage"):
        assert token in out, token


def test_query_region_json_is_parseable_and_one_based():
    r = cli("query_region.py", EXAMPLE, "scaffold_1:2701-2701", "--format", "json", "--classes", "rna_short,rna_long")
    d = json.loads(r.stdout)
    assert d["region"] == {"seqid": "scaffold_1", "start": 2701, "end": 2701, "coordinates": "1-based-closed"}
    assert [(j["start"], j["end"]) for j in d["classes"]["rna_short"]["junctions"]] == [(2701, 3400)]
    assert [m["id"] for m in d["classes"]["rna_long"]["models"]] == ["LR000002"]     # BED model ending 5200, not the ones ending 2700


def test_query_region_gene_and_flank():
    d = json.loads(cli("query_region.py", EXAMPLE, "--gene", "geneA", "--flank", "10", "--format", "json", "--classes", "annotation").stdout)
    assert (d["region"]["start"], d["region"]["end"]) == (991, 2710)
    assert [g["id"] for g in d["classes"]["annotation"]["genes"]] == ["geneA"]


def test_query_region_filters_are_shown_next_to_the_records():
    out = cli("query_region.py", EXAMPLE, "scaffold_1:1-30000", "--classes", "rna_short").stdout
    assert "min_total_unique_reads=1" in out


def test_query_region_no_tabix_gives_identical_json():
    a = cli("query_region.py", EXAMPLE, "scaffold_1:1000-5300", "--format", "json", "--coverage-bins", "0").stdout
    b = cli("query_region.py", EXAMPLE, "scaffold_1:1000-5300", "--format", "json", "--coverage-bins", "0", "--no-tabix").stdout
    assert json.loads(a) == json.loads(b)


def test_query_region_errors():
    assert cli("query_region.py", EXAMPLE, "nosuchscaffold:1-10").returncode == 2
    assert cli("query_region.py", EXAMPLE, "scaffold_1:10-5").returncode == 2
    assert cli("query_region.py", EXAMPLE, "--gene", "nosuchgene").returncode == 2
    assert cli("query_region.py", EXAMPLE, "scaffold_1:1-10", "--classes", "bogus").returncode == 2


def test_gene_lookup_by_name(bundle):
    assert ec.find_gene(bundle.annotation_path, "geneB").id == "geneB"
    with pytest.raises(ec.EvidenceError):
        ec.find_gene(bundle.annotation_path, "nope")


# ---- summarise_gene: the four teaching cases -----------------------------------------------------
def test_geneA_supported_and_bridged_to_geneB(bundle):
    r = report(bundle, "geneA")
    t = r["transcripts"][0]
    assert [i["status"] for i in t["introns"]] == ["supported", "supported"]
    assert set(t["introns"][0]["supported_by"]) == {"rna_short", "rna_long", "proteins"}
    assert t["introns"][0]["rna_short"]["sample_count"] == 8 and t["introns"][0]["rna_short"]["total_unique_reads"] == 225
    # weak alternative acceptor (10 bp upstream) is listed but not dominant
    comp = t["introns"][1]["competing_junctions"][0]
    assert comp["junction"] == [1901, 2490] and comp["type"] == "alt_acceptor" and comp["shift_bp"] == -10 and not comp["dominant"]
    assert not t["introns"][1]["competing_evidence"]
    x = r["cross_locus"]
    assert [n["id"] for n in x["neighbouring_genes"]] == ["geneB"]
    assert x["spanning"]["rna_long"][0]["id"] == "LR000002" and x["spanning"]["rna_long"][0]["novel_introns"] == [[2701, 3400]]
    assert x["spanning"]["proteins"][0]["genes_with_protein_aa_range"] == {"geneA": [1, 200], "geneB": [201, 466]}
    j = x["spanning"]["rna_short_junctions"][0]
    assert (j["junction"], j["from_gene"], j["to_gene"], j["sample_count"], j["total_unique_reads"]) == ([2701, 3400], "geneA", "geneB", 4, 10)
    assert x["spanning"]["predictions"] == ["AUGUSTUS"] and x["single_gene"]["predictions"]["sources_touching_only_one_gene"] == ["Helixer"]
    assert set(x["bridged_genes"]["geneB"]) == {"rna_short", "rna_long", "proteins", "predictions:AUGUSTUS"}
    # balance: most long-read evidence touches only geneA
    assert x["single_gene"]["rna_long"]["reads"] == 75 and x["single_gene"]["rna_long"]["reads_spanning"] == 9


def test_geneB_exon_skipping_competitor_is_not_dominant(bundle):
    r = report(bundle, "geneB")
    intr = r["transcripts"][0]["introns"]
    c = intr[0]["competing_junctions"][0]
    assert c["junction"] == [3601, 4800] and c["type"] == "exon_skipping"
    assert (c["sample_count"], c["total_unique_reads"], c["max_unique_reads"]) == (1, 100, 100)
    assert c["dominant"] is False and intr[0]["competing_evidence"] is False   # one library only, fewer reads than the annotated intron
    assert r["gene"]["annotated_partial"] is True


def test_geneC_minus_strand_and_absence_of_rna_evidence(bundle):
    r = report(bundle, "geneC")
    assert r["gene"]["strand"] == "-"
    i1, i2 = r["transcripts"][0]["introns"]
    assert i1["intron"] == [1401, 2000] and i2["intron"] == [2301, 3200]
    assert i1["rna_short"]["match"] == "exact" and i1["status"] == "supported"
    # intron 2: no RNA-seq junction retained, but long reads and protein match exactly -> still supported, absence recorded, no verdict
    assert i2["rna_short"]["match"] == "not_observed" and "rna_short" in i2["not_observed_in"]
    assert set(i2["supported_by"]) == {"rna_long", "proteins"} and i2["status"] == "supported"
    assert i2["competing_evidence"] is False
    ends = r["transcripts"][0]["terminal_boundaries"]
    assert ends["start"]["transcript_end"] == "3-prime" and ends["end"]["transcript_end"] == "5-prime"   # minus strand
    assert r["transcripts"][0]["full_chain"]["rna_long_exact_chain"] == {"models": 1, "reads": 31}


def test_geneD_prediction_only_with_dominant_competing_acceptor(bundle):
    r = report(bundle, "geneD")
    i = r["transcripts"][0]["introns"][0]
    assert i["status"] == "prediction_only" and i["supported_by"] == []
    assert i["predictions"]["sources"] == ["AUGUSTUS", "Helixer"]
    assert i["rna_short"]["match"] == "not_observed"
    c = i["competing_junctions"][0]
    assert c["junction"] == [7301, 7790] and c["type"] == "alt_acceptor" and c["shift_bp"] == -10
    assert (c["sample_count"], c["total_unique_reads"]) == (6, 62) and c["dominant"] and i["competing_evidence"]
    assert c["fraction_of_annotated"] is None                     # annotated junction has zero reads: no ratio to report


def test_competing_thresholds_are_parameters(bundle):
    r = report(bundle, "geneD", min_comp=7)                       # only 6 samples support the competitor
    assert r["transcripts"][0]["introns"][0]["competing_evidence"] is False
    r = report(bundle, "geneB", ratio=0.3)                        # 100 reads >= 0.3 * 301 but only 1 sample
    assert r["transcripts"][0]["introns"][0]["competing_evidence"] is False
    r = report(bundle, "geneB", ratio=0.3, min_comp=1)
    assert r["transcripts"][0]["introns"][0]["competing_evidence"] is True


def test_near_match_helper():
    assert sg.near((100, 200), [(100, 200), (103, 198), (110, 200)], 5) == [(103, 198)]


def test_coverage_reader_reports_depth_or_says_it_is_unavailable(bundle):
    r = report(bundle, "geneA")
    assert r["coverage_reader"] == "available" or r["coverage_reader"].startswith("unavailable")
    if r["coverage_reader"] == "available":
        cv = r["transcripts"][0]["introns"][0]["coverage"]
        assert cv["intron_mean_depth"] < 20 < min(cv["flanking_exon_mean_depth"])


# ---- rendering ---------------------------------------------------------------------------------
def test_markdown_keeps_support_absence_and_competition_separate():
    out = cli("summarise_gene.py", EXAMPLE, "geneD").stdout
    assert "prediction_only" in out and "COMPETING EVIDENCE" in out
    assert "not observed in the retained junction table" in out
    assert "absence is not proof of error" in out
    out_c = cli("summarise_gene.py", EXAMPLE, "geneC").stdout
    assert "rna_short: not observed" in out_c and "supported" in out_c


def test_markdown_never_issues_a_verdict_about_the_model():
    for gene in ("geneA", "geneB", "geneC", "geneD"):
        out = cli("summarise_gene.py", EXAMPLE, gene).stdout.lower()
        for phrase in ("is wrong", "is incorrect", "should be merged", "must be merged", "is a fusion", "is an error"):
            assert phrase not in out, (gene, phrase)
    assert "not a finding that the models are wrong" in cli("summarise_gene.py", EXAMPLE, "geneA").stdout


def test_json_output_roundtrips():
    d = json.loads(cli("summarise_gene.py", EXAMPLE, "geneA", "--format", "json").stdout)
    assert d["gene"]["id"] == "geneA" and d["coordinates"].startswith("1-based-closed")


def test_unknown_gene_exits_with_error():
    r = cli("summarise_gene.py", EXAMPLE, "nope")
    assert r.returncode == 2 and "not found" in r.stderr


def test_absent_classes_are_stated_not_guessed(tmp_bundle):
    import yaml
    mp = os.path.join(tmp_bundle, "manifest.yaml")
    m = yaml.safe_load(open(mp))
    m["evidence"]["rna_short"] = {"present": False, "reason": "no RNA-seq available"}
    m["evidence"]["coverage"] = {"present": False, "reason": "not retained"}
    yaml.safe_dump(m, open(mp, "w"), sort_keys=False)
    out = cli("summarise_gene.py", tmp_bundle, "geneA").stdout
    assert "rna_short: NOT PRESENT in this bundle" in out
    assert "rna_short: exact" not in out and "rna_short: not observed" not in out       # no claim either way


# ---- competing-junction logic, both strands (the example has no minus-strand competitor) ------------
def J(start, end, strand, total=10, samples=3):
    return dict(start=start, end=end, strand=strand, sample_count=samples, group_count=2, total_unique_reads=total, max_unique_reads=total)


def test_competitor_types_and_shifts_depend_on_strand():
    intron = (2301, 3200)
    junctions = lambda s: [J(2301, 3190, s), J(2311, 3200, s)]
    plus = {tuple(c["junction"]): c for c in sg.competitors(intron, "+", [intron], junctions("+"), 0, 1.0, 2)}
    # + strand: donor is the LOWER coordinate. Sharing the start = sharing the donor = alternative acceptor.
    assert plus[(2301, 3190)]["type"] == "alt_acceptor" and plus[(2301, 3190)]["shift_bp"] == -10
    assert plus[(2311, 3200)]["type"] == "alt_donor" and plus[(2311, 3200)]["shift_bp"] == 10
    minus = {tuple(c["junction"]): c for c in sg.competitors(intron, "-", [intron], junctions("-"), 0, 1.0, 2)}
    # - strand: donor is the HIGHER coordinate. Sharing the end = sharing the donor = alternative acceptor.
    assert minus[(2301, 3190)]["type"] == "alt_donor" and minus[(2311, 3200)]["type"] == "alt_acceptor"


def test_competitors_ignore_other_strand_and_unrelated_junctions():
    intron = (100, 200)
    out = sg.competitors(intron, "+", [intron], [J(100, 190, "-"), J(500, 600, "+"), J(100, 200, "+"), J(101, 190, "+")], 0, 1.0, 2)
    assert out == []                         # wrong strand / no shared site / the intron itself / shares neither site


def test_exon_skipping_is_recognised_from_the_transcripts_other_introns():
    introns = [(100, 200), (300, 400)]
    out = sg.competitors((100, 200), "+", introns, [J(100, 400, "+")], 50, 1.0, 2)
    assert out[0]["type"] == "exon_skipping"
    out = sg.competitors((300, 400), "+", introns, [J(100, 400, "+")], 50, 1.0, 2)
    assert out[0]["type"] == "exon_skipping"


def test_dominance_uses_reads_and_samples():
    intron = (100, 200)
    c = lambda total, samples, annotated: sg.competitors(intron, "+", [intron], [J(100, 190, "+", total, samples)], annotated, 1.0, 2)[0]["dominant"]
    assert c(50, 3, 50) is True and c(49, 3, 50) is False and c(500, 1, 50) is False and c(5, 2, 0) is True


def test_gene_models_can_be_built_from_a_generator():
    genes = ec.build_gene_models(ec.iter_gff3(os.path.join(EXAMPLE, "annotation.gff3.gz")))
    assert sorted(genes) == ["geneA", "geneB", "geneC", "geneD", "geneE"] and genes["geneA"].transcripts[0].introns == [(1301, 1700), (1901, 2500)]


# ---- scoring a competing model that is not in the bundle's annotation ---------------------------------
COMPETING = os.path.join(os.path.dirname(EXAMPLE), "competing_model_AB.gff3")


def test_competing_merged_model_is_scored_against_the_same_evidence(bundle):
    r = sg.analyse(bundle, "geneAB", 5, 50, 1.0, 2, annotation=COMPETING)
    assert r["annotation_used"] == "competing_model_AB.gff3" and r["gene"]["id"] == "geneAB"
    t = r["transcripts"][0]
    assert [i["intron"] for i in t["introns"]] == [[1301, 1700], [1901, 2500], [2701, 3400], [3601, 4000], [4301, 4800]]
    assert all(i["status"] == "supported" for i in t["introns"])
    connecting = t["introns"][2]
    assert connecting["rna_short"]["sample_count"] == 4 and connecting["rna_short"]["total_unique_reads"] == 10
    assert set(connecting["supported_by"]) == {"rna_short", "rna_long", "proteins"}
    # whole-chain support is much thinner than for the two separate genes: this is what a curator has to weigh
    assert t["full_chain"]["rna_long_exact_chain"] == {"models": 1, "reads": 9}
    assert t["full_chain"]["proteins_containing_all_introns"] == ["PA000002"] and t["full_chain"]["predictions_exact_chain"] == ["AUGUSTUS"]
    assert r["cross_locus"]["neighbouring_genes"] == []           # neighbours come from the override file, not the bundle annotation
    sep = report(bundle, "geneA")["transcripts"][0]["full_chain"]["rna_long_exact_chain"]
    assert sep["reads"] == 75 > t["full_chain"]["rna_long_exact_chain"]["reads"]


def test_annotation_override_via_cli():
    out = cli("summarise_gene.py", EXAMPLE, "geneAB", "--annotation", COMPETING).stdout
    assert "gene model taken from competing_model_AB.gff3" in out and "Intron 3: 2701-3400" in out
    assert cli("summarise_gene.py", EXAMPLE, "geneA", "--annotation", COMPETING).returncode == 2       # gene not in that file


# ---- single-exon long-read models -----------------------------------------------------------------------
def test_single_exon_gene_is_supported_by_single_exon_models_not_by_introns(bundle):
    r = report(bundle, "geneE")
    t = r["transcripts"][0]
    assert t["introns"] == [] and r["intron_status_counts"] == {}
    se = t["single_exon_transcript"]
    m = se["single_exon_models"][0]
    assert (m["id"], m["reads"], m["samples"], m["groups"]) == ("LR000005", 20, 3, 3) and m["reciprocal_overlap"] >= 0.95
    assert se["spliced_models_overlapping"]["models"] == 0
    assert r["single_exon_long_reads"]["models"] == 1 and r["single_exon_long_reads"]["reads"] == 20
    assert t["terminal_boundaries"]["start"]["rna_long"]["models_within_tolerance"] == 1      # single-exon models can test a single-exon transcript's ends
    assert r["cross_locus"]["single_gene"]["rna_long"]["reads"] == 0                          # ... but never count as spliced-model evidence
    assert cli("summarise_gene.py", EXAMPLE, "geneE").returncode == 0


def test_spliced_genes_do_not_count_single_exon_models_as_support(bundle):
    r = report(bundle, "geneA")
    assert r["single_exon_long_reads"]["models"] == 0
    assert "single_exon_transcript" not in r["transcripts"][0]
    assert r["cross_locus"]["single_gene"]["rna_long"]["reads"] == 75                          # unchanged by the new models


def test_single_exon_model_over_a_spliced_locus_is_counted_separately(bundle, tmp_path):
    """A single-exon model overlapping a spliced gene is reported as such, not as intron or chain support."""
    import shutil, gzip
    dst = tmp_path / "b"
    shutil.copytree(EXAMPLE, dst)
    p = dst / "evidence" / "long_read_models.bed.gz"
    lines = gzip.open(p, "rt").read().rstrip("\n").split("\n")
    lines.append("scaffold_1\t1100\t1300\tLR999999\t4\t+\t1100\t1300\t0,0,0\t1\t200,\t0,\t4\t2\t2\t.")
    lines.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    (dst / "evidence" / "long_read_models.bed.gz.tbi").unlink()
    b2 = ec.Bundle(str(dst))
    r = sg.analyse(b2, "geneA", 5, 50, 1.0, 2)
    assert r["single_exon_long_reads"] == {**r["single_exon_long_reads"], "models": 1, "reads": 4}
    assert r["transcripts"][0]["introns"][0]["rna_long"]["reads"] == 84                       # intron support unchanged: 75 + 9


def test_reciprocal_overlap_uses_the_longer_of_the_two_lengths(tmp_path):
    """A single-exon model covering half of a single-exon transcript has reciprocal overlap 0.5, not 1.0."""
    import shutil, gzip
    dst = tmp_path / "b"
    shutil.copytree(EXAMPLE, dst)
    p = dst / "evidence" / "long_read_models.bed.gz"
    lines = gzip.open(p, "rt").read().rstrip("\n").split("\n")
    lines = [l for l in lines if "\tLR000005\t" not in l]                                     # drop the full-length geneE model
    lines.append("scaffold_1\t8000\t8450\tLR777777\t6\t+\t8000\t8450\t0,0,0\t1\t450,\t0,\t6\t2\t2\t.")   # 8001-8450: half of geneE (8001-8900)
    lines.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[1])))
    with gzip.open(p, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    (dst / "evidence" / "long_read_models.bed.gz.tbi").unlink()
    r = sg.analyse(ec.Bundle(str(dst)), "geneE", 5, 50, 1.0, 2)
    m = r["transcripts"][0]["single_exon_transcript"]["single_exon_models"][0]
    assert m["id"] == "LR777777" and m["reciprocal_overlap"] == 0.5
