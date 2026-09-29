"""validate_bundle.py, make_manifest.py and the committed example bundle."""
import gzip
import os
import shutil
import subprocess
import sys

import pytest
import yaml

import evidence_common as ec
import make_manifest as mm
import validate_bundle as vb
from conftest import EXAMPLE, ROOT, needs_htslib

pytest.importorskip("jsonschema")


def validate(bundle, capsys):
    rc = vb.main([bundle, "-q"])
    return rc, capsys.readouterr().out


def edit_file(bundle, rel, fn):
    """Rewrite a gzip'd evidence file and drop size/sha256 from the manifest so the edit is not read as tampering."""
    p = os.path.join(bundle, rel)
    opener = gzip.open if rel.endswith(".gz") else open
    with opener(p, "rt") as fh:
        text = fh.read()
    with opener(p, "wt") as fh:
        fh.write(fn(text))
    for ext in (".tbi", ".csi"):
        if os.path.exists(p + ext):
            os.remove(p + ext)
    mp = os.path.join(bundle, "manifest.yaml")
    m = yaml.safe_load(open(mp))
    for ev in m["evidence"].values():
        for o in ev.get("outputs") or []:
            o.pop("sha256", None)
            o.pop("bytes", None)
    m["annotation"].pop("checksum", None)
    m["annotation"].pop("bytes", None)
    yaml.safe_dump(m, open(mp, "w"), sort_keys=False)


def edit_manifest(bundle, fn):
    mp = os.path.join(bundle, "manifest.yaml")
    m = yaml.safe_load(open(mp))
    fn(m)
    yaml.safe_dump(m, open(mp, "w"), sort_keys=False)


def replace_line(text, contains, new):
    lines = text.split("\n")
    i = next(i for i, l in enumerate(lines) if contains in l and not l.startswith("#"))
    lines[i] = new(lines[i]) if callable(new) else new
    return "\n".join(lines)


# ---- the committed example is valid ------------------------------------------------------------
def test_example_bundle_validates(capsys):
    rc, out = validate(EXAMPLE, capsys)
    assert rc == 0, out
    assert "0 error(s), 0 warning(s)" in out


# ---- negative tests: each corruption is caught with a useful message ------------------------------
def test_tampering_is_detected(tmp_bundle, capsys):
    p = os.path.join(tmp_bundle, "evidence", "rna_junctions.tsv.gz")
    with open(p, "ab") as fh:
        fh.write(b"")
    with gzip.open(p, "wt") as fh:
        fh.write("changed\n")
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "sha256 differs" in out


def test_missing_file_is_detected(tmp_bundle, capsys):
    os.remove(os.path.join(tmp_bundle, "evidence", "repeats.bed.gz"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "listed file is missing: evidence/repeats.bed.gz" in out


def test_junction_max_above_total(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/rna_junctions.tsv.gz",
              lambda t: replace_line(t, "3601\t4800", lambda l: "\t".join(l.split("\t")[:8] + ["999"] + l.split("\t")[9:])))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "max_unique_reads > total_unique_reads" in out


def test_unsorted_junction_file(tmp_bundle, capsys):
    def reverse(t):
        head = [l for l in t.split("\n") if l.startswith("#")]
        body = [l for l in t.split("\n") if l and not l.startswith("#")]
        return "\n".join(head + body[::-1]) + "\n"
    edit_file(tmp_bundle, "evidence/rna_junctions.tsv.gz", reverse)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "not sorted" in out


def test_alias_seqid_is_reported_with_the_primary_name(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/rna_junctions.tsv.gz", lambda t: t.replace("scaffold_1\t1301", "lcl|scaffold_1\t1301"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "alias of 'scaffold_1'" in out


def test_coordinate_beyond_sequence_length(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/rna_junctions.tsv.gz", lambda t: t.replace("1301\t1700", "1301\t99999"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "(length 30000) is invalid" in out


def test_intron_chain_must_match_bed12_blocks(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/long_read_models.bed.gz", lambda t: t.replace("1301-1700;1901-2500\n", "1301-1701;1901-2500\n", 1))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "does not match blocks" in out


def test_bed12_blocks_must_end_at_chromend(tmp_bundle, capsys):
    def shrink(t):
        return replace_line(t, "LR000001", lambda l: "\t".join(l.split("\t")[:2] + ["2699"] + l.split("\t")[3:]))
    edit_file(tmp_bundle, "evidence/long_read_models.bed.gz", shrink)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "must start at chromStart" in out


def test_protein_percent_out_of_range(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/protein_alignments.gff3.gz", lambda t: t.replace("pct_identity=91.00", "pct_identity=191.00", 1))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "outside 0-100" in out


def test_match_part_without_parent(tmp_bundle, capsys):
    edit_file(tmp_bundle, "evidence/protein_alignments.gff3.gz", lambda t: t.replace("Parent=PA000004", "Parent=PA999999"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "PA999999" in out


def test_duplicate_prediction_ids(tmp_bundle, capsys):
    # a gene given the ID of an mRNA: same ID on two different feature types
    edit_file(tmp_bundle, "evidence/predictions.gff3.gz", lambda t: t.replace("ID=helixer_g2;", "ID=helixer_g1.t1;", 1))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "duplicate ID helixer_g1.t1" in out


def test_absent_class_needs_a_reason(tmp_bundle, capsys):
    def absent(m):
        m["evidence"]["repeats"] = {"present": False}
    edit_manifest(tmp_bundle, absent)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "reason" in out


def test_present_class_needs_programs_and_summarisation(tmp_bundle, capsys):
    def strip(m):
        del m["evidence"]["proteins"]["programs"]
        del m["evidence"]["proteins"]["summarisation"]
    edit_manifest(tmp_bundle, strip)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "programs" in out and "summarisation" in out


def test_proteins_present_needs_databases(tmp_bundle, capsys):
    edit_manifest(tmp_bundle, lambda m: m["evidence"]["proteins"].pop("databases"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "databases" in out


def test_unlisted_prediction_source_tag_warns(tmp_bundle, capsys):
    def drop_tag(m):
        for p in m["evidence"]["predictions"]["programs"]:
            p.pop("source_tag", None)
        m["evidence"]["predictions"]["programs"][0]["source_tag"] = "AUGUSTUS"
    edit_manifest(tmp_bundle, drop_tag)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 0 and "Helixer" in out and "WARN" in out


def test_n_datasets_mismatch_warns(tmp_bundle, capsys):
    edit_manifest(tmp_bundle, lambda m: m["evidence"]["rna_short"].update(n_datasets=99))
    rc, out = validate(tmp_bundle, capsys)
    assert "n_datasets=99" in out


def test_missing_index_warns(tmp_bundle, capsys):
    os.remove(os.path.join(tmp_bundle, "evidence", "rna_junctions.tsv.gz.tbi"))
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 0 and "no .tbi/.csi index" in out


def test_not_a_bundle(tmp_path):
    with pytest.raises(ec.EvidenceError):
        ec.Bundle(str(tmp_path))


# ---- make_manifest ---------------------------------------------------------------------------
def test_refget_digest_matches_the_ga4gh_specification_example(tmp_path):
    """The refget spec lists SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2 for the sequence ACGT."""
    fa = tmp_path / "x.fa"
    fa.write_text(">one description\nacg\nT\n>two\nTTTT\n")
    recs = list(mm.fasta_records(str(fa)))
    assert recs[0][:2] == ("one", 4) and recs[0][2] == "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"   # lowercase and line breaks normalised
    assert recs[1][:2] == ("two", 4) and len(recs[1][2]) == 32


def test_sequences_command_writes_aliases(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">s1\nACGT\n>s2\nGGGG\n")
    out = tmp_path / "seq.tsv"
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "make_manifest.py"), "sequences", "--fasta", str(fa),
                    "--alias-prefix", "lcl|", "-o", str(out)], check=True, capture_output=True)
    header, rows = ec.read_tsv_table(str(out))
    assert header == ["seqid", "length", "sha512t24u", "aliases"]
    assert [(r["seqid"], r["length"], r["aliases"]) for r in rows] == [("s1", "4", "lcl|s1"), ("s2", "4", "lcl|s2")]


def test_manifest_records_producer_filters_and_file_facts():
    m = yaml.safe_load(open(os.path.join(EXAMPLE, "manifest.yaml")))
    rs = m["evidence"]["rna_short"]
    assert rs["filters"]["min_total_unique_reads"] == 1 and rs["stats"]["aggregate_star_junctions.py"]["input_rows"] == 55
    assert rs["outputs"][0]["records"] == 10 and len(rs["outputs"][0]["sha256"]) == 64
    assert m["assembly"]["total_length"] == 42000 and m["annotation"]["records"] == 34
    assert any(p["name"] == "aggregate_star_junctions.py" for p in rs["programs"])


@needs_htslib
def test_rebuilding_the_example_reproduces_the_committed_bundle(tmp_path):
    out = tmp_path / "rebuilt"
    env = dict(os.environ, PYTHON=sys.executable)
    res = subprocess.run(["bash", os.path.join(ROOT, "examples", "minimal", "build.sh"), str(out)], env=env, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr
    for dirpath, _, files in os.walk(EXAMPLE):
        for f in files:
            rel = os.path.relpath(os.path.join(dirpath, f), EXAMPLE)
            a, b = os.path.join(EXAMPLE, rel), out / rel
            assert b.exists(), rel
            if f.endswith(".tbi") or f.endswith(".bw"):
                continue                              # index/bigWig bytes depend on the tool build; content is compared below
            if f.endswith(".gz"):
                assert gzip.open(a).read() == gzip.open(b).read(), rel
            elif f == "manifest.yaml":
                ma, mb = yaml.safe_load(open(a)), yaml.safe_load(open(b))
                for m_ in (ma, mb):
                    m_["annotation"].pop("bytes"), m_["annotation"].pop("checksum")
                    for ev in m_["evidence"].values():
                        for o in ev["outputs"]:
                            o.pop("bytes"), o.pop("sha256")
                assert ma == mb
            else:
                assert open(a, "rb").read() == open(b, "rb").read(), rel
    if ec.bigwig_backend():
        assert ec.bigwig_intervals(str(out / "evidence" / "coverage.bw"), "scaffold_1", 1, 30000) == \
               ec.bigwig_intervals(os.path.join(EXAMPLE, "evidence", "coverage.bw"), "scaffold_1", 1, 30000)


def test_multi_line_features_may_share_an_id_but_other_duplicates_are_errors(tmp_bundle, capsys):
    """GFF3 allows a CDS split over several lines with one ID (NCBI does this); a different type with the same ID is an error."""
    def share_cds_id(t):
        return t.replace("Parent=helixer_g1.t1", "ID=cds-shared;Parent=helixer_g1.t1")        # all Helixer geneA exons/CDS lines
    edit_file(tmp_bundle, "evidence/predictions.gff3.gz", share_cds_id)
    rc, out = validate(tmp_bundle, capsys)
    assert "duplicate ID cds-shared" in out and rc == 1                                        # exon and CDS lines now share an ID: different types


def test_same_type_same_id_is_accepted(tmp_bundle, capsys):
    def share(t):
        return "\n".join(l.replace("Parent=helixer_g1.t1", "ID=cds-shared;Parent=helixer_g1.t1") if "\tCDS\t" in l and "helixer_g1.t1" in l else l
                         for l in t.split("\n"))
    edit_file(tmp_bundle, "evidence/predictions.gff3.gz", share)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 0 and "duplicate ID" not in out


def test_manifest_template_is_valid_yaml_with_the_required_classes():
    tpl = yaml.safe_load(open(os.path.join(ROOT, "schema", "manifest.template.yaml")))
    assert set(tpl) >= {"spec_version", "bundle", "assembly", "annotation", "evidence"}
    assert set(tpl["evidence"]) >= {"rna_short", "rna_long", "proteins", "predictions"}
    assert tpl["evidence"]["proteins"]["databases"] and tpl["evidence"]["predictions"]["programs"][0]["source_tag"]
    for cls, ev in tpl["evidence"].items():                          # every class either has the full provenance fields or a reason
        assert ("reason" in ev) if not ev["present"] else {"programs", "summarisation", "outputs"} <= set(ev), cls


def test_single_exon_models_are_valid_but_must_have_one_block(tmp_bundle, capsys):
    rc, out = validate(tmp_bundle, capsys)                              # the example contains single-exon models (chain ".")
    assert rc == 0
    def two_blocks_no_chain(t):
        return replace_line(t, "LR000005", lambda l: "\t".join(l.split("\t")[:9] + ["2", "400,400,", "0,500,"] + l.split("\t")[12:]))
    edit_file(tmp_bundle, "evidence/long_read_models.bed.gz", two_blocks_no_chain)
    rc, out = validate(tmp_bundle, capsys)
    assert rc == 1 and "does not match blocks" in out                 # blocks with a gap imply an intron, so the chain "." is wrong
    def touching_blocks(t):                                             # no intron between the blocks, so the chain "." is right, but it is not one block
        return replace_line(t, "LR000005", lambda l: "\t".join(l.split("\t")[:9] + ["2", "450,450,", "0,450,"] + l.split("\t")[12:]))
    shutil.copytree(EXAMPLE, os.path.join(os.path.dirname(tmp_bundle), "b2"))
    b2 = os.path.join(os.path.dirname(tmp_bundle), "b2")
    edit_file(b2, "evidence/long_read_models.bed.gz", touching_blocks)
    rc, out = validate(b2, capsys)
    assert rc == 1 and "must have exactly one block" in out


def test_several_summaries_of_one_class_are_merged_and_counts_add_up():
    ev = {"present": True}
    a = {"tool": "bam_to_junctions.py", "version": "0.1.0", "parameters": {"x": 1}, "filters": {"min_unique": 0}, "stats": {"junctions_seen": 10, "note": "a"}}
    b = {"tool": "bam_to_junctions.py", "version": "0.1.0", "parameters": {"x": 1}, "filters": {"min_unique": 0}, "stats": {"junctions_seen": 5, "note": "b"}}
    c = {"tool": "aggregate_star_junctions.py", "version": "0.1.0", "parameters": {}, "filters": {"min_samples": 2}, "stats": {"output_rows": 7}}
    for s in (a, b, c):
        mm.merge_summary(ev, s)
    assert ev["stats"] == {"bam_to_junctions.py": {"junctions_seen": 15, "note": "b"}, "aggregate_star_junctions.py": {"output_rows": 7}}
    assert [p["name"] for p in ev["programs"]] == ["bam_to_junctions.py", "aggregate_star_junctions.py"]
    assert ev["filters"] == {"min_unique": 0, "min_samples": 2}
