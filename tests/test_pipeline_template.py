"""examples/pipeline_template.sh, driven with the synthetic example's raw inputs: it must reproduce the committed evidence."""
import glob
import gzip
import os
import subprocess
import sys

import pytest
import yaml

from conftest import EXAMPLE, ROOT, needs_htslib

INPUTS = os.path.join(ROOT, "examples", "minimal", "inputs")


@needs_htslib
def test_template_builds_a_valid_bundle_equal_to_the_committed_one(tmp_path):
    prov = yaml.safe_load(open(os.path.join(ROOT, "examples", "minimal", "manifest.template.yaml")))
    prov["evidence"]["coverage"] = {"present": False, "reason": "not built in this test"}          # no UCSC tools needed
    (tmp_path / "prov.yaml").write_text(yaml.safe_dump(prov, sort_keys=False))
    sj = " ".join(f"{os.path.basename(f).split('.')[0]}={f}" for f in sorted(glob.glob(os.path.join(INPUTS, "star", "*.SJ.out.tab"))))
    lb = " ".join(f"{os.path.basename(f).split('.')[0]}={f}" for f in sorted(glob.glob(os.path.join(INPUTS, "long_reads", "*.bed12"))))
    out = tmp_path / "bundle"
    env = dict(os.environ, OUT=str(out), GENOME=os.path.join(INPUTS, "genome.fa"), ANNOTATION=os.path.join(INPUTS, "annotation.gff3"),
               SAMPLES=os.path.join(INPUTS, "samples.tsv"), PROVENANCE=str(tmp_path / "prov.yaml"), SHORT_SJ=sj, LONG_BED12=lb,
               PROTEINS_GFF=os.path.join(INPUTS, "miniprot.gff"), PROTEINS_DB="synthetic proteins", MINIPROT_VERSION="0.18",
               PREDICTIONS_GFF3=os.path.join(INPUTS, "predictions.gff3"), REPEATS_BED=os.path.join(INPUTS, "repeats.bed"),
               CREATED="2026-09-25", STAR_SETTINGS="not applicable (synthetic)", PYTHON=sys.executable)
    r = subprocess.run(["bash", os.path.join(ROOT, "examples", "pipeline_template.sh")], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "0 error(s), 0 warning(s)" in r.stdout
    for rel in ("evidence/rna_junctions.tsv.gz", "evidence/long_read_models.bed.gz", "evidence/protein_alignments.gff3.gz",
                "evidence/predictions.gff3.gz", "evidence/repeats.bed.gz", "annotation.gff3.gz"):
        a, b = os.path.join(EXAMPLE, rel), out / rel
        read = lambda p: (gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)).read()
        assert read(a) == read(b), rel
    # the committed example also demonstrates aliases (lcl|...); the template adds them only with STRIP_PREFIX
    cols = lambda p: [l.split("\t")[:3] for l in open(p) if not l.startswith("#")]
    assert cols(os.path.join(EXAMPLE, "assembly_sequences.tsv")) == cols(out / "assembly_sequences.tsv")


@needs_htslib
def test_template_refuses_to_overwrite_and_reports_missing_variables(tmp_path):
    env = dict(os.environ)
    for k in ("OUT", "GENOME", "ANNOTATION", "PROVENANCE", "SAMPLES"):
        env.pop(k, None)
    r = subprocess.run(["bash", os.path.join(ROOT, "examples", "pipeline_template.sh")], env=env, capture_output=True, text=True)
    assert r.returncode != 0 and "set OUT" in r.stderr
    (tmp_path / "manifest.yaml").write_text("x")
    env.update(OUT=str(tmp_path), GENOME="g", ANNOTATION="a", PROVENANCE="p", SAMPLES="s")
    r = subprocess.run(["bash", os.path.join(ROOT, "examples", "pipeline_template.sh")], env=env, capture_output=True, text=True)
    assert r.returncode != 0 and "exists" in r.stderr
