# Contributing

This is an early proposal. Criticism of the design is as useful as code; open an issue describing a real annotation
question the current bundle cannot answer, a format that does not fit, or a pipeline you would like supported.

## Ground rules

* **Keep it small.** No workflow manager, no server, no custom binary formats. Prefer BED, GFF3, bigWig, tabix and plain TSV.
* **Coordinates are part of the interface.** Every format states whether it is 0- or 1-based and open or closed. Any change that touches
  coordinates needs a test in `tests/test_coordinates.py` and, ideally, a deliberately injected off-by-one that the test catches.
* **Provenance and filters are recorded, never hidden.** A tool that drops records must say so in its output header and its `--summary-json`.
* **No verdicts.** Reporting tools keep *support*, *absence of retained evidence* and *competing evidence* apart (see `docs/LLM_GUIDE.md`).
  Wording that turns a missing record into a claim that a gene model is wrong will be asked to change.
* **Non-model organisms first.** Do not assume human sequence names (`chr1`), gene symbols or GENCODE-style tags.
* **Public data only.** Example data are synthetic; fixtures must come from public sources and be attributed in `tests/data/README.md`.
  Never commit unpublished sequence data or sample metadata.
* **Sources.** Claims about other tools or formats cite the documentation they came from (`docs/RESEARCH_NOTES.md`).

## Development

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-dev.txt
# htslib (bgzip, tabix) is needed for some tests; they are skipped if it is missing
python -m pytest tests -q
bash examples/minimal/build.sh /tmp/rebuilt      # rebuild the example; must reproduce examples/minimal/bundle (a test checks this)
python scripts/validate_bundle.py examples/minimal/bundle
```

Scripts are standard-library Python (≥ 3.9) apart from PyYAML. Keep them small and readable; shared code lives in `scripts/evidence_common.py`.
Changing a file format means changing `docs/SPECIFICATION.md`, the validator, the producer script, the example and the tests together, and bumping
`spec_version` in `evidence_common.py`.
