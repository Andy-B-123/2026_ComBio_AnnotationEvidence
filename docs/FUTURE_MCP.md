# Possible future access layer (MCP or another API)

**Nothing here is implemented.** This is a sketch of how the evidence bundle could later be exposed to agents and
services. The important point is the direction of dependency:

> An MCP server (or REST, or GraphQL) is an **access layer over the evidence files**. It is not the evidence standard.
> The bundle stays complete and useful without it, and any server must be reproducible from the bundle alone.

## Why the files first

* A bundle is a directory of bgzip + tabix files. Any language, any tool, any decade can read it.
* Region queries are already cheap (`tabix`); a server adds convenience, not capability.
* The scripts `query_region.py` and `summarise_gene.py` are the reference behaviour. A server would call the same
  functions in `scripts/evidence_common.py` rather than re-implement coordinate handling.

## Sketch of tools

All coordinates in requests and responses are **1-based, closed**, exactly as in the command-line tools. Every response
carries the relevant `filters` from the manifest so that absence can be interpreted, and every list is capped with a
`truncated: true` flag rather than silently dropped.

| Tool | Arguments | Returns |
|---|---|---|
| `get_bundle_info()` | – | manifest summary: classes present, `n_datasets`, `filters`, `guided_by_annotation`, spec version |
| `get_gene(gene_id)` | gene ID / name | gene, transcripts, exons, CDS, introns, `partial` flags |
| `get_region(seqid, start, end, classes?)` | region | the same as `query_region.py --format json` |
| `get_splice_evidence(seqid, start, end)` | region | junction rows (`sample_count`, `group_count`, unique/multimapping reads, motif) |
| `get_protein_evidence(gene_id)` | gene ID | protein alignments overlapping the gene with per-block protein ranges |
| `get_long_read_evidence(gene_id)` | gene ID | long-read models with intron chains; exact-chain and contains-all-introns counts |
| `summarise_gene_evidence(gene_id)` | gene ID | the JSON of `summarise_gene.py` |
| `compare_models(gene_id_a, gene_id_b)` | two genes | cross-gene evidence (models, alignments, junctions, predictions spanning both) with single-gene counts for balance |

Design rules a server should keep:

1. **Read-only.** The server never edits a bundle or an annotation.
2. **No verdicts.** Responses give support, absence and competing evidence separately, with the numbers; interpretation belongs to the caller
   ([LLM_GUIDE.md](LLM_GUIDE.md)). Any prose the server generates should follow the same three-state vocabulary.
3. **Provenance in every answer.** Bundle id, spec version and the manifest filters travel with the data.
4. **Bounded output.** Region size limits and record caps, with clear signals when a cap is hit.
5. **Stateless.** One bundle per server instance, or the bundle id is an argument; no server-side memory of previous queries.

## Standards it could ride on

* **GA4GH**: `htsget`-style streaming of region slices; `DRS` for distributing bundles; `Service Info` for self-description;
  `refget` digests (already in `assembly_sequences.tsv`) to confirm a client and server agree on the assembly.
  See [RESEARCH_NOTES.md](RESEARCH_NOTES.md), question 10.
* **MCP** for agent access, **OpenAPI** for everyone else. The same function set can back both.

## What would need to exist first

See "Before an MCP implementation is worthwhile" in the README: several real bundles from different pipelines, a
stable spec 1.0, validated behaviour of the reference scripts on real annotations, and evidence that people or agents actually
ask these questions.
