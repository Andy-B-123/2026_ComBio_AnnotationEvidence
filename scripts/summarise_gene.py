#!/usr/bin/env python3
"""Summarise the retained evidence for one annotated gene, intron by intron.

For every intron of every transcript the report keeps three states apart:
  * SUPPORT       an evidence record matches the intron exactly (same first and last intron base);
  * ABSENCE       the evidence class exists in the bundle but has no exactly matching record
                  (this is NOT evidence that the intron is wrong: see the filters printed at the top);
  * COMPETING     a different junction shares exactly one splice site with the intron
                  (alternative acceptor/donor or exon skipping), shown with its support.
It also reports terminal boundaries, mean depth of exons and introns (if a coverage track and a
bigWig reader are available), and evidence that spans this gene and a neighbouring annotated gene
(possible merge / split indicators). Nothing here decides that a model is right or wrong.

Example:
  summarise_gene.py bundle/ geneA
  summarise_gene.py bundle/ geneA --format json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_common as ec  # noqa: E402

__version__ = "0.1.0"


def compatible(item_strand: str, gene_strand: str) -> bool:
    return item_strand in (".", gene_strand) or gene_strand == "."


def near(intron, others, tol):
    """Introns in ``others`` within ``tol`` bases of ``intron`` at both ends but not identical."""
    s, e = intron
    return [o for o in others if o != intron and abs(o[0] - s) <= tol and abs(o[1] - e) <= tol]


def competitors(intron, strand, all_introns, junctions, annotated_total, ratio, min_samples):
    """Junctions sharing exactly one splice site with ``intron``."""
    s, e = intron
    donor_is_start = strand != "-"
    other_starts = {i[0] for i in all_introns if i != intron}
    other_ends = {i[1] for i in all_introns if i != intron}
    out = []
    for j in junctions:
        if not compatible(j["strand"], strand):
            continue
        shares_start, shares_end = j["start"] == s, j["end"] == e
        if shares_start == shares_end:
            continue
        shares_donor = shares_start == donor_is_start
        kind = "alt_acceptor" if shares_donor else "alt_donor"
        differing = j["end"] if shares_start else j["start"]
        if differing in (other_ends if shares_start else other_starts):
            kind = "exon_skipping"
        out.append({"junction": [j["start"], j["end"]], "type": kind,
                    "shift_bp": (j["end"] - e) if shares_start else (j["start"] - s),
                    "sample_count": j["sample_count"], "group_count": j["group_count"],
                    "total_unique_reads": j["total_unique_reads"], "max_unique_reads": j["max_unique_reads"],
                    "fraction_of_annotated": None if annotated_total == 0 else round(j["total_unique_reads"] / annotated_total, 3),
                    "dominant": j["total_unique_reads"] >= ratio * max(annotated_total, 1) and j["sample_count"] >= min_samples})
    return sorted(out, key=lambda c: -c["total_unique_reads"])[:5]


def class_state(b, cls):
    ev = (b.classes().get(cls) or {})
    return "present" if ev.get("present") else "absent"


def genes_in(path: str, seqid: str, start: int, end: int):
    hi = end
    return [g for g in ec.build_gene_models(ec.fetch_gff3(path, seqid, start, end)).values()
            if g.transcripts and ec.overlaps(g.start, g.end, start, hi)]


def analyse(b: ec.Bundle, name: str, tol: int, end_tol: int, ratio: float, min_comp: int, annotation: str | None = None) -> dict:
    """Evidence report for one gene. ``annotation`` overrides the bundle's annotation, so that a candidate
    (competing) gene model in a separate GFF3 can be scored against the same retained evidence."""
    ann_path = annotation or b.annotation_path
    gene = ec.find_gene(ann_path, name)
    seqid, gs, ge, strand = gene.seqid, gene.start, gene.end, gene.strand
    state = {c: class_state(b, c) for c in ("rna_short", "rna_long", "proteins", "predictions", "coverage")}
    junc = [j for j in b.junctions(seqid, gs - 1 if gs > 1 else 1, ge + 1)] if state["rna_short"] == "present" else []
    all_models = [m for m in b.long_read_models(seqid, gs, ge) if compatible(m["strand"], strand)] if state["rna_long"] == "present" else []
    models = [m for m in all_models if m["introns"]]                 # spliced models: introns, chains, cross-gene evidence
    single = [m for m in all_models if not m["introns"]]             # single-exon models: reported separately, never as intron support
    prots = [a for a in b.protein_alignments(seqid, gs, ge) if compatible(a["strand"], strand)] if state["proteins"] == "present" else []
    preds = [g for g in b.prediction_genes(seqid, gs, ge) if compatible(g.strand, strand)] if state["predictions"] == "present" else []
    lo = min([gs] + [m["start"] for m in all_models] + [a["start"] for a in prots] + [g.start for g in preds])
    hi = max([ge] + [m["end"] for m in all_models] + [a["end"] for a in prots] + [g.end for g in preds])
    neighbours = [g for g in genes_in(ann_path, seqid, lo, hi) if compatible(g.strand, strand)]
    cov = b.coverage_intervals(seqid, lo, hi) if state["coverage"] == "present" else None

    all_introns = sorted({i for t in gene.transcripts for i in t.introns})
    report = {"tool": f"summarise_gene.py {__version__}", "bundle": (b.manifest.get("bundle") or {}).get("id"),
              "annotation_used": os.path.basename(ann_path) if annotation else "bundle annotation",
              "coordinates": "1-based-closed (intron = first..last intron base)",
              "gene": {"id": gene.id, "names": sorted(gene.names), "seqid": seqid, "start": gs, "end": ge, "strand": strand,
                       "annotated_partial": gene.attrs.get("partial") == "true"},
              "parameters": {"near_tolerance_bp": tol, "end_tolerance_bp": end_tol, "competing_min_ratio": ratio,
                             "competing_min_samples": min_comp},
              "classes": state, "filters": {c: b.filters(c) for c in ("rna_short", "rna_long", "proteins") if state[c] == "present"},
              "n_datasets": {c: (b.classes().get(c) or {}).get("n_datasets") for c in ("rna_short", "rna_long", "proteins", "predictions")},
              "guided_by_annotation": (b.classes().get("rna_short") or {}).get("guided_by_annotation"),
              "coverage_reader": "available" if (cov is not None) else ("unavailable (install pyBigWig or UCSC bigWigToBedGraph)" if state["coverage"] == "present" else "no coverage track"),
              "transcripts": []}

    for t in gene.transcripts:
        tr = {"id": t.id, "exons": t.exons, "cds": t.cds, "introns": [], "partial": t.partial}
        for k, (s, e) in enumerate(t.introns, 1):
            rec = {"index": k, "intron": [s, e], "length": e - s + 1}
            sup, absent = [], []
            # short-read RNA junctions
            if state["rna_short"] == "present":
                hit = [j for j in junc if (j["start"], j["end"]) == (s, e) and compatible(j["strand"], strand)]
                if hit:
                    j = hit[0]
                    rec["rna_short"] = {"match": "exact", **{k2: j[k2] for k2 in ("sample_count", "group_count", "total_unique_reads", "max_unique_reads",
                                                                                   "median_unique_reads", "total_multimapping_reads", "max_overhang", "motif")}}
                    sup.append("rna_short")
                else:
                    rec["rna_short"] = {"match": "not_observed"}
                    absent.append("rna_short")
                ann_total = rec["rna_short"].get("total_unique_reads", 0)
                rec["competing_junctions"] = competitors((s, e), strand, t.introns, junc, ann_total, ratio, min_comp)
            # long reads
            if state["rna_long"] == "present":
                exact = [m for m in models if (s, e) in m["introns"]]
                nr = [m["id"] for m in models if near((s, e), m["introns"], tol) and (s, e) not in m["introns"]]
                if exact:
                    rec["rna_long"] = {"match": "exact", "models": len(exact), "reads": sum(m.get("read_count") or 0 for m in exact),
                                       "max_sample_count": max(m.get("sample_count") or 0 for m in exact), "model_ids": [m["id"] for m in exact][:6]}
                    sup.append("rna_long")
                else:
                    rec["rna_long"] = {"match": "not_observed", "near_match_models": nr[:6]}
                    absent.append("rna_long")
            # proteins
            if state["proteins"] == "present":
                exact = [a for a in prots if (s, e) in a["introns"]]
                if exact:
                    ids = [a["attrs"].get("pct_identity") for a in exact]
                    best = max((float(x) for x in ids if x not in (None, ".")), default=None)
                    rec["proteins"] = {"match": "exact", "alignments": len(exact), "alignment_ids": [a["id"] for a in exact][:6],
                                       "proteins": sorted({a["protein"] for a in exact})[:6], "best_pct_identity": best}
                    sup.append("proteins")
                else:
                    nr = [a["id"] for a in prots if near((s, e), a["introns"], tol)]
                    rec["proteins"] = {"match": "not_observed", "near_match_alignments": nr[:6]}
                    absent.append("proteins")
            # predictions
            pred_sources = []
            if state["predictions"] == "present":
                pred_sources = sorted({g.source for g in preds for pt in g.transcripts if (s, e) in pt.introns})
                near_src = sorted({g.source for g in preds for pt in g.transcripts if near((s, e), pt.introns, tol)})
                rec["predictions"] = {"match": "exact" if pred_sources else "not_observed", "sources": pred_sources, "near_match_sources": near_src}
            rec["supported_by"] = sup
            rec["not_observed_in"] = absent
            rec["status"] = "supported" if sup else ("prediction_only" if pred_sources else "no_retained_evidence")
            rec["competing_evidence"] = any(c["dominant"] for c in rec.get("competing_junctions", []))
            if cov is not None:
                left = ec.mean_depth(cov, max(lo, t.exons[k - 1][0]), t.exons[k - 1][1])
                right = ec.mean_depth(cov, t.exons[k][0], t.exons[k][1])
                rec["coverage"] = {"intron_mean_depth": round(ec.mean_depth(cov, s, e), 2), "flanking_exon_mean_depth": [round(left, 2), round(right, 2)]}
            tr["introns"].append(rec)

        # full-chain and terminal-boundary support
        chain = t.introns
        tr["full_chain"] = {}
        if state["rna_long"] == "present" and chain:
            tr["full_chain"]["rna_long_exact_chain"] = {"models": sum(1 for m in models if m["introns"] == chain),
                                                        "reads": sum(m.get("read_count") or 0 for m in models if m["introns"] == chain)}
            tr["full_chain"]["rna_long_models_containing_all_introns"] = {
                "models": sum(1 for m in models if all(i in m["introns"] for i in chain)),
                "reads": sum(m.get("read_count") or 0 for m in models if all(i in m["introns"] for i in chain))}
        if state["proteins"] == "present" and chain:
            tr["full_chain"]["proteins_containing_all_introns"] = [a["id"] for a in prots if all(i in a["introns"] for i in chain)][:6]
        if state["predictions"] == "present" and chain:
            tr["full_chain"]["predictions_exact_chain"] = sorted({g.source for g in preds for pt in g.transcripts if pt.introns == chain})
        cds_lo = t.cds[0][0] if t.cds else t.start
        cds_hi = t.cds[-1][1] if t.cds else t.end
        tr["terminal_boundaries"] = {}
        for label, coord, cds_coord in (("start", t.start, cds_lo), ("end", t.end, cds_hi)):
            b_rec = {"coordinate": coord, "cds_coordinate": cds_coord,
                     "transcript_end": "5-prime" if (label == "start") == (strand != "-") else "3-prime"}
            if state["rna_long"] == "present":
                pool = models if chain else all_models                # a single-exon transcript can be tested by single-exon models
                near_m = [m for m in pool if abs((m["start"] if label == "start" else m["end"]) - coord) <= end_tol]
                b_rec["rna_long"] = {"models_within_tolerance": len(near_m), "reads": sum(m.get("read_count") or 0 for m in near_m),
                                     "closest_delta_bp": min((((m["start"] if label == "start" else m["end"]) - coord) for m in pool),
                                                              key=abs, default=None)}
            if state["proteins"] == "present":
                near_a = [a for a in prots if abs((a["start"] if label == "start" else a["end"]) - cds_coord) <= end_tol]
                b_rec["proteins"] = {"alignments_within_tolerance_of_cds_end": len(near_a)}
            if state["predictions"] == "present":
                b_rec["predictions"] = sorted({g.source for g in preds for pt in g.transcripts
                                               if abs((pt.start if label == "start" else pt.end) - coord) <= end_tol})
            tr["terminal_boundaries"][label] = b_rec
        if state["rna_long"] == "present" and not chain:
            tl = t.end - t.start + 1
            def recip(m):
                ovl = min(t.end, m["end"]) - max(t.start, m["start"]) + 1
                return 0.0 if ovl <= 0 else round(ovl / max(tl, m["end"] - m["start"] + 1), 3)
            ranked = sorted(all_models, key=lambda m: (-(recip(m) if not m["introns"] else 0), -(m.get("read_count") or 0)))
            tr["single_exon_transcript"] = {
                "single_exon_models": [{"id": m["id"], "start": m["start"], "end": m["end"], "reads": m.get("read_count"),
                                        "samples": m.get("sample_count"), "groups": m.get("group_count"), "reciprocal_overlap": recip(m)}
                                       for m in ranked if not m["introns"] and recip(m) > 0][:5],
                "spliced_models_overlapping": {"models": sum(1 for m in models if recip(m) > 0),
                                               "reads": sum(m.get("read_count") or 0 for m in models if recip(m) > 0),
                                               "ids": [m["id"] for m in models if recip(m) > 0][:5]}}
        if cov is not None:
            tr["coverage"] = {"exons": [[s, e, round(ec.mean_depth(cov, s, e), 2)] for s, e in t.exons],
                              "introns": [[s, e, round(ec.mean_depth(cov, s, e), 2)] for s, e in t.introns]}
        report["transcripts"].append(tr)

    report["single_exon_long_reads"] = ({"models": len(single), "reads": sum(m.get("read_count") or 0 for m in single),
                                          "max_sample_count": max((m.get("sample_count") or 0 for m in single), default=0),
                                          "ids": [m["id"] for m in single][:8],
                                          "note": "unspliced or single-exon reads: counted separately; they may be fragments, unspliced RNA or single-exon genes"}
                                         if state["rna_long"] == "present" else None)
    report["cross_locus"] = cross_locus(gene, neighbours, junc, models, prots, preds)
    counts = {}
    for tr in report["transcripts"]:
        for rec in tr["introns"]:
            counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            if rec["competing_evidence"]:
                counts["with_competing_evidence"] = counts.get("with_competing_evidence", 0) + 1
    report["intron_status_counts"] = counts
    return report


def cross_locus(gene, neighbours, junc, models, prots, preds):
    """Evidence touching this gene and another annotated gene on the same strand."""
    others = [g for g in neighbours if g.id != gene.id]

    def hits(exons):
        out = []
        for g in neighbours:
            gex = [x for t in g.transcripts for x in t.exons]
            if any(ec.overlaps(a, b, c, d) for a, b in exons for c, d in gex):
                out.append(g)
        return out

    res = {"neighbouring_genes": [{"id": g.id, "start": g.start, "end": g.end, "annotated_partial": g.attrs.get("partial") == "true"} for g in others],
           "spanning": {"rna_long": [], "proteins": [], "predictions": [], "rna_short_junctions": []}, "single_gene": {}}
    annotated_introns = {i for g in neighbours for t in g.transcripts for i in t.introns}
    only, multi = 0, 0
    reads_only, reads_multi = 0, 0
    for m in models:
        h = hits(m["exons"])
        if len(h) >= 2:
            multi += 1
            reads_multi += m.get("read_count") or 0
            res["spanning"]["rna_long"].append({"id": m["id"], "genes": [g.id for g in h], "reads": m.get("read_count"),
                                               "samples": m.get("sample_count"), "groups": m.get("group_count"),
                                               "novel_introns": [list(i) for i in m["introns"] if i not in annotated_introns]})
        else:
            only += 1
            reads_only += m.get("read_count") or 0
    res["single_gene"]["rna_long"] = {"models_touching_only_one_gene": only, "reads": reads_only,
                                      "models_spanning_multiple_genes": multi, "reads_spanning": reads_multi}
    only, multi = 0, 0
    for a in prots:
        h = hits(a["exons"])
        if len(h) >= 2:
            multi += 1
            per_gene = {}
            for g in h:
                gex = [x for t in g.transcripts for x in t.exons]
                tr_ = [(bl["target_start"], bl["target_end"]) for bl in a["blocks"] if any(ec.overlaps(bl["start"], bl["end"], c, d) for c, d in gex)]
                per_gene[g.id] = [min(x[0] for x in tr_), max(x[1] for x in tr_)] if tr_ else None
            res["spanning"]["proteins"].append({"id": a["id"], "protein": a["protein"], "protein_length": a["attrs"].get("query_length"),
                                               "pct_identity": a["attrs"].get("pct_identity"), "genes_with_protein_aa_range": per_gene})
        else:
            only += 1
    res["single_gene"]["proteins"] = {"alignments_touching_only_one_gene": only, "alignments_spanning_multiple_genes": multi}
    only, multi = [], []
    for g in preds:
        (multi if len(hits([x for t in g.transcripts for x in t.exons])) >= 2 else only).append(g.source)
    res["single_gene"]["predictions"] = {"sources_touching_only_one_gene": sorted(set(only)), "sources_spanning_multiple_genes": sorted(set(multi))}
    res["spanning"]["predictions"] = sorted(set(multi))
    ends, starts = {}, {}
    for g in neighbours:
        for t in g.transcripts:
            for x in t.exons:
                ends.setdefault(x[1], set()).add(g.id)
                starts.setdefault(x[0], set()).add(g.id)
    for j in junc:
        left, right = ends.get(j["start"] - 1, set()), starts.get(j["end"] + 1, set())
        pairs = [(a, b) for a in left for b in right if a != b]
        if pairs:
            res["spanning"]["rna_short_junctions"].append({"junction": [j["start"], j["end"]], "from_gene": pairs[0][0], "to_gene": pairs[0][1],
                                                          "sample_count": j["sample_count"], "group_count": j["group_count"],
                                                          "total_unique_reads": j["total_unique_reads"], "max_unique_reads": j["max_unique_reads"]})
    bridged = {}
    for cls, items in (("rna_long", res["spanning"]["rna_long"]), ("proteins", res["spanning"]["proteins"]), ("rna_short", res["spanning"]["rna_short_junctions"])):
        for it in items:
            gs_ = it.get("genes") or list(it.get("genes_with_protein_aa_range", {})) or [it.get("from_gene"), it.get("to_gene")]
            for gid in gs_:
                if gid and gid != gene.id:
                    bridged.setdefault(gid, set()).add(cls)
    for src in res["spanning"]["predictions"]:
        for g in others:
            bridged.setdefault(g.id, set()).add(f"predictions:{src}")
    res["bridged_genes"] = {k: sorted(v) for k, v in bridged.items()}
    return res


# --------------------------------------------------------------------------------------
def fmt_filters(f):
    return "; ".join(f"{k}={v}" for k, v in f.items()) if f else "none recorded"


def render_md(r: dict) -> str:
    g = r["gene"]
    L = [f"# Evidence summary for {g['id']}", "",
         f"Bundle `{r['bundle']}`" + (f" (gene model taken from {r['annotation_used']}, not the bundle annotation)" if r["annotation_used"] != "bundle annotation" else "")
         + f". Gene {g['id']} at {g['seqid']}:{g['start']}-{g['end']} ({g['strand']})"
         + (" - annotated as partial" if g["annotated_partial"] else "") + ".",
         f"All coordinates are 1-based, closed; an intron is given as its first and last intron base.", ""]
    L.append("## What absence can and cannot mean here")
    for c in ("rna_short", "rna_long", "proteins", "predictions", "coverage"):
        if r["classes"][c] == "absent":
            L.append(f"- {c}: NOT PRESENT in this bundle, so nothing below says anything about it.")
    for c, f in r["filters"].items():
        n = r["n_datasets"].get(c)
        L.append(f"- {c}: {n} dataset(s); filters applied before the records were kept: {fmt_filters(f)}.")
    if r.get("guided_by_annotation"):
        L.append("- rna_short: mapping was GUIDED by an annotation, so junction support may be partly circular.")
    L.append("- A junction that is not listed may still have been seen and filtered out, or the tissue may not have been sampled.")
    L.append("")
    for t in r["transcripts"]:
        L.append(f"## Transcript {t['id']}: {len(t['exons'])} exon(s), {len(t['introns'])} intron(s)" + (" (annotated partial)" if t["partial"] else ""))
        L.append(f"Exons: {spans(t['exons'])}. CDS: {spans(t['cds'])}.")
        L.append("")
        for rec in t["introns"]:
            flag = "  + COMPETING EVIDENCE" if rec["competing_evidence"] else ""
            L.append(f"### Intron {rec['index']}: {rec['intron'][0]}-{rec['intron'][1]} ({rec['length']} bp)  ->  {rec['status']}{flag}")
            if "rna_short" in rec:
                x = rec["rna_short"]
                if x["match"] == "exact":
                    gc = "?" if x["group_count"] is None else x["group_count"]
                    L.append(f"- rna_short: exact match. {x['sample_count']} sample(s) / {gc} group(s); unique reads total {x['total_unique_reads']}, "
                             f"max in one sample {x['max_unique_reads']}, median {x['median_unique_reads']:g}; multimapping {x['total_multimapping_reads']}; "
                             f"max overhang {x['max_overhang']}; motif {x['motif']}.")
                else:
                    L.append("- rna_short: not observed in the retained junction table.")
            if "rna_long" in rec:
                x = rec["rna_long"]
                if x["match"] == "exact":
                    L.append(f"- rna_long: exact match in {x['models']} model(s) with {x['reads']} read(s) (largest sample count {x['max_sample_count']}): {', '.join(x['model_ids'])}.")
                else:
                    L.append("- rna_long: not observed in any retained model" + (f" (near matches within tolerance: {', '.join(x['near_match_models'])})" if x["near_match_models"] else "") + ".")
            if "proteins" in rec:
                x = rec["proteins"]
                if x["match"] == "exact":
                    bi = "" if x["best_pct_identity"] is None else f", best identity {x['best_pct_identity']:g}%"
                    L.append(f"- proteins: exact match in {x['alignments']} alignment(s) ({', '.join(x['proteins'])}){bi}.")
                else:
                    L.append("- proteins: not observed at this exact intron" + (f" (near matches: {', '.join(x['near_match_alignments'])})" if x["near_match_alignments"] else "") + ".")
            if "predictions" in rec:
                x = rec["predictions"]
                L.append(f"- predictions: exact match by {', '.join(x['sources'])}." if x["sources"] else
                         "- predictions: no exact match" + (f" (near: {', '.join(x['near_match_sources'])})" if x["near_match_sources"] else "") + ".")
            for c in rec.get("competing_junctions", []):
                fr = "" if c["fraction_of_annotated"] is None else f", {c['fraction_of_annotated']:.1%} of the annotated junction's unique reads"
                gc = "?" if c["group_count"] is None else c["group_count"]
                L.append(f"- competing junction {c['junction'][0]}-{c['junction'][1]} ({c['type']}, shift {c['shift_bp']:+d} bp): {c['sample_count']} sample(s)/{gc} group(s), "
                         f"{c['total_unique_reads']} unique reads (max {c['max_unique_reads']}){fr}" + ("  [DOMINANT: at least as well supported as the annotated intron]" if c["dominant"] else "") + ".")
            if "coverage" in rec:
                cv = rec["coverage"]
                L.append(f"- coverage: intron mean depth {cv['intron_mean_depth']}; flanking exons {cv['flanking_exon_mean_depth'][0]} and {cv['flanking_exon_mean_depth'][1]}.")
            L.append("")
        fc = t["full_chain"]
        if fc:
            L.append("### Whole-transcript intron chain")
            if "rna_long_exact_chain" in fc:
                a, b_ = fc["rna_long_exact_chain"], fc["rna_long_models_containing_all_introns"]
                L.append(f"- rna_long: {a['models']} model(s) ({a['reads']} reads) have exactly this intron chain; {b_['models']} model(s) ({b_['reads']} reads) contain all of its introns (chain may be longer).")
            if "proteins_containing_all_introns" in fc:
                L.append(f"- proteins containing all introns: {', '.join(fc['proteins_containing_all_introns']) or 'none'}.")
            if "predictions_exact_chain" in fc:
                L.append(f"- predictions with exactly this chain: {', '.join(fc['predictions_exact_chain']) or 'none'}.")
            L.append("")
        if "single_exon_transcript" in t:
            se = t["single_exon_transcript"]
            L.append("### Single-exon transcript: long-read evidence (no intron to test; overlap = shared bases / longer of the two)")
            if se["single_exon_models"]:
                for m in se["single_exon_models"]:
                    L.append(f"- single-exon model {m['id']} {m['start']}-{m['end']}: reciprocal overlap {m['reciprocal_overlap']:.2f}; "
                             f"{m['reads']} read(s), {m['samples']} sample(s), {m['groups']} group(s).")
            else:
                L.append("- no single-exon long-read model overlaps this transcript.")
            sp = se["spliced_models_overlapping"]
            L.append(f"- spliced long-read models overlapping this transcript: {sp['models']} ({sp['reads']} reads)" + (f": {', '.join(sp['ids'])}. A spliced structure at a locus annotated without introns is competing evidence." if sp["models"] else "."))
            L.append("")
        L.append("### Terminal boundaries (junction data cannot test these; coordinates are genomic, left = lower)")
        for lab, x in t["terminal_boundaries"].items():
            bits = []
            if "rna_long" in x:
                bits.append(f"rna_long: {x['rna_long']['models_within_tolerance']} model(s) / {x['rna_long']['reads']} read(s) end within {r['parameters']['end_tolerance_bp']} bp (closest delta {x['rna_long']['closest_delta_bp']} bp)")
            if "proteins" in x:
                bits.append(f"proteins: {x['proteins']['alignments_within_tolerance_of_cds_end']} alignment(s) end within tolerance of the CDS boundary {x['cds_coordinate']}")
            if "predictions" in x:
                bits.append("predictions: " + (", ".join(x["predictions"]) or "none") + " within tolerance")
            what = "; ".join(bits) if bits else "no long-read, protein or prediction class in this bundle can test this end"
            L.append(f"- genomic {'left' if lab == 'start' else 'right'} end {x['coordinate']} (the {x['transcript_end'].replace('-', ' ')} end of this transcript): {what}.")
        L.append("")
        if "coverage" in t:
            L.append("### Coverage (mean depth, quantised sum over libraries)")
            L.append("- exons: " + ", ".join(f"{s}-{e}: {d}" for s, e, d in t["coverage"]["exons"]))
            if t["coverage"]["introns"]:
                L.append("- introns: " + ", ".join(f"{s}-{e}: {d}" for s, e, d in t["coverage"]["introns"]))
            L.append("")
    se = r.get("single_exon_long_reads")
    if se:
        L.append("## Single-exon long-read models overlapping this gene")
        if se["models"]:
            L.append(f"- {se['models']} model(s), {se['reads']} read(s) (largest sample count {se['max_sample_count']}): {', '.join(se['ids'])}. {se['note']}.")
        else:
            L.append("- none overlap this gene.")
        L.append("")
    x = r["cross_locus"]
    L.append("## Evidence spanning this gene and a neighbouring annotated gene")
    if not x["neighbouring_genes"]:
        L.append("- No other annotated gene on the same strand lies within the span of the evidence that overlaps this gene.")
    for n in x["neighbouring_genes"]:
        L.append(f"- neighbouring gene {n['id']} {n['start']}-{n['end']}" + (" (annotated partial)" if n["annotated_partial"] else ""))
    sl = x["single_gene"]
    if "rna_long" in sl:
        L.append(f"- rna_long (spliced models only): {sl['rna_long']['models_touching_only_one_gene']} model(s) ({sl['rna_long']['reads']} reads) touch only one gene; "
                 f"{sl['rna_long']['models_spanning_multiple_genes']} model(s) ({sl['rna_long']['reads_spanning']} reads) span several genes.")
    for m in x["spanning"]["rna_long"]:
        L.append(f"  - {m['id']}: genes {', '.join(m['genes'])}; {m['reads']} read(s), {m['samples']} sample(s), {m['groups']} group(s); introns not in the annotation: "
                 + (", ".join(f"{a}-{b}" for a, b in m["novel_introns"]) or "none"))
    if "proteins" in sl:
        L.append(f"- proteins: {sl['proteins']['alignments_touching_only_one_gene']} alignment(s) touch only one gene; {sl['proteins']['alignments_spanning_multiple_genes']} span several.")
    for a in x["spanning"]["proteins"]:
        L.append(f"  - {a['id']} ({a['protein']}, {a['protein_length']} aa, identity {a['pct_identity']}%): " + "; ".join(
            f"{gid}: aa {'-'.join(map(str, rng)) if rng else '?'}" for gid, rng in a["genes_with_protein_aa_range"].items()))
    for j in x["spanning"]["rna_short_junctions"]:
        L.append(f"- rna_short junction {j['junction'][0]}-{j['junction'][1]} joins an exon end of {j['from_gene']} to an exon start of {j['to_gene']}: "
                 f"{j['sample_count']} sample(s), {j['group_count']} group(s), {j['total_unique_reads']} unique reads (max {j['max_unique_reads']}).")
    if "predictions" in sl:
        L.append(f"- predictions: spanning several genes: {', '.join(sl['predictions']['sources_spanning_multiple_genes']) or 'none'}; "
                 f"touching only one: {', '.join(sl['predictions']['sources_touching_only_one_gene']) or 'none'}.")
    if x["bridged_genes"]:
        L.append("- classes that bridge to a neighbouring gene: " + "; ".join(f"{k}: {', '.join(v)}" for k, v in x["bridged_genes"].items()) + ".")
        L.append("  This is a REASON TO LOOK CLOSER (one gene wrongly split? read-through transcription? chimeric alignment?), not a finding that the models are wrong. "
                 "Compare the bridging evidence with the evidence that supports each gene on its own.")
    L.append("")
    L.append("## Summary")
    L.append("- intron status counts: " + (", ".join(f"{k}={v}" for k, v in sorted(r["intron_status_counts"].items())) if r["intron_status_counts"] else "none (no introns)") + ".")
    L.append("- statuses: `supported` = an exact match in rna_short, rna_long or proteins; `prediction_only` = only ab initio predictions match exactly; "
             "`no_retained_evidence` = no exact match in any class (absence is not proof of error).")
    return "\n".join(L) + "\n"


def spans(iv):
    return ", ".join(f"{s}-{e}" for s, e in iv) if iv else "-"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n\n", 1)[1])
    ap.add_argument("bundle", help="bundle directory")
    ap.add_argument("gene", help="gene ID, Name or locus_tag as in annotation.gff3.gz")
    ap.add_argument("--format", choices=("md", "json"), default="md")
    ap.add_argument("--annotation", help="GFF3 (plain or bgzip) holding the gene model to examine instead of the bundle's annotation, e.g. a "
                    "competing model proposed by a curator; it is scored against the same retained evidence")
    ap.add_argument("--near-tolerance", type=int, default=5, help="bases within which an intron counts as a near match (default 5)")
    ap.add_argument("--end-tolerance", type=int, default=50, help="bases within which a transcript end counts as supported (default 50)")
    ap.add_argument("--competing-ratio", type=float, default=1.0,
                    help="a competing junction is DOMINANT if its unique reads >= ratio x the annotated junction's (default 1.0)")
    ap.add_argument("--min-competing-samples", type=int, default=2, help="samples a competing junction needs to be DOMINANT (default 2)")
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args(argv)
    b = ec.Bundle(args.bundle)
    r = analyse(b, args.gene, args.near_tolerance, args.end_tolerance, args.competing_ratio, args.min_competing_samples, args.annotation)
    print(json.dumps(r, indent=1) if args.format == "json" else render_md(r))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ec.EvidenceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
