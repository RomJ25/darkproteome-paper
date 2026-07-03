#!/usr/bin/env python3
"""
class_decoy_ledger.py — emit a per-class decoy ledger ("antigen datasheet")
from target-decoy proteomics outputs you ALREADY have, making D_N a free output.

The single missing statistic for re-verifiable non-canonical (cryptic) antigen claims is the
per-class accepted DECOY count D_N (see "Canonical self in cryptic cancer-epitope catalogues
and the class-decoy ledger needed to verify non-canonical antigen claims", Rom Jan).
This tool reads existing outputs — a mokapot/Percolator PSM table, or an MSFragger pepXML
intermediate (which RETAINS decoys + class-labelable accessions) — and emits, per class:
accepted targets T, accepted decoys D, and the class-FDR estimate (D+1)/T at a chosen
threshold. The ledger is written as JSON + TSV; with --emit-diagfdr the tool also writes the
per-PSM diagFDR contract (id, is_decoy, q, pep, run, score), so a pepXML or PSM table can feed
diagFDR (which does not read pepXML) rather than rival it.

Stdlib only (no install needed). MIT.

Examples
--------
# From an MSFragger pepXML (no q-values -> target-decoy FDR sweep on hyperscore):
python3 class_decoy_ledger.py --pepxml run.pepXML --alpha 0.03 \
        --canonical-fasta swissprot_human.fasta --out ledger

# From a mokapot / Percolator / diagFDR-style PSM table (has q-values + decoy labels):
python3 class_decoy_ledger.py --psms mokapot.psms.tsv --id-col PSMId \
        --decoy-col is_decoy --q-col mokapot_qvalue --accession-col Proteins \
        --alpha 0.01 --unit psm --out ledger
"""
import argparse, csv, json, hashlib, re, os
import xml.etree.ElementTree as ET

DEFAULT_DECOY_PREFIX = "rev_"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def classify(accession, class_map):
    """De-decoyed accession -> class label. Default heuristic matches the manuscript's
    convention (UniProt = canonical; ENSP*/-Mut = variant; everything else = non-canonical)."""
    if class_map:
        for pat, cls in class_map:
            if pat in accession or re.search(pat, accession):
                return cls
        return "other"
    if accession.startswith(("sp|", "tr|")):
        return "canonical"
    if accession.startswith("ENSP") or "-Mut" in accession:
        return "variant"
    return "noncanonical"


def load_pepxml(path, decoy_prefix, score_name="hyperscore"):
    """pepXML -> top-hit PSMs. `score_name` is the search_score to read (MSFragger:
    hyperscore; other engines name theirs differently -- pass --pepxml-score-name and the
    matching --score-ascending if lower-is-better, e.g. an e-value/expect score)."""
    loc = lambda t: t.rsplit("}", 1)[-1]
    rows, idx = [], 0
    for _, el in ET.iterparse(path, events=("end",)):
        if loc(el.tag) == "search_hit" and el.get("hit_rank") == "1":
            prot = el.get("protein", "")
            is_dec = prot.startswith(decoy_prefix)
            base = prot[len(decoy_prefix):] if is_dec else prot
            score = None
            for ss in el:
                if loc(ss.tag) == "search_score" and ss.get("name") == score_name:
                    try:
                        score = float(ss.get("value"))
                    except ValueError:
                        pass
                    break
            if score is not None:
                idx += 1
                rows.append({"id": f"psm{idx}", "peptide": el.get("peptide", ""),
                             "score": score, "is_decoy": is_dec, "accession": base})
            el.clear()
    return rows


def load_psms_tsv(path, a):
    rows = []
    with open(path, newline="") as f:
        for i, row in enumerate(csv.DictReader(f, delimiter="\t")):
            acc = row.get(a.accession_col, "") if a.accession_col else ""
            dec = str(row.get(a.decoy_col, "")).strip().lower() in ("1", "true", "t", "yes", "decoy")
            if a.decoy_prefix and acc.startswith(a.decoy_prefix):
                dec = True
                acc = acc[len(a.decoy_prefix):]
            d = {"id": row.get(a.id_col, f"psm{i}"), "is_decoy": dec, "accession": acc,
                 "peptide": row.get(a.peptide_col, "") if a.peptide_col else ""}
            if a.q_col and row.get(a.q_col):
                d["q"] = float(row[a.q_col])
            if a.score_col and row.get(a.score_col):
                d["score"] = float(row[a.score_col])
            if a.class_col and row.get(a.class_col):
                d["class"] = row[a.class_col]
            if a.run_col and row.get(a.run_col):
                d["run"] = row[a.run_col]
            if a.pep_col and row.get(a.pep_col):
                d["pep"] = row[a.pep_col]
            rows.append(d)
    return rows


def accept_by_q(psms, alpha):
    return [p for p in psms if p.get("q") is not None and p["q"] <= alpha]


def dedupe_by_peptide(accepted):
    """Collapse accepted PSM rows to one row per (class, peptide, is_decoy) — for
    --unit peptide/unique-peptide, so a peptide backed by several PSMs is counted once."""
    seen = {}
    for p in accepted:
        key = (p["class"], p.get("peptide", ""), p["is_decoy"])
        if key not in seen:
            seen[key] = p
    return list(seen.values())


def multiplicity_strata(accepted):
    """Per-class PSM-replication-depth buckets (n=1 / n=2 / n>=3) -> {T, D}, for
    --stratify-multiplicity. Verified empirically on PXD055609 T1-T5 (pooled): singleton
    (n=1) PSMs carry a materially higher class-FDR than multi-PSM peptides in EVERY class
    (not cryptic-specific), while canonical vs. noncanonical PSM-multiplicity distributions
    are themselves statistically indistinguishable."""
    counts = {}
    for p in accepted:
        key = (p["class"], p.get("peptide", ""), p["is_decoy"])
        counts[key] = counts.get(key, 0) + 1
    strata = {}
    for (cls, pep, is_dec), n in counts.items():
        bucket = "n=1" if n == 1 else ("n=2" if n == 2 else "n>=3")
        c = strata.setdefault(cls, {}).setdefault(bucket, {"T": 0, "D": 0})
        c["D" if is_dec else "T"] += 1
    return strata


def assign_qvalues_sweep(psms, ascending):
    """Target-decoy q-value sweep: attach the monotone q-value p['q'] to every scored PSM.
    The accepted set {p : p['q'] <= alpha} is then the q-value rule, identical to a direct sweep."""
    ps = [p for p in psms if p.get("score") is not None]
    ps.sort(key=lambda p: p["score"], reverse=not ascending)
    cumT, cumD, fdr = 0, 0, []
    for p in ps:
        if p["is_decoy"]:
            cumD += 1
        else:
            cumT += 1
        fdr.append(cumD / max(cumT, 1))
    m = 1.0
    for i in range(len(ps) - 1, -1, -1):
        m = min(m, fdr[i]); ps[i]["q"] = m


def main():
    ap = argparse.ArgumentParser(description="Per-class decoy ledger (D_N) from target-decoy outputs.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pepxml", help="MSFragger pepXML (intermediate; retains decoys)")
    src.add_argument("--psms", help="PSM table TSV (mokapot/Percolator/diagFDR shape)")
    ap.add_argument("--id-col", default="PSMId"); ap.add_argument("--decoy-col", default="is_decoy")
    ap.add_argument("--q-col"); ap.add_argument("--score-col")
    ap.add_argument("--accession-col"); ap.add_argument("--peptide-col"); ap.add_argument("--class-col")
    ap.add_argument("--decoy-prefix", default=DEFAULT_DECOY_PREFIX)
    ap.add_argument("--class-map", help="TSV: pattern<TAB>class (substring/regex -> class)")
    ap.add_argument("--pepxml-score-name", default="hyperscore",
                    help="search_score name to read from --pepxml (MSFragger: hyperscore; other "
                         "search engines name theirs differently -- pair with --score-ascending "
                         "if it's lower-is-better, e.g. an e-value/expect score)")
    ap.add_argument("--score-ascending", action="store_true", help="lower score = better (e.g. e-value)")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--unit", default="psm", choices=["psm", "peptide", "unique-peptide"],
                    help="psm = count every accepted PSM; peptide/unique-peptide = dedupe to one "
                         "row per (class, peptide, decoy-status) first, then count")
    ap.add_argument("--stratify-multiplicity", action="store_true",
                    help="also report per-class class-FDR stratified by PSM-replication-depth "
                         "(n=1/n=2/n>=3 PSMs per peptide) — singleton IDs carry materially higher "
                         "FDR in every class (PXD055609 T1-T5)")
    ap.add_argument("--canonical-fasta", help="SwissProt FASTA -> flag canonical-self (non-novelty axis)")
    ap.add_argument("--emit-diagfdr", metavar="PATH",
                    help="also write the per-PSM diagFDR contract TSV (id, is_decoy, q, pep, run, score)")
    ap.add_argument("--run-col", help="PSM-table column naming the run/file (for --emit-diagfdr)")
    ap.add_argument("--run", help="single run label for --emit-diagfdr (default: analysis_id)")
    ap.add_argument("--pep-col", help="PSM-table column with the posterior error probability (for --emit-diagfdr)")
    ap.add_argument("--out", required=True, help="output prefix (.json + .tsv)")
    a = ap.parse_args()

    class_map = None
    if a.class_map:
        class_map = [tuple(l.rstrip("\n").split("\t")[:2]) for l in open(a.class_map)
                     if l.strip() and not l.startswith("#")]

    if a.pepxml:
        psms = load_pepxml(a.pepxml, a.decoy_prefix, a.pepxml_score_name)
        if not psms:
            ap.error(f"0 PSMs parsed from {a.pepxml}: no search_hit had a search_score "
                     f"name={a.pepxml_score_name!r}. Check --pepxml-score-name (MSFragger "
                     "uses 'hyperscore'; other search engines use a different name).")
        ascending = a.score_ascending
        software, use_q = "MSFragger(pepXML)" if a.pepxml_score_name == "hyperscore" else "pepXML", False
        fdr_basis, inpath = f"target-decoy sweep on {a.pepxml_score_name}", a.pepxml
    else:
        if not a.accession_col:
            ap.error("--accession-col is required with --psms (without it every PSM silently "
                     "classifies as 'noncanonical' -- see classify()).")
        if not a.peptide_col and (a.unit != "psm" or a.stratify_multiplicity):
            ap.error("--peptide-col is required with --psms when --unit is peptide/unique-peptide "
                     "or --stratify-multiplicity is set (without it every PSM's peptide is '', and "
                     "the per-peptide dedup silently collapses everything into one row).")
        psms, ascending = load_psms_tsv(a.psms, a), a.score_ascending
        use_q = a.q_col is not None
        if use_q and not any(p.get("q") is not None for p in psms):
            ap.error(f"--q-col {a.q_col!r} matched no values in {a.psms} -- check the column "
                     "name (the ledger would otherwise silently come out empty).")
        software = "PSM-table"
        fdr_basis = "reported q-values" if use_q else f"target-decoy sweep on {a.score_col}"
        inpath = a.psms
    for p in psms:
        p.setdefault("class", classify(p["accession"], class_map))

    if not use_q:                       # no reported q-value -> derive q from the target-decoy sweep
        assign_qvalues_sweep(psms, ascending)
    accepted = accept_by_q(psms, a.alpha)
    counted = accepted if a.unit == "psm" else dedupe_by_peptide(accepted)

    classes = {}
    for p in counted:
        c = classes.setdefault(p["class"], {"T": 0, "D": 0})
        c["D" if p["is_decoy"] else "T"] += 1
    gT = sum(c["T"] for c in classes.values()); gD = sum(c["D"] for c in classes.values())

    canon = None
    if a.canonical_fasta:
        seqs, cur = [], []
        for line in open(a.canonical_fasta):
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else: cur.append(line.strip())
        if cur: seqs.append("".join(cur))
        big = "\n".join(s.upper() for s in seqs)
        canon = {}
        for cls in classes:
            if cls in ("canonical", "variant"):
                continue
            peps = {p["peptide"].upper() for p in accepted
                    if p["class"] == cls and not p["is_decoy"] and p["peptide"]}
            self_ = sum(1 for pp in peps if pp in big)
            canon[cls] = {"n_target_peptides": len(peps), "canonical_self": self_,
                          "canonical_self_frac": round(self_ / max(len(peps), 1), 4)}

    rows = [{"class": cls, "T_class": c["T"], "D_class": c["D"],
             "class_fdr_hat": round((c["D"] + 1) / max(c["T"], 1), 5),
             "f": round(c["T"] / max(gT, 1), 5)} for cls, c in sorted(classes.items())]
    meta = {"analysis_id": os.path.basename(inpath), "software": software, "unit": a.unit,
            "alpha": a.alpha, "fdr_convention": "(D+1)/T", "fdr_basis": fdr_basis,
            "global_T": gT, "global_D": gD, "global_fdr_hat": round((gD + 1) / max(gT, 1), 5),
            "input_sha256": {os.path.basename(inpath): sha256(inpath)}, "decoy_prefix": a.decoy_prefix}
    out = {"meta": meta, "ledger": rows}
    if canon:
        out["canonical_self"] = canon

    strata = multiplicity_strata(accepted) if a.stratify_multiplicity else None
    if strata:
        out["multiplicity_strata"] = strata

    with open(a.out + ".json", "w") as f:
        json.dump(out, f, indent=2)
    with open(a.out + ".tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["analysis_id", "software", "unit", "alpha", "fdr_convention", "class",
                    "T_class", "D_class", "class_fdr_hat", "f", "global_T", "global_D"])
        for r in rows:
            w.writerow([meta["analysis_id"], software, a.unit, a.alpha, "(D+1)/T", r["class"],
                        r["T_class"], r["D_class"], r["class_fdr_hat"], r["f"], gT, gD])

    n_emit = 0
    if a.emit_diagfdr:                  # per-PSM diagFDR contract (so a pepXML / PSM table can feed diagFDR)
        run_default = a.run or meta["analysis_id"]
        with open(a.emit_diagfdr, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["id", "is_decoy", "q", "pep", "run", "score"])
            for p in psms:
                if p.get("q") is None:
                    continue
                w.writerow([p["id"], int(bool(p["is_decoy"])), p["q"],
                            p.get("pep", ""), p.get("run", run_default), p.get("score", "")])
                n_emit += 1

    print(f"# class-decoy ledger  ({fdr_basis}; unit={a.unit}; alpha={a.alpha})")
    print(f"  GLOBAL: T={gT} D={gD}  FDR_hat=(D+1)/T={(gD + 1) / max(gT, 1):.2%}")
    for r in rows:
        print(f"  {r['class']:12s} T={r['T_class']:6d} D={r['D_class']:5d}  "
              f"class-FDR (D+1)/T = {r['class_fdr_hat']:.2%}   (f={r['f']:.1%})")
    if canon:
        for cls, cc in canon.items():
            print(f"  [non-novelty] {cls}: {cc['canonical_self']}/{cc['n_target_peptides']} "
                  f"canonical-self = {cc['canonical_self_frac']:.1%}")
    if strata:
        print("  [multiplicity] class-FDR by PSM-replication-depth (n PSMs/peptide):")
        for cls in sorted(strata):
            for b in ("n=1", "n=2", "n>=3"):
                c = strata[cls].get(b, {"T": 0, "D": 0})
                print(f"    {cls:12s} {b:6s} T={c['T']:5d} D={c['D']:4d}  "
                      f"class-FDR (D+1)/T = {(c['D'] + 1) / max(c['T'], 1):.2%}")
    msg = f"wrote {a.out}.json + {a.out}.tsv"
    if a.emit_diagfdr:
        msg += f" + {a.emit_diagfdr} ({n_emit} PSMs; diagFDR contract)"
    print(msg)


if __name__ == "__main__":
    main()
