"""Reference model / evidence-control layer -- the capstone of the as-reported audit.

Turns the 306k as-reported claim catalog into the reusable referee artifact:
  (1) a per-peptide SCORED TABLE  -> data/claim_catalog_scored.csv
  (2) the EVIDENCE-SURVIVORSHIP CURVE printed for non-canonical claims AND the canonical
      cancer-testis control (MAGE/SSX), run through the SAME filters.

Scores strictly AS-REPORTED + our computed provenance/specificity features (canonical-self,
normal-overlap, critical-organ). It does NOT recompute FDR from spectra -- that would be a
separate re-analysis of raw spectra (the natural next study; not required for these results,
which stand on the published record and the identifiability theory -- Methods).
The honest ceiling: this completes the *reporting* referee; the immunogenicity cliff at the end
is the point, not a bug.

    python3 src/darkproteome/reference_model.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- centralized data paths
import deepen_specificity as d  # canonical_hits, load, load_hla_ligand_atlas, tier_of, NORMAL

csv.field_size_limit(10_000_000)
CAT = os.path.join(paths.REPO, "data", "claim_catalog_scaled.csv")
OUT = os.path.join(paths.REPO, "data", "claim_catalog_scored.csv")

NR = {"not reported", "not stated", "", "insufficient-info"}
REUSE_FIELDS = ["hla_allele", "reported_fdr", "periodicity_pct", "tumor_specificity_basis",
                "orf_id_or_locus", "underlying_dataset_accession", "validation_level",
                "citation_doi_pmid"]
TCELL = {"T-cell-validated", "in-vivo"}


def present(v):
    return (v or "").strip().lower() not in NR


def aggregate():
    nc, ctrl = {}, {}
    with open(CAT, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            seq = (r.get("peptide_sequence") or "").strip().upper()
            if not seq or not seq.isalpha():
                continue
            tbl = ctrl if r.get("_canonical") == "yes" else nc
            rec = tbl.setdefault(seq, {"allele": False, "source_rigor": False,
                                       "specificity": False, "immuno": False, "reuse": 0})
            if present(r.get("hla_allele")):
                rec["allele"] = True
            if present(r.get("periodicity_pct")) or present(r.get("reported_fdr")):
                rec["source_rigor"] = True
            if present(r.get("tumor_specificity_basis")):
                rec["specificity"] = True
            if (r.get("validation_level") or "").strip() in TCELL:
                rec["immuno"] = True
            rec["reuse"] = max(rec["reuse"], sum(present(r.get(f)) for f in REUSE_FIELDS))
    return nc, ctrl


def survivor(seqs, flags, canon, normal_ref):
    not_canon = {p for p in seqs if p not in canon}
    restricted = {p for p in not_canon if p not in normal_ref}
    allele = {p for p in restricted if flags[p]["allele"]}
    subst = {p for p in allele if flags[p]["source_rigor"]}
    immuno = {p for p in subst if flags[p]["immuno"]}
    return [("all unique peptides", len(seqs)),
            ("not canonical-self (provenance-clean)", len(not_canon)),
            ("absent from normal (tumor-restricted)", len(restricted)),
            ("allele-assigned (actionable presentation)", len(allele)),
            ("source-translation substantiated", len(subst)),
            ("human T-cell validated (antigen)", len(immuno))]


def tier_of_pep(p, flags, canon, normal_ref):
    if p in canon:
        return "0_canonical-self"
    if p in normal_ref:
        return "1_normal-present"
    if not flags[p]["allele"]:
        return "2_provenance-clean"
    if not flags[p]["source_rigor"]:
        return "3_actionable-presentation"
    if not flags[p]["immuno"]:
        return "4_substantiated"
    return "5_validated-antigen"


def curve(title, rows):
    print(f"\n{title}")
    base = rows[0][1] or 1
    for label, n in rows:
        bar = "#" * round(40 * n / base)
        print(f"  {n:>8,}  {100*n/base:5.1f}%  {bar}  {label}")


def main():
    print("loading normal references ...")
    normal_d = d.load(d.NORMAL)                  # IEAtlas normal (cryptic-space)
    normal_seqs = set(normal_d)
    hla = d.load_hla_ligand_atlas()              # canonical-space normal
    hla_seqs = set(hla)

    print("aggregating catalog per unique peptide ...")
    nc, ctrl = aggregate()
    nc_seqs, ctrl_seqs = set(nc), set(ctrl)
    print(f"  non-canonical unique peptides: {len(nc_seqs):,}; canonical-control: {len(ctrl_seqs):,}")

    print("computing canonical-self (SwissProt substring) for both sets ...")
    canon_nc = d.canonical_hits(nc_seqs)
    canon_ctrl = d.canonical_hits(ctrl_seqs)

    # non-canonical normal ref = IEAtlas cryptic-space (HLA Ligand Atlas can't contain cryptic peptides);
    # canonical-control normal ref = both (canonical peptides CAN appear in HLA Ligand Atlas).
    curve("EVIDENCE-SURVIVORSHIP — non-canonical claims (n=%d unique peptides)" % len(nc_seqs),
          survivor(nc_seqs, nc, canon_nc, normal_seqs))
    curve("CONTROL — canonical cancer-testis antigens MAGE/SSX (same filters)",
          survivor(ctrl_seqs, ctrl, canon_ctrl, normal_seqs | hla_seqs))

    # reusability score distribution (0..8 fields present)
    import collections
    dist = collections.Counter(nc[p]["reuse"] for p in nc_seqs)
    print("\nREUSABILITY score (per-peptide, of 8 key fields present):")
    for k in range(9):
        if dist.get(k):
            print(f"  {k}/8 fields: {dist[k]:>8,} ({100*dist[k]/len(nc_seqs):4.1f}%)")
    mean_reuse = sum(nc[p]["reuse"] for p in nc_seqs) / max(1, len(nc_seqs))
    print(f"  mean reusability: {mean_reuse:.2f}/8")

    # write the scored table (the reusable artifact)
    print(f"\nwriting scored per-peptide table -> {OUT}")
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["peptide", "canonical_self", "normal_present", "normal_critical_organ",
                    "allele_assigned", "source_rigor_reported", "specificity_reported",
                    "immunogenicity_evidence", "reusability_score_8", "tier"])
        for p in sorted(nc_seqs):
            norm = p in normal_seqs
            crit = norm and any(d.tier_of(t) == "CRITICAL" for t in normal_d.get(p, ()))
            w.writerow([p, int(p in canon_nc), int(norm), int(crit),
                        int(nc[p]["allele"]), int(nc[p]["source_rigor"]),
                        int(nc[p]["specificity"]), int(nc[p]["immuno"]),
                        nc[p]["reuse"], tier_of_pep(p, nc, canon_nc, normal_seqs)])
    print("done.")


if __name__ == "__main__":
    main()
