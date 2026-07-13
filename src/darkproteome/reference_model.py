"""Reference model / evidence-control layer -- the capstone of the as-reported audit.

Turns the 306k as-reported claim catalog into the reusable audit artifact:
  (1) a per-peptide SCORED TABLE  -> data/claim_catalog_scored.csv
  (2) the EVIDENCE-SURVIVORSHIP CURVE printed for non-canonical claims AND the canonical
      cancer-testis control (MAGE/SSX), run through the SAME filters.

Scores strictly AS-REPORTED + our computed provenance/specificity features (canonical-self,
normal-overlap, critical-organ). It does NOT recompute FDR from spectra -- that would be a
separate re-analysis of raw spectra (the natural next study; not required for these results,
which stand on the published record and the identifiability theory -- Methods).
The honest ceiling: this completes the *reporting* audit; the immunogenicity cliff at the end
is the point, not a bug.

    python3 src/darkproteome/reference_model.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- centralized data paths
import evidence_dimensions as ed  # noqa: E402  -- THE scorer; this module must not score
import deepen_specificity as d  # canonical_hits, load, load_hla_ligand_atlas, tier_of, NORMAL

csv.field_size_limit(10_000_000)
CAT = os.path.join(paths.REPO, "data", "claim_catalog_scaled.csv")
OUT = os.path.join(paths.REPO, "data", "claim_catalog_scored.csv")

NR = {"not reported", "not stated", "", "insufficient-info"}
REUSE_FIELDS = ["hla_allele", "reported_fdr", "periodicity_pct", "tumor_specificity_modality",
                "orf_id_or_locus", "underlying_dataset_accession", "validation_level",
                "citation_doi_pmid"]


def present(v):
    return (v or "").strip().lower() not in NR


# REPORTING-COMPLETENESS vs EVIDENCE-STRENGTH. Until this module scored the first and
# CALLED it the second: `source_rigor` was `present(periodicity_pct) or present(reported_fdr)` --
# a test that the FIELD IS POPULATED -- and the survivorship curve then labelled that stage
# "source-translation substantiated". An FDR of 1.0 counted exactly like an FDR of 0.001. That is
# the very error this paper indicts, committed by this paper's own code, and it was a THIRD
# scorer disagreeing with the authoritative one and with `consensus_bar.py`.
#
# Both quantities are legitimate and we keep both -- but they are now named apart and never
# substituted for one another:
#   *_reported  -- is the field populated? REPORTING COMPLETENESS. NOT the reviewer's `A`
#                  estimand: `A_i` is source-attribution RESOLUTION (does the record leave only
#                  the nominated source compatible?), a different object entirely. Conflating the
#                  two was the error in the first replacement table -- do not reintroduce it.
#   *_strict    -- does it clear the bar? Delegated ENTIRELY to evidence_dimensions.py.

def aggregate():
    nc, ctrl = {}, {}
    with open(CAT, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            seq = (r.get("peptide_sequence") or "").strip().upper()
            if not seq or not seq.isalpha():
                continue
            dims = ed.score_all(r)               # THE authoritative scorer
            tbl = ctrl if r.get("_canonical") == "yes" else nc
            rec = tbl.setdefault(seq, {
                "allele": False, "source_rigor_reported": False,
                "specificity_reported": False, "immuno_reported": False, "reuse": 0,
                "source_strict": False, "hla_strict": False,
                "specificity_strict": False, "immuno_strict": False,
                "source_decidable": False, "specificity_decidable": False,
            })
            # -- reporting completeness (a claim about the RECORD, never about the evidence) --
            if present(r.get("hla_allele")):
                rec["allele"] = True
            if present(r.get("periodicity_pct")) or present(r.get("reported_fdr")):
                rec["source_rigor_reported"] = True
            if present(r.get("tumor_specificity_modality")):
                rec["specificity_reported"] = True
            if present(r.get("validation_level")):
                rec["immuno_reported"] = True
            # -- evidence strength (evidence_dimensions.py, and only it) --
            if dims["source_translation"]["outcome"] == ed.SUPPORTS:
                rec["source_strict"] = True
            if dims["hla_elution"]["outcome"] == ed.SUPPORTS:
                rec["hla_strict"] = True
            if dims["normal_presentation"]["adjudicable"]:
                rec["specificity_strict"] = True
            if dims["human_tcell_assay"]["outcome"] == ed.SUPPORTS:
                rec["immuno_strict"] = True
            if dims["source_translation"]["adjudicable"]:
                rec["source_decidable"] = True
            if dims["normal_presentation"]["adjudicable"]:
                rec["specificity_decidable"] = True
            rec["reuse"] = max(rec["reuse"], sum(present(r.get(f)) for f in REUSE_FIELDS))
    return nc, ctrl


def survivor_reporting(seqs, flags, canon, normal_ref):
    """What the RECORD contains. Every stage is 'is this reported?', never 'is this true?'."""
    not_canon = {p for p in seqs if p not in canon}
    restricted = {p for p in not_canon if p not in normal_ref}
    allele = {p for p in restricted if flags[p]["allele"]}
    subst = {p for p in allele if flags[p]["source_rigor_reported"]}
    immuno = {p for p in subst if flags[p]["immuno_reported"]}
    return [("all unique peptides", len(seqs)),
            ("no canonical-sequence match (provenance-clean)", len(not_canon)),
            ("not seen in a normal ligandome", len(restricted)),
            ("an HLA allele is REPORTED", len(allele)),
            ("a translation statistic is REPORTED (any value)", len(subst)),
            ("an immunogenicity level is REPORTED (any level)", len(immuno))]


def survivor_evidence(seqs, flags, canon, normal_ref):
    """What the record SUBSTANTIATES: adjudicable AND outcome=supports, per evidence_dimensions."""
    not_canon = {p for p in seqs if p not in canon}
    restricted = {p for p in not_canon if p not in normal_ref}
    hla = {p for p in restricted if flags[p]["hla_strict"]}
    subst = {p for p in hla if flags[p]["source_strict"]}
    immuno = {p for p in subst if flags[p]["immuno_strict"]}
    return [("all unique peptides", len(seqs)),
            ("no canonical-sequence match (provenance-clean)", len(not_canon)),
            ("not seen in a normal ligandome", len(restricted)),
            ("HLA elution: adjudicable + supports", len(hla)),
            ("source translation: adjudicable + supports", len(subst)),
            ("human T-cell assay: adjudicable + POSITIVE", len(immuno))]


def tier_of_pep(p, flags, canon, normal_ref):
    if p in canon:
        return "0_canonical-self"
    if p in normal_ref:
        return "1_normal-present"
    if not flags[p]["allele"]:
        return "2_provenance-clean"
    if not flags[p]["source_rigor_reported"]:
        return "3_actionable-presentation"
    if not flags[p]["immuno_reported"]:
        return "4_translation-stat-reported"   # was "4_substantiated" -- it never was
    return "5_immunogenicity-reported"         # was "5_validated-antigen" -- nor was it


def curve(title, rows):
    print(f"\n{title}")
    base = rows[0][1] or 1
    for label, n in rows:
        bar = "#" * round(40 * n / base)
        print(f"  {n:>8,}  {100*n/base:5.1f}%  {bar}  {label}")


def main():
    paths.require(paths.SPROT, paths.IEATLAS_NORMAL, paths.HLALA)
    if not os.path.exists(CAT):
        sys.exit(f"need {CAT} (regenerate via ingest_atlases.py)")
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
    #
    # TWO curves, because they answer two different questions and were previously conflated.
    # The reporting curve is about the RECORD; the evidence curve is about the CLAIMS. Quoting
    # a number off the first while describing the second is what produced "0 of 306,844".
    curve("REPORTING CASCADE — what the record CONTAINS (non-canonical, n=%d unique peptides)"
          % len(nc_seqs), survivor_reporting(nc_seqs, nc, canon_nc, normal_seqs))
    curve("CONTROL, reporting cascade — curated cancer-testis antigens (same filters)",
          survivor_reporting(ctrl_seqs, ctrl, canon_ctrl, normal_seqs | hla_seqs))

    curve("EVIDENCE CASCADE — what the record SUBSTANTIATES (adjudicable + supports only)",
          survivor_evidence(nc_seqs, nc, canon_nc, normal_seqs))
    curve("CONTROL, evidence cascade — curated cancer-testis antigens (same filters)",
          survivor_evidence(ctrl_seqs, ctrl, canon_ctrl, normal_seqs | hla_seqs))

    # A zero in the evidence cascade means nothing unless the axis was DECIDABLE. Say so, loudly.
    n_src_dec = sum(1 for p in nc_seqs if nc[p]["source_decidable"])
    n_spec_dec = sum(1 for p in nc_seqs if nc[p]["specificity_decidable"])
    # ADJUDICABILITY -- deliberately NOT called `A`. `A_i` is source-attribution RESOLUTION
    # (does the record leave only the nominated source compatible?); this is whether the record
    # carries enough claim-linked information to apply a reporting criterion at all. Conflating
    # the two was the error in the first replacement table.
    print("\nADJUDICABILITY of the evidence cascade (NOT the `A` estimand -- see evidence_dimensions):")
    print(f"  source translation  adjudicable for {n_src_dec:,}/{len(nc_seqs):,} peptides"
          + ("   <<< its zero says nothing about the claims" if not n_src_dec else ""))
    print(f"  normal presentation adjudicable for {n_spec_dec:,}/{len(nc_seqs):,} peptides"
          + ("   <<< likewise" if not n_spec_dec else ""))
    if not n_src_dec:
        print("  -> Do NOT read the source-translation drop as attrition. No claim in this corpus")
        print("     reports the statistic the criterion needs, so none COULD satisfy it. That is a")
        print("     measurement of the REUSABLE RECORD, not a biological success rate.")

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
        # `*_reported` = the field is populated (reusability / `A`). `*_strict` = it clears the
        # evidence_dimensions bar. Different columns because they are different things.
        w.writerow(["peptide", "canonical_self", "normal_present", "normal_critical_organ",
                    "allele_reported", "source_rigor_reported", "specificity_reported",
                    "immunogenicity_reported",
                    "source_orf_strict", "hla_presentation_strict",
                    "tumor_specificity_strict", "immunogenicity_strict",
                    "source_orf_decidable", "tumor_specificity_decidable",
                    "reusability_score_8", "tier"])
        for p in sorted(nc_seqs):
            norm = p in normal_seqs
            crit = norm and any(d.tier_of(t) == "CRITICAL" for t in normal_d.get(p, ()))
            f = nc[p]
            w.writerow([p, int(p in canon_nc), int(norm), int(crit),
                        int(f["allele"]), int(f["source_rigor_reported"]),
                        int(f["specificity_reported"]), int(f["immuno_reported"]),
                        int(f["source_strict"]), int(f["hla_strict"]),
                        int(f["specificity_strict"]), int(f["immuno_strict"]),
                        int(f["source_decidable"]), int(f["specificity_decidable"]),
                        f["reuse"], tier_of_pep(p, nc, canon_nc, normal_seqs)])
    print("done.")


if __name__ == "__main__":
    main()
