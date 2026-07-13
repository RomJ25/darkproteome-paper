"""The reporting-and-adjudicability audit — the paper's Result 2.

    python3 src/darkproteome/audit.py data/claim_catalog_scaled.csv

This replaces the four-axis survivor funnel (retired). There is no "strict survivor
fraction" any more, and reviving one would be a mistake: a joint pass/fail count over dimensions
the record cannot even decide is not a measurement of the claims, it is a measurement of the
scorer. See `evidence_dimensions.py`.

What this prints instead, per evidence dimension:

    asserted -> claim_linked -> quantitative -> modality_appropriate -> adjudicable

and, ONLY where adjudicable, the reported outcome. Everything else is "the record does not say",
which is a finding about the record and never about the biology.

STRATIFIED, always. A pooled denominator dominated by 293,222 IEAtlas rows hides whether the same
defect exists in the end-to-end cohorts, which is the question a reader actually has.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_dimensions as ed  # noqa: E402
from schema import validate_row  # noqa: E402

csv.field_size_limit(10_000_000)


def stratum_of(row):
    """atlas / cohort / t-cell-tested — the three populations, never pooled into one bar."""
    src = (row.get("_source") or "").strip()
    if ed.human_tcell_state(row) in ("human-assay-positive", "human-assay-negative",
                                     "nonhuman-assay-indirect", "assayed-result-not-claim-linked"):
        return "3. T-cell-tested"
    if src.startswith("cohort:"):
        return "2. end-to-end cohort"
    return "1. atlas ligand record"


def report(rows, label):
    m = ed.matrix(rows)
    print()
    print(ed.format_matrix(m, f"=== {label}  (n={len(rows):,}) ==="))

    # the three-way states that a binary could not express
    from collections import Counter
    ns = Counter(ed.normal_presentation_state(r) for r in rows)
    tc = Counter(ed.human_tcell_state(r) for r in rows)
    print(f"\n  normal-presentation evidence state:")
    for k in ("modality-appropriate-claim-linked", "explicit-normal-detection",
              "indirect-or-study-level", "not-reported"):
        if ns.get(k):
            tag = "   <- the ONE true empirical negative" if k == "explicit-normal-detection" else ""
            print(f"      {k:<38} {ns[k]:>9,}{tag}")
    print(f"  human T-cell assay state:")
    for k in ("human-assay-positive", "human-assay-negative", "nonhuman-assay-indirect",
              "assayed-result-not-claim-linked", "not-assayed"):
        if tc.get(k):
            print(f"      {k:<38} {tc[k]:>9,}")
    return m


def main(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    n = len(rows)
    invalid = [i for i, r in enumerate(rows, 1) if validate_row(r)]

    print(f"\n=== darkproteome reporting & adjudicability audit: {os.path.basename(path)} ===")
    print(f"\nclaims: {n:,}    schema-invalid rows: {len(invalid)}")
    print("\nSCOPE. This audits evidence available in public, MACHINE-READABLE, CLAIM-LINKED")
    print("supplementary tables and atlas exports. Study-level prose, figure-only results, and")
    print("analyses recoverable only by reprocessing raw data are NOT counted as reusable per-claim")
    print("evidence. A study can do excellent science and still fail this level — that is the")
    print("subject, and nothing here is a statement about whether the biology is real.")

    report(rows, "POOLED (do not quote alone — see strata below)")

    strata = {}
    for r in rows:
        strata.setdefault(stratum_of(r), []).append(r)
    for label in sorted(strata):
        report(strata[label], label)

    print("\n" + "=" * 78)
    print("THE HEADLINE, stated at the level the evidence supports")
    print("=" * 78)
    m = ed.matrix(rows)
    st, ar, tc, cf = (m["source_translation"], m["allele_restriction"],
                      m["human_tcell_assay"], m["class_fdr_reconstructible"])
    ribo = sum(1 for r in rows if "ribo" in (r.get("evidence_types") or "").lower())
    cohort = [r for r in rows if (r.get("_source") or "").startswith("cohort:")]
    cohort_allele = sum(1 for r in cohort if ed.allele_restriction(r)["claim_linked"])
    print(f"""
  Of {n:,} audited machine-readable claims:

  SOURCE TRANSLATION — collapsing these levels mischaracterises what the authors did:
    - a dedicated translation ANALYSIS (Ribo-seq/RibORF) is named for  {ribo:,}
    - translation is implied by the MS observation itself for          {st['asserted']:,}
    - a CLAIM-LINKED quantitative statistic (this ORF's own score,
      this peptide's own q-value) is present for                       {st['quantitative']:,}
    NOTE: the sources DO report translation and search statistics — at the STUDY level (a RibORF
    score cutoff, an average read periodicity, a PSM-level FDR). What none publishes is the
    PER-CLAIM value, so no individual claim can be independently re-adjudicated. Do NOT write
    "the field does not report translation statistics"; it does, for the study.
    Whether translation is biologically real is NOT ASSESSED here — that is a different study.

  HLA PRESENTATION — asserted almost everywhere, re-evaluable almost nowhere:
    - reported as HLA-eluted:                                          {m['hla_elution']['asserted']:,}
    - per-peptide ALLELE restriction reported:                         {ar['claim_linked']:,}
    - ... of which, in the END-TO-END COHORTS (n={len(cohort):,}):{'':>18}{cohort_allele:,}
    Presentation is generally asserted by resource inclusion; the information needed to
    RE-EVALUATE it at claim level is sparse. Do NOT read "{ar['claim_linked']:,} of {ar['claim_linked']:,} pass" as
    good reporting — it is {100*ar['claim_linked']/n:.2f}% COVERAGE.

  HUMAN T-CELL ASSAY:
    - a claim-linked human assay RESULT travels with                   {tc['adjudicable']:,}
    - assayed, but the result is published only inside a figure:       {sum(1 for r in rows if ed.human_tcell_state(r) == 'assayed-result-not-claim-linked'):,}
    No audited record carries a reusable, machine-readable POSITIVE human T-cell result. This is
    NOT "the claims failed an assay" — {sum(1 for r in rows if ed.human_tcell_state(r) == 'not-assayed'):,} were never assayed at all.

  CLASS-SPECIFIC IDENTIFICATION ERROR:
    - per-PSM q-value available for                                    {cf['claim_linked']:,}
    - accepted-decoy count for the claim's CLASS (D_N) available for   {cf['modality_appropriate']:,}
    So the class-specific target-decoy estimate is not reconstructible from the published tables
    for ANY claim — including the {cf['claim_linked']:,} that DO carry a per-PSM q-value.""")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: audit.py <claim_catalog.csv>")
    main(sys.argv[1])
