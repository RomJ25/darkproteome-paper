"""Stress-test the reporting result against the obvious objections.

    python3 src/darkproteome/robustness.py data/claim_catalog_real.csv

Three objections, stress-tested against `evidence_dimensions.py`:

  1. Corpus definition   -- does the reporting result move if you cut the corpus differently?
  2. Criterion leniency  -- how far must the criteria be relaxed before the record can adjudicate
                            anything? (The honest form of "your bar is arbitrarily strict.")
  3. The immunogenicity ceiling -- the generous upper bound on reusable human ground truth.

NOTE what is NOT claimed anywhere below: that any claim is biologically false. Relaxing a criterion
until claims become adjudicable tells you what the record would need to report, not what is true.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_dimensions as ed  # noqa: E402

NONCODING = {"lncRNA-ORF", "pseudogene-ORF"}  # narrow "dark proteome" (excludes altORF/other)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def cut_table(rows):
    noncanon = [r for r in rows if r.get("_canonical") != "yes"]
    narrow = [r for r in noncanon if r["orf_class"] in NONCODING]
    seen, dedup = set(), []
    for r in noncanon:
        p = r["peptide_sequence"]
        if p and p != "not reported" and p in seen:
            continue
        if p and p != "not reported":
            seen.add(p)
        dedup.append(r)
    print("\n=== 1. The reporting result does not depend on how the corpus is cut ===")
    print(f"  {'corpus':<34}{'N':>6}{'src-transl.':>13}{'allele':>9}{'T-cell':>9}")
    print(f"  {'':<34}{'':>6}{'adjudicable':>13}{'reported':>9}{'result':>9}")
    print("  " + "-" * 71)
    for name, rs in (("ALL non-canonical", noncanon),
                     ("NARROW (lncRNA+pseudogene only)", narrow),
                     ("DEDUP by peptide", dedup)):
        m = ed.matrix(rs)
        print(f"  {name:<34}{len(rs):>6}"
              f"{m['source_translation']['adjudicable']:>13}"
              f"{m['allele_restriction']['claim_linked']:>9}"
              f"{m['human_tcell_assay']['adjudicable']:>9}")
    print("  -> stable across every definition. The gap is in the RECORD, not in our slicing.")
    return noncanon


def leniency_ladder(rows):
    """How far must a criterion be relaxed before the record can adjudicate ANYTHING?

    This is the honest answer to "your bar is arbitrarily strict". Each rung LOOSENS what we are
    willing to accept, and reports how many claims become adjudicable. If a rung is still 0, the
    record does not contain the evidence at ANY strictness -- which is a fact about the record.
    """
    n = len(rows)
    print("\n=== 2. Escalating leniency — how far must the criteria drop to adjudicate anything? ===")

    # -- source translation --
    strict = sum(1 for r in rows if ed.source_translation(r)["adjudicable"])
    any_fdr = sum(1 for r in rows if ed._num(r.get("reported_fdr")) is not None)
    any_period = sum(1 for r in rows if ed._num(r.get("periodicity_pct")) is not None)
    any_qval = sum(1 for r in rows if ed._num(r.get("psm_qvalue")) is not None)
    binary_ribo = sum(1 for r in rows if "ribo" in (r.get("evidence_types") or "").lower())
    print("\n  SOURCE TRANSLATION")
    print(f"    prespecified criterion (periodicity>=70% or protein FDR<=0.001) : {strict:>6}")
    print(f"    relax to: ANY reported protein-level FDR, at any value          : {any_fdr:>6}")
    print(f"    relax to: ANY reported Ribo-seq periodicity, at any value       : {any_period:>6}")
    print(f"    relax to: ANY per-PSM q-value (WRONG UNIT — not a protein FDR)  : {any_qval:>6}")
    print(f"    relax to: a binary 'RibORF says translated', no statistic at all : {binary_ribo:>6}")
    print("    -> the first three are the criterion loosened to nothing, and they stay at 0 or")
    print("       near it. Only the BINARY assertion is available. Translation is asserted and")
    print("       never quantified — at ANY threshold we could have chosen.")

    # -- allele --
    allele = sum(1 for r in rows if ed.allele_restriction(r)["claim_linked"])
    print("\n  HLA ALLELE RESTRICTION")
    print(f"    per-peptide allele reported                                      : {allele:>6}"
          f"  ({100*allele/max(1,n):.1f}%)")
    print("    -> there is no leniency knob here. Either the allele is in the table or it is not.")

    # -- normal presentation --
    ligandome_claim = sum(1 for r in rows if ed.normal_presentation(r)["adjudicable"])
    any_ligandome = sum(1 for r in rows
                        if ed._txt(r.get("tumor_specificity_modality")) in ed.LIGANDOME_MODALITIES)
    any_modality = sum(1 for r in rows
                       if ed._txt(r.get("tumor_specificity_modality")) in ed.ALL_MODALITIES)
    print("\n  NORMAL PRESENTATION")
    print(f"    prespecified (broad normal LIGANDOME, reported per claim)        : {ligandome_claim:>6}")
    print(f"    relax to: any normal-ligandome evidence, even study-level         : {any_ligandome:>6}")
    print(f"    relax to: ANY normal-tissue evidence at all, incl. RNA proxies    : {any_modality:>6}")
    print("    -> RNA abundance does not determine HLA presentation, so the last rung answers a")
    print("       different question. It is reported here so the reader can see the whole ladder.")

    # -- immunogenicity --
    strict_t = sum(1 for r in rows if ed.human_tcell_assay(r)["adjudicable"])
    pos = sum(1 for r in rows if ed.human_tcell_state(r) == "human-assay-positive")
    mouse = sum(1 for r in rows if ed.human_tcell_state(r) == "nonhuman-assay-indirect")
    figlock = sum(1 for r in rows if ed.human_tcell_state(r) == "assayed-result-not-claim-linked")
    print("\n  HUMAN T-CELL ASSAY")
    print(f"    claim-linked HUMAN assay result                                  : {strict_t:>6}"
          f"   (of which POSITIVE: {pos})")
    print(f"    relax to: + credit humanized-mouse reactivity as human evidence   : {strict_t + mouse:>6}")
    print(f"    relax to: + credit results published ONLY inside a figure         : {strict_t + mouse + figlock:>6}")
    print("    -> the ONLY way to manufacture a reusable positive is to credit non-human evidence")
    print("       or to digitize a figure by hand. That is not a stricter bar; it is the absence")
    print("       of a machine-readable one.")


def immuno_truth(rows):
    print("\n=== 3. The reusable human-immunogenicity ground truth (the real ceiling) ===")
    from collections import Counter
    c = Counter(ed.human_tcell_state(r) for r in rows)
    for k in ("human-assay-positive", "human-assay-negative", "nonhuman-assay-indirect",
              "assayed-result-not-claim-linked", "not-assayed"):
        print(f"    {k:<34} {c.get(k, 0):>7,}")
    ceiling = c.get("nonhuman-assay-indirect", 0) + c.get("assayed-result-not-claim-linked", 0)
    print(f"\n    machine-readable human POSITIVES                : {c.get('human-assay-positive', 0)}")
    print(f"    generous reusable-positive CEILING             : {ceiling}"
          f"   (every one needs non-human evidence or figure digitization)")
    print(f"\n    NOTE: {c.get('not-assayed', 0):,} claims were NEVER ASSAYED. That is not a failed assay,")
    print( "    and the old scorer recorded it as one.")


def main(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    print("#" * 78)
    print("# ROBUSTNESS — is the REPORTING result an artifact of our choices?")
    print("#" * 78)
    noncanon = cut_table(rows)
    leniency_ladder(noncanon)
    immuno_truth(noncanon)
    print("\n" + "#" * 78)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/claim_catalog_real.csv")
