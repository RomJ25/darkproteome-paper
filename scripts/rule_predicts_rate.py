"""THE NATURAL EXPERIMENT: the retention rule predicts the canonical-overlap rate.

Three public non-canonical immunopeptidome resources. Same domain, same measurement (exact
substring of a reviewed canonical human protein, unique peptide sequences). They differ in ONE
design decision -- how they decide what counts as "non-canonical" -- and the canonical-overlap rate
tracks that decision across ~2,000x.

  EXCLUSION rule ("drop it if it ALSO maps to canonical")
    CrypticProteinDB  "BLASTP was used to eliminate all proteins with alignment to canonical
                       proteins from GENCODE (v.39) with an E value less than 0.01."
                       + PSM FDR <0.01, PeptideProphet 1%
    Raja ovarian      "peptides mapping to 'protein_coding,' 'IG_C_gene,' and 'IG_V_gene' or HLA
                       transcript biotypes were EXCLUDED to derive the cryptic peptides."
                       + 3% PSM FDR (Percolator)

  RETENTION rule ("keep it if it DERIVES from a non-coding locus")
    IEAtlas           "Files were searched against both our integrated benchmarked ncORF library
                       and the canonical human proteome ... A peptide spectrum match FDR of 0.05
                       was used, and NO PROTEIN FDR WAS SET ... Only epitopes derived from
                       non-coding regions were RETAINED."

A peptide sequence compatible with BOTH a canonical protein and an ncORF is DROPPED by an exclusion
rule and KEPT by a retention rule. The rules are not equivalent, and the difference is measurable.

WHAT THIS DOES AND DOES NOT SHOW
  DOES: the design decision has a large, measurable consequence, and two of the three resources
        already make the safer choice. This is a constructive finding, not an accusation: the field
        demonstrably knows how to do this.
  DOES NOT: establish that any individual IEAtlas peptide is canonically derived. An identical
        amino-acid sequence is ONE sequence compatible with MULTIPLE source loci -- MS identifies
        the sequence, never the locus. What a retention rule does is admit that ambiguity into the
        catalogue WITHOUT RESOLVING IT.

    python3 scripts/rule_predicts_rate.py
"""

import sys as _sys
print("""
================================================================================
  RETRACTED -- DO NOT QUOTE THIS SCRIPT'S OUTPUT
================================================================================
  This script argued that the false-discovery rate could not EXPLAIN the 56.3% canonical overlap,
  on the grounds that a composition-matched SHUFFLE places chance canonical overlap near 0.1%.

  That is the wrong null object. A false target PSM is not an arbitrary shuffled amino-acid string;
  it is an accepted, incorrect candidate drawn from the ACTUAL SEARCH DATABASE. Its class
  composition is not described by a shuffle. The argument is withdrawn.

  The paper needs no such argument. Source ambiguity is present even when the sequence is CORRECTLY
  identified: FDR concerns whether the spectrum was assigned to the right sequence; canonical
  overlap concerns whether that correctly-identified sequence determines a source. Different
  objects. See manuscript_v2.md R5.

  Superseded by: no replacement -- the claim is retracted, not re-derived

  The code below is kept for the record, not for use. It is not run.
================================================================================
""", file=_sys.stderr)
_sys.exit(1)

import csv
import os
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "manuscript"))
csv.field_size_limit(10_000_000)

TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")

RULES = {
    "CrypticProteinDB": ("EXCLUSION", "BLASTP eliminates any protein aligning to GENCODE canonical "
                                      "(E<0.01)", "<1% PSM + 1% PeptideProphet"),
    "Raja ovarian":     ("EXCLUSION", "peptides mapping to protein_coding / IG_C / IG_V / HLA "
                                      "biotypes are EXCLUDED", "3% PSM (Percolator)"),
    "IEAtlas":          ("RETENTION", "only epitopes DERIVED FROM non-coding regions are RETAINED",
                                      "5% PSM, NO protein FDR"),
}


def main():
    for p in (TIER1, SCALED, SCORED):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    # canonical_self per unique peptide, from the scored artifact
    selfmap = {}
    with open(SCORED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            selfmap[r["peptide"]] = int(r["canonical_self"])

    # atlas peptides by source
    seen = defaultdict(set)
    with open(SCALED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            src = (r.get("_source") or "").strip()
            p = (r.get("peptide_sequence") or "").strip().upper()
            if p and p.isalpha() and not src.startswith("cohort:"):
                seen[src].add(p)

    rows = []
    # CrypticProteinDB = the union of its two files (one database, one bar)
    cpdb = seen["CrypticProteinDB-immuno"] | seen["CrypticProteinDB-epitopes"]
    for label, peps in (("CrypticProteinDB", cpdb), ("IEAtlas", seen["IEAtlas"])):
        n = sum(1 for p in peps if p in selfmap)
        k = sum(1 for p in peps if selfmap.get(p))
        rows.append((label, k, n))

    # Raja, from the tier-1 artifact (same measurement)
    raja = [r for r in csv.DictReader(open(TIER1)) if r["cohort"].startswith("Raja")]
    rows.append(("Raja ovarian", sum(int(r["canonical_self_exact"]) for r in raja), len(raja)))

    order = {"CrypticProteinDB": 0, "Raja ovarian": 1, "IEAtlas": 2}
    rows.sort(key=lambda r: order[r[0]])

    print("=" * 96)
    print("THE RULE PREDICTS THE RATE")
    print("  measurement: exact substring of a reviewed canonical human protein (Swiss-Prot),")
    print("  unique peptide sequences. Identical for all three resources.")
    print("=" * 96)
    print(f"  {'resource':<18}{'rule':<11}{'pooled FDR':<26}{'canonical overlap':>22}")
    print("  " + "-" * 92)
    for label, k, n in rows:
        rule, _desc, fdr = RULES[label]
        pct = 100 * k / n if n else 0
        pstr = f"{k:,}/{n:,} = {pct:.3f}%" if pct < 1 else f"{k:,}/{n:,} = {pct:.1f}%"
        print(f"  {label:<18}{rule:<11}{fdr:<26}{pstr:>22}")
    print()
    for label, _k, _n in rows:
        rule, desc, _f = RULES[label]
        print(f"  {label:<18} [{rule}] {desc}")

    lo = min((100 * k / n) for _l, k, n in rows if n)
    hi = max((100 * k / n) for _l, k, n in rows if n)
    print()
    print(f"  spread: {lo:.3f}%  ->  {hi:.1f}%   ({hi/lo:,.0f}x)")
    print()
    print("  A peptide sequence compatible with BOTH a canonical protein and an ncORF is DROPPED by")
    print("  an exclusion rule and KEPT by a retention rule. The two rules are not equivalent, and")
    print("  the difference is measurable at ~2,000x.")
    print()
    print("  CONSTRUCTIVE, NOT ACCUSATORY: two of the three resources already make the safer choice.")
    print("  The field demonstrably knows how to do this. Nothing here shows that any individual")
    print("  IEAtlas peptide is canonically derived -- MS identifies the SEQUENCE, never the LOCUS.")
    print("  What a retention rule does is admit the ambiguity into the catalogue WITHOUT RESOLVING")
    print("  IT.")


def confounds():
    """A referee will attack the comparison on two grounds. Both fail."""
    import re
    PSEUDO = re.compile(r"^[A-Z0-9\-]+P\d+$")
    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

    print("\n" + "=" * 96)
    print("CONFOUND 1 — 'it is the ORF-CLASS COMPOSITION, not the rule'")
    print("=" * 96)
    print("  IEAtlas's library explicitly includes pseudogenes, the highest-overlap class. If the")
    print("  exclusion-rule resources simply lack pseudogenes, the library -- not the rule -- would")
    print("  explain everything. So HOLD THE CLASS FIXED.\n")

    ie = os.path.join(REPO, "data", "external", "atlases",
                      "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
    pseudo, other = set(), set()
    with open(ie, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 3 or not r[0]:
                continue
            seq = r[0].strip().upper()
            if not seq.isalpha():
                continue
            gene = (r[2] or "").split("_")[0].upper()
            (pseudo if PSEUDO.match(gene) else other).add(seq)

    def rate(s):
        n = sum(1 for p in s if p in selfmap)
        k = sum(1 for p in s if selfmap.get(p))
        return k, n, (100 * k / n if n else 0)

    raja_ps = [r for r in csv.DictReader(open(TIER1, newline="", encoding="utf-8"))
               if r["cohort"].startswith("Raja") and r["orf_class"] == "pseudogene-ORF"]
    rk = sum(int(r["canonical_self_exact"]) for r in raja_ps)

    kp, np_, pp = rate(pseudo)
    ko, no_, po = rate(other)
    print(f"  {'PSEUDOGENE-derived peptides':<40}{'canonical overlap':>26}")
    print("  " + "-" * 66)
    print(f"  {'IEAtlas          [RETENTION]':<40}{kp:>10,}/{np_:<8,} = {pp:5.1f}%")
    print(f"  {'Raja ovarian     [EXCLUSION]':<40}{rk:>10,}/{len(raja_ps):<8,} = "
          f"{100*rk/len(raja_ps):5.1f}%")
    print()
    print("  Same ORF class. Opposite rule. 60% vs 0%.")
    print()
    print(f"  And the overlap is NOT a pseudogene artifact -- it is PERVASIVE across IEAtlas:")
    print(f"    IEAtlas, pseudogene-symbol ORFs : {kp:>7,}/{np_:<8,} = {pp:5.1f}%")
    print(f"    IEAtlas, ALL OTHER ORF types    : {ko:>7,}/{no_:<8,} = {po:5.1f}%")
    print("  Non-pseudogene ORFs sit within a few points of the pseudogenes, so class composition")
    print("  explains almost none of the 56.3%. CONFOUND 1 FAILS.")

    print("\n" + "=" * 96)
    print("CONFOUND 2 — 'it is the FDR, not the rule'")
    print("=" * 96)
    print("  The three resources also differ in pooled FDR (1%, 3%, 5%). But a looser FDR admits")
    print("  more FALSE identifications, and a false ID is not preferentially a canonical substring")
    print("  -- a composition-matched shuffle puts chance canonical overlap near 0.1%. No FDR can")
    print("  manufacture a 56% canonical-substring rate.")
    print()
    print("  And the two exclusion-rule resources differ in FDR by 3x (1% -> 3%) while their overlap")
    print("  moves 0.026% -> 0.168% (~6x). The RULE moves it to 56.3% (~2,000x). The rule dominates")
    print("  the FDR by two to three orders of magnitude. CONFOUND 2 FAILS.")


if __name__ == "__main__":
    main()
    confounds()
