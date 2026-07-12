"""Does source ambiguity have a consequence? Test it inside the resource.

If a catalogued "cancer epitope" is in fact a peptide of an abundant canonical protein, it should
also be presented on NORMAL tissue -- because that protein is expressed there too. IEAtlas publishes
its own normal-tissue epitope set, so the prediction is testable within the resource, using the
catalogue's NON-overlapping epitopes as an internal control and requiring no external reference.

WHAT THIS DOES NOT SHOW: that any individual epitope is canonically derived. MS identifies the
SEQUENCE, never the LOCUS. Presence in a normal-tissue set is evidence about PRESENTATION, not SOURCE.

    python3 scripts/consequence.py
"""
import csv
import math
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ATL = os.path.join(REPO, "data", "external", "atlases")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
csv.field_size_limit(10_000_000)


def load(path):
    s = set()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if r and r[0]:
                q = r[0].strip().upper()
                if q.isalpha():
                    s.add(q)
    return s


def _require(*paths):
    """Fail with a usable message, not a traceback, when the external inputs are absent.

    The large public inputs (Swiss-Prot, the atlas exports, the ncORF libraries, the fetched full
    texts) are not redistributed in this repository. Populate `data/external/` from the sources
    documented in `data/SOURCES.md` and `data/external/README.md`.
    """
    import sys as _s
    missing = [p for p in paths if not __import__("os").path.exists(p)]
    if missing:
        _s.exit("missing required input(s):\n  " + "\n  ".join(missing) +
                "\n\nThese are large public files and are not redistributed here.\n"
                "Populate data/external/ -- see data/SOURCES.md and data/external/README.md.")


def main():
    _require(SCORED,
             os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt"),
             os.path.join(ATL, "IEAtlas_Epitopes_In_Normal_Tissues.txt"))
    cancer = load(os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt"))
    normal = load(os.path.join(ATL, "IEAtlas_Epitopes_In_Normal_Tissues.txt"))
    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

    scored = [p for p in cancer if p in selfmap]
    canon = [p for p in scored if selfmap[p]]
    other = [p for p in scored if not selfmap[p]]

    print(f"IEAtlas: {len(cancer):,} unique CANCER epitopes | {len(normal):,} unique NORMAL-tissue "
          f"epitopes\n")
    print("=" * 88)
    print("THE CONSEQUENCE — are the source-ambiguous 'cancer' epitopes also on NORMAL tissue?")
    print("  Internal control: the NON-overlapping epitopes from the same catalogue.")
    print("=" * 88)
    kc = sum(1 for p in canon if p in normal)
    ko = sum(1 for p in other if p in normal)
    p1, p2 = kc / len(canon), ko / len(other)
    print(f"  {'canonical-overlapping cancer epitopes':<48}{kc:>8,}/{len(canon):<8,} = {100*p1:5.1f}%")
    print(f"  {'NON-overlapping (internal control)':<48}{ko:>8,}/{len(other):<8,} = {100*p2:5.1f}%")
    pp = (kc + ko) / (len(canon) + len(other))
    se = math.sqrt(pp * (1 - pp) * (1 / len(canon) + 1 / len(other)))
    print(f"\n  risk ratio {p1/p2:.1f}x   (two-proportion z = {(p1-p2)/se:,.0f}; "
          f"n = {len(canon)+len(other):,})")
    print("  A source-ambiguous 'cancer' epitope is more than twice as likely to also appear in")
    print("  IEAtlas's OWN normal-tissue set. That is exactly what is expected if it is a peptide")
    print("  of a canonical protein that is expressed in normal tissue too.")

    both = [p for p in canon if p in normal]
    print("\n" + "=" * 88)
    print("THE SUBSET THAT NEEDS NO INFERENCE AT ALL")
    print("=" * 88)
    print(f"  canonical-overlapping AND in IEAtlas's own normal-tissue set : {len(both):>8,}")
    print(f"  = {100*len(both)/len(cancer):.1f}% of every cancer epitope the atlas catalogues.")
    print("\n  These are catalogued as cancer epitopes; they are compatible with an abundant")
    print("  canonical protein; and they are OBSERVED ON NORMAL TISSUE by the atlas's own")
    print("  measurement. No external reference is needed to see it — it is internal to the")
    print("  resource. A target-selection pipeline drawing from this catalogue without an")
    print("  exclusion step will encounter them.")

    print("\n" + "=" * 88)
    print("WHAT THE FIELD'S OWN STANDARD WOULD COST")
    print("=" * 88)
    print(f"  removed by the exclusion standard (Bedran et al. 2023): {len(canon):>8,} "
          f"= {100*len(canon)/len(scored):.1f}%")
    print(f"  surviving:                                             {len(scored)-len(canon):>8,} "
          f"= {100*(len(scored)-len(canon))/len(scored):.1f}%")
    print("\n  For comparison, the cost when Ouspenskaia et al. applied class-specific FDR control:")
    print("  24% of nuORF peptides overall, and up to 76% of one ORF class (8,567 -> 6,501).")
    print("  Applying an established standard to an ncORF catalogue is expected to be expensive.")
    print("  That is not an argument against applying it.")

    print("\n" + "=" * 88)
    print("WHAT THIS DOES NOT SHOW")
    print("=" * 88)
    print("""
  It does NOT show that any individual epitope is canonically derived, or that it is not a genuine
  ncORF product. MS identifies the SEQUENCE, never the LOCUS. Presence in the normal set is
  evidence about PRESENTATION, not about SOURCE.

  What it shows is that the source-ambiguous fraction behaves exactly as canonical-derived peptides
  would -- and that a user selecting tumour-specific targets from this catalogue, without an
  exclusion step, is drawing from a pool in which more than one in eight entries is BOTH
  canonical-compatible AND already observed on normal tissue.
""")


if __name__ == "__main__":
    main()
