"""The cross-catalogue comparison, restricted to what we can actually control.

WHAT WAS WRONG (external review, upheld). The draft compared IEAtlas's 56.3% against the 1.4-5%
range published by Bedran et al. for four other catalogues, and reported the ratio as "11-40x" --
including in the title ("an order of magnitude"). That is NOT a controlled comparison. Those four
rates were produced with a different reference, normalization, deduplication and peptide unit. A
fold-change across pipelines is arithmetic, not measurement.

Two honest options: remeasure all four under our pipeline, or stop doing fold arithmetic. We do not
hold the four catalogues' peptide lists, so we take the second -- and salvage the comparison that IS
controlled: the catalogues we reprocessed ourselves, end to end, against the same reference, same
exact-substring criterion, same unique-sequence unit.

Also per review: exact canonical-substring probability is strongly length-dependent, so a resource
rich in 8-mers is not comparable to one dominated by 10-11mers. Every rate below is therefore
reported BOTH crude AND directly standardized to a common length distribution (IEAtlas's), so the
contrast cannot be a composition artifact.

    python3 scripts/cross_catalogue.py
"""
import csv
import os
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
csv.field_size_limit(10_000_000)


def main():
    for p in (SCALED, SCORED, TIER1):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

    cats = defaultdict(set)
    for r in csv.DictReader(open(SCALED, newline="", encoding="utf-8")):
        s = (r.get("_source") or "").strip()
        p = (r.get("peptide_sequence") or "").strip().upper()
        if p and p.isalpha() and not s.startswith("cohort:"):
            cats[s].add(p)
    crypt = cats["CrypticProteinDB-immuno"] | cats["CrypticProteinDB-epitopes"]

    groups = {"IEAtlas": cats["IEAtlas"], "CrypticProteinDB": crypt}

    # Raja lives in the non-novelty floor table, with its own canonical_self_exact column.
    raja = [r for r in csv.DictReader(open(TIER1, newline="", encoding="utf-8"))
            if r["cohort"].startswith("Raja")]

    # Reference length distribution = IEAtlas's (the resource under audit).
    ref = defaultdict(int)
    for p in groups["IEAtlas"]:
        if p in selfmap:
            ref[len(p)] += 1
    tot = sum(ref.values())
    weights = {L: c / tot for L, c in ref.items()}

    print("SAME PIPELINE: same reference, same exact-substring criterion, same unique-sequence unit.")
    print("(Bedran et al.'s 1.4-5% values are NOT here. They were produced by a different pipeline and")
    print(" are reported in the paper as CONTEXT ONLY -- no fold arithmetic against them.)\n")

    print(f"  {'catalogue':<20}{'n':>10}{'crude':>10}{'length-standardized':>22}")
    print("  " + "-" * 62)

    rows = []
    for name, peps in groups.items():
        sc = [p for p in peps if p in selfmap]
        if not sc:
            continue
        k = sum(selfmap[p] for p in sc)
        crude = 100 * k / len(sc)
        num = defaultdict(lambda: [0, 0])
        for p in sc:
            num[len(p)][0] += selfmap[p]
            num[len(p)][1] += 1
        std, wsum = 0.0, 0.0
        for L, w in weights.items():
            kk, nn = num[L]
            if nn >= 20:
                std += w * kk / nn
                wsum += w
        std = 100 * std / wsum if wsum else float("nan")
        rows.append((name, len(sc), crude, std, wsum, k))
        # With a handful of events, a standardized rate is falsely precise. Say so instead of
        # printing 0.000%.
        cell = f"{std:.3f}%" if k >= 20 else f"n/a ({k} event{'s' if k != 1 else ''})"
        print(f"  {name:<20}{len(sc):>10,}{crude:>9.3f}%{cell:>22}")

    rk = sum(int(r["canonical_self_exact"]) for r in raja)
    print(f"  {'Raja (cohort)':<20}{len(raja):>10,}{100*rk/len(raja):>9.3f}%{'n/a (no len col)':>22}")

    print("\n  coverage of the standardization: "
          + ", ".join(f"{n} {100*w:.0f}%" for n, _, _, _, w, _ in rows))

    print("\n" + "=" * 84)
    print("WHAT THIS LICENSES")
    print("=" * 84)
    ie = next(r for r in rows if r[0] == "IEAtlas")
    cp = next(r for r in rows if r[0] == "CrypticProteinDB")
    print(f"""
  Under one pipeline, IEAtlas is at {ie[2]:.1f}% and CrypticProteinDB at {cp[2]:.3f}%
  ({cp[5]} overlapping sequence of {cp[1]:,}). IEAtlas's rate is not a length-composition artifact:
  standardizing to a common length distribution leaves it unchanged ({ie[2]:.1f}% -> {ie[3]:.1f}%).
  CrypticProteinDB's rate is NOT standardized -- with {cp[5]} event it cannot be, and a standardized
  figure there would be falsely precise. The contrast rests on the raw counts, which is enough:
  {ie[5]:,}/{ie[1]:,} against {cp[5]}/{cp[1]:,}.

  SAY: "Among the catalogues we reprocessed under a single pipeline, IEAtlas's canonical-sequence
  overlap is orders of magnitude higher than CrypticProteinDB's, and IEAtlas's rate is not explained
  by peptide-length composition."

  DO NOT SAY: "11-40x higher than other catalogues." Those four rates come from a different
  pipeline; the ratio is not measured. Cite them as published context and leave the arithmetic out.

  The title must not claim a fold-change we did not measure.
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
