"""Measure the UNION of IEAtlas's ncORF search libraries -- and bound what the missing one can do.

WHAT THIS FIXES. library_ambiguity.py measured nuORFdb v1.2 alone: 34.1% of its distinct 9-mers are
also 9-mers of the reviewed canonical human proteome. We then claimed 34.1% was a LOWER BOUND on the
ambiguity of IEAtlas's COMBINED search library (RPFdb v2.0 + nuORFdb v1.2 + TransLnc). A reviewer
destroyed that claim and was RIGHT:

    |(A u B) n C| / |A u B|  is NOT monotone in adding B.

If B contributes mostly non-canonical k-mers, the union rate FALLS below A's rate. A rate over a union
is not bounded below by the rate over one member. The claim was retracted; the paper says the combined
library's ambiguity is UNKNOWN.

This script closes as much of that gap as is honestly possible, in two moves.

  MOVE 1 -- MEASURE THE UNION DIRECTLY, for the two of the three sources we actually hold, in the exact
  versions IEAtlas used:
      (a) nuORFdb v1.2 alone      -- reproduces 34.1%; this VALIDATES the pipeline. If it does not
                                     reproduce, the pipeline is wrong and nothing below is usable.
      (b) TransLnc alone          -- new.
      (c) nuORFdb u TransLnc      -- new. THE KEY NUMBER: 2 of IEAtlas's 3 sources.
  This is the EMPIRICAL test of the reviewer's monotonicity objection. Adding a real second library
  either raises or lowers the rate. We report the direction plainly whichever way it comes out. If it
  FALLS, that CONFIRMS the reviewer with data, and we say so.

  MOVE 2 -- BOUND THE THIRD SOURCE INSTEAD OF GUESSING IT. Let U = nuORFdb u TransLnc, with
  u = |U| distinct 9-mers and h = |U n C| canonical ones. Let RPFdb contribute m NOVEL 9-mers (not
  already in U), of which x are canonical. Then

      combined rate  =  (h + x) / (u + m),     0 <= x <= min(m, |C| - h)

  x <= m is trivial. x <= |C| - h is NOT: a novel canonical 9-mer must be drawn from the canonical
  9-mer universe C, and h of those are already in U. That second constraint is what makes the upper
  bound SHARP -- and it yields a hard, distribution-free CEILING on the full three-source library's
  ambiguity that holds no matter what RPFdb contains, obtained by maximising over m as well as x.

  There is NO corresponding positive floor: h/(u+m) -> 0 as m grows without bound. That absence IS the
  reviewer's objection, stated quantitatively. We do not have a lower bound and we do not claim one.

WHY RPFdb IS NOT HERE (and is not faked). RPFdb v2.0 is genuinely unavailable in usable form:
  * the live RPFdb site now serves v3.0 only; IEAtlas integrated v2.0;
  * RPFdb distributes RibORF output -- genomic COORDINATES, not amino-acid sequences.
Reconstructing an amino-acid library from it would need genome extraction + frame/start-codon
translation choices with many free parameters, and the result would be OUR library, not IEAtlas's.
Approximating it and calling it RPFdb would manufacture the very number this script exists to bound.
So m is left as the single free unknown, and the answer is reported as an interval in m.

REFERENCES USED. Primary: the modern reviewed human proteome (swissprot_human.fasta), so the number is
directly comparable to library_ambiguity.py's 34.1%. Secondary: SwissProt 2022_01 human (OX=9606), the
release IEAtlas actually searched -- the era-correct arm, per era_correct_reference.py.

CONVENTIONS (identical to library_ambiguity.py, deliberately):
  * K = 9, a representative HLA-I ligand length.
  * The FULL library is measured. NO SAMPLING. A sampled estimate is biased UPWARD here and the bias is
    large: sampling 4,000 of nuORFdb's 229,251 ORFs gives 43.6%, 20,000 gives 40.7%, the full library
    gives 34.1%. Small samples contain fewer distinct ncORF-specific k-mers, so the canonical-shared
    ones are over-weighted. Never sample this measurement.
  * A k-mer counts as canonical iff it occurs as an exact substring of some reviewed human protein.
    k-mers are not allowed to span two proteins.

    python3 scripts/library_union.py
"""
import gzip
import json
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")

SPROT_NOW = os.path.join(EXT, "swissprot_human.fasta")            # modern reviewed human proteome
SPROT_2022 = os.path.join(EXT, "uniprot_sprot.fasta.gz")          # SwissProt 2022_01 (IEAtlas era)

NUORFDB = os.path.join(EXT, "nuorfdb", "PA_nuORFdb_v1.2_protein.fasta")
TRANSLNC = os.path.join(EXT, "translnc", "lncRNA_peptide_AA_seq.fasta")

ART = os.path.join(REPO, "data", "derived_library_union.json")

K = 9
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# The one source we do NOT hold, and why. Kept as data so the artifact records it too.
RPFDB_MISSING = (
    "RPFdb v2.0 -- the third source IEAtlas integrates -- is NOT measured here and is NOT approximated. "
    "The live RPFdb site serves v3.0 only (IEAtlas used v2.0), and RPFdb distributes RibORF genomic "
    "COORDINATES, not amino-acid sequences; translating them back would require many free parameters "
    "and would produce our library, not IEAtlas's. It is treated as genuinely unavailable, and its "
    "contribution is carried through as the single free unknown m (novel 9-mers it adds)."
)


def read_fasta(path):
    """Same reader as library_ambiguity.py."""
    seqs, cur = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).upper())
    return seqs


def read_fasta_human_2022(path):
    """SwissProt 2022_01, human subset by OX=9606 -- same selector as era_correct_reference.py."""
    seqs, keep, cur = [], False, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    seqs.append("".join(cur).upper())
                keep = "OX=9606 " in line or line.rstrip().endswith("OX=9606")
                cur = []
            elif keep:
                cur.append(line.strip())
    if keep and cur:
        seqs.append("".join(cur).upper())
    return seqs


def kmers_of(seqs, k=K):
    out = set()
    for s in seqs:
        for i in range(len(s) - k + 1):
            out.add(s[i:i + k])
    return out


def pct(a, b):
    return 100.0 * a / b if b else 0.0


def measure(name, canon_seqs, NU, TL, both):
    """All three rates + the union, against one canonical reference. Returns a dict."""
    print(f"  indexing canonical {K}-mers for {name} ...", flush=True)
    C = kmers_of(canon_seqs)
    print(f"    {len(C):,} distinct canonical {K}-mers "
          f"({len(canon_seqs):,} proteins, {sum(map(len, canon_seqs)):,} residues)", flush=True)

    nu_hit = len(NU & C)
    tl_hit = len(TL & C)
    both_hit = len(both & C)

    # Inclusion-exclusion -- exact, and it never materialises the union set.
    u = len(NU) + len(TL) - len(both)
    h = nu_hit + tl_hit - both_hit

    res = {
        "reference": name,
        "canonical_proteins": len(canon_seqs),
        "canonical_kmers": len(C),
        "nuorfdb": {"kmers": len(NU), "canonical": nu_hit, "pct": round(pct(nu_hit, len(NU)), 2)},
        "translnc": {"kmers": len(TL), "canonical": tl_hit, "pct": round(pct(tl_hit, len(TL)), 2)},
        "shared_nuorfdb_translnc": {"kmers": len(both), "canonical": both_hit,
                                    "pct": round(pct(both_hit, len(both)), 2)},
        "union": {"kmers": u, "canonical": h, "pct": round(pct(h, u), 2)},
        "union_effect_pp": round(pct(h, u) - pct(nu_hit, len(NU)), 2),
    }
    del C
    return res


def interval_table(u, h, canon_kmers):
    """The sharp achievable interval for the 3-source library as a function of the single unknown m.

    combined = (h + x) / (u + m),  with  0 <= x <= min(m, canon_kmers - h).

    lo(m) = h / (u + m)                                  [RPFdb adds NO canonical 9-mers]
    hi(m) = (h + min(m, canon_kmers - h)) / (u + m)      [RPFdb adds as many as the canonical universe
                                                          still allows -- x cannot exceed the canonical
                                                          9-mers not already covered by U]
    """
    headroom = canon_kmers - h                       # canonical 9-mers U does not already contain
    grid = [
        ("0  (RPFdb adds nothing new)", 0),
        ("0.10 x u", int(0.10 * u)),
        ("0.25 x u", int(0.25 * u)),
        ("0.50 x u", int(0.50 * u)),
        # m = headroom is where the MAX column peaks -- the ceiling. Shown so the non-monotone MAX
        # column is legible: it climbs to this row, then falls.
        ("= headroom  <<< MAX peaks here", headroom),
        ("1.00 x u", int(1.00 * u)),
        ("2.00 x u", int(2.00 * u)),
        ("5.00 x u", int(5.00 * u)),
        ("10.0 x u", int(10.00 * u)),
    ]
    grid.sort(key=lambda t: t[1])
    rows = []
    for label, m in grid:
        lo = pct(h, u + m)
        hi = pct(h + min(m, headroom), u + m)
        rows.append({"label": label, "m": m, "lo_pct": round(lo, 2), "hi_pct": round(hi, 2)})

    # The ceiling over ALL m: hi(m) rises while m <= headroom, then falls. So it peaks at m = headroom.
    m_star = headroom
    ceiling = pct(canon_kmers, u + headroom)         # = (h + headroom) / (u + headroom)
    return rows, headroom, m_star, ceiling


def main():
    for p in (SPROT_NOW, NUORFDB, TRANSLNC):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    print("=" * 96)
    print("THE UNION OF IEAtlas's ncORF SEARCH LIBRARIES -- measured, not assumed")
    print("=" * 96)
    print(f"""
  IEAtlas integrated THREE sources: RPFdb v2.0 + nuORFdb v1.2 + TransLnc.
  We hold TWO of the three, in the versions IEAtlas used. We measure their union directly, in full,
  with no sampling. RPFdb is unavailable and is carried as an explicit unknown, not approximated.
""")

    print("  building library k-mer sets (full, no sampling) ...", flush=True)
    nu_orfs = read_fasta(NUORFDB)
    NU = kmers_of(nu_orfs)
    print(f"    nuORFdb v1.2 : {len(nu_orfs):,} ORFs      -> {len(NU):,} distinct {K}-mers", flush=True)
    del nu_orfs

    tl_pep = read_fasta(TRANSLNC)
    TL = kmers_of(tl_pep)
    print(f"    TransLnc     : {len(tl_pep):,} peptides -> {len(TL):,} distinct {K}-mers", flush=True)
    del tl_pep

    both = NU & TL
    print(f"    shared       : {len(both):,} {K}-mers occur in BOTH libraries "
          f"({pct(len(both), len(NU) + len(TL) - len(both)):.1f}% of their union)\n", flush=True)

    # k-mers containing a non-standard residue can never match canonical and would deflate the rates.
    # Count them so the reader can see the effect is nil rather than take our word for it.
    nu_bad = sum(1 for km in NU if not set(km) <= STD_AA)
    tl_bad = sum(1 for km in TL if not set(km) <= STD_AA)
    print(f"  sanity: {K}-mers containing a non-standard residue (X/U/...): "
          f"nuORFdb {nu_bad:,}, TransLnc {tl_bad:,} -- negligible, left in (no filtering).\n")

    refs = [("SwissProt modern (human)", read_fasta(SPROT_NOW))]
    if os.path.exists(SPROT_2022):
        refs.append(("SwissProt 2022_01 (human, OX=9606) [era-correct]",
                     read_fasta_human_2022(SPROT_2022)))
    else:
        print(f"  NOTE: {SPROT_2022} absent -- era-correct arm skipped.\n")

    results = []
    for name, seqs in refs:
        results.append(measure(name, seqs, NU, TL, both))
        print()

    primary = results[0]

    # ---- pipeline validation gate -------------------------------------------------------------
    nu_pct = primary["nuorfdb"]["pct"]
    print("=" * 96)
    print("PIPELINE VALIDATION GATE")
    print("=" * 96)
    ok = abs(nu_pct - 34.1) <= 0.5
    print(f"  nuORFdb v1.2 alone, modern reference : {nu_pct:.1f}%   "
          f"(library_ambiguity.py reported 34.1%)")
    if not ok:
        print("\n  *** GATE FAILED. The pipeline does not reproduce the published 34.1%. Every number")
        print("      below is therefore untrustworthy and MUST NOT be used. Stop and fix the pipeline.\n")
        json.dump({"gate": "FAILED", "nuorfdb_pct": nu_pct, "expected": 34.1},
                  open(ART, "w"), indent=2)
        return 1
    print("  GATE PASSED -- the pipeline reproduces the published figure. Proceed.\n")

    # ---- the three rates ----------------------------------------------------------------------
    print("=" * 96)
    print(f"THE THREE RATES -- fraction of each library's distinct {K}-mers that are ALSO canonical")
    print("=" * 96)
    for r in results:
        print(f"\n  reference: {r['reference']}  ({r['canonical_kmers']:,} canonical {K}-mers)")
        print(f"    {'library':<44}{'distinct 9-mers':>16}{'also canonical':>16}{'rate':>9}")
        print("    " + "-" * 85)
        for key, label in (("nuorfdb", "(a) nuORFdb v1.2  [as IEAtlas used]"),
                           ("translnc", "(b) TransLnc      [as IEAtlas used]"),
                           ("union", "(c) nuORFdb UNION TransLnc  <<< KEY")):
            d = r[key]
            print(f"    {label:<44}{d['kmers']:>16,}{d['canonical']:>16,}{d['pct']:>8.1f}%")
        print(f"    {'    (the ' + str(r['shared_nuorfdb_translnc']['kmers']) + ' shared 9-mers, for reference)':<44}"
              f"{'':>16}{r['shared_nuorfdb_translnc']['canonical']:>16,}"
              f"{r['shared_nuorfdb_translnc']['pct']:>8.1f}%")

    # ---- the monotonicity verdict -------------------------------------------------------------
    delta = primary["union_effect_pp"]
    u = primary["union"]["kmers"]
    h = primary["union"]["canonical"]
    print("\n" + "=" * 96)
    print("THE REVIEWER'S MONOTONICITY OBJECTION, TESTED EMPIRICALLY")
    print("=" * 96)
    print(f"""
  The objection: |(A u B) n C| / |A u B| is not monotone in adding B, so nuORFdb's 34.1% cannot be a
  lower bound on the combined library. That is a statement about what CAN happen. Here is what DOES
  happen when a real second IEAtlas source is added.

    nuORFdb alone            {primary['nuorfdb']['pct']:>6.1f}%
    nuORFdb u TransLnc       {primary['union']['pct']:>6.1f}%
    ------------------------------------
    effect of adding TransLnc {delta:>+6.1f} pp
""")
    if delta < 0:
        print(f"""  DIRECTION: the union rate FALLS. The reviewer is CONFIRMED, empirically and not just in
  principle. Adding a genuine second IEAtlas source drops the rate by {abs(delta):.1f} pp, because
  TransLnc contributes overwhelmingly non-canonical {K}-mers ({primary['translnc']['pct']:.1f}% canonical, versus
  {primary['nuorfdb']['pct']:.1f}% for nuORFdb) and enlarges the denominator faster than the numerator. The retracted
  "lower bound" was not merely unproven -- it is FALSE in the direction the reviewer named, on the very
  data it was asserted over. The retraction stands and is now evidenced.""")
    elif delta > 0:
        print(f"""  DIRECTION: the union rate RISES, by {delta:.1f} pp. This is informative but it does NOT rescue the
  retracted claim. A single library that happens to raise the rate does not make the union operator
  monotone; the reviewer's objection was about what the operator GUARANTEES, and it still guarantees
  nothing. RPFdb -- the source we do not hold -- can still drive the combined rate anywhere in the
  interval below. 34.1% remains NOT a lower bound.""")
    else:
        print("""  DIRECTION: the rate is unchanged to the reported precision. Uninformative either way: it neither
  confirms nor rescues anything, and 34.1% remains NOT a lower bound.""")

    # ---- the sharp interval -------------------------------------------------------------------
    rows, headroom, m_star, ceiling = interval_table(u, h, primary["canonical_kmers"])
    print("\n" + "=" * 96)
    print("THE FULL 3-SOURCE LIBRARY: A SHARP INTERVAL IN THE SINGLE UNKNOWN")
    print("=" * 96)
    print(f"""
  U = nuORFdb u TransLnc.   u = |U| = {u:,} distinct {K}-mers.   h = |U n C| = {h:,}.
  Let RPFdb v2.0 add m NOVEL {K}-mers (not already in U), of which x are canonical.

      combined rate = (h + x) / (u + m),    0 <= x <= min(m, |C| - h)

  x <= |C| - h because a novel canonical {K}-mer must come from the canonical {K}-mer universe
  (|C| = {primary['canonical_kmers']:,}) and h of those are already inside U. Canonical headroom = {headroom:,}.
""")
    print(f"    {'m (novel 9-mers RPFdb adds)':<32}{'m':>14}{'MIN rate':>12}{'MAX rate':>12}")
    print("    " + "-" * 70)
    for r in rows:
        print(f"    {r['label']:<32}{r['m']:>14,}{r['lo_pct']:>11.1f}%{r['hi_pct']:>11.1f}%")
    print(f"""
  THE ONLY HARD STATEMENT WE CAN MAKE. Maximising over BOTH x and m: hi(m) climbs while m <= headroom
  and falls after, so it peaks at m = {m_star:,}, giving

      CEILING: the 3-source library's {K}-mer canonical-ambiguity rate CANNOT EXCEED {ceiling:.1f}%,
      whatever RPFdb v2.0 contains. This is distribution-free and needs no assumption about RPFdb.

  AND THE STATEMENT WE CANNOT MAKE. There is NO positive floor: lo(m) = h/(u+m) -> 0 as m grows. If
  RPFdb is large and overwhelmingly non-canonical, the combined rate can be driven arbitrarily close to
  zero. Our measurement does NOT determine the combined rate, and no measurement of a SUBSET of the
  libraries ever could. That is exactly the reviewer's point, now quantified.

  The union rate {primary['union']['pct']:.1f}% is the combined rate ONLY in the m = 0 corner -- i.e. only if RPFdb v2.0
  contributed no {K}-mer that nuORFdb and TransLnc had not already contributed. We have no evidence for
  that and do not assume it.
""")

    # ---- what remains unknown -----------------------------------------------------------------
    print("=" * 96)
    print("A COMPARISON THIS SCRIPT REFUSES TO MAKE")
    print("=" * 96)
    print(f"""
  It is tempting to notice that IEAtlas's published catalogue sits at 56.3% canonical overlap, which is
  ABOVE the {ceiling:.1f}% ceiling just derived, and to conclude that the output is "more canonical-ambiguous
  than any possible version of its own search library."

  DO NOT. The two numbers are in DIFFERENT UNITS and the comparison is arithmetic, not measurement:
    * {ceiling:.1f}% is over DISTINCT {K}-mers of the SEARCH SPACE (every {K}-mer the library could in principle
      yield, each counted once).
    * 56.3% is over DISTINCT CATALOGUED PEPTIDES of the OUTPUT, at their native lengths (8-12+), after
      the search, the FDR filter, and the atlas's own dedup.
  Different objects, different denominators, different lengths. This project has already retracted one
  cross-pipeline ratio ("11-40x outlier") for exactly this error. The ceiling constrains the LIBRARY.
  It says nothing directly about the OUTPUT rate, and no claim here rests on comparing them.
""")

    print("=" * 96)
    print("WHAT REMAINS UNKNOWN -- and why we did not paper over it")
    print("=" * 96)
    print(f"\n  {RPFDB_MISSING}\n")
    print(f"""  Concretely, to close the interval a future pass would need EITHER
    (i)  the RPFdb v2.0 amino-acid library as IEAtlas built it (archived copy, or from the authors), OR
    (ii) IEAtlas's integrated search FASTA itself -- which would make all three sources moot and give
         the combined rate directly, with no interval at all.
  Until one of those exists, the combined library's ambiguity is bounded above by {ceiling:.1f}% and is
  otherwise UNKNOWN. The paper should say precisely that, and no more.
""")

    payload = {
        "k": K,
        "sampling": "none -- full libraries measured; sampling is biased upward here",
        "libraries_measured": {
            "nuORFdb_v1.2": os.path.relpath(NUORFDB, REPO),
            "TransLnc": os.path.relpath(TRANSLNC, REPO),
        },
        "rpfdb_v2.0": {"measured": False, "reason": RPFDB_MISSING},
        "nonstandard_kmers": {"nuorfdb": nu_bad, "translnc": tl_bad},
        "by_reference": results,
        "primary_reference": primary["reference"],
        "monotonicity_test": {
            "nuorfdb_alone_pct": primary["nuorfdb"]["pct"],
            "union_pct": primary["union"]["pct"],
            "effect_of_adding_translnc_pp": delta,
            "direction": "DOWN" if delta < 0 else ("UP" if delta > 0 else "FLAT"),
            "verdict": ("reviewer CONFIRMED empirically: adding a real second IEAtlas source LOWERS the "
                        "union rate; 34.1% is not a lower bound" if delta < 0 else
                        "union rate rises here, but this does NOT establish monotonicity and 34.1% is "
                        "still not a lower bound"),
        },
        "interval": {
            "u_union_kmers": u,
            "h_union_canonical": h,
            "canonical_kmers": primary["canonical_kmers"],
            "canonical_headroom": headroom,
            "formula": "(h + x) / (u + m), 0 <= x <= min(m, |C| - h)",
            "rows": rows,
            "ceiling_pct": round(ceiling, 2),
            "ceiling_attained_at_m": m_star,
            "floor_pct": None,
            "floor_note": "no positive floor exists: h/(u+m) -> 0 as m grows without bound",
            "combined_rate_determined": False,
        },
    }
    json.dump(payload, open(ART, "w"), indent=2)
    print(f"  wrote {os.path.relpath(ART, REPO)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
