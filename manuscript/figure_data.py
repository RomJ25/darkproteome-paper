"""Every number in every figure, COMPUTED from the live artifacts. Nothing typed in.

    python3 manuscript/figure_data.py          # print the manifest
    from figure_data import fig1, fig2         # use it

WHY THIS FILE EXISTS. `make_figures.py` used to hardcode its data:

    vals   = [0.2, 0.0, 0.0, 37.1, 56.3, 0.0]              # Fig 1a
    counts = [180650, 82453, 75471, 80, 0, 0]              # Fig 2

Both were stale and one was never right at all:
  * 37.1% was the pseudogene rate computed from a claim table that kept only the FIRST peptide
    of each multi-peptide gene row. The true peptide-level rate is 57.7%. The figure kept
    printing 37.1% long after the code stopped producing it.
  * CrypticProteinDB was a literal `0.0`, rendering as a zero-height "0%" bar. Its real rate is
    ~0.026% -- small, but NOT zero, and the difference is the entire point of showing it (the
    overlap is AVOIDABLE, and one atlas nearly avoids it). Its denominator is computed here as
    the UNION of the database's two files, which resolves a 3,774-vs-3,810 discrepancy that had
    been open since.
  * Fig 2's fourth stage was labelled "source-translation substantiated" and drew a count of 0 as
    if it were attrition. Nothing in the corpus reports a translation statistic, so that zero was
    structural -- the funnel was drawing a property of the scorer as a property of the claims.

A figure that carries its own copy of the numbers is a fourth scorer. This module is the fix:
the figures import their data, so a figure CANNOT disagree with the analysis that produced it.

`manifest()` prints every value with its provenance so the caption can be checked against it.
"""
import csv
import os
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(REPO, "src", "darkproteome"))
import evidence_dimensions as ed  # noqa: E402

csv.field_size_limit(10_000_000)
TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")

RAJA = "Raja ovarian (ads7405)"
HCC = "HCC / Camarena-Albà (adn3628)"


def _need(path, how):
    if not os.path.exists(path):
        sys.exit(f"missing {path}\n  regenerate with: {how}")


def fig1():
    """Panel (a): exact canonical-substring rate, by class x cohort, + the atlases.

    Every cell is a live numerator/denominator from tier1_nonnovelty.py's own output. The atlas
    rates are computed from the scaled catalog. Nothing is typed.
    """
    _need(TIER1, "python3 src/darkproteome/tier1_nonnovelty.py")
    _need(SCALED, "python3 src/darkproteome/ingest_atlases.py")

    tot, hit = defaultdict(int), defaultdict(int)
    with open(TIER1, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            k = (r["cohort"], r["orf_class"])
            tot[k] += 1
            hit[k] += int(r["canonical_self_exact"])

    bars = []
    for cohort, short in ((RAJA, "Raja"), (HCC, "HCC")):
        for cls, lab in (("altORF", "altORF"), ("lncRNA-ORF", "lncRNA"),
                         ("pseudogene-ORF", "pseudo."), ("other", "other")):
            n = tot[(cohort, cls)]
            if not n:
                continue
            k = hit[(cohort, cls)]
            bars.append({"label": f"{lab}\n({short})", "k": k, "n": n,
                         "pct": 100 * k / n, "cohort": short, "cls": cls})

    # ATLAS bars. These MUST be computed, not typed: the old figure hardcoded CrypticProteinDB as
    # a literal `0.0`, which rendered as a zero-height "0%" bar. Its true rate is small but NOT
    # zero, and the difference is the entire reason the bar is in the figure -- it shows the
    # overlap is AVOIDABLE. `canonical_self` per peptide comes from the scored table.
    _need(SCORED, "python3 src/darkproteome/reference_model.py")
    selfmap = {}
    with open(SCORED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            selfmap[r["peptide"]] = int(r["canonical_self"])

    seen = defaultdict(set)
    with open(SCALED, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            src = (r.get("_source") or "").strip()
            p = (r.get("peptide_sequence") or "").strip().upper()
            if not p or not p.isalpha() or src.startswith("cohort:"):
                continue
            seen[src].add(p)

    # CrypticProteinDB ships TWO files and the "CrypticProteinDB" bar means the whole database.
    # This denominator is exactly where a long-standing discrepancy lived: the manuscript said
    # 1/3,774 while a check computed 3,810 and the difference was never resolved.
    # It is simply immuno (3,774 unique peptides) + epitopes (36) = 3,810. Stating the union
    # explicitly is the fix; leaving it implicit is how the two numbers drifted for nine days.
    ATLAS_GROUPS = {
        "IEAtlas\n(atlas)": ("IEAtlas",),
        "Cryptic-\nProteinDB": ("CrypticProteinDB-immuno", "CrypticProteinDB-epitopes"),
    }
    for lab, srcs in ATLAS_GROUPS.items():
        peps = set().union(*(seen.get(s, set()) for s in srcs))
        n = sum(1 for p in peps if p in selfmap)      # only peptides the scorer actually saw
        k = sum(1 for p in peps if selfmap.get(p))
        if not n:
            continue
        bars.append({"label": lab, "k": k, "n": n, "pct": 100 * k / n,
                     "cohort": "atlas", "cls": "+".join(srcs)})

    return {"bars": bars, "atlas_peptides": {k: len(v) for k, v in seen.items()}}


def fig2():
    """The reporting-and-adjudicability matrix — REPLACES the survivorship funnel.

    The funnel is dead twice over: it drew a structurally-impossible dimension as attrition, and
    the survivor model it depicted no longer exists. This returns the matrix `audit.py` prints,
    pooled AND by stratum, because a pooled denominator dominated by 293k atlas rows hides the
    cohorts entirely.
    """
    _need(SCALED, "python3 src/darkproteome/ingest_atlases.py")
    rows = []
    with open(SCALED, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    def stratum(r):
        if ed.human_tcell_state(r) != "not-assayed":
            return "T-cell-tested"
        return "cohort" if (r.get("_source") or "").startswith("cohort:") else "atlas"

    strata = defaultdict(list)
    for r in rows:
        strata[stratum(r)].append(r)

    out = {"pooled": ed.matrix(rows), "n_pooled": len(rows), "strata": {}}
    for k in ("atlas", "cohort", "T-cell-tested"):
        out["strata"][k] = {"matrix": ed.matrix(strata[k]), "n": len(strata[k])}

    # the T-cell states, which no binary could express
    from collections import Counter
    out["tcell_states"] = dict(Counter(ed.human_tcell_state(r) for r in rows))
    out["normal_states"] = dict(Counter(ed.normal_presentation_state(r) for r in rows))
    return out


def manifest():
    print("=" * 78)
    print("FIGURE DATA MANIFEST — every value computed live; check captions against this")
    print("=" * 78)

    f1 = fig1()
    print("\nFIG 1a — exact canonical-substring rate (class x cohort)")
    print(f"  source: data/primary_tier1_nonnovelty.csv")
    for b in f1["bars"]:
        print(f"    {b['cohort']:<5} {b['cls']:<45} {b['k']:>5}/{b['n']:<6} = {b['pct']:6.3f}%")
    print("  atlas unique peptides:", f1["atlas_peptides"])
    print("  NOTE the pseudogene(HCC) cell: the figure used to print 37.1% here, from a claim")
    print("       table that kept only the FIRST peptide of each gene row. It is 57.7%.")

    f2 = fig2()
    print(f"\nFIG 2 — reporting & adjudicability matrix (n={f2['n_pooled']:,})")
    print("  source: data/claim_catalog_scaled.csv via evidence_dimensions.py")
    print(f"  {'dimension':<27}{'asserted':>10}{'claim-lnk':>11}{'quant':>8}"
          f"{'modality':>10}{'ADJUDIC.':>10}{'supports':>10}")
    for k in ed.DIMENSIONS:
        d = f2["pooled"][k]
        print(f"  {k:<27}{d['asserted']:>10,}{d['claim_linked']:>11,}{d['quantitative']:>8,}"
              f"{d['modality_appropriate']:>10,}{d['adjudicable']:>10,}{d[ed.SUPPORTS]:>10,}")
    print("\n  by stratum (allele reported / T-cell adjudicable):")
    for k, v in f2["strata"].items():
        m = v["matrix"]
        print(f"    {k:<16} n={v['n']:>7,}   allele={m['allele_restriction']['claim_linked']:>6,}"
              f"   T-cell={m['human_tcell_assay']['adjudicable']:>3,}")
    print(f"\n  human T-cell states: {f2['tcell_states']}")
    print(f"  normal-presentation states: {f2['normal_states']}")
    print("\n  The funnel this replaces drew stage 'source-translation substantiated' = 0 as")
    print("  ATTRITION. It was structural: nothing in the corpus reports a translation")
    print("  statistic, so the bar was a property of the scorer, not of the claims.")


if __name__ == "__main__":
    manifest()
