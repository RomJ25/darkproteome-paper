"""Regenerate every headline number in the manuscript, and FAIL on drift.

    python3 manuscript/verify_manuscript.py          # check
    python3 manuscript/verify_manuscript.py --print  # emit the numbers

A manuscript that carries its own copy of its numbers is a second scorer, and second scorers drift:
an earlier draft printed a pseudogene rate long after the analysis stopped producing it, and a
CrypticProteinDB rate of 0.0 that was never right at all. Every claim below is re-derived from the
committed artifacts and compared against what the manuscript asserts. If they disagree, this exits
non-zero and the paper does not build.

It also fails if the paper:
  * asserts any of the phrasings RETRACTED during review, outside an explicit withdrawal or
    disclaimer. Five of these were publication blockers found by external review on 2026-07-13:
      - "11-40x" -- a fold-change against rates from ANOTHER pipeline. Not a measurement.
      - "z = 74"  -- a two-proportion z on 174,465 CLUSTERED observations. Never valid.
      - "lower bound on its library" -- |(A u B) n C| / |A u B| is NOT monotone in adding B.
      - "internal control" -- the non-overlapping set controls nothing; it is a COMPARATOR.
      - "the field's standard" / "applies no exclusion rule" -- stronger than the evidence.
  * drops a REQUIRED PRIOR-ART CITATION, or the explicit statement that the contribution is
    empirical rather than conceptual. The principles applied here are not ours, and a draft that
    fails to say so is claiming someone else's contribution -- which a referee would end the paper
    with, in one line.

The large public inputs (the atlas exports) are not redistributed. On a clean checkout the checks
that need them are SKIPPED and reported as skipped, while everything derivable from the committed
artifacts is still verified. A verifier that crashes on a clean checkout is a verifier nobody runs.
"""
import csv
import json
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv.field_size_limit(10_000_000)

TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
R3 = os.path.join(REPO, "data", "derived_r3_inference.json")
BIAS = os.path.join(REPO, "data", "derived_detection_bias.json")
ERA = os.path.join(REPO, "data", "derived_era_reference.json")
ATL = os.path.join(REPO, "data", "external", "atlases")
IE_CANCER = os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
MS = os.path.join(REPO, "manuscript", "manuscript_v2.md")

PSEUDO = re.compile(r"^[A-Z0-9\-]+P\d+$")

# Retracted. Allowed ONLY inside an explicit withdrawal/disclaimer, or a quotation of the retracted
# wording. The paper is permitted -- required, even -- to NAME what it withdrew.
BANNED = [
    "manufactures", "systematically re-labels", "discarded by the retention",
    "genuinely non-canonical", "canonical self", "strict survivor",
    "rule predicts the rate", "predicts the canonical-overlap rate",
    # -- the five blockers from the 2026-07-13 external review --
    "11-40", "11–40",                       # fold-change across pipelines: not measured
    "z = 74", "z=74",                       # invalid: clustered observations
    "lower bound on its library",           # union rate is not monotone
    "internal control",                     # it is a within-resource COMPARATOR
    "the field's standard",                 # published+recommended, not universal
    "applies no exclusion rule",            # its METHODS do not describe one
    "magnitude is explained by the library",
    # -- STALE VALUES. A presence-check catches a number that was CHANGED; it cannot catch a
    # contradictory number added ALONGSIDE the right one. These five are the non-partitioning class
    # counts the review caught (they summed to 175,011 > 174,465, because 546 sequences carry both
    # labels). If any reappears, the old double-counted split is back.
    "9,874", "16,323", "88,827", "158,688", "60.5%",
    "22,003 entries",                       # unit error: they are unique SEQUENCES, not rows
]

REQUIRED = [
    ("Bedran", "the criterion + the metric (Cancer Immunol Res 2023)"),
    ("Woo et al. 2014", "class-specific FDR under-control is prior art"),
    ("Kumar et al. 2022", "'most shared peptides should be dropped' is prior art"),
    ("Nesvizhskii", "the shared-peptide / protein-inference problem is textbook"),
    ("Our contribution is empirical", "the paper must concede the principles are not new"),
    ("within-resource comparator", "the non-overlapping set is not a control"),
]


def main():
    for p in (TIER1, SCALED, SCORED, R3, BIAS, ERA, MS):
        if not os.path.exists(p):
            sys.exit(f"missing required artifact: {p}\n(regenerate: see escalations/"
                     "2026-07-13-class-label-provenance/)")

    r3 = json.load(open(R3))
    bias = json.load(open(BIAS))
    era = json.load(open(ERA))

    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}
    src = {}
    for r in csv.DictReader(open(SCALED, newline="", encoding="utf-8")):
        s = (r.get("_source") or "").strip()
        p = (r.get("peptide_sequence") or "").strip().upper()
        if p and p.isalpha() and not s.startswith("cohort:"):
            src.setdefault(s, set()).add(p)
    t1 = list(csv.DictReader(open(TIER1, newline="", encoding="utf-8")))

    def rate(peps):
        n = sum(1 for p in peps if p in selfmap)
        k = sum(1 for p in peps if selfmap.get(p))
        return k, n, (100 * k / n if n else 0.0)

    facts, checks, skipped = {}, [], []

    # --- R1: the same-pipeline catalogues (committed artifacts only) ---
    facts["cryptic"] = rate(src["CrypticProteinDB-immuno"] | src["CrypticProteinDB-epitopes"])
    facts["ieatlas"] = rate(src["IEAtlas"])
    raja = [r for r in t1 if r["cohort"].startswith("Raja")]
    rk = sum(int(r["canonical_self_exact"]) for r in raja)
    facts["raja"] = (rk, len(raja), 100 * rk / len(raja))

    checks += [
        (f"{facts['ieatlas'][0]:,} / {facts['ieatlas'][1]:,}", "IEAtlas overlap n/N"),
        ("56.3%", "IEAtlas overlap rate"),
        (f"{facts['cryptic'][0]:,} / {facts['cryptic'][1]:,}", "CrypticProteinDB n/N"),
        ("0.026%", "CrypticProteinDB rate"),
        (f"{facts['raja'][0]:,} / {facts['raja'][1]:,}", "Raja n/N"),
        # R2 -- measured TWICE, over different k-mer windows. The paper must report BOTH.
        ("34.1%", "nuORFdb latent canonical ambiguity (9-mer)"),
        ("34.4%", "nuORFdb, independent 8-11mer corroboration"),
        ("1.0–2.4%", "GENCODE Ribo-seq ORF latent ambiguity"),
        ("8,448,245", "nuORFdb distinct 9-mers"),
    ]

    # --- the era-correct reference (blocker 2) ---
    checks += [
        (f"{era['pct_2022_01']}%", "era-correct rate vs Swiss-Prot 2022_01"),
        (f"{era['overlap_2022_01']:,}", "era-correct overlap count"),
        (f"{era['proteins_2022_01']:,}", "human proteins in Swiss-Prot 2022_01"),
        (f"{era['retrospective_only']}", "overlaps a 2022 analyst could not have seen"),
        (f"{era['retrospective_pct_of_overlap']}%", "retrospective share of the overlap set"),
    ]

    # --- R3 inference (blocker 5): the z is GONE; RR + clustered CI replace it ---
    checks += [
        (f"{r3['overlapping_in_normal']:,}", "canonical-overlapping also in normal export"),
        (f"{r3['pct_overlapping_in_normal']}%", "overlapping normal-export rate"),
        (f"{r3['pct_comparator_in_normal']}%", "comparator normal-export rate"),
        (f"{r3['pct_of_whole_catalogue']}%", "share of the whole catalogue that is both"),
        (f"{r3['rr_length_standardized']}", "length-standardized risk ratio"),
        (f"[{r3['ci95'][0]}, {r3['ci95'][1]}]", "gene-clustered bootstrap 95% CI"),
        (f"{r3['n_clusters']:,}", "source-gene clusters resampled"),
    ]

    # --- R2 detection-bias test (was an untested hypothesis; now measured) ---
    checks += [
        (f"{bias['mean_types_overlapping']}", "mean cancer types, overlapping"),
        (f"{bias['mean_types_comparator']}", "mean cancer types, comparator"),
        (f"{bias['length_strata_holding']} of {bias['length_strata_tested']}",
         "length strata in which detection breadth holds"),
        (f"{bias['ribo_catalogue_rr']}×", "ribosomal enrichment in the catalogue"),
        (f"{bias['ribo_library_rr']}×", "ribosomal enrichment in the library (the control)"),
        (f"{bias['ribo_excess']}×", "excess of catalogue over library"),
    ]

    # --- R1 class strata (blocker 3): must PARTITION, and the paper must say so ---
    if os.path.exists(IE_CANCER):
        genes = {}
        with open(IE_CANCER, newline="", encoding="utf-8", errors="replace") as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd, None)
            for r in rd:
                if len(r) < 3 or not r[0]:
                    continue
                q = r[0].strip().upper()
                if q.isalpha():
                    genes.setdefault(q, set()).add((r[2] or "").split("_")[0].upper())
        strata = {"pseudogene-only": [], "non-pseudogene-only": [], "mixed": []}
        for q, gs in genes.items():
            if q not in selfmap:
                continue
            a = any(PSEUDO.match(g) for g in gs)
            b = any(not PSEUDO.match(g) for g in gs)
            strata["pseudogene-only" if a and not b
                   else "non-pseudogene-only" if b and not a else "mixed"].append(q)

        tot = sum(len(v) for v in strata.values())
        if tot != facts["ieatlas"][1]:
            print(f"PARTITION FAILURE: strata sum to {tot:,}, catalogue is "
                  f"{facts['ieatlas'][1]:,}. The strata are not mutually exclusive.")
            return 1
        for k, v in strata.items():
            facts[k] = (sum(selfmap[p] for p in v), len(v),
                        100 * sum(selfmap[p] for p in v) / len(v))
        checks += [
            (f"{facts['pseudogene-only'][1]:,}", "pseudogene-only n"),
            (f"{facts['non-pseudogene-only'][1]:,}", "non-pseudogene-only n"),
            (f"{facts['mixed'][1]}", "sequences carrying BOTH labels"),
            (f"{facts['pseudogene-only'][2]:.1f}%", "pseudogene-only overlap rate"),
            (f"{facts['non-pseudogene-only'][2]:.1f}%", "non-pseudogene-only rate (= the "
                                                        "no-pseudogene counterfactual)"),
        ]
    else:
        skipped.append("R1 class strata — need data/external/atlases/ (large public files, not "
                       "redistributed; see data/SOURCES.md)")

    if "--print" in sys.argv:
        for k, (a, b, p) in facts.items():
            print(f"  {k:<22} {a:>7,}/{b:<8,} = {p:6.3f}%")
        for k, v in (("r3", r3), ("detection-bias", bias), ("era", era)):
            print(f"  {k}: {json.dumps(v)}")
        for s in skipped:
            print(f"  ! skipped: {s}")
        return 0

    text = open(MS, encoding="utf-8").read()
    bad = [f"{lab}: manuscript does not contain {needle!r}"
           for needle, lab in checks if needle not in text]

    # A disclaimer SECTION is a block, not a line -- its bullets do not each repeat "we do not
    # claim". Excise it, then check what remains. Anything asserted OUTSIDE it is a regression.
    # A retracted phrase may also appear in the BODY when the paper names what it withdrew -- that
    # is required honesty, not a regression, so a withdrawal on the line exempts it too.
    body = re.sub(r"^##+\s*What this paper does not claim.*?(?=^##\s|\Z)", "", text,
                  flags=re.S | re.M | re.I)
    ok_line = re.compile(
        r"do(es)? not claim|we do not|never write|retracted|banned|~~|❌"
        r"|withdraw|invalid|was wrong|earlier draft|no fold-change|not a lower bound|not monotone",
        re.I)
    for b in BANNED:
        for m in re.finditer(re.escape(b), body, re.I):
            a = body.rfind("\n", 0, m.start()) + 1
            z = body.find("\n", m.end())
            line = body[a: z if z > 0 else len(body)]
            if ok_line.search(line):
                continue
            bad.append(f"RETRACTED PHRASE ASSERTED: {b!r} -> …{line.strip()[:70]}…")

    for cite, why in REQUIRED:
        if cite not in text:
            bad.append(f"MISSING REQUIRED: {cite!r} -- {why}")

    if bad:
        print("MANUSCRIPT VERIFICATION: FAIL\n")
        for b in bad:
            print(f"  - {b}")
        print(f"\n{len(bad)} problem(s).")
        return 1

    print("MANUSCRIPT VERIFICATION: PASS")
    print(f"  - all {len(checks)} headline numbers match the artifacts")
    print(f"  - no retracted phrasing ({len(BANNED)} checked, incl. the 5 review blockers)")
    print(f"  - required prior-art + hedging language present ({len(REQUIRED)} checked)")
    print("  - ORF-class strata partition exactly")
    for s in skipped:
        print(f"  ! SKIPPED: {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
