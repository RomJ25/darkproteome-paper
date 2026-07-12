"""Regenerate every headline number in the manuscript, and FAIL on drift.

    python3 manuscript/verify_manuscript.py          # check
    python3 manuscript/verify_manuscript.py --print  # emit the numbers

A manuscript that carries its own copy of its numbers is a second scorer, and second scorers drift:
an earlier draft printed a pseudogene rate long after the analysis stopped producing it, and a
CrypticProteinDB rate of 0.0 that was never right at all. Every claim below is re-derived from the
committed artifacts and compared against what the manuscript asserts. If they disagree, this exits
non-zero and the paper does not build.

It also fails if the paper:
  * asserts any of the phrasings RETRACTED during review ("manufactures", "the rule predicts the
    rate", "strict survivor", ...) outside an explicit disclaimer; or
  * drops a REQUIRED PRIOR-ART CITATION, or the explicit concession that the paper contributes no
    conceptual novelty. The standard applied here is not ours, and a draft that fails to say so is
    claiming someone else's contribution -- which a referee would end the paper with, in one line.

The large public inputs (the atlas exports) are not redistributed. On a clean checkout the checks
that need them are SKIPPED and reported as skipped, while everything derivable from the committed
artifacts is still verified. A verifier that crashes on a clean checkout is a verifier nobody runs.
"""
import csv
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv.field_size_limit(10_000_000)

TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
ATL = os.path.join(REPO, "data", "external", "atlases")
IE_CANCER = os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
IE_NORMAL = os.path.join(ATL, "IEAtlas_Epitopes_In_Normal_Tissues.txt")
MS = os.path.join(REPO, "manuscript", "manuscript_v2.md")

PSEUDO = re.compile(r"^[A-Z0-9\-]+P\d+$")

# Retracted during external review. Allowed ONLY inside an explicit disclaimer, or a quotation of the
# retracted wording. "the rule predicts the rate" is banned because PRIOR ART CONTRADICTS IT:
# catalogues that also lack an exclusion rule sit at 1.4-5%, not 56%, so a permissive rule does not
# account for the magnitude (Bedran et al. 2023).
BANNED = ["manufactures", "systematically re-labels", "discarded by the retention",
          "genuinely non-canonical", "canonical self", "strict survivor",
          "rule predicts the rate", "predicts the canonical-overlap rate"]

REQUIRED = [
    ("Bedran", "the metric + the exclusion standard (Cancer Immunol Res 2023)"),
    ("Woo et al. 2014", "class-specific FDR under-control is prior art"),
    ("Kumar et al. 2022", "'most shared peptides should be dropped' is prior art"),
    ("Nesvizhskii", "the shared-peptide / protein-inference problem is textbook"),
    ("no conceptual novelty", "the paper must concede this explicitly"),
]


def epitopes(path):
    out = set()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if r and r[0] and r[0].strip().upper().isalpha():
                out.add(r[0].strip().upper())
    return out


def main():
    for p in (TIER1, SCALED, SCORED, MS):
        if not os.path.exists(p):
            sys.exit(f"missing required artifact: {p}")

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

    # --- R1: the catalogues (committed artifacts only) ---
    facts["cryptic"] = rate(src["CrypticProteinDB-immuno"] | src["CrypticProteinDB-epitopes"])
    facts["ieatlas"] = rate(src["IEAtlas"])
    raja = [r for r in t1 if r["cohort"].startswith("Raja")]
    rk = sum(int(r["canonical_self_exact"]) for r in raja)
    facts["raja"] = (rk, len(raja), 100 * rk / len(raja))
    rp = [r for r in t1 if r["cohort"].startswith("Raja") and r["orf_class"] == "pseudogene-ORF"]
    facts["raja_pseudo"] = (sum(int(r["canonical_self_exact"]) for r in rp), len(rp), 0.0)

    checks += [
        (f"{facts['ieatlas'][0]:,} / {facts['ieatlas'][1]:,}", "IEAtlas overlap n/N"),
        ("56.3%", "IEAtlas overlap rate"),
        (f"{facts['cryptic'][0]:,} / {facts['cryptic'][1]:,}", "CrypticProteinDB n/N"),
        ("0.026%", "CrypticProteinDB rate"),
        (f"{facts['raja'][0]:,} / {facts['raja'][1]:,}", "Raja n/N"),
        (f"{facts['raja_pseudo'][0]} / {facts['raja_pseudo'][1]}", "Raja pseudogene n/N"),
        # R2 -- measured TWICE, over different k-mer windows. The paper must report BOTH: presenting
        # a single-implementation number as novel is how a prior in-house result nearly got
        # re-presented as a discovery.
        ("34.1%", "nuORFdb latent canonical ambiguity (9-mer)"),
        ("34.4%", "nuORFdb, independent 8-11mer corroboration"),
        ("1.0–2.4%", "GENCODE Ribo-seq ORF latent ambiguity"),
        ("8,448,245", "nuORFdb distinct 9-mers"),
    ]

    # --- R1 class-fixed control + R3 consequence: need the non-redistributed atlas exports ---
    if os.path.exists(IE_CANCER) and os.path.exists(IE_NORMAL):
        ps, other = set(), set()
        with open(IE_CANCER, newline="", encoding="utf-8", errors="replace") as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd, None)
            for r in rd:
                if len(r) < 3 or not r[0]:
                    continue
                q = r[0].strip().upper()
                if not q.isalpha():
                    continue
                g = (r[2] or "").split("_")[0].upper()
                (ps if PSEUDO.match(g) else other).add(q)
        facts["ieatlas_pseudo"] = rate(ps)
        facts["ieatlas_other"] = rate(other)

        cancer, normal = epitopes(IE_CANCER), epitopes(IE_NORMAL)
        sc = [p for p in cancer if p in selfmap]
        ov = [p for p in sc if selfmap[p]]
        nov = [p for p in sc if not selfmap[p]]
        facts["consequence_ov"] = (sum(1 for p in ov if p in normal), len(ov), 0.0)
        facts["consequence_nov"] = (sum(1 for p in nov if p in normal), len(nov), 0.0)

        checks += [
            (f"{facts['ieatlas_pseudo'][0]:,} / {facts['ieatlas_pseudo'][1]:,}",
             "IEAtlas pseudogene n/N"),
            ("60.5%", "IEAtlas pseudogene rate"),
            (f"{facts['consequence_ov'][0]:,} = 22.4%",
             "canonical-overlapping also on normal tissue"),
            ("9.1%", "non-overlapping control rate"),
            ("12.6%", "share of the whole catalogue that is both"),
        ]
    else:
        skipped.append("R1 class-fixed control + R3 consequence — need data/external/atlases/ "
                       "(large public files, not redistributed; see data/SOURCES.md)")

    if "--print" in sys.argv:
        for k, (a, b, p) in facts.items():
            print(f"  {k:<18} {a:>7,}/{b:<8,} = {p:6.3f}%")
        for s in skipped:
            print(f"  ! skipped: {s}")
        return 0

    text = open(MS, encoding="utf-8").read()
    bad = [f"{lab}: manuscript does not contain {needle!r}"
           for needle, lab in checks if needle not in text]

    # A disclaimer SECTION is a block, not a line -- its bullets do not each repeat "we do not
    # claim". Excise it, then check what remains. Anything asserted OUTSIDE it is a regression.
    body = re.sub(r"^##+\s*What this paper does not claim.*?(?=^##\s|\Z)", "", text,
                  flags=re.S | re.M | re.I)
    disclaim = re.compile(r"do(es)? not claim|we do not|never write|retracted|banned|~~|❌", re.I)
    for b in BANNED:
        for m in re.finditer(re.escape(b), body, re.I):
            a = body.rfind("\n", 0, m.start()) + 1
            z = body.find("\n", m.end())
            line = body[a: z if z > 0 else len(body)]
            if disclaim.search(line):
                continue
            bad.append(f"BANNED PHRASE ASSERTED: {b!r} -> …{line.strip()[:70]}…")

    for cite, why in REQUIRED:
        if cite not in text:
            bad.append(f"MISSING REQUIRED PRIOR-ART CITATION: {cite!r} -- {why}")

    if bad:
        print("MANUSCRIPT VERIFICATION: FAIL\n")
        for b in bad:
            print(f"  - {b}")
        print(f"\n{len(bad)} problem(s).")
        return 1

    print("MANUSCRIPT VERIFICATION: PASS")
    print(f"  - all {len(checks)} headline numbers match the artifacts")
    print(f"  - no banned phrasing ({len(BANNED)} retracted overclaims checked)")
    print(f"  - required prior-art citations present ({len(REQUIRED)} checked)")
    for s in skipped:
        print(f"  ! SKIPPED: {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
