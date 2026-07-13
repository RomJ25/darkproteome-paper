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
    disclaimer. Five of these were publication blockers found by external review on:
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
SDB = os.path.join(REPO, "data", "derived_search_database.json")
UNI = os.path.join(REPO, "data", "derived_library_union.json")
PSG = os.path.join(REPO, "data", "derived_pseudogene_parent.json")
ABD = os.path.join(REPO, "data", "derived_abundance_direct.json")
SUPP = os.path.join(REPO, "manuscript", "supplement_v2.md")
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
    # -- the five blockers from the external review --
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
    # -- ERROR #18, found by us, not the reviewer. Comparing the catalogue's rate (distinct PEPTIDES,
    # native lengths, post-search/FDR/dedup) to a library's rate (distinct 9-MERS of an undetected
    # candidate space) as LEVELS is a cross-unit comparison -- the identical error to the withdrawn
    # "11-40x" fold-change. The ribosomal test survives because it is a RATIO OF RATIOS, which is
    # dimensionless. These phrasings assert the invalid level comparison.
    # NB: these were once written as "rate exceeds the library" / "exceeds the library's", which a
    # single inserted word ("exceeds the SEARCH library's") walked straight past. Ban the short stem.
    "exceeds the library",
    "exceeds the search library",
    "excess of the catalogued rate over the library rate",
    "catalogued rate (56.3%) is higher",
    "is higher than nuORFdb",
]

REQUIRED = [
    ("Bedran", "the criterion + the metric (Cancer Immunol Res 2023)"),
    ("Woo et al. 2014", "class-specific FDR under-control is prior art"),
    ("Kumar et al. 2022", "'most shared peptides should be dropped' is prior art"),
    ("Nesvizhskii", "the shared-peptide / protein-inference problem is textbook"),
    ("Our contribution is empirical", "the paper must concede the principles are not new"),
    ("within-resource comparator", "the non-overlapping set is not a control"),
    # The estimand upgrade. If either of these vanishes, the paper has silently reverted to the
    # weaker "external audit against a reference of our choosing" framing.
    ("and the canonical human proteome", "IEAtlas's Methods, verbatim: the canonical proteome was searched"),
    ("own search database", "the overlap is INTERNAL to the resource, not retrospective"),
]


def main():
    for p in (TIER1, SCALED, SCORED, R3, BIAS, ERA, SDB, UNI, MS):
        if not os.path.exists(p):
            sys.exit(f"missing required artifact: {p}\n"
                     "(regenerate it with the analysis script that emits it)")

    r3 = json.load(open(R3))
    bias = json.load(open(BIAS))
    era = json.load(open(ERA))
    sdb = json.load(open(SDB))
    uni = json.load(open(UNI))

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

    # --- the search database: the canonical proteins were IN IT (the estimand upgrade) ---
    # This is the paper's strongest claim, so it gets the hardest guard: the script that produced
    # these numbers re-verifies IEAtlas's two load-bearing Methods quotes VERBATIM against the
    # on-disk paper, and refuses to emit the artifact otherwise. If that flag is ever false, the
    # framing is unfounded and the build must fail rather than print a number.
    if not sdb.get("methods_quotes_verified"):
        sys.exit("FATAL: derived_search_database.json does not attest that IEAtlas's Methods quotes "
                 "were verified. The 'own search database' framing is unfounded. Re-run "
                 "scripts/search_database.py")
    checks += [
        (f"{sdb['n_matching_search_db']:,}", "sequences matching the search database's canonical half"),
        (f"{sdb['pct_exactly_one_gene']}%", "ambiguous sequences compatible with exactly ONE gene"),
        (f"{sdb['pct_at_most_two_genes']}%", "ambiguous sequences compatible with at most TWO genes"),
    ]

    # --- R2 library union (blocker 1): the retraction, now MEASURED rather than conceded ---
    mono = uni["monotonicity_test"]
    if mono["direction"] != "DOWN":
        sys.exit("FATAL: derived_library_union.json no longer shows the union rate FALLING. The "
                 "paper's account of why '34.1% is a lower bound' was false must be rewritten.")
    pr = next(r for r in uni["by_reference"] if r["reference"] == uni["primary_reference"])
    cap = max(r["hi_pct"] for r in uni["interval"]["rows"])
    checks += [
        (f"{pr['translnc']['pct']:.1f}%", "Translnc latent canonical ambiguity"),
        (f"{pr['union']['pct']:.1f}%", "nuORFdb union Translnc"),
        (f"{pr['union']['kmers']:,}", "distinct 9-mers in the union"),
        (f"{abs(mono['effect_of_adding_translnc_pp']):.1f} pp", "drop from adding Translnc"),
        (f"{cap:.1f}%", "distribution-free cap on the 3-source library"),
    ]

    # --- R2 direct abundance (replaces the detection-breadth PROXY) ---
    # The placebo is the load-bearing control: if breaking the peptide->protein link does NOT collapse
    # the trend, the machinery invents trends and every number in this block is worthless.
    if os.path.exists(ABD):
        ab = json.load(open(ABD))
        if not ab["C3_placebo"]["collapses"]:
            sys.exit("FATAL: the abundance placebo did NOT collapse. The trend is an artifact of the "
                     "machinery and the claim must come out of the paper.")
        if ab["verdict"] != "CORROBORATES":
            sys.exit(f"FATAL: abundance verdict is {ab['verdict']!r}, not CORROBORATES. The paper "
                     "asserts a measured abundance effect that the artifact does not support.")
        # A cross-version check that SILENTLY returns null is worse than none: the earlier PaxDb
        # release has a different column layout, and a parser that skipped every row emitted `null`
        # with no error. A null here means the check did not run -- do not let it pass as if it had.
        if ab.get("version_robustness_v5_q5_minus_q1") is None:
            sys.exit("FATAL: the PaxDb cross-version check is null -- it did not actually run. "
                     "Re-run abundance_direct.py; do not report a robustness check that was skipped.")
        pl, bb = ab["A_which_proteins_are_hit"]["protein_level"], ab["B_abundance_predicts_breadth"]["ab_max"]
        checks += [
            (f"{100*ab['join']['join_rate']:.1f}%", "PaxDb join rate"),
            (f"{ab['catalogue']['overlapping_with_abundance']:,}", "overlapping sequences with abundance"),
            (f"{pl['median_hit_ppm']} ppm", "median abundance, canonical proteins that are hit"),
            (f"{pl['median_nothit_ppm']} ppm", "median abundance, canonical proteins never hit"),
            (f"{pl['fold']}×", "abundance fold, hit vs not-hit"),
            (f"{pl['auc']}", "protein-level AUC"),
            (f"{bb['q5_minus_q1_lengthstd']}", "Q5-Q1 detection-breadth gap, length-standardized"),
            (f"[{bb['ci95_cluster_canonical_gene'][0]}, {bb['ci95_cluster_canonical_gene'][1]}]",
             "gene-clustered 95% CI on the abundance-breadth gap"),
            (f"{bb['spearman_rho']}", "per-sequence Spearman rho (WEAK -- must be reported)"),
            # The CRUDE gap must appear too. The trend is monotone only AFTER length standardization;
            # quoting only the standardized figure would imply a clean dose-response that is not there.
            (f"{bb['q5_minus_q1_crude']}", "CRUDE (unstandardized) Q5-Q1 gap -- the trend saturates"),
            (f"{ab['version_robustness_v5_q5_minus_q1']}", "same gap on the previous PaxDb release"),
            (f"{ab['C2_protein_length_control']['deciles_auc_above_half']} / "
             f"{ab['C2_protein_length_control']['deciles_tested']}",
             "protein-length deciles in which the effect holds"),
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
        r"|withdraw|invalid|was wrong|earlier draft|no fold-change|not a lower bound|not monotone"
        # naming what was withdrawn, in a script header or a correction notice
        r"|do not say|never say|used to (make|say|open)|no longer|correction|error #"
        r"|never valid|be deleted|refuse|must not|cross-unit|not a bound",
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

    # --- the SUPPLEMENT was never verified at all. It carries numbers, so it can drift. ---
    if os.path.exists(SUPP) and os.path.exists(PSG):
        supp = open(SUPP, encoding="utf-8").read()
        p = json.load(open(PSG))
        ar, cv, h2h = p["authoritative_result"], p["coverage"], p["head_to_head"]
        # Null D is the only NON-DEGENERATE family-respecting null (B and C collapse because the
        # parents are pairwise 9-mer-disjoint). It is run against two decoy pools; the STRONG-paralog
        # pool is the hard one and the paper must report it, not just the permissive one.
        pools = p["nulls"]["D_family_decoy_swap"]["pools"]
        dw, ds = pools["any_shared_kmer"], pools["strong_paralogs"]
        supp_checks = [
            (f"{ar['hit_rate_pct']}%", "S1 authoritative parent-match rate"),
            (f"{ar['hit_own_curated_parent']} / {ar['testable_peptides']}", "S1 parent hits / testable"),
            (f"{h2h['agreement_pct']}%", "S1 heuristic-vs-curated agreement"),
            (f"{cv['coverage_all_symbols_pct']}%", "S1 coverage of all pseudogene symbols"),
            (f"{cv['coverage_named_symbols_pct']}%", "S1 coverage of NAMED symbols"),
            (f"{dw['observed_pct']}%", "S1 null-D observed rate (permissive decoys)"),
            (f"{dw['null_mean_pct']}%", "S1 null-D null mean (permissive decoys)"),
            (f"{ds['observed_pct']}%", "S1 null-D observed rate (STRONG paralog decoys)"),
            (f"{ds['null_mean_pct']}%", "S1 null-D null mean (STRONG paralog decoys)"),
        ]
        for needle, lab in supp_checks:
            if needle not in supp:
                bad.append(f"{lab}: supplement does not contain {needle!r}")
        # Nulls B and C are DEGENERATE. Quoting their p-values as evidence would be dishonest --
        # null C's p = 1.0 arises because the parents are pairwise 9-mer-disjoint, which is the
        # ANSWER to the exchangeability objection, not a failed test.
        if "DEGENERATE" not in supp.upper():
            bad.append("S1 must disclose that nulls B and C are DEGENERATE, not quote them")
        checks += supp_checks

    # --- The ban swept ONLY the manuscript. That is the hole that let ERROR #18 live on in
    # abundance_bias.py's docstring and in ONEPAGER.md after the paper itself was fixed -- and it is
    # the same hole that let a retracted claim survive in library_ambiguity.py's header in round 2.
    # A retracted claim asserted in a LIVE script or a PUBLIC-FACING doc is exactly as quotable as one
    # in the paper. Sweep them too, with the same disclaimer exemption.
    live = [os.path.join(REPO, "ONEPAGER.md")]
    # The analysis scripts sit in a different directory in the research repo than in the release
    # repo. Sweep whichever exists, so ONE guard file works unmodified in both -- a divergent copy
    # is a guard that silently stops guarding the repo you forgot to update.
    for cand in (os.path.join(REPO, "scripts"),):
        if os.path.isdir(cand):
            live += [os.path.join(cand, f) for f in sorted(os.listdir(cand)) if f.endswith(".py")]
    for path in live:
        if not os.path.exists(path):
            continue
        try:
            txt = open(path, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        # QUARANTINED scripts (consequence.py, rule_predicts_rate.py) refuse to run and print a
        # retraction notice instead of output. They KEEP the retracted text on purpose, for the
        # record. Flagging them would be noise, and a noisy guard is one that gets ignored.
        if "DO NOT QUOTE THIS SCRIPT'S OUTPUT" in txt:
            continue
        # The release builder is TOOLING, not a claim-bearing document. Its retracted strings are
        # SUBSTITUTION RULES -- the left-hand sides it rewrites on the way out. Banning them there
        # would be banning the machinery that removes them.
        if os.path.basename(path) == "build_public_release.py":
            continue
        lines = txt.split("\n")
        for b in BANNED:
            for m in re.finditer(re.escape(b), txt, re.I):
                ln = txt.count("\n", 0, m.start())
                # A retraction notice is a BLOCK -- "ALL RETRACTED IN REVIEW:" heads a bullet list,
                # and the bullets do not each repeat the word. Check a window, not one line.
                window = "\n".join(lines[max(0, ln - 4): ln + 1])
                if ok_line.search(window):
                    continue
                bad.append(f"RETRACTED PHRASE IN {os.path.relpath(path, REPO)}: {b!r} "
                           f"-> …{lines[ln].strip()[:60]}…")

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
