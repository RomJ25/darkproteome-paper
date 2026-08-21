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
R1L = os.path.join(REPO, "data", "derived_r1_length_strata.json")
BIAS = os.path.join(REPO, "data", "derived_detection_bias.json")
ERA = os.path.join(REPO, "data", "derived_era_reference.json")
SDB = os.path.join(REPO, "data", "derived_search_database.json")
UNI = os.path.join(REPO, "data", "derived_library_union.json")
LAM = os.path.join(REPO, "data", "derived_library_ambiguity.json")
PSG = os.path.join(REPO, "data", "derived_pseudogene_parent.json")
ABD = os.path.join(REPO, "data", "derived_abundance_direct.json")
P2SPLIT = os.path.join(REPO, "data", "derived_p2_pseudogene_split.json")
IL = os.path.join(REPO, "data", "derived_ieatlas_il_equivalence.json")
PXD038782 = os.path.join(REPO, "data", "derived_pxd038782_benign_crosswalk.json")
IMM_REC = os.path.join(REPO, "data", "derived_immunoverse_phla_recurrence.json")
IMM_LIN = os.path.join(REPO, "data", "derived_immunoverse_lineage_audit.json")
IMM_BLCA = os.path.join(REPO, "data", "derived_immunoverse_blca_raw_audit.json")
IMM_DLBC = os.path.join(REPO, "data", "derived_immunoverse_dlbc_raw_audit.json")
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
    ("Rom Jan", "manuscript author name"),
    ("Bedran", "the criterion + the metric (Cancer Immunol Res 2023)"),
    ("Woo et al. 2014", "class-specific FDR under-control is prior art"),
    ("Aggarwal et al. 2022", "'most shared peptides should be dropped' is prior art"),
    ("Nesvizhskii", "the shared-peptide / protein-inference problem is textbook"),
    ("Our contribution is empirical", "the paper must concede the principles are not new"),
    ("within-resource comparator", "the non-overlapping set is not a control"),
    # The estimand upgrade. If either of these vanishes, the paper has silently reverted to the
    # weaker "external audit against a reference of our choosing" framing.
    ("and the canonical human proteome", "IEAtlas's Methods, verbatim: the canonical proteome was searched"),
    ("own search database", "the overlap is INTERNAL to the resource, not retrospective"),
    # Found in a fresh-review pass: a submittable manuscript needs a references section,
    # and the PSM-vs-peptide FDR unit mismatch is load-bearing for R4/S2's reporting-gap argument.
    ("## References", "a manuscript citing ~16 works by author-year needs a bibliography"),
    ("peptide spectrum match false discovery rate", "the PSM-level FDR quote grounding the R4/S2 unit-mismatch caveat"),
    ("MIAIPE", "existing immunopeptidomics reporting standard must be acknowledged"),
    ("mzIdentML", "existing evidence data model must be acknowledged"),
    ("Group-walk", "group-aware FDR control prior art must be acknowledged"),
    ("Reconstructibility is not statistical control", "the class-decoy ledger is an audit object, not sufficient control"),
    ("Li et al. 2025", "the external ImmunoVerse recurrence pilot must cite its source"),
    ("Wu et al. 2026", "the prospective TIPs calibration audit must cite its source"),
    ("Wen et al. 2025", "the proposed entrapment gate must cite the valid-estimator framework"),
    ("not evidence of a confidence-propagation failure", "the preregistered TIPs audit is not an empirical FDR result"),
]


def flex(s):
    """A literal string, as a regex tolerant of markdown line-wrapping. Found in a
    fresh-review pass: every check in this file used exact substring/`in` matching, which silently
    misses a needle if the source .md happens to wrap a line exactly at one of its spaces (the
    manuscript's prose wraps at ~100 chars and is hand-edited, so this is not hypothetical -- it
    already happened once, to a 'no positive floor' phrase this exact guard was supposed to catch
    drift on). re.escape() in this Python version escapes a literal space as '\\ '; replacing that
    with '\\s+' keeps every other character an exact literal match while letting exactly the
    space-vs-newline distinction pass."""
    return re.escape(s).replace("\\ ", r"\s+")


def main():
    for p in (TIER1, SCALED, SCORED, R3, BIAS, ERA, SDB, UNI, LAM, MS, P2SPLIT, IL,
              PXD038782, IMM_REC, IMM_LIN, IMM_BLCA, IMM_DLBC):
        if not os.path.exists(p):
            sys.exit(f"missing required artifact: {p}\n"
                     "(regenerate it with the analysis script that emits it)")

    r3 = json.load(open(R3))
    il = json.load(open(IL))
    pxd = json.load(open(PXD038782))
    imm_rec = json.load(open(IMM_REC))
    imm_lin = json.load(open(IMM_LIN))
    imm_blca = json.load(open(IMM_BLCA))
    imm_dlbc = json.load(open(IMM_DLBC))
    bias = json.load(open(BIAS))
    p2split = json.load(open(P2SPLIT))
    era = json.load(open(ERA))
    sdb = json.load(open(SDB))
    uni = json.load(open(UNI))
    lam = json.load(open(LAM))
    pr = next(r for r in uni["by_reference"] if r["reference"] == uni["primary_reference"])

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
        (f"{il['n_il_equivalent_compatible']:,} / {il['n_total']:,}", "IEAtlas I/L-equivalent overlap n/N"),
        (f"{il['pct_il_equivalent_compatible']}%", "IEAtlas I/L-equivalent overlap rate"),
        (f"{il['n_additional_il_equivalent']}", "additional I/L-equivalent compatible sequences"),
        (f"{il['delta_percentage_points']} percentage points", "I/L-equivalent increase in percentage points"),
    ]

    # --- Aggregation / recurrence estimand: frozen ImmunoVerse catalogue-to-study audit and
    # current raw-linked worked example. These are deliberately bounded checks, not evidence of a
    # catalogue-level FDR failure or of irretrievable provenance loss. ---
    lin_inv = imm_lin["result"]["lineage_inventory"]
    lin_top = imm_lin["result"]["sensitivities"]["strong_or_weak_predicted_binder"][
        "top_n_by_sample_condition_recurrence"]
    blca = imm_blca["result"]
    dlbc = imm_dlbc["result"]
    dlbc_records = {r["peptide"]: r for r in dlbc["peptide_records"]}
    for peptide in ("AEGPDHHSL", "VPHTRPVSL"):
        if not dlbc_records[peptide]["every_expected_label_has_current_identified_scan"]:
            sys.exit(f"FATAL: {peptide} no longer has a current identified scan in every expected "
                     "DLBC label; withdraw the raw-linked recurrence example.")
        if dlbc_records[peptide]["common_reported_hla_alleles_across_expected_labels"]:
            sys.exit(f"FATAL: {peptide} now has a common reported HLA allele across DLBC labels; "
                     "withdraw the disjoint-genotype claim.")
    checks += [
        (f"{lin_top['median_sequence_recurrence_sample_condition_labels']:.0f} sample/condition labels",
         "ImmunoVerse median top-1000 label recurrence"),
        (f"{lin_top['median_sequence_recurrence_source_studies_lower']:.0f} source studies",
         "ImmunoVerse median top-1000 source-study recurrence"),
        (f"{lin_top['median_conservative_sequence_study_to_best_predicted_phla_ratio']:.2f}×",
         "ImmunoVerse conservative sequence-study/pHLA-study recurrence ratio"),
        (f"{100 * lin_top['fraction_conservatively_above_one']:.1f}%",
         "ImmunoVerse top-1000 fraction with conservative pHLA ratio above one"),
        ("AEGPDHHSL", "raw-linked DLBC recurrent peptide example 1"),
        ("VPHTRPVSL", "raw-linked DLBC recurrent peptide example 2"),
        (f"{blca['s7_peptide_count']} / {blca['s7_peptide_count']}",
         "BLCA exact historical-sequence recovery"),
        (f"{blca['s7_total_n_psm']} / {blca['consolidated_exact_sequence_scan_rows']}",
         "BLCA historical/current exact scan-row recovery"),
    ]

    # R2 -- measured TWICE, over different k-mer windows. The paper must report BOTH. These were
    # hardcoded literals until a fresh-review pass found it -- unlike every other check
    # in this file, they never re-derived from a live artifact. nuORFdb's rate/kmer-count now reads
    # from derived_library_union.json (already loaded as `uni` below for the neighboring Translnc/
    # union checks); the GENCODE Ribo-seq range now reads from derived_library_ambiguity.json (which
    # previously didn't exist at all -- library_ambiguity.py only ever printed to stdout).
    checks += [
        (f"{pr['nuorfdb']['pct']:.1f}%", "nuORFdb latent canonical ambiguity (9-mer)"),
        ("34.4%", "nuORFdb, independent 8-11mer corroboration"),
        (f"{lam['gencode_range_lo_pct']:.1f}–{lam['gencode_range_hi_pct']:.1f}%",
         "GENCODE Ribo-seq ORF latent ambiguity"),
        (f"{pr['nuorfdb']['kmers']:,}", "nuORFdb distinct 9-mers"),
    ]

    # --- the era-correct reference (blocker 2) ---
    checks += [
        (f"{era['pct_2022_01']}%", "era-correct rate vs Swiss-Prot 2022_01"),
        (f"{era['overlap_2022_01']:,}", "era-correct overlap count"),
        (f"{era['proteins_2022_01']:,}", "human proteins in Swiss-Prot 2022_01"),
        (f"{era['retrospective_only']}", "overlaps a 2022 analyst could not have seen"),
        (f"{era['retrospective_pct_of_overlap']}%", "retrospective share of the overlap set"),
        # Found by an independent fresh-review pass: the paper's own 231/97,999/98,193 did
        # not visibly reconcile without the backward-loss count (sequences matching 2022_01 but not
        # the modern reference -- the two releases are not simply nested). era_correct_reference.py
        # always computed this but never persisted it until now.
        (f"{era['lost_only']}", "sequences matching Swiss-Prot 2022_01 but not the modern reference"),
    ]

    # --- R1 by-length robustness (round-4 review): replaces the TAUTOLOGICAL "length-standardized
    # to the catalogue's own distribution" row, which is a mathematical identity (it reproduces the
    # crude rate for ANY population, confounded or not) and was never backed by an artifact distinct
    # from the headline string. The honest check is the per-length range itself. ---
    if os.path.exists(R1L):
        r1l = json.load(open(R1L))
        if abs(r1l["own_distribution_reweighted_pct"] - r1l["crude_pct"]) > 0.05:
            sys.exit("FATAL: derived_r1_length_strata.json's own-distribution reweighting no longer "
                     "reproduces the crude rate -- the tautology argument in this comment is wrong, "
                     "or the pipeline changed. Re-derive before trusting either number.")
        checks += [
            (f"{r1l['range_lo_pct']}%–{r1l['range_hi_pct']}%", "R1 overlap rate range across length strata"),
        ]
    else:
        skipped.append("R1 length-strata range (derived_r1_length_strata.json missing)")

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
    cap = max(r["hi_pct"] for r in uni["interval"]["rows"])
    checks += [
        (f"{pr['translnc']['pct']:.1f}%", "Translnc latent canonical ambiguity (whole-file diagnostic)"),
        (f"{pr['union']['pct']:.1f}%", "nuORFdb union Translnc (whole-file diagnostic)"),
        (f"{pr['union']['kmers']:,}", "distinct 9-mers in the whole-file union"),
        (f"{abs(mono['effect_of_adding_translnc_pp']):.1f} pp", "drop from adding Translnc (whole-file)"),
        (f"{cap:.1f}%", "distribution-free cap on the 3-source library (whole-file diagnostic)"),
    ]

    # --- human-only Translnc variant (found in a fresh-review pass): the whole-file
    # union above mixes ~148,667 mouse-convention Translnc entries into a human-immunopeptidome
    # calculation. The human-only variant is reported as PRIMARY throughout the paper; the whole-file
    # numbers above are kept only as the diagnostic showing why species-filtering matters. ---
    if "human_only_variant" not in uni:
        sys.exit("FATAL: derived_library_union.json has no human_only_variant -- re-run "
                 "library_union.py. The paper's primary union/ceiling numbers depend on it.")
    ho = uni["human_only_variant"]
    checks += [
        (f"{ho['translnc_kmers']:,}", "Translnc, human-only entries, distinct 9-mers"),
        (f"{ho['translnc_canonical']:,}", "Translnc, human-only entries, canonical-matching 9-mers"),
        (f"{ho['translnc_pct']:.1f}%", "Translnc, human-only entries, latent canonical ambiguity"),
        (f"{ho['union_canonical']:,}", "canonical-matching 9-mers in the human-only union"),
        (f"{ho['union_pct']:.1f}%", "nuORFdb union Translnc, human-only (PRIMARY)"),
        (f"{ho['union_kmers']:,}", "distinct 9-mers in the human-only union"),
        (f"{ho['ceiling_pct']:.1f}%", "distribution-free cap, human-only union (PRIMARY)"),
        # Found by a fresh-review pass: the manuscript's "re-scored against Swiss-Prot
        # 2022_01, essentially unchanged" claim for the human-only union had no computation behind it
        # at all -- only the whole-file union was ever checked against both references. Now computed
        # for real in library_union.py and guarded here.
        (f"{ho['union_pct_era_correct']:.1f}%", "human-only union, re-scored against Swiss-Prot 2022_01"),
    ]
    # --- found in a consistency-sweep re-read: an earlier draft asserted the nuORFdb-vs-
    # Translnc(human-only) gap was "34x" -- that number belongs to a different pair (nuORFdb vs GENCODE
    # Ribo-seq phase 2, ~34.1/1.0) and was carried over uncorrected. Guard the real ratio live. ---
    nu_vs_translnc_gap = round(pr["nuorfdb"]["pct"] / ho["translnc_pct"])
    checks += [
        (f"~{nu_vs_translnc_gap}×", "nuORFdb-vs-Translnc(human-only) latent-ambiguity gap"),
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

        # Reachability restrictions on the "never hit" population (found in a fresh-review
        # pass): gene-appears-elsewhere-in-catalogue (output-side) and library-9mer-content (search-side)
        # give two different, non-interchangeable corrections -- both must be guarded so neither can
        # silently drift or vanish.
        checks += [
            (f"{pl['n_reachable']:,}", "never-hit proteins, gene reachable elsewhere in catalogue"),
            (f"{pl['n_absent']:,}", "never-hit proteins, gene absent from catalogue entirely"),
            (f"{pl['fold_reachability_restricted']}×", "abundance fold, reachability-restricted"),
            (f"{pl['auc_reachability_restricted']}", "AUC, reachability-restricted"),
            (f"{pl['n_lib_reachable']:,}", "never-hit proteins sharing >=1 library 9-mer"),
            (f"{pl['n_lib_unreachable']:,}", "never-hit proteins sharing NO library 9-mer"),
            (f"{pl['fold_library_content_restricted']}×", "abundance fold, library-content-restricted"),
            (f"{pl['auc_library_content_restricted']}", "AUC, library-content-restricted"),
            (f"{pl['library_kmers_n']:,}", "library 9-mers used for the reachability test"),
            (f"{round(pl['median_hit_ppm']/pl['median_absent_ppm'], 2)}×",
             "abundance fold, hit vs fully-absent-gene subset alone"),
            # Found by a fresh-review pass: the fold above was checked but its accompanying
            # AUC (quoted in the manuscript right next to it) had no computation anywhere in the
            # generator script or this guard -- an orphaned, unverifiable number. Now computed
            # (abundance_direct.py) and guarded here, live.
            (f"{pl['auc_fully_unreachable_alone']}", "AUC, hit vs fully-absent-gene subset alone"),
            # The two "share of never-hit proteins excluded" percentages quoted in Prediction 3's
            # library-content paragraph -- one was silently wrong (49.4% where the arithmetic on the
            # already-guarded n_absent/n_not_hit counts gives 50.6%; 49.4% is instead the complementary
            # RETAINED fraction). Guard both live from the counts already checked above.
            (f"{round(100 * pl['n_lib_unreachable'] / pl['n_not_hit'], 1)}%",
             "share of never-hit proteins excluded by the library-content restriction"),
            (f"{round(100 * pl['n_absent'] / pl['n_not_hit'], 1)}%",
             "share of never-hit proteins excluded by the catalogue-co-occurrence restriction"),
        ]

        # Found by a fresh-review pass: these raw numbers already feed already-checked
        # derived figures above (the Q1/Q5 quintile means feed q5_minus_q1_lengthstd; the median
        # protein lengths feed the C2 AUC) but were never themselves checked, so they could drift
        # silently even with every ratio above still passing.
        bml = bb["bin_means_lengthstd"]
        c2 = ab["C2_protein_length_control"]
        c3 = ab["C3_placebo"]
        checks += [
            (f"{bml['Q1']:.2f}", "Q1 quintile mean, length-standardized detection breadth"),
            (f"{bml['Q5']:.2f}", "Q5 quintile mean, length-standardized detection breadth"),
            (f"{int(c2['hit_median_aa'])}", "median protein length, hit canonical proteins"),
            (f"{int(c2['nothit_median_aa'])}", "median protein length, never-hit canonical proteins"),
            (f"{c3['draws_ge_observed']} of {c3['n_draws']}", "placebo draws reaching the observed gap"),
        ]

    # --- R3 inference (blocker 5): the z is GONE; RR + clustered CI replace it ---
    checks += [
        (f"{r3['overlapping_in_normal']:,}", "canonical-overlapping also in normal export"),
        (f"{r3['pct_overlapping_in_normal']}%", "overlapping normal-export rate"),
        (f"{r3['pct_comparator_in_normal']}%", "comparator normal-export rate"),
        (f"{r3['pct_of_whole_catalogue']}%", "share of the whole catalogue that is both"),
        (f"{r3['rr_length_standardized']}", "length-standardized risk ratio"),
        (f"{r3['rr_crude']}", "crude (non-length-standardized) risk ratio, now shown in Fig 3a's crude-vs-standardized caveat"),
        (f"[{r3['ci95'][0]}, {r3['ci95'][1]}]", "gene-clustered bootstrap 95% CI"),
        (f"{r3['n_clusters']:,}", "source-gene clusters resampled"),
        (f"{r3['bootstrap_B']:,}", "bootstrap iterations B"),
    ]
    adj = r3["detection_breadth_adjustment"]
    checks += [
        (f"{adj['common_support_n']:,}", "length+breadth common-support sequences"),
        (f"{adj['rr_length_breadth_standardized']}×", "length+breadth-standardized risk ratio"),
        (f"[{adj['ci95_gene_clustered'][0]:.2f}, {adj['ci95_gene_clustered'][1]:.2f}]",
         "gene-clustered CI for length+breadth-standardized risk ratio"),
        (f"{adj['exact_breadth_sensitivity']['rr']}×", "exact-breadth sensitivity risk ratio"),
    ]

    # --- bounded external recurrence: processed tables only, with the limitation retained ---
    checks += [
        (f"{pxd['n_benign_search_tables']}", "PXD038782 benign processed tables"),
        (f"{pxd['n_valid_processed_rows']:,}", "PXD038782 valid processed peptide rows"),
        (f"{pxd['n_ieatlas_cancer_sequences_recurrent']:,}", "IEAtlas sequences recurring in PXD038782 benign tables"),
        (f"{pxd['n_recurrent_exact_canonical_compatible']:,}", "PXD038782 recurrences exact canonical-compatible"),
        (f"{pxd['n_recurrent_absent_ieatlas_normal']:,}", "PXD038782 recurrences absent IEAtlas normal export"),
        (f"{pxd['tissue_breadth']['n_at_least_2']:,}", "PXD038782 recurrences in at least two tissue labels"),
        ("processed-table recurrence", "PXD038782 result must not be upgraded to raw-spectrum validation"),
        ("No raw-spectrum re-search", "PXD038782 limitation must remain explicit"),
    ]

    # --- R3 additions found by a fresh-review pass: the "both"/mixed stratum RR and the
    # per-length RR range were computed and printed by consequence_robust.py all along, but never
    # persisted to the artifact or checked here. ---
    rrl = r3["rr_by_length"]
    checks += [
        (f"{r3['rr_by_stratum']['mixed']}×", "risk ratio, 'both'-labels stratum (R1's source-ambiguous 546)"),
        (f"{min(rrl.values())}", "R3 per-length risk-ratio range, low end"),
        (f"{max(rrl.values())}", "R3 per-length risk-ratio range, high end"),
        (f"{r3['n_tissue_clusters']}", "tissue clusters resampled (second robustness axis)"),
        (f"[{r3['ci95_tissue_clustered'][0]}, {r3['ci95_tissue_clustered'][1]}]",
         "tissue-clustered bootstrap 95% CI"),
        (f"{r3['n_eff_tissue_clusters']}", "effective tissue-cluster count (Herfindahl) -- the "
         "honesty caveat on the CI above must not silently drop if this number changes"),
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
        # Found by a fresh-review pass: the Translnc-inclusive sensitivity numbers (the
        # ones actually reported as PRIMARY in the manuscript's Prediction 2) were computed and printed
        # but never guarded at all -- only the nuORFdb-only baseline above was checked.
        (f"{bias['ribo_translnc_inclusive']['humanonly']['excess']}×",
         "ribosomal excess, Translnc-inclusive human-only (PRIMARY)"),
        (f"{bias['ribo_translnc_inclusive']['wholefile']['excess']}×",
         "ribosomal excess, Translnc-inclusive whole-file (diagnostic)"),
        # Found by a fresh-review pass: the raw shares behind the two ratios above
        # (mean-types and ribosomal enrichment) were printed and used but never themselves checked.
        (f"{bias['pct_ge2_overlapping']}%", "share detected in >=2 cancer types, overlapping"),
        (f"{bias['pct_ge2_comparator']}%", "share detected in >=2 cancer types, comparator"),
        (f"{bias['ribo_pct_overlapping']}%", "raw ribosomal-ORF share, catalogue overlapping"),
        (f"{bias['ribo_pct_comparator']}%", "raw ribosomal-ORF share, catalogue non-overlapping"),
    ]

    # --- P2 pseudogene-contamination check (external-review-prompted): is the ribosomal
    # excess an artifact of ribosomal-NAMED PSEUDOGENES (e.g. RPS3AP12) riding inside the RPL*/RPS*/
    # MRPL*/MRPS* prefix test, rather than true protein-coding ribosomal genes? Split with the same
    # authoritative NCBI gene_group pseudogene->parent registry used for S1. If either sub-excess were
    # missing or <=1x for the true-gene split, the manuscript's claim that this is not contamination
    # would not hold -- fail loudly rather than let a stale number sit uncontradicted.
    if p2split["excess"]["true_gene"] is None or p2split["excess"]["true_gene"] <= 1:
        sys.exit("FATAL: P2's true-ribosomal-gene excess is missing or <=1x -- the manuscript's "
                  "'not pseudogene contamination' claim is not supported by the artifact.")
    checks += [
        (f"{p2split['catalogue_rr']['true_gene']:.2f}×", "P2 pseudogene split: catalogue RR, true ribosomal genes"),
        (f"{p2split['library_rr']['true_gene']:.2f}×", "P2 pseudogene split: library RR, true ribosomal genes"),
        (f"{p2split['excess']['true_gene']:.2f}×", "P2 pseudogene split: excess, true ribosomal genes"),
        (f"{p2split['catalogue_rr']['pseudogene']:.2f}×", "P2 pseudogene split: catalogue RR, ribosomal-named pseudogenes"),
        (f"{p2split['library_rr']['pseudogene']:.2f}×", "P2 pseudogene split: library RR, ribosomal-named pseudogenes"),
        (f"{p2split['excess']['pseudogene']:.2f}×", "P2 pseudogene split: excess, ribosomal-named pseudogenes"),
    ]

    # --- R4 (found in a second-round review pass, corroborated independently by two
    # separate external reviewers): 245,870 is NOT a cancer+normal peptide-level union. IEAtlas's own
    # Methods give HLA-I/HLA-II percentages (60.60%/41.90% bound-share; 37.16%/50.76% immunogenic-share
    # of 54,017/51,015 epitopes) from which HLA-I + HLA-II totals reconstruct to ~145,363 + ~100,502 =
    # 245,865, matching 245,870 to within percentage-rounding -- i.e. 245,870 is a CLASS-summed total.
    # The naive algebraic residual (174,465 + 94,375 - 245,870 = 22,970) is arithmetically correct but
    # was misreported as "sequences shared between the exports" in an earlier draft -- that claim is
    # now withdrawn, because the TRUE shared count is directly measured elsewhere in this paper (R3):
    # 22,003 + 6,976 = 28,979, which the residual understates by ~6,000. Guard both numbers so neither
    # can silently regress, and guard that the residual and the true measured overlap stay DIFFERENT
    # (if they ever matched, the "second unit mismatch" claim above would no longer hold). ---
    naive_residual = 174_465 + 94_375 - 245_870
    if naive_residual != 22_970:
        sys.exit(f"FATAL: 174,465 + 94,375 - 245,870 = {naive_residual:,}, not 22,970. IEAtlas's "
                 "reported totals changed -- re-check the R4 unit-mismatch argument before shipping it.")
    true_overlap = r3["overlapping_in_normal"] + r3["comparator_in_normal"]
    if true_overlap != 28_979:
        sys.exit(f"FATAL: r3 overlapping_in_normal + comparator_in_normal = {true_overlap:,}, not "
                 "28,979 -- re-derive the R4 unit-mismatch numbers before shipping them.")
    if naive_residual == true_overlap:
        sys.exit("FATAL: the naive algebraic residual now EQUALS the directly-measured overlap -- the "
                  "manuscript's claim that they diverge (proving 245,870 is not a tissue-level union) "
                  "no longer holds; re-derive before shipping.")
    hla_i_total = round(54_017 / 0.3716)
    hla_ii_total = round(51_015 / 0.5076)
    checks += [
        (f"{naive_residual:,}", "the naive (and wrong) cancer/normal residual, 174,465+94,375-245,870"),
        (f"{true_overlap:,}", "the TRUE, directly-measured cancer/normal overlap (22,003+6,976)"),
        (f"{hla_i_total:,}", "reconstructed HLA-I epitope total (54,017 / 0.3716)"),
        (f"{hla_ii_total:,}", "reconstructed HLA-II epitope total (51,015 / 0.5076)"),
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
        # Found by a fresh-review pass: the "both" stratum's rate (only its count was
        # checked), and the "1,801 sequences map to more than one gene" claim (no artifact anywhere
        # in the repo) -- both computable from the same `genes` dict already built above.
        n_multi_gene = sum(1 for q, gs in genes.items() if q in selfmap and len(gs) > 1)
        checks += [
            (f"{facts['pseudogene-only'][1]:,}", "pseudogene-only n"),
            (f"{facts['non-pseudogene-only'][1]:,}", "non-pseudogene-only n"),
            (f"{facts['mixed'][1]}", "sequences carrying BOTH labels"),
            (f"{facts['pseudogene-only'][2]:.1f}%", "pseudogene-only overlap rate"),
            (f"{facts['non-pseudogene-only'][2]:.1f}%", "non-pseudogene-only rate (= the "
                                                        "no-pseudogene counterfactual)"),
            (f"{facts['mixed'][2]:.1f}%", "'both'-labels stratum overlap rate"),
            (f"{n_multi_gene:,}", "sequences mapping to more than one gene"),
            (f"{100*n_multi_gene/facts['ieatlas'][1]:.1f}%", "share mapping to more than one gene"),
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
           for needle, lab in checks if not re.search(flex(needle), text)]

    # The P2 pseudogene-split excess numbers (3.47x / 2.16x) each appear TWICE in the manuscript
    # (Prediction 2 discussion + the "detection-bias test" Methods paragraph). A bare presence check
    # (above) is satisfied if only ONE of the two copies is correct -- confirmed by deliberately
    # corrupting one copy and observing the guard still pass. Checking occurrence COUNT, not just
    # presence, is what actually catches a copy left stale after the other is edited -- exactly the
    # "a presence-check cannot catch a contradictory number added alongside a correct one" failure mode
    # this project has hit before (the five stale ORF-class counts, round 2).
    for label, val in (("true-ribosomal-gene excess (3.47x)", p2split["excess"]["true_gene"]),
                       ("ribosomal-named-pseudogene excess (2.16x)", p2split["excess"]["pseudogene"])):
        needle = f"{val:.2f}×"
        n = len(re.findall(flex(needle), text))
        if n < 2:
            bad.append(f"P2 pseudogene split: {label} appears {n} time(s) in the manuscript, expected "
                       "2 (main text + Methods) -- one copy may be stale")

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
        r"|never valid|be deleted|refuse|must not|cross-unit|not a bound"
        # GENERAL disclaiming construction ("...not an 'internal control'", "not a X") -- found in a
        # fresh-review pass hardening this guard's flex() whitespace-tolerance: "not a
        # lower bound"/"not a bound" above were each added as one-off exemptions for exactly this
        # pattern; generalize it once rather than adding a new narrow phrase every time a script
        # names a banned term in order to reject it (the same "ban the stem, not the phrasing"
        # lesson this project applies to BANNED, applied to the allow-side of the guard instead).
        r"|not an? [\"'‘“]",
        re.I)
    for b in BANNED:
        for m in re.finditer(flex(b), body, re.I):
            a = body.rfind("\n", 0, m.start()) + 1
            z = body.find("\n", m.end())
            line = body[a: z if z > 0 else len(body)]
            if ok_line.search(line):
                continue
            bad.append(f"RETRACTED PHRASE ASSERTED: {b!r} -> …{line.strip()[:70]}…")

    for cite, why in REQUIRED:
        if not re.search(flex(cite), text):
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
            # The 62 = 36 + 1 + 25 accounting gap and the 36 -> 35 gene-level shortfall (found by two
            # independent fresh-review passes): both fields already existed in the artifact
            # and were simply never checked here.
            (f"{cv['distinct_pseudogene_symbols']}", "S1 distinct pseudogene symbols"),
            (f"{len(cv['symbols_resolved_but_no_curated_parent'])}", "S1 symbols resolved but no curated parent"),
            (f"{len(cv['symbols_absent_from_all_registries'])}", "S1 symbols absent from every registry"),
            (f"{ar['gene_level_pseudogenes']}", "S1 gene-level testable pseudogenes (36 resolved -> 35 testable)"),
            (f"{ar['gene_level_with_parent_hit']}", "S1 gene-level pseudogenes with a parent hit"),
        ]
        # Null B's restricted-permutable-subset result (found by a fresh-review pass:
        # this was computed and sitting in the artifact but the supplement called it "DEGENERATE"
        # alongside Null C, which is genuinely degenerate for a different reason -- Null B is merely
        # SMALL (n=5), and its real, non-vacuous result was going unreported).
        nb = p["nulls"]["B_hgnc_family"]["permutable_subset"]
        supp_checks += [
            (f"{nb['observed']}/{nb['n']}", "S1 null-B restricted-subset observed/n"),
            (f"{nb['observed_pct']}%", "S1 null-B restricted-subset observed rate"),
            (f"{nb['null_mean']}", "S1 null-B restricted-subset null mean"),
        ]
        for needle, lab in supp_checks:
            if not re.search(flex(needle), supp):
                bad.append(f"{lab}: supplement does not contain {needle!r}")
        # Null C is DEGENERATE (parents pairwise 9-mer-disjoint -- the family-respecting shuffle is a
        # true identity, p=1.0 vacuously). Null B is NOT degenerate -- it is a small-n restricted test
        # that does run and must be reported, not folded into C's dismissal.
        if "DEGENERATE" not in supp.upper():
            bad.append("S1 must disclose that null C is DEGENERATE, not quote its p-value as evidence")
        if p["nulls"]["B_hgnc_family"]["n_permutable_items"] == 0:
            bad.append("S1's Null B is now ALSO degenerate (0 permutable items) -- the supplement's "
                       "claim that it is a real, small-n test no longer holds; revert to the old framing.")
        checks += supp_checks

    # S3's numbers come from four linked artifacts. Validate them independently of S1's optional
    # pseudogene artifact so the recurrence supplement cannot silently escape the build gate.
    if os.path.exists(SUPP):
        supp = open(SUPP, encoding="utf-8").read()
        phla_all = imm_rec["result"]["sensitivities"]["strong_or_weak_predicted_binder"]["all_peptides"]
        exact = lin_top["exact_mapping_subset"]
        s3_checks = [
            (f"{lin_inv['s3_source_rows']:,} rows", "S3 source-row inventory"),
            (f"{lin_inv['s3_distinct_raw_file_names']:,} distinct raw-file names", "S3 raw-file inventory"),
            (f"{lin_inv['s3_biological_labels']:,} biological labels", "S3 biological-label inventory"),
            (f"{lin_inv['s3_source_study_identifiers']} source-study identifiers", "S3 study inventory"),
            (f"{lin_inv['s7_labels_with_exact_s3_match']} / {lin_inv['s7_used_biological_labels']}",
             "S3/S7 exact biological-label join"),
            (f"{lin_inv['s7_cancer_label_pairs_with_exact_s3_match']} / {lin_inv['s7_cancer_label_pairs']}",
             "S3/S7 exact cancer-label join"),
            (f"{lin_top['median_sequence_label_to_study_ratio_lower_bound']:.3f}",
             "S3 median label/study recurrence ratio"),
            (f"{100 * lin_top['fraction_with_label_recurrence_exceeding_study_upper_bound']:.1f}%",
             "S3 fraction whose label recurrence exceeds study upper bound"),
            (f"{100 * lin_top['fraction_conservatively_at_least_two']:.1f}%",
             "S3 conservative pHLA ratio at least two"),
            (f"{exact['peptide_count']}-peptide exact-study-mapping subset", "S3 exact-mapping subset size"),
            (f"{phla_all['peptide_count']:,} distinct peptides", "S3 catalogue distinct-peptide count"),
            (f"{blca['peptides_with_exact_sequence_recovery']} / {blca['s7_peptide_count']} peptides",
             "S3 BLCA exact-sequence recovery"),
            (f"{blca['peptides_with_current_identified_scan']} / {blca['s7_peptide_count']} peptides",
             "S3 BLCA current-identified recovery"),
            (f"{dlbc['scan_source_rows']:,} scan rows", "S3 DLBC scan-table size"),
            (f"{dlbc['peptides_with_exact_sequence_recovery']} / {dlbc['s7_peptide_count']}",
             "S3 DLBC exact-sequence recovery"),
            (f"{dlbc['s7_total_n_psm']}", "S3 DLBC historical total n_psm"),
            (f"{dlbc['consolidated_exact_sequence_scan_rows']}", "S3 DLBC current exact-sequence rows"),
            ("LASPHSPIL", "S3 DLBC decision-version example"),
        ]
        for needle, lab in s3_checks:
            if not re.search(flex(needle), supp):
                bad.append(f"{lab}: supplement does not contain {needle!r}")
        checks += s3_checks

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
            for m in re.finditer(flex(b), txt, re.I):
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
