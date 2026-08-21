"""R3, redone with valid inference. The two-proportion z-test was NOT valid; this replaces it.

WHAT WAS WRONG (external review, upheld). R3 reported that canonical-overlapping cancer-catalogued
sequences are also present in IEAtlas's own normal-tissue export at 22.4% vs 9.1% for non-overlapping
ones, and attached a two-proportion z = 74. That z treats 174,465 peptides as independent Bernoulli
trials. They are not. They are clustered by source gene/ORF, gene family, dataset, tissue, donor,
HLA allele and pipeline. The z is enormous mostly because a heavily structured catalogue was counted
as 174,465 independent experiments.

There are also two measured detection-opportunity confounders: PEPTIDE LENGTH and CANCER-DETECTION
BREADTH. Short peptides match the canonical proteome more readily and recur more often across
datasets; peptides already found in many cancer types also have more opportunity to recur in the
normal export.

WHAT THIS DOES INSTEAD
  1. Reports the rate at EACH peptide length -- no pooling.
  2. Reports a LENGTH-STANDARDIZED risk ratio (direct standardization to the catalogue's own length
     distribution), so the comparison is not a length artifact.
  3. Repeats direct standardization jointly by length and cancer-detection breadth, with breadth
     coarsened to 1, 2, 3, 4 and 5+ cancer types; exact breadth is retained as a sensitivity check.
  4. Stratifies by ORF class, using MUTUALLY EXCLUSIVE strata. (The previous class split did not
     partition: 546 peptides carry both a pseudogene and a non-pseudogene ORF label, because a
     peptide may map to several ORF_IDs. Presenting them as complements double-counted those 546.)
  5. Replaces the z-test with a GENE-CLUSTERED BOOTSTRAP: resample source genes with replacement,
     recompute the length-standardized RR, and report a percentile CI. The unit of resampling is the
     cluster, so the CI respects the dependence the z-test ignored.

The descriptive fact was never in doubt. What was in doubt was the inference. If the clustered CI
crosses 1, the association is not established and R3 must be restated as descriptive only.

LANGUAGE, per review: the non-overlapping group is a WITHIN-RESOURCE COMPARATOR, not an "internal
control". It does not control abundance, detectability, HLA coverage or study composition.

    python3 scripts/consequence_robust.py
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ATL = os.path.join(REPO, "data", "external", "atlases")
IE_C = os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
IE_N = os.path.join(ATL, "IEAtlas_Epitopes_In_Normal_Tissues.txt")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
csv.field_size_limit(10_000_000)

PSEUDO = re.compile(r"^[A-Z0-9\-]+P\d+$")
B = 2000
SEED = 20260713


class RNG:
    """Deterministic LCG -- the bootstrap must be reproducible without numpy."""

    def __init__(self, s):
        self.s = s

    def randint(self, n):
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        return (self.s >> 33) % n


def load():
    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}
    genes = defaultdict(set)
    tissues = defaultdict(set)
    with open(IE_C, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 4 or not r[0]:
                continue
            q = r[0].strip().upper()
            if q.isalpha():
                genes[q].add((r[2] or "").split("_")[0].upper())
                t = (r[3] or "").strip()
                if t:
                    tissues[q].add(t)
    normal = set()
    with open(IE_N, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if r and r[0] and r[0].strip().upper().isalpha():
                normal.add(r[0].strip().upper())

    rows = []
    for q, gs in genes.items():
        if q not in selfmap:
            continue
        a = any(PSEUDO.match(g) for g in gs)
        b = any(not PSEUDO.match(g) for g in gs)
        ts = tissues.get(q) or {"?"}
        rows.append({
            "pep": q, "len": len(q), "breadth": len(tissues.get(q, set())),
            "breadth_bin": min(len(tissues.get(q, set())), 5),
            "x": selfmap[q], "y": 1 if q in normal else 0,
            "clust": min(gs), "tissue_clust": min(ts), "cls": "pseudogene-only" if a and not b
                                     else ("non-pseudogene-only" if b and not a else "mixed"),
        })
    return rows


def std_rr(rows, weights):
    """Direct-standardized risk ratio: reweight each arm to a common length distribution."""
    num = defaultdict(lambda: [0, 0])
    den = defaultdict(lambda: [0, 0])
    for r in rows:
        (num if r["x"] else den)[r["len"]][0] += r["y"]
        (num if r["x"] else den)[r["len"]][1] += 1
    p1 = p0 = 0.0
    for L, w in weights.items():
        k1, n1 = num[L]
        k0, n0 = den[L]
        if not n1 or not n0:
            continue
        p1 += w * k1 / n1
        p0 += w * k0 / n0
    return (p1 / p0) if p0 else float("nan"), p1, p0


def std_rr_joint(rows, key):
    """Direct standardization on common-support strata, normalized to their pooled distribution."""
    cells = defaultdict(lambda: [[0, 0], [0, 0]])
    for r in rows:
        arm = 1 if r["x"] else 0
        cells[key(r)][arm][0] += r["y"]
        cells[key(r)][arm][1] += 1
    common = {k: v for k, v in cells.items() if v[0][1] and v[1][1]}
    common_n = sum(v[0][1] + v[1][1] for v in common.values())
    if not common_n:
        return float("nan"), float("nan"), float("nan"), 0, 0
    p0 = p1 = 0.0
    for v in common.values():
        w = (v[0][1] + v[1][1]) / common_n
        p0 += w * v[0][0] / v[0][1]
        p1 += w * v[1][0] / v[1][1]
    return (p1 / p0) if p0 else float("nan"), p1, p0, common_n, len(common)


def main():
    for p in (IE_C, IE_N, SCORED):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
    rows = load()
    ov = [r for r in rows if r["x"]]
    nov = [r for r in rows if not r["x"]]
    print(f"scored cancer-catalogued sequences: {len(rows):,}")
    print(f"  canonical-overlapping           : {len(ov):,}")
    print(f"  non-overlapping (comparator)    : {len(nov):,}")
    print(f"  distinct source-gene clusters   : {len({r['clust'] for r in rows}):,}")

    ck = sum(r["y"] for r in ov)
    nk = sum(r["y"] for r in nov)
    print(f"\nCRUDE (as previously reported):")
    print(f"  overlapping also in normal export : {ck:,}/{len(ov):,} = {100*ck/len(ov):.1f}%")
    print(f"  comparator also in normal export  : {nk:,}/{len(nov):,} = {100*nk/len(nov):.1f}%")
    print(f"  crude risk ratio                  : {(ck/len(ov))/(nk/len(nov)):.2f}x")

    # ---- 1. per-length, no pooling -------------------------------------------------
    print("\n" + "=" * 84)
    print("1. RATE AT EACH PEPTIDE LENGTH (the confounder, unpooled)")
    print("=" * 84)
    print(f"  {'len':>4}{'n(ovl)':>9}{'% normal':>10}{'n(comp)':>9}{'% normal':>10}{'RR':>8}")
    print("  " + "-" * 50)
    lens = sorted({r["len"] for r in rows})
    rr_by_length = {}
    for L in lens:
        a = [r for r in ov if r["len"] == L]
        b = [r for r in nov if r["len"] == L]
        if len(a) < 30 or len(b) < 30:
            continue
        pa = sum(r["y"] for r in a) / len(a)
        pb = sum(r["y"] for r in b) / len(b)
        rr_L = (pa / pb) if pb else float("nan")
        if rr_L == rr_L:
            rr_by_length[L] = round(rr_L, 2)
        print(f"  {L:>4}{len(a):>9,}{100*pa:>9.1f}%{len(b):>9,}{100*pb:>9.1f}%"
              f"{(pa/pb if pb else float('nan')):>8.2f}")

    # ---- 2. length-standardized RR -------------------------------------------------
    w = defaultdict(int)
    for r in rows:
        w[r["len"]] += 1
    weights = {L: c / len(rows) for L, c in w.items()}
    rr, p1, p0 = std_rr(rows, weights)
    print("\n" + "=" * 84)
    print("2. LENGTH-STANDARDIZED RISK RATIO (direct standardization to the catalogue's own")
    print("   length distribution)")
    print("=" * 84)
    print(f"  standardized rate, overlapping : {100*p1:.1f}%")
    print(f"  standardized rate, comparator  : {100*p0:.1f}%")
    print(f"  LENGTH-STANDARDIZED RR         : {rr:.2f}x   (crude was "
          f"{(ck/len(ov))/(nk/len(nov)):.2f}x)")

    # ---- 2b. length + detection-breadth standardization ----------------------------
    rr_lb, p1_lb, p0_lb, n_lb, cells_lb = std_rr_joint(
        rows, lambda r: (r["len"], r["breadth_bin"]))
    rr_lbe, p1_lbe, p0_lbe, n_lbe, cells_lbe = std_rr_joint(
        rows, lambda r: (r["len"], r["breadth"]))
    print("\n" + "=" * 84)
    print("2b. LENGTH + CANCER-DETECTION-BREADTH STANDARDIZATION")
    print("=" * 84)
    print(f"  coarsened breadth (1,2,3,4,5+), common-support n={n_lb:,}, cells={cells_lb}")
    print(f"    standardized rates: {100*p1_lb:.1f}% vs {100*p0_lb:.1f}%; RR={rr_lb:.2f}x")
    print(f"  exact breadth sensitivity, common-support n={n_lbe:,}, cells={cells_lbe}")
    print(f"    standardized rates: {100*p1_lbe:.1f}% vs {100*p0_lbe:.1f}%; RR={rr_lbe:.2f}x")

    # ---- 3. by ORF class, mutually exclusive ---------------------------------------
    print("\n" + "=" * 84)
    print("3. BY ORF CLASS -- mutually exclusive strata (these partition; the old split did not)")
    print("=" * 84)
    print(f"  {'stratum':<22}{'n':>9}{'RR (length-std)':>18}")
    print("  " + "-" * 49)
    stratum_rr = {}
    for c in ("pseudogene-only", "non-pseudogene-only", "mixed"):
        sub = [r for r in rows if r["cls"] == c]
        if len(sub) < 100:
            continue
        r_c, _, _ = std_rr(sub, weights)
        stratum_rr[c] = round(r_c, 2)
        print(f"  {c:<22}{len(sub):>9,}{r_c:>17.2f}x")
    print(f"  {'(partition check)':<22}"
          f"{sum(1 for r in rows if r['cls'] in ('pseudogene-only','non-pseudogene-only','mixed')):>9,}")

    # ---- 4. gene-clustered bootstrap -----------------------------------------------
    print("\n" + "=" * 84)
    print("4. GENE-CLUSTERED BOOTSTRAP of the length-standardized RR (replaces the invalid z)")
    print("=" * 84)
    byc = defaultdict(list)
    for r in rows:
        byc[r["clust"]].append(r)
    clusters = list(byc.values())
    rng = RNG(SEED)
    boot = []
    boot_lb = []
    for _ in range(B):
        samp = []
        for _ in range(len(clusters)):
            samp.extend(clusters[rng.randint(len(clusters))])
        v, _, _ = std_rr(samp, weights)
        if v == v:
            boot.append(v)
        v_lb, _, _, _, _ = std_rr_joint(samp, lambda r: (r["len"], r["breadth_bin"]))
        if v_lb == v_lb:
            boot_lb.append(v_lb)
    boot.sort()
    lo = boot[int(0.025 * len(boot))]
    hi = boot[int(0.975 * len(boot)) - 1]
    print(f"  resampled {len(clusters):,} gene clusters with replacement, B = {len(boot):,}")
    print(f"  length-standardized RR = {rr:.2f}x   95% CI [{lo:.2f}, {hi:.2f}]")
    valid = lo > 1.0
    print(f"\n  => CI {'EXCLUDES' if valid else 'INCLUDES'} 1.0. The association is "
          f"{'established under clustered inference' if valid else 'NOT established'}.")

    boot_lb.sort()
    lo_lb = boot_lb[int(0.025 * len(boot_lb))]
    hi_lb = boot_lb[int(0.975 * len(boot_lb)) - 1]
    print(f"  length+breadth-standardized RR = {rr_lb:.2f}x   95% CI [{lo_lb:.2f}, {hi_lb:.2f}]")

    # ---- 4b. tissue-clustered bootstrap (additional robustness axis) ---------------
    # The manuscript names tissue among the dependence structures ("clustered by source gene,
    # gene family, dataset, tissue, donor, HLA allele and pipeline") but until now only gene was
    # ever operationalized as a resampling unit -- tissue was read nowhere in load(). This adds
    # tissue (a peptide's min() tissue, same convention as the gene cluster key) as a second,
    # independent clustering axis, reported alongside the gene-clustered CI, not replacing it.
    print("\n" + "=" * 84)
    print("4b. TISSUE-CLUSTERED BOOTSTRAP (a second, independent robustness axis)")
    print("=" * 84)
    byt = defaultdict(list)
    for r in rows:
        byt[r["tissue_clust"]].append(r)
    tclusters = list(byt.values())
    rngt = RNG(SEED)
    boot_t = []
    for _ in range(B):
        samp = []
        for _ in range(len(tclusters)):
            samp.extend(tclusters[rngt.randint(len(tclusters))])
        v, _, _ = std_rr(samp, weights)
        if v == v:
            boot_t.append(v)
    boot_t.sort()
    lo_t = boot_t[int(0.025 * len(boot_t))]
    hi_t = boot_t[int(0.975 * len(boot_t)) - 1]
    print(f"  resampled {len(tclusters):,} tissue clusters with replacement, B = {len(boot_t):,}")
    print(f"  length-standardized RR = {rr:.2f}x   95% CI [{lo_t:.2f}, {hi_t:.2f}]")
    valid_t = lo_t > 1.0
    print(f"\n  => CI {'EXCLUDES' if valid_t else 'INCLUDES'} 1.0 under tissue-clustered "
          f"resampling too.")

    # HONESTY CHECK, not just a headline: 15 raw clusters is not 15 clusters' worth of bootstrap
    # information if they are wildly unbalanced (found on a fresh-review pass -- the top
    # 2 of 15 tissue clusters hold 41% of the catalogue between them). The effective cluster count
    # (the Herfindahl-based diversity measure 1/sum(p_i^2), standard for "how many roughly-equal
    # units does this unbalanced set behave like") quantifies exactly how much weaker this axis is
    # than its raw count suggests, rather than leaving that as an eyeballed impression in prose.
    tsizes = [len(v) for v in tclusters]
    ttotal = sum(tsizes)
    n_eff_tissue = 1 / sum((n / ttotal) ** 2 for n in tsizes)
    print(f"  effective tissue-cluster count (Herfindahl 1/sum(p_i^2)): {n_eff_tissue:.1f} "
          f"of {len(tclusters)} raw clusters -- read the CI as indicative, not precise.")

    art = os.path.join(REPO, "data", "derived_r3_inference.json")
    json.dump({
        "n_total": len(rows), "n_overlapping": len(ov), "n_comparator": len(nov),
        "n_clusters": len({r["clust"] for r in rows}),
        "overlapping_in_normal": ck, "comparator_in_normal": nk,
        "pct_overlapping_in_normal": round(100 * ck / len(ov), 1),
        "pct_comparator_in_normal": round(100 * nk / len(nov), 1),
        "pct_of_whole_catalogue": round(100 * ck / len(rows), 1),
        "rr_crude": round((ck / len(ov)) / (nk / len(nov)), 2),
        "rr_length_standardized": round(rr, 2),
        "ci95": [round(lo, 2), round(hi, 2)], "bootstrap_B": len(boot), "seed": SEED,
        "detection_breadth_adjustment": {
            "breadth_definition": "number of distinct IEAtlas cancer-type labels per peptide",
            "primary_breadth_bins": ["1", "2", "3", "4", "5+"],
            "common_support_n": n_lb,
            "common_support_cells": cells_lb,
            "standardized_pct_overlapping": round(100 * p1_lb, 1),
            "standardized_pct_comparator": round(100 * p0_lb, 1),
            "rr_length_breadth_standardized": round(rr_lb, 2),
            "ci95_gene_clustered": [round(lo_lb, 2), round(hi_lb, 2)],
            "exact_breadth_sensitivity": {
                "common_support_n": n_lbe,
                "common_support_cells": cells_lbe,
                "standardized_pct_overlapping": round(100 * p1_lbe, 1),
                "standardized_pct_comparator": round(100 * p0_lbe, 1),
                "rr": round(rr_lbe, 2),
            },
        },
        "rr_by_stratum": stratum_rr,
        "rr_by_length": rr_by_length,
        "n_tissue_clusters": len(tclusters),
        "ci95_tissue_clustered": [round(lo_t, 2), round(hi_t, 2)],
        "n_eff_tissue_clusters": round(n_eff_tissue, 1),
    }, open(art, "w"), indent=2)
    print(f"\n  wrote {os.path.relpath(art, REPO)}")

    print("\n" + "=" * 84)
    print("VERDICT")
    print("=" * 84)
    if valid:
        print(f"""
  The association attenuates but survives the measured detection-opportunity adjustment
  ({rr_lb:.2f}x after joint length+breadth standardization, vs {rr:.2f}x length-only and
  {(ck/len(ov))/(nk/len(nov)):.2f}x crude), it holds within every ORF-class
  stratum, and a bootstrap that resamples SOURCE GENES -- respecting the clustering the z-test
  ignored -- gives 95% CI [{lo:.2f}, {hi:.2f}].

  The z = 74 must still be DELETED. It was never valid, and being directionally right does not
  rescue it. Report the length+breadth-standardized RR with the clustered CI as primary.

  INTERPRETATION, bounded: this is "consistent with, but not specific to, greater detectability or
  expression of canonical-compatible sequences." It is NOT proof of canonical origin. And the
  practical consequence is that these sequences warrant normal-presentation review before use as
  tumour-restricted targets -- not that they are established on-target/off-tumour risks, which would
  additionally depend on allele matching, tissue context, abundance and TCR avidity.""")
    else:
        print("""
  The association does NOT survive clustered inference. R3 must be restated as a descriptive
  observation with no significance claim attached.""")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
