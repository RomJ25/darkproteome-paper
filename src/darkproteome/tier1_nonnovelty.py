"""The non-novelty floor -- the UNCONDITIONAL dark-proteome indictment.

A 'cryptic' peptide that is an exact substring of the canonical proteome is
MS-UNFALSIFIABLE as non-canonical: identical sequence -> identical spectrum, so
by parsimony it must be attributed canonically. This floor needs NO decoys, NO
global FDR, NO null model -- pure sequence string-matching -- so it survives even
a DISHONEST global FDR (see the manuscript's Methods for the full identifiability
derivation).

Two floors per cohort:
  * EXACT canonical-substring          -> conservative floor
  * I/L-collapsed canonical-substring  -> MS-REALISTIC floor. Leu/Ile are isobaric
    (113.084 Da) and indistinguishable on MS. (We collapse ONLY I/L -- on the
    Orbitrap these studies used, K/Q is mass-resolvable, so we do NOT over-collapse.)

Empirically, careful PRIMARY studies (Raja ovarian; HCC) come back ~0% -- they already
exclude canonical matches -- while the catalog-wide rate is 54.4%, driven by AGGREGATOR
atlases (IEAtlas/CPDB). So this floor indicts sloppy aggregator DATABASES; primary studies
still require the manuscript's budget bound to rule out a subtler per-study FDR failure.

    python3 src/darkproteome/tier1_nonnovelty.py

Public data only; stdlib + the local SwissProt FASTA used by reference_model.py.
"""
import csv
import math
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- centralized data paths
import deepen_specificity as d  # SPROT path, load() for IEAtlas normal

random.seed(0)
csv.field_size_limit(10_000_000)

REAL = os.path.join(paths.REPO, "data", "claim_catalog_real.csv")
SCORED = os.path.join(paths.REPO, "data", "claim_catalog_scored.csv")
OUT = os.path.join(paths.REPO, "data", "primary_tier1_nonnovelty.csv")

# label -> citation DOI (the two careful PRIMARY immunopeptidomics studies)
COHORTS = {
    "Raja ovarian (ads7405)": "10.1126/sciadv.ads7405",
    "HCC / Camarena-Albà (adn3628)": "10.1126/sciadv.adn3628",
}


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / den
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def canonical_hits_both(query, sprot=d.SPROT):
    """One streaming SwissProt pass. Returns (exact_hits, il_hits): subsets of
    `query` occurring as an exact / I-L-collapsed substring of any protein.
    exact_hits is a subset of il_hits by construction."""
    if not query:
        return set(), set()
    il_to_orig = {}
    for p in query:
        il_to_orig.setdefault(p.replace("I", "L"), []).append(p)
    lengths = sorted({len(p) for p in query})
    exact, il, prot = set(), set(), []

    def scan(seq):
        sil = seq.replace("I", "L")
        n = len(seq)
        for L in lengths:
            for i in range(n - L + 1):
                if seq[i:i + L] in query:
                    exact.add(seq[i:i + L])
                o = il_to_orig.get(sil[i:i + L])
                if o:
                    il.update(o)

    with open(sprot, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if prot:
                    scan("".join(prot))
                    prot = []
            else:
                prot.append(line.strip())
        if prot:
            scan("".join(prot))
    return exact, il


def load_cohorts(real=REAL):
    """label -> (set seqs with REAL sequences, n_gene_level_dropped, dict seq->orf_class)."""
    seqs = {lab: set() for lab in COHORTS}
    dropped = {lab: 0 for lab in COHORTS}
    klass = {lab: {} for lab in COHORTS}
    doi2lab = {doi: lab for lab, doi in COHORTS.items()}
    with open(real, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("_canonical") == "yes":
                continue
            lab = doi2lab.get(r.get("citation_doi_pmid"))
            if not lab:
                continue
            s = (r.get("peptide_sequence") or "").strip().upper()
            if s and s.isalpha():
                seqs[lab].add(s)
                klass[lab].setdefault(s, r.get("orf_class"))
            else:
                dropped[lab] += 1
    return {lab: (seqs[lab], dropped[lab], klass[lab]) for lab in COHORTS}


def band(seqs, lo=8, hi=11):
    return {p for p in seqs if lo <= len(p) <= hi}


def audit(label, seqs, normal, seq2class):
    n = len(seqs)
    print(f"\n=== {label}  (N={n} unique peptide-level cryptic seqs) ===")
    print(f"  length dist: {dict(sorted(Counter(len(p) for p in seqs).items()))}")
    exact, il = canonical_hits_both(seqs)
    # composition-matched null. `seqs` is a set -- its iteration order depends on
    # per-process string hash randomization, so random.seed(0) alone does NOT make this
    # shuffle reproducible across runs unless the loop order is pinned first.
    sh = set()
    for p in sorted(seqs):
        cs = list(p)
        random.shuffle(cs)
        sh.add("".join(cs))
    sh_exact, _ = canonical_hits_both(sh)
    null_rate = len(sh_exact) / max(1, len(sh))

    ke, ki = len(exact), len(il)
    pe, le, ue = wilson(ke, n)
    pi, li, ui = wilson(ki, n)
    enr = (pe / null_rate) if null_rate > 0 else float("inf")
    print(f"  EXACT canonical-substring (conservative floor): {ke}/{n} = {100*pe:.1f}%  [95% CI {100*le:.1f}-{100*ue:.1f}]")
    print(f"  I/L-collapsed (MS-realistic floor):             {ki}/{n} = {100*pi:.1f}%  [95% CI {100*li:.1f}-{100*ui:.1f}]")
    # class-resolved: the floor is concentrated in pseudogene-ORFs (parent-gene substrings)
    cls_tot, cls_self = Counter(), Counter()
    for p in seqs:
        c = seq2class.get(p, "?")
        cls_tot[c] += 1
        if p in exact:
            cls_self[c] += 1
    print("  canonical-self BY orf_class (the decisive cut):")
    for c in sorted(cls_tot, key=lambda c: -cls_tot[c]):
        print(f"      {c:16s}: {cls_self[c]:3d}/{cls_tot[c]:4d} = {100*cls_self[c]/max(1,cls_tot[c]):5.1f}%")
    # length-controlled: 8-11mers only (matches the HLA-I regime where chance is non-trivial)
    b = band(seqs)
    be = len(b & exact)
    print(f"  ...restricted to 8-11mers (chance-meaningful): {be}/{len(b)} = "
          f"{100*be/max(1,len(b)):.1f}%   (length-controlled check)")
    print(f"  [null] composition-matched shuffle rate: {100*null_rate:.1f}%  -> enrichment {enr:.1f}x")
    if normal is not None:
        kn = len(seqs & normal)
        print(f"  also exact-listed on normal tissue (IEAtlas, secondary): {kn}/{n} = {100*kn/max(1,n):.1f}%")
    return dict(label=label, n=n, exact=ke, il=ki, pe=pe, pi=pi, exact_set=exact, il_set=il,
                normal=(len(seqs & normal) if normal is not None else None))


def main():
    cohorts = load_cohorts()
    print("PRIMARY-STUDY non-novelty audit (peptide-level cryptic sequences):")
    for lab, (seqs, dropped, _k) in cohorts.items():
        print(f"  {lab}: {len(seqs)} peptide-level seqs  (+{dropped} gene/locus-level rows with no peptide, dropped)")

    try:
        normal = set(d.load(d.NORMAL))
    except Exception:
        normal = None

    results = {}
    for lab, (seqs, _dropped, klass) in cohorts.items():
        results[lab] = audit(lab, seqs, normal, klass)

    # self-consistency + catalog-wide 54% reproduction
    if os.path.exists(SCORED):
        stored, tot_self, tot = {}, 0, 0
        with open(SCORED, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                stored[r["peptide"]] = r["canonical_self"] == "1"
                tot += 1
                tot_self += (r["canonical_self"] == "1")
        allseqs = [(lab, p) for lab, (seqs, _d, _k) in cohorts.items() for p in seqs]
        cov = sum(1 for lab, p in allseqs if p in stored)
        mism = sum(1 for lab, p in allseqs
                   if p in stored and stored[p] != (p in results[lab]["exact_set"]))
        print(f"\nself-consistency vs claim_catalog_scored.csv: {cov}/{len(allseqs)} present, "
              f"{mism} mismatches (expect 0)")
        print(f"catalog-wide reproduction: canonical_self = {tot_self:,}/{tot:,} = "
              f"{100*tot_self/max(1,tot):.1f}%  (the 54% = AGGREGATOR-driven; cf. primary studies above)")

    # per-peptide artifact
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "peptide", "len", "orf_class",
                    "canonical_self_exact", "canonical_self_IL", "normal_listed"])
        for lab, (seqs, _d, klass) in cohorts.items():
            es, ils = results[lab]["exact_set"], results[lab]["il_set"]
            for p in sorted(seqs):
                w.writerow([lab, p, len(p), klass.get(p, "?"), int(p in es), int(p in ils),
                            int(normal is not None and p in normal)])
    print(f"\nwrote per-peptide table -> {OUT}")

    print("\n=== CLASS-RESOLVED VERDICT (not a clean primary-vs-aggregator dichotomy) ===")
    print("  altORF / lncRNA-ORF:  ~0% canonical-self in BOTH primary studies -> genuinely non-canonical.")
    print("  pseudogene-ORF:       Raja 0/98 vs HCC 43/116 = 37% -> pseudogene antigens are intrinsically")
    print("                        MS-unfalsifiable vs their PARENT gene (e.g. RPS3AP12 -> RPS3A); also a")
    print("                        tumor-specificity risk. A careful study (Raja) excludes them; HCC did not.")
    print("  AGGREGATOR atlases:   54.4% canonical-self (indiscriminate) -- the worst.")
    print("  -> This is a CLASS-RESOLVED floor: it indicts pseudogene-ORF claims + aggregators;")
    print("     altORF/lncRNA primary claims pass. Per-study FDR still needs the manuscript's budget bound.")


if __name__ == "__main__":
    main()
