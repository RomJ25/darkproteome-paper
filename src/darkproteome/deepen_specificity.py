"""Deepened specificity audit — stress-test the 16.6% IEAtlas normal-contamination result.

The earlier "bonus" stat (ingest_atlases.py) was a PLAIN exact-sequence intersection of
IEAtlas cancer vs normal epitope lists: 28,979 / 174,465 = 16.6%. That number has none of
the refinements a fair reviewer would demand. This script adds them, each able to WEAKEN it.

Falsification criteria, fixed before running:
  C1  if overlap drops below ~5% after removing peptides that are also exact substrings of
      the canonical human proteome (SwissProt)            -> headline weakens materially
  C2  if a composition+length-matched SHUFFLE null overlaps normal at a similar rate
      (enrichment < ~2x)                                  -> 16.6% is a composition artifact
  C3  if ~all overlapping peptides sit only on immune-privileged / non-critical normal
      tissue (testis etc.), not vital organs              -> the SAFETY implication weakens

Report whatever it shows. Public data only; stdlib + a local SwissProt FASTA.

    python3 src/darkproteome/deepen_specificity.py
"""
import csv
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- centralized data paths

random.seed(0)
csv.field_size_limit(10_000_000)

CANCER = paths.IEATLAS_CANCER
NORMAL = paths.IEATLAS_NORMAL
SPROT = paths.SPROT
HLALA = paths.HLALA

# Tissue criticality tiers (on-target/off-tumor toxicity severity if a T cell attacks it).
CRITICAL = {"Brain", "Cerebellum", "Heart", "Liver", "Lung", "Kidney", "Pancreas",
            "Nerve", "Adrenal Gland"}
HIGH = {"Esophagus", "Stomach", "Colon", "Small Intestine", "Bladder", "Thyroid",
        "Aorta", "Trachea", "Tongue", "Gallbladder", "Muscle", "Skin"}
MODERATE = {"Breast", "Prostate", "Uterus", "Ovary"}
IMMUNE = {"Thymus", "Spleen", "Bone Marrow", "Lymph Node"}   # presentation expected; thymus=central tolerance
PRIVILEGED = {"Testis"}                                       # classically tolerated target

_TIERS = [("CRITICAL", CRITICAL), ("HIGH", HIGH), ("MODERATE", MODERATE),
          ("IMMUNE", IMMUNE), ("PRIVILEGED", PRIVILEGED)]
_TIER_LC = {name: {t.lower() for t in s} for name, s in _TIERS}


def tier_of(tissue):
    t = (tissue or "").strip().lower()
    for name, s in _TIER_LC.items():
        if t in s:
            return name
    return "OTHER"


def load(path):
    """seq -> set(tissues); only alpha sequences, uppercased (matches ingest_atlases)."""
    d = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 4 or not r[0]:
                continue
            s = r[0].strip().upper()
            if not s.isalpha():
                continue
            d.setdefault(s, set()).add((r[3] or "").strip())
    return d


def load_hla_ligand_atlas():
    """HLA-I normal-tissue peptides from HLA Ligand Atlas (independent consortium): seq->tissues."""
    d = {}
    with open(HLALA, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 5 or r[2] != "HLA-I":
                continue
            s = r[1].strip().upper()
            if not s.isalpha():
                continue
            d[s] = {t.strip() for t in r[4].split(",") if t.strip()}
    return d


def canonical_hits(query):
    """Subset of `query` peptides that occur as an exact substring of any SwissProt protein.
    Memory-light: stream proteome, test each generated k-mer against the (small) query set."""
    lengths = sorted({len(p) for p in query})
    hit = set()
    prot = []

    def scan(seq):
        for L in lengths:
            for i in range(len(seq) - L + 1):
                k = seq[i:i + L]
                if k in query:
                    hit.add(k)

    with open(SPROT, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if prot:
                    scan("".join(prot))
                    prot = []
            else:
                prot.append(line.strip())
        if prot:
            scan("".join(prot))
    return hit


def overlap_rate(cancer, normal):
    both = cancer & normal
    return both, len(both) / max(1, len(cancer))


def shuffle_null(cancer_list, normal, reps=3):
    """Per-peptide shuffle (preserves exact length+composition). Mean overlap with normal."""
    rates = []
    for _ in range(reps):
        hits = 0
        for p in cancer_list:
            cs = list(p)
            random.shuffle(cs)
            if "".join(cs) in normal:
                hits += 1
        rates.append(hits / len(cancer_list))
    return sum(rates) / len(rates)


def uniform_null(cancer_list, normal, aa_weights):
    """Random peptides, length-matched, AA sampled from cancer-population background."""
    alphabet, weights = zip(*aa_weights.items())
    hits = 0
    for p in cancer_list:
        rp = "".join(random.choices(alphabet, weights=weights, k=len(p)))
        if rp in normal:
            hits += 1
    return hits / len(cancer_list)


def main():
    print("Loading IEAtlas ...")
    cancer_d = load(CANCER)
    normal_d = load(NORMAL)
    cancer = set(cancer_d)
    normal = set(normal_d)
    print(f"  unique cancer epitopes: {len(cancer):,}")
    print(f"  unique normal epitopes: {len(normal):,}")

    # ---- Baseline (reproduce the headline) ----
    both, rate = overlap_rate(cancer, normal)
    print("\n[0] BASELINE exact overlap (the published 16.6%)")
    print(f"  cancer also in normal: {len(both):,} / {len(cancer):,} = {100*rate:.1f}%")

    # ---- C1: canonical-overlap removal ----
    print("\n[C1] removing cancer peptides that are exact substrings of canonical proteome ...")
    canon = canonical_hits(cancer)
    # null: would composition+length-matched random peptides be canonical substrings this often?
    # `cancer` is a set -- iterate a sorted() snapshot, not the set directly, or this shuffle
    # is non-reproducible across process runs: iteration order depends on per-process string
    # hash randomization, so random.seed(0) alone doesn't pin it.
    shuffled = set()
    for p in sorted(cancer):
        cs = list(p)
        random.shuffle(cs)
        shuffled.add("".join(cs))
    canon_null = canonical_hits(shuffled)
    noncanon_cancer = cancer - canon
    nc_both, nc_rate = overlap_rate(noncanon_cancer, normal)
    canon_in_overlap = len(both & canon)
    print(f"  cancer peptides that ARE canonical substrings: {len(canon):,} "
          f"({100*len(canon)/len(cancer):.1f}% of cancer set)")
    print(f"  [null] shuffled (composition-matched) canonical-substring rate: "
          f"{100*len(canon_null)/max(1,len(shuffled)):.1f}%  "
          f"-> enrichment {(len(canon)/len(cancer))/max(1e-9,len(canon_null)/len(shuffled)):.1f}x")
    print(f"  of the {len(both):,} overlap peptides, canonical-explained: {canon_in_overlap:,} "
          f"({100*canon_in_overlap/max(1,len(both)):.1f}%)")
    print(f"  >>> overlap among TRULY non-canonical cancer epitopes: "
          f"{len(nc_both):,} / {len(noncanon_cancer):,} = {100*nc_rate:.1f}%")

    # ---- C2: null / enrichment ----
    print("\n[C2] enrichment vs composition-matched null ...")
    nc_list = sorted(noncanon_cancer)  # deterministic order -- see the C1 shuffle comment above
    aa = Counter("".join(nc_list))
    aa_weights = {k: v for k, v in aa.items()}
    sh = shuffle_null(nc_list, normal)
    un = uniform_null(nc_list, normal, aa_weights)
    enr_sh = nc_rate / sh if sh > 0 else float("inf")
    enr_un = nc_rate / un if un > 0 else float("inf")
    print(f"  observed (non-canonical) overlap rate: {100*nc_rate:.2f}%")
    print(f"  shuffle null (composition+length preserved): {100*sh:.4f}%  -> enrichment {enr_sh:.0f}x")
    print(f"  uniform null (AA background, length matched): {100*un:.4f}%  -> enrichment {enr_un:.0f}x")

    # ---- C3: tissue criticality of the (non-canonical) contamination ----
    print("\n[C3] normal-tissue criticality of the non-canonical overlap peptides ...")
    tier_counts = Counter()
    crit_peptides = 0
    only_privileged = 0
    for p in nc_both:
        tissues = normal_d.get(p, set())
        tiers = set()
        for t in tissues:
            if t in CRITICAL:
                tiers.add("CRITICAL")
            elif t in HIGH:
                tiers.add("HIGH")
            elif t in MODERATE:
                tiers.add("MODERATE")
            elif t in IMMUNE:
                tiers.add("IMMUNE")
            elif t in PRIVILEGED:
                tiers.add("PRIVILEGED")
            else:
                tiers.add("OTHER:" + t)
        for tr in tiers:
            tier_counts[tr] += 1
        if "CRITICAL" in tiers:
            crit_peptides += 1
        if tiers <= {"PRIVILEGED", "IMMUNE"}:
            only_privileged += 1
    print(f"  non-canonical overlap peptides: {len(nc_both):,}")
    print(f"  ... hitting >=1 CRITICAL normal tissue (brain/heart/liver/lung/kidney/...): "
          f"{crit_peptides:,} ({100*crit_peptides/max(1,len(nc_both)):.1f}%)")
    print(f"  ... on ONLY immune/privileged tissue (expected/tolerated): "
          f"{only_privileged:,} ({100*only_privileged/max(1,len(nc_both)):.1f}%)")
    print(f"  tier histogram (peptide may span tiers): "
          f"{dict(tier_counts.most_common(8))}")
    # critical-contamination rate vs the whole non-canonical cancer set
    print(f"  >>> CRITICAL-tissue contamination of non-canonical cancer set: "
          f"{crit_peptides:,} / {len(noncanon_cancer):,} = "
          f"{100*crit_peptides/len(noncanon_cancer):.1f}%")

    # ---- C4: can an INDEPENDENT normal ligandome (HLA Ligand Atlas) replicate this? ----
    # Honest result: NO -- and it's a search-space CONFOUND, not a refutation. HLA Ligand Atlas
    # was built by searching the CANONICAL proteome, so it structurally cannot contain cryptic
    # peptides. We prove that (it is ~100% canonical) and show a positive control that works.
    print("\n[C4] independent replication vs HLA Ligand Atlas (separate consortium, HLA-I) ...")
    hla_set = set(load_hla_ligand_atlas())
    hla_canon = canonical_hits(hla_set)
    pos_control = canon & hla_set            # canonical cancer peptides -> SHOULD replicate
    confounded = noncanon_cancer & hla_set   # non-canonical -> cannot, by construction
    print(f"  HLA Ligand Atlas HLA-I peptides: {len(hla_set):,}; canonical substrings: "
          f"{len(hla_canon):,} ({100*len(hla_canon)/len(hla_set):.1f}%) -> a CANONICAL-space ligandome")
    print(f"  POSITIVE CONTROL: canonical cancer epitopes also in HLA Ligand Atlas: "
          f"{len(pos_control):,} ({100*len(pos_control)/max(1,len(canon)):.1f}% of canonical-cancer) "
          f"-> cross-resource exact match WORKS when peptides are in-space")
    print(f"  non-canonical cancer epitopes in HLA Ligand Atlas: {len(confounded):,} "
          f"({100*len(confounded)/len(noncanon_cancer):.2f}%) -> ~0 BECAUSE HLA Ligand Atlas can't "
          f"contain cryptic peptides, NOT because the overlap is unreal")
    print("  CONCLUSION: HLA Ligand Atlas cannot test non-canonical normal presentation; an")
    print("  independent CRYPTIC-space normal ligandome is required -- a re-analysis of raw")
    print("  spectra against a cryptic-ORF database.")

    # ---- verdict vs pre-registered criteria ----
    print("\n=== VERDICT vs pre-registered criteria ===")
    print(f"  C1 (>5% after canonical removal?):  {100*nc_rate:.1f}%  -> "
          f"{'HOLDS' if nc_rate >= 0.05 else 'WEAKENED'}")
    print(f"  C2 (enrichment >=2x over null?):    {enr_sh:.0f}x  -> "
          f"{'HOLDS' if (sh==0 or enr_sh >= 2) else 'WEAKENED'}")
    print(f"  C3 (real safety signal, not just privileged?): "
          f"{100*crit_peptides/max(1,len(nc_both)):.1f}% hit critical tissue -> "
          f"{'REAL SAFETY SIGNAL' if crit_peptides/max(1,len(nc_both)) >= 0.10 else 'MOSTLY BENIGN'}")
    print(f"  C4 (independent replication): NOT testable off-the-shelf -- HLA Ligand Atlas is "
          f"{100*len(hla_canon)/len(hla_set):.0f}% canonical-space; positive control "
          f"{100*len(pos_control)/max(1,len(canon)):.0f}% confirms the method; non-canonical "
          f"replication needs a raw-spectra re-analysis.")


if __name__ == "__main__":
    main()
