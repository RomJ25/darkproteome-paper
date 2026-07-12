"""Re-derive the headline against the reference IEAtlas ACTUALLY SEARCHED (SwissProt 2022_01).

REVIEWER OBJECTION (upheld). Our 56.3% canonical-overlap rate was computed against a MODERN human
SwissProt. IEAtlas was built against SwissProt as of FEBRUARY 2022. Sequence novelty is
reference-relative -- N_i(R) -- so a peptide we call "canonical" may have entered the reference AFTER
IEAtlas was built. If enough of the 56.3% is retrospective, the audit is anachronistic: we would be
faulting IEAtlas for failing to exclude sequences that were not yet in anyone's canonical proteome.

This settles it by rebuilding the reference at the era IEAtlas used and re-deriving the number.

    data/external/uniprot_sprot.fasta.gz -- 566,996 entries, the exact size of SwissProt release
    2022_01 (23 Feb 2022). Human subset taken by OX=9606.

The audit only stands if the rate is essentially UNCHANGED. If it collapses, the finding is an
artifact of our reference choice and the paper must say so.

    python3 scripts/era_correct_reference.py
"""
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
SPROT_2022 = os.path.join(EXT, "uniprot_sprot.fasta.gz")     # release 2022_01
SPROT_NOW = os.path.join(EXT, "swissprot_human.fasta")       # modern, used for the published 56.3%
IE = os.path.join(EXT, "atlases", "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
csv.field_size_limit(10_000_000)

PSEUDO = re.compile(r"^[A-Z0-9\-]+P\d+$")


def human_2022():
    seqs, keep, cur = [], False, []
    with gzip.open(SPROT_2022, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    seqs.append("".join(cur))
                keep = "OX=9606 " in line or line.rstrip().endswith("OX=9606")
                cur = []
            elif keep:
                cur.append(line.strip())
    if keep and cur:
        seqs.append("".join(cur))
    return seqs


def modern():
    seqs, cur = [], []
    with open(SPROT_NOW, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def main():
    for p in (SPROT_2022, SPROT_NOW, IE):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    print("Building the two references ...")
    old, new = human_2022(), modern()
    print(f"  SwissProt 2022_01, human (OX=9606) : {len(old):,} proteins   [what IEAtlas searched]")
    print(f"  SwissProt modern, human            : {len(new):,} proteins   [what we scored against]")

    # A peptide is "canonical" if it occurs as a substring of ANY canonical protein -- the same
    # exact-substring criterion used throughout, applied to each reference.
    def substr_index(seqs, k):
        idx = set()
        for s in seqs:
            for i in range(len(s) - k + 1):
                idx.add(s[i:i + k])
        return idx

    peps, genes = set(), {}
    with open(IE, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 3 or not r[0]:
                continue
            q = r[0].strip().upper()
            if q.isalpha():
                peps.add(q)
                genes.setdefault(q, (r[2] or "").split("_")[0].upper())
    print(f"  IEAtlas distinct cancer epitopes   : {len(peps):,}\n")

    # Index by every peptide length present, so the test is an exact substring test, not a 9-mer proxy.
    bylen = defaultdict(set)
    for q in peps:
        bylen[len(q)].add(q)

    res = {}
    for tag, seqs in (("2022_01", old), ("modern", new)):
        hit = set()
        for L, group in bylen.items():
            idx = substr_index(seqs, L)
            hit |= {q for q in group if q in idx}
        res[tag] = hit
        print(f"  canonical-overlapping vs {tag:<8}: {len(hit):>7,} / {len(peps):,} = "
              f"{100*len(hit)/len(peps):.1f}%")

    o, n = res["2022_01"], res["modern"]
    only_new = n - o
    print("\n" + "=" * 88)
    print("RETROSPECTIVE FRACTION -- overlaps that exist ONLY because the reference grew since 2022")
    print("=" * 88)
    print(f"  canonical in modern but NOT in 2022_01: {len(only_new):,} "
          f"({100*len(only_new)/len(n) if n else 0:.2f}% of the modern overlap set)")
    print(f"  canonical in 2022_01 but NOT in modern: {len(o - n):,}")

    delta = 100 * len(n) / len(peps) - 100 * len(o) / len(peps)
    print(f"\n  rate shift: {100*len(n)/len(peps):.1f}% (modern) -> {100*len(o)/len(peps):.1f}% "
          f"(era-correct)   [{delta:+.1f} pp]")

    ps = {q for q in peps if PSEUDO.match(genes.get(q, ""))}
    if ps:
        print(f"\n  pseudogene-ORF epitopes, era-correct: "
              f"{len(ps & o):,}/{len(ps):,} = {100*len(ps & o)/len(ps):.1f}%")

    art = os.path.join(REPO, "data", "derived_era_reference.json")
    json.dump({
        "n_peptides": len(peps),
        "proteins_2022_01": len(old), "proteins_modern": len(new),
        "overlap_2022_01": len(o), "overlap_modern": len(n),
        "pct_2022_01": round(100 * len(o) / len(peps), 1),
        "pct_modern": round(100 * len(n) / len(peps), 1),
        "retrospective_only": len(only_new),
        "retrospective_pct_of_overlap": round(100 * len(only_new) / len(n), 2),
        "shift_pp": round(delta, 1),
    }, open(art, "w"), indent=2)
    print(f"\n  wrote {os.path.relpath(art, REPO)}")

    print("\n" + "=" * 88)
    print("VERDICT")
    print("=" * 88)
    if abs(delta) < 2.0:
        print(f"""
  The objection does NOT bite. Against the reference IEAtlas actually searched, the canonical-overlap
  rate is {100*len(o)/len(peps):.1f}% -- a shift of {delta:+.1f} percentage points. The overlap was
  almost entirely knowable in February 2022; it is not an artifact of a reference that grew
  afterwards. The audit is not anachronistic, and the outlier status (vs 1.4-5%) is untouched.

  Report the era-correct number as the primary robustness check, and keep saying that novelty is
  reference-relative -- because it is, even though here it does not change the answer.""")
    else:
        print(f"""
  The objection BITES. The rate moves {delta:+.1f} pp when scored against the era-correct reference,
  so a material part of the 56.3% is retrospective. The headline must be restated against SwissProt
  2022_01 and the anachronism disclosed.""")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
