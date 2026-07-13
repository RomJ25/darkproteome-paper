"""Are the ncORF CLASS LABELS trustworthy enough to build a class-structured claim on?

The reviewer's objection, and it is serious: the class-structured result (pseudogene-derived
claims carry a 57.7% canonical-substring rate; altORF/lncRNA claims carry ~0%) is now a PRINCIPAL
CONCLUSION -- and the class labels come from source-supplied biotype strings. That is the same
kind of annotation-derived field that produced the cancer-testis-control fiasco, where a
`gene_type == protein_coding` test was described in prose as "a curated panel of known-real
antigens" and was nothing of the sort.

So: can the class label be corroborated INDEPENDENTLY of the annotation that supplied it?

For pseudogene-derived claims there is a test that does not rely on the biotype label at all. A
processed pseudogene is a retro-copy of a parent gene, so a peptide encoded by a pseudogene ORF can
be an exact substring of the PARENT PROTEIN. If canonical-sequence compatibility is structured by
that homology, then for a pseudogene-labelled claim whose peptide matches the canonical proteome,
the matched canonical protein should be the pseudogene's OWN PARENT -- not an unrelated protein.

  RPS3AP12  (pseudogene)  ->  parent RPS3A  ->  peptide should match RPS3A, specifically.

WHAT THIS DOES AND DOES NOT SHOW -- read before quoting any number below.

  IT DOES show that canonical-sequence compatibility is CONCENTRATED IN THE ANNOTATED PARENT: the
  ambiguity sits exactly where sequence homology predicts it.

  IT DOES NOT resolve provenance. `DEVAFRKF` is encoded by BOTH RPS3AP12 and RPS3A. MS identifies
  the peptide SEQUENCE; it never chose a locus. A "parent hit" is therefore not evidence that the
  canonical parent produced the peptide, nor that the pseudogene did not. It is evidence that the
  two cannot be told apart by sequence.

  THE NULL MUST RESPECT THE CONDITIONING. A uniform "1 gene in ~20,400" null is INVALID and is not
  used here: it ignores the selection (we look only at peptides already known to match SOMETHING
  canonical), protein lengths, paralogy, shared k-mers, and -- decisively -- the fact that the
  pseudogene's parent annotation is ITSELF derived from homology. Instead we permute only the
  pseudogene->parent PAIRING while holding every peptide's canonical hit set FIXED. That breaks the
  pseudogene-parent link and changes nothing else.

Nothing here reads `orf_class`; the parent is derived from the gene SYMBOL by stripping the
pseudogene suffix, and the match is computed against SwissProt directly.

    python3 scripts/parent_gene_test.py
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SPROT = os.path.join(REPO, "data", "external", "swissprot_human.fasta")
TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
REAL = os.path.join(REPO, "data", "claim_catalog_real.csv")

csv.field_size_limit(10_000_000)

# A pseudogene symbol is conventionally PARENT + "P" + optional number:  RPS3AP12 -> RPS3A
# Also handles the bare form (RPS3AP -> RPS3A).
PSEUDO_RE = re.compile(r"^(?P<parent>[A-Z0-9\-]+?)P\d*$")


def load_swissprot():
    """gene symbol -> set of protein sequences, and a flat list for scanning."""
    by_gene = defaultdict(list)
    seqs = []
    gene = None
    cur = []
    with open(SPROT, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur and gene:
                    s = "".join(cur)
                    by_gene[gene].append(s)
                    seqs.append((gene, s))
                m = re.search(r"\bGN=([^\s]+)", line)
                gene = m.group(1).upper() if m else None
                cur = []
            else:
                cur.append(line.strip())
    if cur and gene:
        s = "".join(cur)
        by_gene[gene].append(s)
        seqs.append((gene, s))
    return by_gene, seqs


def putative_parent(symbol):
    """Derive the parent gene symbol from a pseudogene symbol. None if not derivable."""
    s = (symbol or "").strip().upper()
    if not s or not s[0].isalpha():
        return None                      # ENSG-style / clone-style names carry no parent info
    m = PSEUDO_RE.match(s)
    if not m:
        return None
    p = m.group("parent")
    return p if len(p) >= 2 else None


def main():
    for p in (SPROT, TIER1, REAL):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    by_gene, seqs = load_swissprot()
    print(f"SwissProt: {len(seqs):,} reviewed human proteins, {len(by_gene):,} gene symbols\n")

    # peptide -> its source gene symbol, for HCC pseudogene-labelled claims
    pep_gene = {}
    with open(REAL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("orf_class") != "pseudogene-ORF":
                continue
            p = (r.get("peptide_sequence") or "").strip().upper()
            if p and p.isalpha():
                pep_gene.setdefault(p, r.get("orf_id_or_locus") or "")

    # which of them are canonical substrings (from the non-novelty floor artifact)
    canon = set()
    with open(TIER1, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["orf_class"] == "pseudogene-ORF" and r["canonical_self_exact"] == "1":
                canon.add(r["peptide"])

    tested = [(p, pep_gene[p]) for p in sorted(canon) if p in pep_gene]
    print("=" * 78)
    print("THE PARENT-GENE TEST")
    print("=" * 78)
    print(f"  pseudogene-labelled peptides that ARE canonical substrings : {len(tested)}")

    derivable = [(p, g, putative_parent(g)) for p, g in tested]
    have_parent = [(p, g, par) for p, g, par in derivable if par]
    no_parent = [(p, g) for p, g, par in derivable if not par]
    print(f"    of which the symbol yields a putative parent            : {len(have_parent)}")
    print(f"    symbol carries NO parent information (clone/ENSG-style)  : {len(no_parent)}")
    if no_parent:
        print(f"      e.g. {[g for _, g in no_parent[:6]]}")

    if not have_parent:
        print("\n  Cannot run the test: no pseudogene symbol in the matched set names a parent.")
        return

    # CONTROL, and it decides how the negatives may be read: is the derived parent even IN
    # SwissProt? `FAM8A6P` yields a putative parent `FAM8A6`, which is not a reviewed protein --
    # the real parent is FAM8A1. A peptide that fails to land in a gene THAT DOES NOT EXIST in the
    # reference is UNTESTABLE, not a refutation of the label. Scoring it as a miss would silently
    # invent evidence against the class exactly the way a structural zero invents evidence for it.
    # SYMBOL DRIFT is real and it makes the naive test UNDERCOUNT parent hits:
    #   PHB -> PHB1,  FAM8A6 -> FAM8A1,  MKRN4 -> MKRN1
    # FAM8A6P's peptide landing in FAM8A1 is parent retention under a RENAMED symbol, not a miss.
    # A hit gene counts as the parent if either symbol is a prefix of the other and the remainder is
    # purely numeric (FAM8A6/FAM8A1 share FAM8A; PHB/PHB1 share PHB) -- deliberately conservative:
    # it requires a shared root plus a numeric suffix, not mere similarity.
    #
    # ORDER MATTERS HERE, and getting it wrong silently zeroed this bucket: an earlier version
    # binned "parent absent from SwissProt" as UNTESTABLE *before* testing for drift, so exactly the
    # cases drift exists to rescue never reached the check. Resolve the hit FIRST; only call a claim
    # untestable if it neither hits the parent, nor a drift-variant, nor has a parent in the
    # reference at all.
    def same_gene(a, b):
        if a == b:
            return True
        root = os.path.commonprefix([a, b])
        return len(root) >= 3 and (a[len(root):].isdigit() or b[len(root):].isdigit())

    hit_parent, hit_drift, hit_other, untestable = [], [], [], []
    for pep, gsym, parent in have_parent:
        genes = {g for g, s in seqs if pep in s}
        if parent in genes:
            hit_parent.append((pep, gsym, parent, genes))
        elif any(same_gene(parent, g) for g in genes):
            hit_drift.append((pep, gsym, parent, genes))
        elif parent not in by_gene:
            untestable.append((pep, gsym, parent))   # parent isn't in the reference at all
        else:
            hit_other.append((pep, gsym, parent, genes))
    testable = hit_parent + hit_drift + hit_other

    print(f"\n  CONTROL — is the derived parent a reviewed protein at all?")
    print(f"    parent EXISTS in SwissProt (the test can run)  : {len(testable):>4}")
    print(f"    parent ABSENT from SwissProt (UNTESTABLE)      : {len(untestable):>4}")
    if untestable:
        print(f"      e.g. {[(g, p) for _, g, p in untestable[:5]]}")
        print(f"      -> these are NOT misses. A peptide cannot land in a protein that the")
        print(f"         reference does not contain. Excluded from the denominator.")

    n = len(testable)
    if not n:
        print("\n  No testable claims.")
        return
    par_all = len(hit_parent) + len(hit_drift)
    print(f"\n  Of the {n} TESTABLE peptides, the canonical protein containing the peptide is:")
    print(f"    the pseudogene's OWN PARENT, exact symbol : {len(hit_parent):>4}  "
          f"({100*len(hit_parent)/n:.1f}%)")
    print(f"    the parent under a DRIFTED symbol         : {len(hit_drift):>4}  "
          f"({100*len(hit_drift)/n:.1f}%)   e.g. FAM8A6P -> FAM8A1")
    print(f"    ------------------------------------------------------")
    print(f"    PARENT RETENTION, total                   : {par_all:>4}  "
          f"({100*par_all/n:.1f}%)")
    print(f"    some genuinely OTHER canonical gene       : {len(hit_other):>4}  "
          f"({100*len(hit_other)/n:.1f}%)")

    # Are the "other" hits just short-peptide coincidences? Length is the discriminator.
    import statistics as st
    lp = [len(p) for p, _, _, _ in hit_parent + hit_drift]
    lo = [len(p) for p, _, _, _ in hit_other]
    if lp and lo:
        print(f"\n  Peptide LENGTH — a short peptide can land in an unrelated protein by chance,")
        print(f"  so length is the discriminator between retention and coincidence:")
        print(f"    parent-hits : n={len(lp):>3}  median {st.median(lp):.0f}aa  "
              f"({sum(1 for x in lp if x <= 9)} are <=9aa)")
        print(f"    other-hits  : n={len(lo):>3}  median {st.median(lo):.0f}aa  "
              f"({sum(1 for x in lo if x <= 9)} are <=9aa)")

    # GENE-level: does each pseudogene retain its parent at all?
    by_pseudo = defaultdict(lambda: {"parent": 0, "other": 0})
    for _, g, p, _ in hit_parent + hit_drift:
        by_pseudo[g]["parent"] += 1
    for _, g, p, _ in hit_other:
        by_pseudo[g]["other"] += 1
    genes_any = sum(1 for v in by_pseudo.values() if v["parent"])
    print(f"\n  GENE-level (the claim the paper actually makes is about the CLASS, not the peptide):")
    print(f"    pseudogenes with >=1 canonical-matching peptide      : {len(by_pseudo):>4}")
    print(f"    ... of which >=1 peptide lands in their OWN parent   : {genes_any:>4}  "
          f"({100*genes_any/max(1,len(by_pseudo)):.1f}%)")

    print("\n  examples where the peptide lands in the pseudogene's own parent:")
    for pep, gsym, parent, genes in hit_parent[:10]:
        others = sorted(genes - {parent})
        extra = f"   (also: {others[:3]})" if others else ""
        print(f"    {gsym:<14} -> parent {parent:<10} peptide {pep:<14} FOUND IN {parent}{extra}")

    if hit_other:
        print("\n  examples where it lands somewhere else (the label's weak cases):")
        for pep, gsym, parent, genes in hit_other[:8]:
            print(f"    {gsym:<14} -> expected {parent:<10} but found in {sorted(genes)[:3]}")

    # --- IS A "PARENT HIT" EVEN A UNIQUE ASSIGNMENT? (it usually is not) ---
    only_par = sum(1 for p, g, par, gs in hit_parent + hit_drift if len(gs) == 1)
    multi = len(hit_parent) + len(hit_drift) - only_par
    tot_par = len(hit_parent) + len(hit_drift)
    print("\n" + "=" * 78)
    print("A 'PARENT HIT' IS NOT A UNIQUE SOURCE ASSIGNMENT")
    print("=" * 78)
    print(f"  of the {tot_par} peptides that land in their own parent:")
    print(f"    match ONLY the parent gene                : {only_par:>4}  ({100*only_par/tot_par:.1f}%)")
    print(f"    match the parent AND further canonical genes: {multi:>4}  ({100*multi/tot_par:.1f}%)")
    print("  So even where the parent is hit, the SEQUENCE is compatible with the parent, the")
    print("  pseudogene ORF, and often further proteins besides. MS identifies the sequence; it")
    print("  never chose a locus.")

    # --- CONDITIONAL PERMUTATION NULL ---
    import random
    parents = [par for _, _, par, _ in hit_parent + hit_drift + hit_other]
    hitsets = [gs for _, _, _, gs in hit_parent + hit_drift + hit_other]
    obs = tot_par
    random.seed(0)
    null = []
    for _ in range(10000):
        sh = parents[:]
        random.shuffle(sh)
        null.append(sum(1 for gs, par in zip(hitsets, sh)
                        if par in gs or any(same_gene(par, g) for g in gs)))
    mean = sum(null) / len(null)
    ge = sum(1 for x in null if x >= obs)
    print("\n" + "=" * 78)
    print("THE NULL — conditional permutation (hit sets FIXED, parent pairing shuffled)")
    print("=" * 78)
    print(f"  observed  : {obs}/{n}")
    print(f"  null mean : {mean:.1f}/{n} = {100*mean/n:.1f}%     null max: {max(null)}")
    print(f"  p         : {(ge+1)/(len(null)+1):.1e}   (10,000 permutations; {ge} >= observed)")
    print("  This null preserves the selection, every peptide's canonical hit set, and the pool of")
    print("  parents. It breaks ONLY the pseudogene-parent pairing. The naive 1/20,213 null is")
    print("  DELETED -- it ignored all of the above, and the parent annotation is itself homology-")
    print("  derived, so a uniform null was never the right comparison.")

    print("\n" + "=" * 78)
    print("WHAT THIS DOES AND DOES NOT LICENSE")
    print("=" * 78)
    print("  DOES: canonical-sequence compatibility is CONCENTRATED IN THE ANNOTATED PARENT --")
    print("        the ambiguity sits exactly where sequence homology predicts. The pseudogene")
    print("        class label is not arbitrary: it is corroborated by sequence structure.")
    print("  DOES NOT: resolve provenance. It does NOT show the canonical parent produced the")
    print("        peptide, nor that the pseudogene did not. It shows the two CANNOT BE TOLD")
    print("        APART BY SEQUENCE. It also says nothing about lncRNA-ORF or altORF labels, for")
    print("        which no analogous test exists, nor about the clone-style symbols excluded above.")
    print(f"  THE {len(hit_other)} NON-PARENT MATCHES ARE AN UNRESOLVED RESIDUAL. We do not explain them.")
    print("        (Out-of-frame retro-copy ORFs and incidental short-peptide matches are both")
    print("        plausible; parent-hits and other-hits share a median length of 9aa, so the one")
    print("        test we ran does not separate them. Reported, not explained.)")


if __name__ == "__main__":
    main()
