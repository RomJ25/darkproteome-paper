"""Measure the latent canonical ambiguity of the ncORF LIBRARY, not the filter.

IEAtlas sits at 56.3% canonical-sequence overlap (98,193 / 174,465 unique cancer-catalogued
sequences). ORF-class composition does not account for it: with no pseudogene ORFs at all the rate
would be 55.8%, a shift of 0.5 pp.

So: how much of the ambiguity is already present in the LIBRARY that was searched, before any filter?

NOTE -- THREE CLAIMS THIS SCRIPT'S HEADER USED TO MAKE, ALL RETRACTED IN REVIEW:
  * "IEAtlas is an 11-40x outlier" vs the catalogues Bedran et al. audited. Their 1.4-5% came from a
    DIFFERENT pipeline, reference, dedup and peptide unit; a ratio across pipelines is arithmetic,
    not measurement. Cited as context only, with no fold arithmetic. (See cross_catalogue.py.)
  * "the FDR is ruled out, because chance canonical overlap is ~0.1% under a shuffle." Wrong null
    object: a false target PSM is not a shuffled string, it is an incorrect candidate drawn from the
    ACTUAL SEARCH DATABASE. (See the retraction notice in rule_predicts_rate.py.)
  * "IEAtlas's non-pseudogene ORFs are themselves at 56.0%" -- computed over a denominator that did
    not partition (546 peptides carry both a pseudogene and a non-pseudogene ORF label and were
    double-counted). The correct mutually-exclusive figure is 55.8%.

AND WHAT THIS SCRIPT'S RESULT DOES *NOT* LICENSE: nuORFdb's 34.1% is NOT a lower bound on IEAtlas's
combined library (nuORFdb + RPFdb + Translnc). |(A u B) n C| / |A u B| is not monotone in adding B --
if B contributes mostly non-canonical k-mers the combined rate FALLS. The combined proportion is
UNKNOWN.

THE HYPOTHESIS. An ncORF library whose entries are themselves substrings of canonical proteins has a
peptide space that is canonical-ambiguous BY CONSTRUCTION. Any peptide drawn from such an ORF is
compatible with both sources, and no filter downstream can undo that -- only an exclusion rule at the
peptide level can.

THE TEST. For each library we hold, what fraction of its HLA-I-length peptide space (9-mers) also
occurs in the reviewed canonical human proteome? This is the same measurement applied to the
catalogues, but applied to the SEARCH SPACE instead of the OUTPUT.

  nuORFdb                 -- one of the three sources IEAtlas integrates (RPFdb + nuORFdb + Translnc)
  GENCODE Ribo-seq ORFs   -- the community phase-1 / phase-2 ncORF sets

A CRITICAL CONTROL ON OUR OWN HYPOTHESIS: Ouspenskaia et al. ALSO searched nuORFdb, and their
published catalogue is at 3%. So if nuORFdb's peptide space is heavily canonical-ambiguous, the
library CANNOT by itself explain IEAtlas -- and the difference must lie in what each pipeline did with
the ambiguity. That is a real finding either way, and it is the honest way to run this test: the
result is informative whichever direction it goes.

    python3 scripts/library_ambiguity.py
"""
import os
import random
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
SPROT = os.path.join(EXT, "swissprot_human.fasta")

LIBS = [
    ("nuORFdb v1.2", os.path.join(EXT, "nuorfdb", "PA_nuORFdb_v1.2_protein.fasta"),
     "integrated by IEAtlas; ALSO searched by Ouspenskaia (published catalogue: 3%)"),
    ("GENCODE Ribo-seq ORFs (phase 1)",
     os.path.join(EXT, "gencode_riboseq_orfs", "phase1_Ribo-seq_ORFs.all.faa"),
     "community ncORF set"),
    ("GENCODE Ribo-seq ORFs (phase 2)",
     os.path.join(EXT, "gencode_riboseq_orfs", "phase2_Ribo-seq_ORFs.comprehensive.faa"),
     "community ncORF set, comprehensive"),
]

K = 9            # a representative HLA-I ligand length
# NO SAMPLING. A sampled estimate is BIASED UPWARD here and the bias is large: sampling 4,000 of
# nuORFdb's 229,251 ORFs gives 43.6%, 20,000 gives 40.7%, and the FULL library gives 34.1%. Small
# samples contain fewer distinct ncORF-specific k-mers, so the canonical-shared ones are
# over-weighted. Every library below is measured in full.


def read_fasta(path):
    seqs, cur = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).upper())
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).upper())
    return seqs


def main():
    if not os.path.exists(SPROT):
        sys.exit(f"missing {SPROT}")
    canon = read_fasta(SPROT)
    blob = "\x00".join(canon)          # \x00 prevents matches spanning two proteins
    print(f"canonical reference: {len(canon):,} reviewed human proteins, "
          f"{sum(map(len, canon)):,} residues")
    # Substring-searching an 11 MB blob per k-mer is O(n*m) and does not finish. Build the set of
    # canonical k-mers ONCE, then every lookup is O(1). ~11M k-mers; a few seconds and ~1 GB.
    print(f"indexing canonical {K}-mers ...", flush=True)
    CANON_K = set()
    for s_ in canon:
        for i in range(len(s_) - K + 1):
            CANON_K.add(s_[i:i + K])
    print(f"  {len(CANON_K):,} distinct canonical {K}-mers\n")

    print("=" * 94)
    print(f"HOW MUCH OF EACH ncORF LIBRARY'S {K}-mer PEPTIDE SPACE IS ALREADY CANONICAL?")
    print("=" * 94)
    print(f"  {'library':<34}{'ORFs':>9}{f'distinct {K}-mers':>13}{'also canonical':>12}{'':>8}")
    print("  " + "-" * 82)

    results = []
    for name, path, note in LIBS:
        if not os.path.exists(path):
            print(f"  {name:<34}  MISSING: {path}")
            continue
        orfs = read_fasta(path)
        kmers = set()
        for s in orfs:
            for i in range(len(s) - K + 1):
                kmers.add(s[i:i + K])
        hit = len(kmers & CANON_K)
        pct = 100 * hit / len(kmers) if kmers else 0
        results.append((name, len(orfs), len(orfs), len(kmers), hit, pct, note))
        print(f"  {name:<34}{len(orfs):>9,}{len(kmers):>13,}{hit:>12,} ={pct:6.1f}%")

    print()
    for name, n_orf, n_s, n_k, hit, pct, note in results:
        print(f"  {name:<34} {note}")

    # Whole-ORF containment: the strongest form of the problem.
    print("\n" + "=" * 94)
    print("STRONGER: is the ENTIRE ncORF protein a substring of a canonical protein?")
    print("  (If so, EVERY peptide it can yield is canonical-ambiguous by construction.)")
    print("=" * 94)
    for name, path, _note in LIBS:
        if not os.path.exists(path):
            continue
        orfs = [s for s in read_fasta(path) if len(s) >= K]
        rng = random.Random(0)
        # whole-ORF containment IS sampled -- it is an O(n*m) substring scan and unbiased under
        # random sampling (unlike the k-mer union, which is not).
        samp = orfs if len(orfs) <= 4000 else rng.sample(orfs, 4000)
        contained = sum(1 for s in samp if s in blob)
        print(f"  {name:<34}{contained:>7,}/{len(samp):<7,} = "
              f"{100*contained/max(1,len(samp)):5.1f}%  of ORFs lie wholly inside a canonical "
              f"protein (sampled)")

    print("\n" + "=" * 94)
    print("HOW TO READ THIS")
    print("=" * 94)
    print("""
  If a library's peptide space is heavily canonical-ambiguous, then the ambiguity is created at
  LIBRARY-CONSTRUCTION time, not at filtering time -- and no downstream rule other than a
  peptide-level exclusion can remove it.

  BUT: Ouspenskaia et al. also searched nuORFdb, and their published catalogue sits at 3%. So a
  high nuORFdb ambiguity would NOT by itself explain IEAtlas's 56.3%. It would instead show that
  the ambiguity is LATENT IN THE LIBRARY and that pipelines differ in whether they RESOLVE it --
  which is a sharper claim, and one that points at the exclusion rule as the operative difference
  after all.

  A LOW library ambiguity would REFUTE the library hypothesis outright, and we would report that.
""")


if __name__ == "__main__":
    main()
