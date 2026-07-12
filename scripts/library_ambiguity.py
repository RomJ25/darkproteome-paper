"""Latent canonical ambiguity of an ncORF library.

A catalogue of "non-canonical" peptides inherits the properties of the library it was searched
against. If a library's peptide space overlaps the canonical proteome, then peptides drawn from it
are compatible with BOTH a canonical protein and an ncORF -- and no downstream step other than a
peptide-level exclusion can undo that.

This measures, for each ncORF library, what fraction of its HLA-I-length peptide space (distinct
9-mers) also occurs in the reviewed canonical human proteome. It is the catalogue metric applied to
the SEARCH SPACE instead of the OUTPUT.

  NO SAMPLING. A sampled estimate is biased UPWARD here and the bias is large: 4,000 of nuORFdb's
  229,251 ORFs gives 43.6%, 20,000 gives 40.7%, and the full library gives 34.1%. Small samples hold
  fewer distinct ncORF-specific k-mers, so canonical-shared ones are over-weighted. Whole-ORF
  containment IS sampled -- that statistic is an O(n*m) substring scan and is unbiased under random
  sampling, unlike a k-mer union.

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


def _require(*paths):
    """Fail with a usable message, not a traceback, when the external inputs are absent.

    The large public inputs (Swiss-Prot, the atlas exports, the ncORF libraries, the fetched full
    texts) are not redistributed in this repository. Populate `data/external/` from the sources
    documented in `data/SOURCES.md` and `data/external/README.md`.
    """
    import sys as _s
    missing = [p for p in paths if not __import__("os").path.exists(p)]
    if missing:
        _s.exit("missing required input(s):\n  " + "\n  ".join(missing) +
                "\n\nThese are large public files and are not redistributed here.\n"
                "Populate data/external/ -- see data/SOURCES.md and data/external/README.md.")


def main():
    _require(SPROT)
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
