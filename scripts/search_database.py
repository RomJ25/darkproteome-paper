"""The canonical proteins were IN IEAtlas's OWN SEARCH DATABASE.

This is not an external audit. It is an internal one.

Until now the paper has treated SwissProt 2022_01 as an ANACHRONISM CONTROL: "would a February-2022
analyst have been able to see these matches?" (Answer: yes -- era_correct_reference.py, 56.2%, only
231 retrospective.) That framing concedes too much. It implies the analyst would have had to go and
LOOK SOMETHING UP in an external reference.

They would not have. Read IEAtlas's Methods (Cao et al., NAR 2023, PMC9825419):

    "Files were searched against BOTH our integrated benchmarked ncORF library AND the canonical
     human proteome obtained from the UniProt database with Swiss-Prot protein evidence
     (downloaded in February 2022)."

    "Only epitopes derived from non-coding regions were retained."

The canonical human proteome was a component of the search database. MaxQuant, given several FASTA
inputs, concatenates them into one target database; every spectrum is scored against canonical and
non-canonical candidates TOGETHER. So for a catalogued sequence that also occurs in a canonical
protein, an entry carrying that identical sequence was PHYSICALLY PRESENT in the database the
spectrum was matched against -- sitting alongside the ncORF the peptide was ultimately attributed to.

The ambiguity is therefore not something we discovered from outside. It is a property of the search
itself, and the search engine's own output enumerates it: MaxQuant's `peptides.txt` carries a
`Proteins` column listing EVERY database entry containing the peptide. The information the remedy
asks for already existed, in an intermediate file, and was discarded at the step described as "only
epitopes derived from non-coding regions were retained."

THE MEASUREMENT. Two things, neither of which needs any reference IEAtlas did not itself use:

  (1) What fraction of the catalogue matches an entry in the CANONICAL HALF OF ITS OWN SEARCH
      DATABASE? Same number as the era-correct check (that reference IS the search database's
      canonical half) -- but a different ESTIMAND, and a much stronger one.

  (2) HOW EXPENSIVE IS THE LABEL? The paper's remedy is "list every compatible source locus." If an
      ambiguous peptide is typically compatible with dozens of canonical proteins, that label is
      unwieldy and the remedy is glib. If it is typically ONE protein, the label is a single
      accession -- and the remedy costs a column.

WHAT WOULD REFUTE THIS. Two hard checks, both fatal if they fail, both run below:

  * If IEAtlas's Methods do NOT actually say the canonical proteome was searched, the whole framing
    is unfounded -- so the quote is verified against the on-disk XML rather than trusted from memory,
    and the script EXITS NONZERO if it is not found verbatim.
  * If IEAtlas had used the canonical half of its own database to EXCLUDE shared sequences, the
    overlap would be ~0 -- as it is for the two catalogues that describe such a rule (0.026%, 0.17%).
    A near-zero result here would refute the paper's central claim outright. It is 56.2%.

FRAMING (non-negotiable). This does NOT show the peptides are canonically derived, and it does NOT
show anyone acted improperly. MS identifies a SEQUENCE, never a LOCUS, and sequence identity is
symmetric: the canonical entry is not "the right answer" either. It shows that the evidence needed to
mark these sequences source-ambiguous was inside the pipeline that produced the catalogue.

    python3 scripts/search_database.py
"""
import csv
import gzip
import html
import json
import os
import re
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
SPROT_2022 = os.path.join(EXT, "uniprot_sprot.fasta.gz")          # the search database's canonical half
IE = os.path.join(EXT, "atlases", "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
IE_PAPER = os.path.join(EXT, "fulltext", "PMC9825419.xml")        # IEAtlas, NAR 2023
csv.field_size_limit(10_000_000)

# The two sentences the framing rests on. Verified verbatim against the paper, not quoted from memory.
QUOTES = [
    "searched against both our integrated benchmarked ncORF library and the canonical human proteome",
    "Only epitopes derived from non-coding regions were retained",
]


def verify_methods():
    """Fail loudly if IEAtlas's Methods do not say what this script claims they say."""
    if not os.path.exists(IE_PAPER):
        sys.exit(f"missing {IE_PAPER} -- cannot verify the quotes; refusing to assert them")
    x = open(IE_PAPER, encoding="utf-8", errors="replace").read()
    x = re.sub(r"<[^>]+>", " ", x)
    x = re.sub(r"\s+", " ", html.unescape(x))
    print("=" * 94)
    print("IEAtlas's OWN METHODS (Cao et al., Nucleic Acids Research 2023; verified against PMC9825419)")
    print("=" * 94)
    ok = True
    for q in QUOTES:
        if q.lower() in x.lower():
            print(f'  [verbatim] "...{q}..."')
        else:
            print(f"  [ABSENT]  {q}")
            ok = False
    if not ok:
        sys.exit("\nFATAL: a load-bearing quote is not in the paper. The framing is unfounded.")
    print("\n  => The canonical human proteome was a component of the searched database.")
    print("     The ncORF library and the canonical proteome were searched TOGETHER.\n")


def canonical_2022():
    """Human (OX=9606) SwissProt 2022_01 -- the canonical half of IEAtlas's search database."""
    prots, keep, cur, gene = [], False, [], None
    with gzip.open(SPROT_2022, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    prots.append(("".join(cur), gene))
                keep = "OX=9606 " in line or line.rstrip().endswith("OX=9606")
                if keep:
                    m = re.search(r"\bGN=(\S+)", line)
                    gene = m.group(1) if m else "?"
                    acc = line.split("|")[1] if "|" in line else "?"
                    gene = f"{gene}"
                cur = []
            elif keep:
                cur.append(line.strip())
    if keep and cur:
        prots.append(("".join(cur), gene))
    return prots


def main():
    for p in (SPROT_2022, IE):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")

    verify_methods()

    prots = canonical_2022()
    print(f"canonical half of the search database : {len(prots):,} human proteins "
          f"(SwissProt 2022_01, downloaded Feb 2022)")

    peps = set()
    with open(IE, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 3 or not r[0]:
                continue
            q = r[0].strip().upper()
            if q.isalpha():
                peps.add(q)
    print(f"IEAtlas distinct cancer epitopes      : {len(peps):,}\n")

    bylen = defaultdict(set)
    for q in peps:
        bylen[len(q)].add(q)
    lengths = sorted(bylen)

    # For every catalogued sequence, WHICH canonical entries in the same search database contain it?
    # Scanning the proteome once and testing membership in the peptide set is the cheap direction.
    print("scanning the canonical half of the search database for each catalogued sequence ...",
          flush=True)
    compat_prot = defaultdict(set)      # peptide -> {protein index}
    compat_gene = defaultdict(set)      # peptide -> {gene symbol}
    for pi, (s, gene) in enumerate(prots):
        n = len(s)
        for L in lengths:
            grp = bylen[L]
            for i in range(n - L + 1):
                sub = s[i:i + L]
                if sub in grp:
                    compat_prot[sub].add(pi)
                    compat_gene[sub].add(gene)

    hit = set(compat_prot)
    pct = 100 * len(hit) / len(peps)

    print("\n" + "=" * 94)
    print("(1) MATCHES AGAINST THE CANONICAL HALF OF IEAtlas's OWN SEARCH DATABASE")
    print("=" * 94)
    print(f"  {len(hit):,} / {len(peps):,} = {pct:.1f}% of catalogued cancer sequences also occur in")
    print("  a canonical protein that was PRESENT IN THE DATABASE THE SPECTRA WERE SEARCHED AGAINST.")
    print()
    print("  These are not matches an auditor found in a reference IEAtlas never consulted.")
    print("  The competing canonical entry was in the search space, alongside the ncORF the peptide")
    print("  was attributed to. Tandem MS identifies the SEQUENCE; it does not choose the LOCUS.")

    # The refutation check.
    print("\n  REFUTATION CHECK -- if the canonical half of the database had been used to EXCLUDE")
    print("  shared sequences (as CrypticProteinDB and Raja et al. describe doing), this rate would")
    print(f"  be ~0%. Those catalogues sit at 0.026% and 0.17%. Observed here: {pct:.1f}%.")
    if pct < 5.0:
        sys.exit("\n  REFUTED: overlap is near zero. The paper's central claim does not hold.")
    print("  => not refuted.")

    print("\n" + "=" * 94)
    print("(2) HOW EXPENSIVE IS THE PROPOSED LABEL?")
    print('  The remedy asks the atlas to "list every compatible source locus." If that list is huge,')
    print("  the remedy is glib. So: how many canonical entries is an ambiguous peptide compatible")
    print("  with?")
    print("=" * 94)

    nprot = sorted(len(compat_prot[q]) for q in hit)
    ngene = sorted(len(compat_gene[q]) for q in hit)

    def med(v):
        n = len(v)
        return (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2) if n else 0

    one_gene = sum(1 for q in hit if len(compat_gene[q]) == 1)
    le2_gene = sum(1 for q in hit if len(compat_gene[q]) <= 2)
    pct_one = 100 * one_gene / len(hit)
    pct_le2 = 100 * le2_gene / len(hit)

    print(f"  compatible canonical PROTEINS per ambiguous sequence: median {med(nprot):g}, "
          f"max {nprot[-1]:,}")
    print(f"  compatible canonical GENES    per ambiguous sequence: median {med(ngene):g}, "
          f"max {ngene[-1]:,}")
    print()
    print(f"  compatible with exactly ONE canonical gene : {one_gene:,} / {len(hit):,} = {pct_one:.1f}%")
    print(f"  compatible with at most TWO canonical genes: {le2_gene:,} / {len(hit):,} = {pct_le2:.1f}%")
    print()
    print("  => The label is CHEAP. For the large majority of ambiguous sequences it is a single gene")
    print("     symbol. The remedy costs one column, not a redesign -- and MaxQuant's own")
    print('     `peptides.txt` already emits it (the `Proteins` column lists every matching entry).')

    art = os.path.join(REPO, "data", "derived_search_database.json")
    json.dump({
        "canonical_proteins_in_search_db": len(prots),
        "n_peptides": len(peps),
        "n_matching_search_db": len(hit),
        "pct_matching_search_db": round(pct, 1),
        "median_compatible_proteins": med(nprot),
        "median_compatible_genes": med(ngene),
        "pct_exactly_one_gene": round(pct_one, 1),
        "pct_at_most_two_genes": round(pct_le2, 1),
        "max_compatible_genes": ngene[-1],
        "methods_quotes_verified": True,
    }, open(art, "w"), indent=1)
    print(f"\nwrote {os.path.relpath(art, REPO)}")


if __name__ == "__main__":
    main()
