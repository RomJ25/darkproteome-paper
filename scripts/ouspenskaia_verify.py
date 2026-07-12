"""Verify the positive control AT SOURCE, and establish the contrast with IEAtlas.

Ouspenskaia et al. 2022 (Nat Biotechnol) is the paper that makes this project's theoretical claim
EMPIRICAL. It is a PMC author manuscript: the web page serves a JavaScript shell to `curl`, which is
why a naive fetch returns 23 words. E-utilities serves the full text.

    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC10198624&rettype=xml"

WHAT IT ESTABLISHES -- and the load-bearing number is one that no summary quotes:

    "While global FDR was set to 1%, FDR for nuORF peptides was 4.6% overall, and as high as 14%
     for 3' dORFs. We devised a group-based filtering approach to reduce the nuORF FDR rate to 1%
     across different types of nuORFs. This approach REMOVED 24% of nuORF peptides overall, and up
     to 76% of peptides assigned to 3' overlap dORFs, retaining 6,501 high confidence (FDR<1%)
     peptides from 3,261 nuORFs."

Class-specific FDR control DELETES a quarter of the nuORF peptides, and three quarters of one class.
That is the predicted inflation, measured in the literature, by authors who then CORRECTED it.
Ouspenskaia is a POSITIVE EXEMPLAR, not a failure case, and the 4.6% is a PRE-correction diagnostic.
Citing it as the FDR of their final catalogue would misrepresent a paper that did the right thing.

    python3 scripts/ouspenskaia_verify.py
"""
import hashlib
import html
import os
import re
import sys
import textwrap

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FT = os.path.join(REPO, "data", "external", "fulltext")
OUSP = os.path.join(FT, "ousp_efetch.xml")      # PMC10198624 via E-utilities
IEATLAS = os.path.join(FT, "PMC9825419.xml")    # IEAtlas, gkac776


def text_of(p):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ",
                  open(p, encoding="utf-8", errors="replace").read())))


def main():
    for p in (OUSP, IEATLAS):
        if not os.path.exists(p):
            sys.exit(f"missing {p} -- fetch the full texts first (see docstring)")

    o, ie = text_of(OUSP), text_of(IEATLAS)
    print("=== SOURCES ===")
    for name, p in (("Ouspenskaia (PMC10198624)", OUSP), ("IEAtlas (PMC9825419)", IEATLAS)):
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
        print(f"  {name:<28} sha256={h}…  {len(text_of(p).split()):,} words")

    print("\n" + "=" * 80)
    print("1. THE POSITIVE CONTROL, VERBATIM")
    print("=" * 80)
    quotes = [
        "While global FDR was set to 1%",
        "much higher FDR for nuORFs",
        "removed 24% of nuORF peptides",
        "6,501",
    ]
    for q in quotes:
        m = re.search(re.escape(q), o, re.I)
        if not m:
            print(f"  *** NOT FOUND: {q!r} -- DO NOT CITE ***")
            continue
        a, b = max(0, m.start() - 90), min(len(o), m.end() + 300)
        print(textwrap.fill(f"  …{o[a:b]}…", 96, subsequent_indent="   "))
        print()

    print("=" * 80)
    print("2. THE CONTRAST — did IEAtlas apply class-specific FDR control?")
    print("=" * 80)
    pats = [r"group[- ]based", r"group[- ]specific", r"class[- ]specific", r"nuORF FDR",
            r"separate.{0,30}FDR", r"per[- ](type|class|group)"]
    hits = [m.group(0) for p in pats for m in re.finditer(p, ie, re.I)]
    print(f"  mentions of group/class/type-specific FDR control in IEAtlas: {len(hits)}")
    if not hits:
        print("  => NONE. Exhaustive search of the full text finds no group-, class-, or")
        print("     type-specific FDR control anywhere in IEAtlas.")
    print()
    for label, pat, txt in (("IEAtlas pooled FDR", r"peptide spectrum match false discovery rate[^.]*", ie),
                            ("IEAtlas retention", r"Only epitopes derived[^.]*", ie),
                            ("IEAtlas database", r"searched against both[^.]*", ie)):
        m = re.search(pat, txt, re.I)
        if m:
            print(textwrap.fill(f"  [{label}] {m.group(0)}.", 96, subsequent_indent="      "))

    print("\n" + "=" * 80)
    print("3. THE CONTRAST, TABULATED")
    print("=" * 80)
    print("""
                            Ouspenskaia 2022          IEAtlas 2023
    search database         annotated + nuORF         canonical + ncORF   (both COMBINED)
    pooled FDR              1%                        5%    <- 5x more permissive
    protein-level FDR       --                        "no protein FDR was set"
    class-specific control  YES (group-based)         NONE FOUND
    effect of that control  removed 24% of nuORF
                            peptides; up to 76% of
                            3' overlap dORFs          n/a
    retained                6,501 at <1% class FDR    "only epitopes derived from
                                                       non-coding regions"

  The failure mode Ouspenskaia DIAGNOSED AND CORRECTED in 2022 is unaddressed in the atlas
  published the following year, at a pooled threshold five times more permissive.

  This is a documented methodological difference between two public resources, both of which
  describe their own procedures accurately. It is NOT an accusation of misconduct.
""")

    print("=" * 80)
    print("4. TWO DISTINCT PROBLEMS — DO NOT CONFLATE (this would be the next error)")
    print("=" * 80)
    print("""
  A. SOURCE AMBIGUITY (our sequence measurement)
     The peptide is CORRECTLY identified as a sequence; its SOURCE LOCUS is unresolved,
     because the same sequence is encoded by a canonical protein AND an ncORF.
     NOT an FDR problem -- no identification is wrong.

  B. CLASS-SPECIFIC FDR UNDER-CONTROL (Ouspenskaia's measurement)
     A pooled threshold under-controls the minority class, so a fraction of ncORF peptide
     IDENTIFICATIONS ARE WRONG -- the spectrum was not that peptide.

  Both inflate the ncORF catalogue, by DIFFERENT routes. A peptide can be a genuine
  identification with an ambiguous source (A), a false identification (B), or both.
  Merging them is a category error.
""")


if __name__ == "__main__":
    main()
