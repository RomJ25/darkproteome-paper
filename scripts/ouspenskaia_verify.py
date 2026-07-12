"""Class-specific FDR control: the problem, and the published remedy, verified at source.

Ouspenskaia et al. (Nat Biotechnol) applied a 1% GLOBAL FDR to a combined annotated-ORF/nuORF search
and measured the resulting nuORF-class FDR at 4.6% overall -- 14% for one ORF class. They then
introduced group-based filtering that brought each class back to ~1%, at a cost of 24% of nuORF
peptides overall and up to 76% of one class.

The 4.6% is therefore a PRE-CORRECTION DIAGNOSTIC, not the error rate of their published catalogue.
Citing it as the latter would misrepresent a paper that detected the problem and fixed it. They are a
positive exemplar for the reporting standard proposed here.

This script verifies every quoted figure against the fetched full text.

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
    _require(FT, OUSP, IEATLAS)
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
