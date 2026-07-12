"""What the primary Methods sections actually report.

An audit of the reporting record is only as good as its reading of that record. Supplementary tables
do not carry a study's search parameters; Methods sections do. This script extracts and verifies,
from the fetched full texts, what each audited source reports about its FDR handling, its translation
statistics, and its search database.

The distinction that matters throughout: a STUDY-LEVEL threshold (a global FDR, a score cutoff, an
average periodicity) is not a CLAIM-LINKED statistic (this ORF's own score, this peptide's own
q-value). Sources report the former. None publishes the latter, so no individual claim can be
independently re-adjudicated from the reported record.

    python3 scripts/methods_audit.py
"""
import hashlib
import html
import os
import re
import sys
import textwrap

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FT = os.path.join(REPO, "data", "external", "fulltext")
PAPERS = {"HCC (adn3628)": "PMC11235171.xml", "ovarian (ads7405)": "PMC11837991.xml"}


def text_of(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", s)))


def show(t, kw, before=180, after=220, limit=3):
    out = []
    for m in re.finditer(kw, t, re.I):
        a, b = max(0, m.start() - before), min(len(t), m.end() + after)
        out.append(t[a:b])
        if len(out) >= limit:
            break
    return out


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
    _require(FT)
    if not os.path.isdir(FT):
        sys.exit(f"missing {FT} -- fetch the open-access full texts first")

    texts = {}
    print("=== SOURCES ===")
    for name, fn in PAPERS.items():
        p = os.path.join(FT, fn)
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        texts[name] = text_of(p)
        print(f"  {name:<20} {fn}  sha256={h[:16]}…  {len(texts[name].split()):,} words")

    print("\n" + "=" * 84)
    print("WHAT THE METHODS REPORT — and what we wrongly said they do not")
    print("=" * 84)

    print("""
  OUR CLAIM (in code and in the draft): `reported_fdr` and `periodicity_pct` are populated on
  0 of 307,318 rows, therefore "no source reports a protein-level FDR per claim" and
  "translation is ASSERTED and never QUANTIFIED".

  THE SECOND HALF OF THAT IS FALSE. Both papers report translation and search statistics. They
  report them at the STUDY level, not per claim -- which is a different, narrower, and still
  damning finding. But as written, our claim mischaracterises what the authors did.
""")

    facts = [
        ("HCC (adn3628)", "FDR threshold set at 5%",
         "their OWN immunopeptidomics search (Comet + Percolator via MHCquant) was filtered at "
         "5% FDR. Not 1%, not 3% -- those numbers in the text refer to OTHER studies they "
         "re-analysed."),
        ("HCC (adn3628)", "RibORF score of at least 0.5",
         "translation was called by a RibORF score threshold of >=0.5 with >=5 footprints. A "
         "THRESHOLD, applied per ORF -- but the per-ORF SCORE is never published."),
        ("HCC (adn3628)", "average periodicity",
         "an AVERAGE periodicity over 0.5 (0.58 / 0.51 for 28/29-bp reads) is reported for the "
         "READ POOL. This is a study-level property of the Ribo-seq library, not a per-ORF "
         "statistic, and cannot adjudicate any individual claim."),
        ("HCC (adn3628)", "103,789 sequences",
         "THE SEARCH DATABASE INCLUDED Swiss-Prot/TrEMBL WITH ISOFORMS (103,789 sequences, "
         "downloaded 21 Apr 2023) alongside 5,021 ncORF sequences."),
        ("HCC (adn3628)", "uniquely matching peptides",
         "the Methods state that ONLY UNIQUELY MATCHING PEPTIDES were considered."),
        ("ovarian (ads7405)", "3% at the PSM level",
         "FDR controlled by Percolator, filtered to 3% at the PSM level for both search engines."),
    ]
    for paper, kw, gloss in facts:
        hits = show(texts[paper], re.escape(kw), limit=1)
        print(f"\n  [{paper}]  \"{kw}\"")
        print(textwrap.fill(gloss, 96, initial_indent="    => ", subsequent_indent="       "))
        for h in hits:
            print(textwrap.fill(f"…{h}…", 96, initial_indent="    | ", subsequent_indent="    | "))

    print("\n" + "=" * 84)
    print("THE CORRECTED CLAIM")
    print("=" * 84)
    print("""
  WRONG : "The papers do not report translation statistics or an FDR."
  RIGHT : "Both papers report STUDY-LEVEL search and translation thresholds (HCC: 5% FDR,
           RibORF score >=0.5, >=5 footprints; ovarian: 3% PSM FDR). NEITHER publishes the
           PER-CLAIM statistic -- no ncORF's own RibORF score, and no peptide's own q-value --
           so no individual claim can be independently re-adjudicated from the reported record."

  The reporting-ladder result is UNCHANGED and is now correctly grounded: `asserted` is satisfied
  (the analysis is named and thresholded), `claim_linked` fails for the statistic (only a binary
  pass/fail reaches the reader), and so `quantitative` and `adjudicable` are 0. The LADDER was
  right. The PROSE was wrong.
""")

    print("=" * 84)
    print("A NEW FINDING, AND IT IS SHARPER THAN ANYTHING WE HAD")
    print("=" * 84)
    print("""
  The HCC Methods state BOTH of the following:

    (a) the search database comprised Swiss-Prot/TrEMBL INCLUDING ISOFORMS -- 103,789 sequences --
        alongside 5,021 predicted ncORF sequences; and
    (b) "We only considered uniquely matching peptides."

  A peptide that is an exact substring of a reviewed Swiss-Prot protein is NOT uniquely matching
  against a database that contains that protein.

  Yet 213 of 369 (57.7%) of the reported HCC pseudogene-ORF peptides ARE exact substrings of
  reviewed Swiss-Prot proteins -- and for 33 of 34 such pseudogenes, at least one peptide lands in
  the pseudogene's OWN PARENT gene.

  These two facts are in TENSION. We cannot say from the published record which database entry any
  peptide was assigned to, because the per-PSM protein assignment is not published (sheet S26 gives
  43 PSMs of the many thousands reported). So we do not claim a pipeline error, and we should not:

    THE FINDING IS THE IRRESOLVABILITY. The paper states a uniqueness requirement against a
    database that contains the canonical proteins these peptides match, and the reported record
    does not let a reader determine how that requirement was satisfied for them.

  This is a far stronger, far more specific claim than "the field reports poorly" -- and it is
  available ONLY because we read the Methods. It was invisible in the supplementary tables.
""")

    print("=" * 84)
    print("WHAT THIS SETTLES ABOUT THE AUDIT'S SCOPE")
    print("=" * 84)
    print("""
  A supplements-only audit is NOT SAFE. It produced a false claim ("the papers do not report an
  FDR") that we were one step from publishing, in a paper whose entire authority rests on
  characterising the reporting record accurately. One such error would be fatal.

  Both papers were open access and were fetched in a single call. There is no defence for not
  having read them. Every audited primary source must have its Methods read before publication.
""")


if __name__ == "__main__":
    main()
