# Revision notice — the manuscript in this directory is superseded

`manuscript.md` and `manuscript.html` in this directory are the version that went to external
review. **No `manuscript.pdf` is published**: a PDF cannot be updated once downloaded, and a PDF is
what gets cited — an unmarked PDF is exactly how a retracted result outlives its retraction. The
sources remain, so it is regenerable by anyone who wants it. **Two of their central claims have been withdrawn.** They are retained
unaltered so the correction is inspectable, not to be cited.

## 1. "0 of 306,844 claims pass all four evidence axes" — WITHDRAWN (true by construction)

The scorer could record a source-ORF pass only via a Ribo-seq **periodicity percentage** or a
protein-level **FDR**. Neither field is populated on a single claim in the corpus. So no claim
*could* pass that axis, so none could pass all four, and the zero followed from the scoring code
before any datum was read. It measured the scorer, not the literature.

## 2. Unassayed claims were scored as immunogenicity FAILURES — WITHDRAWN

A peptide observed by MS with no T-cell assay is not a negative immunogenicity result. It was
never assayed. Nearly the entire corpus is in that state.

Both errors are the same mistake — **absence of evidence recorded as evidence of absence** —
which is precisely the failure this project exists to criticise.

## 3. The pseudogene canonical-substring rate was 37.1%; it is 57.7%

The old figure came from a claim table that kept only the **first** peptide of each multi-peptide
gene row, so a gene-level rate was reported as a peptide-level one. Correcting the claim unit
moves the rate *against* the paper's own convenience: the overlap is larger, not smaller.

## What stands

The sequence-overlap results are unaffected and reproduce: the 56.3% atlas canonical-substring
rate, the class-resolved floors, and the identifiability argument. Note that the overlap is
**sequence non-uniqueness** — the reported record cannot distinguish the nominated ncORF from a
canonical source by sequence alone — and never a claim that the biology is absent.

## What is authoritative now

`src/` — and in particular the reporting-and-adjudicability matrix produced by
`src/darkproteome/audit.py`. Run `src/darkproteome/scoring_conformance.py` before quoting any
number from this repository.

---

## The typeset PDF of the superseded paper has been removed from this repository (2026-07-13)

`manuscript/manuscript.pdf` and `manuscript/tex/manuscript.pdf` rendered the **withdrawn** paper and
**carried no supersession notice inside them**. `manuscript.md` and `manuscript.html` do carry one; a
PDF cannot, and a PDF is the artifact that actually gets downloaded, shared and cited. Leaving an
unmarked render of a withdrawn result in a public repository is how a retracted number outlives its
retraction.

The Markdown and LaTeX sources remain, under their supersession notices, so the record is intact and
the PDF is regenerable by anyone who wants it. **The current paper is `manuscript_v2.md`.**
