"""RETIRED. Use `evidence_dimensions.py`.

This module scored four "independent axes" as a pass/fail survivor funnel. That framing produced
two defects of the same shape — a scorer deciding something the reported record never spoke to —
and both of them reached the manuscript:

  1. STRUCTURAL ZERO. `score_source_orf` could strict-pass only via a Ribo-seq periodicity % or a
     protein-level FDR. NEITHER FIELD IS POPULATED ON ANY OF 307,318 ROWS, so no claim could pass,
     so "0 of 306,844 claims pass all four axes" — the abstract — was true before any datum was
     read. It measured the scorer, not the literature.

  2. UNASSAYED READ AS FAILED. `score_immunogenicity` returned `fail` for `MS-presented`. An
     MS-presented peptide with no T-cell assay is not a negative immunogenicity result; it is NOT
     ASSAYED. 307,274 of 307,318 rows are MS-presented, so this reported adjudicability as 307,278
     when the number of claims carrying a human T-cell result is 2.

Both are absence of evidence recorded as evidence of absence. `evidence_dimensions.py` cannot
express that error, because "the record does not say" is a first-class state there.

Also retired with it: the conjunction in `score_hla_presentation` (eluted AND allele AND length),
which let near-universal ELUTION mask the fact that an ALLELE is reported for 275 of 307,318
claims. Presentation and allele restriction are now separate dimensions.

Importing this raises. There is no compatibility shim on purpose: a silent fallback is how the
first three scorers came to coexist and disagree.
"""

raise ImportError(
    "axes.py is RETIRED -- it scored unassayed claims as failures and had a "
    "structurally unpassable source axis. Use evidence_dimensions.py; see this file's docstring."
)
