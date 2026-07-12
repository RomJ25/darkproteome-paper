# Supplement

Companion to `manuscript_v2.md`. Two items sit here rather than in the main text, for two different
reasons: **S1** is a real result that is **not yet headline-grade**, and says so; **S2** is a
derivation of a result that is **not ours** (Woo et al. 2014) and is reproduced only because the main
text depends on its exact form.

---

## S1. Where it can be interrogated, the ambiguity is structured by homology

**This quantifies a known phenomenon.** That a processed pseudogene shares sequence with its parent
gene, and that the resulting peptides are not attributable to one locus, is the classical
shared-peptide / protein-inference problem (Nesvizhskii & Aebersold). We measure it; we do not
discover it.

**Why this is in the Supplement.** Parentage here is **derived from the gene symbol** (strip a
trailing `P` + digits, allowing documented root-renames). That is not an authoritative relationship.
Gene-family naming — `FAM*`, zinc fingers, ribosomal proteins, and duplicated genes generally — makes
string-derived parentage hazardous, and an external reviewer was right that this is too thin to carry
a headline biological mechanism. **Promoting this result requires an authoritative, versioned
pseudogene–parent annotation and manual adjudication of all 34 gene-level cases.** Until that exists,
it is reported as suggestive and nothing more.

### The measurement

A processed pseudogene is a retro-copy of a parent gene, so a peptide encoded by a pseudogene ORF can
be an exact substring of the **parent protein**.

Of the pseudogene-labelled IEAtlas peptides that are canonical substrings, **132 are testable** (the
rest carry clone-style symbols from which no parent is derivable, or a derived parent absent from the
reference). Of those:

| | |
|---|---:|
| land in the pseudogene's **own parent gene** | **79 / 132 = 59.8%** |
| at gene level, pseudogenes with ≥1 peptide in their own parent | **33 / 34 = 97.1%** |
| match the parent **alone** | 52 / 79 = 65.8% |
| match the parent **and further canonical genes** | 27 / 79 = 34.2% |
| land in a canonical gene that is **not** the derived parent | **53 / 132 = 40.2%** |

Peptides whose symbol yields no parent, or whose derived parent is absent from the reference, are
**untestable and excluded from the denominator** — scoring them as misses would invent evidence
against the class.

### The null

A uniform "one gene in ~20,400" null is **invalid**, and an earlier draft that used it was wrong:
pseudogene parent annotations are themselves homology-derived, and we select peptides already known to
match *something* canonical. Both facts destroy the uniformity assumption.

We instead hold each peptide's canonical hit set **fixed** — preserving the selection, the protein
lengths, the paralogy and the k-mer structure — and permute **only** the pseudogene→parent pairing
(10,000 permutations, seed 0).

    observed   79 / 132
    null mean   4.9 / 132 = 3.7%    (null max 14; 0 permutations ≥ observed)
    p < 1e-4

**This null is itself imperfect**, and we flag it rather than hide it: it assumes parent labels are
exchangeable across the pool, which they are not — gene families, parent lengths and paralogy
structure all differ. A stronger test would permute within biologically meaningful strata, or test
against an authoritative candidate-parent set. This is a further reason S1 is not a headline.

### What it does and does not license

**It does** show that canonical-sequence compatibility is *concentrated in the annotated parent* — the
ambiguity sits exactly where sequence homology predicts, so the pseudogene class label is corroborated
by sequence structure rather than being arbitrary.

**It does not resolve provenance.** `DEVAFRKF` is encoded by both *RPS3AP12* and *RPS3A*; MS
identifies the sequence, not the locus. A parent hit is rarely even a *unique* assignment: 34.2% of
parent hits also match further canonical genes — more ambiguity, not less.

**The 53 non-parent matches (40.2%) are an unresolved residual.** Out-of-frame retro-copy ORFs and
incidental short-peptide matches are both plausible; parent-hits and non-parent hits share a median
length of 9 aa, so length does not separate them. **We report this and do not explain it.**

Across other classes we report only descriptive heterogeneity in canonical overlap. The retro-copy
mechanism is claimed for **processed pseudogenes only** — the one class where the label is corroborated
by sequence. Low overlap in lncRNA-ORF (0.5%) and altORF (0.2%) classes argues against one specific
failure mode (wholesale pseudogene-like contamination of those classes) and nothing more; it does not
validate those annotations.

---

## S2. Set-identification of a class-specific FDR (Woo et al. 2014; reproduced, not original)

Let a target–decoy procedure accept a set of identifications at a reported pooled false-discovery rate
*q*. Partition the accepted set into a non-canonical class *N* and its complement, and let *f* be the
fraction of accepted identifications in class *N*.

Write *FDR_N* for the false-discovery rate **within class *N***, and *FDR_C* for that in the
complement. The pooled rate is the *f*-weighted mixture:

    q = f · FDR_N + (1 − f) · FDR_C

with the only constraints being that both class-specific rates are probabilities:

    0 ≤ FDR_N ≤ 1,   0 ≤ FDR_C ≤ 1

Solving for *FDR_N* and imposing those bounds on *FDR_C* gives

    FDR_N = ( q − (1 − f) · FDR_C ) / f

which is decreasing in *FDR_C*. Substituting the two extremes *FDR_C* = 1 and *FDR_C* = 0, and
intersecting with [0, 1]:

    Θ_N(q, f) = [ max(0, (q − (1 − f)) / f),  min(1, q / f) ]

**Sharpness.** Every point in Θ_N is attained by some admissible *FDR_C* ∈ [0, 1], and no point outside
it is: the endpoints correspond exactly to *FDR_C* = 1 and *FDR_C* = 0. So **given *q* and *f* alone,
Θ_N is the sharp identified set** — the data-generating process is not pinned down more tightly by that
information.

**What this does *not* say.** Θ_N is sharp *with respect to q and f*. It is **not** the claim that no
further information could narrow it. Calibrated class-specific posterior error probabilities,
entrapment measurements, or a validated mixture model would each add information and could tighten the
interval. The per-class accepted decoy count *D_N* is simply the **cheapest sufficient** such object,
and one the pipeline already computes: with *T_N* and *D_N* and the stated threshold, unit and
convention, the selected class-specific target–decoy estimate becomes **reconstructible**. *D_N* does
not identify the true class-specific false-discovery *proportion*, and we do not claim it does.

**Why it bites here.** IEAtlas reports *q* = 0.05 and a non-canonical count (245,870), but **no
canonical count**. *f* is therefore unknown, and Θ_N(0.05, *f*) is **unconstrained** — for small *f*
the upper endpoint min(1, *q*/*f*) reaches 1. The interval cannot be evaluated at all from what the
resource publishes. That is the reporting gap, and it is closed by one small table.

**Worked illustration.** At *q* = 0.05:

| *f* (non-canonical share of accepted IDs) | Θ_N(0.05, *f*) |
|---:|---|
| 0.50 | [0.00, 0.10] |
| 0.20 | [0.00, 0.25] |
| 0.05 | [0.00, 1.00] |
| 0.01 | [0.00, 1.00] |

The scarcer the class, the less a pooled threshold says about it — which is the entire point of Woo et
al. 2014, restated here only because the resource under audit does not report the quantity that would
resolve it.
