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

**Why this is still in the Supplement.** An external reviewer raised two objections: parentage was
**derived from the gene symbol** (strip a trailing `P` + digits), which is not authoritative; and the
permutation null assumed an exchangeability that **gene families violate**. Both have now been
addressed — with an authoritative annotation and a family-respecting null, below — and the result
survives both. It nonetheless stays here, for two reasons the reader should weigh: the testable set is
**133 peptides** out of 174,465, and the authoritative mapping is **not independent** of the symbol
heuristic it replaces (see the caveat below). It is a mechanistic vignette that corroborates the main
result; it is not itself a headline measurement.

### The authoritative parent mapping

Parentage is now taken from **NCBI Gene `gene_group`** — the curated *"Related pseudogene"* relation
(build 2026-07-13; 13,437 human pseudogenes with a curated parent; no multi-parent cases), with symbols
resolved through `Homo_sapiens.gene_info`, HGNC (`hgnc_complete_set`, build 2026-07-10) including previous
and alias symbols, and GENCODE v26 as a clone-name→ENSG bridge. (pseudogene.org's psiDR and PseudoPipe
were downloaded and **not used**: they are built on GENCODE v7/hg19 (2012) and Ensembl 90 (2017)
respectively, too stale to adjudicate current symbols.)

**Coverage, stated plainly.** Of the 62 distinct pseudogene symbols carrying a canonical-substring
peptide, **36 resolve to a curated parent (58.1%)**. The shortfall is **not** annotation failure: the
25 unresolved symbols are clone-style identifiers (`AC005262.1`, `AL158050.1`, …) that are absent from
**every** registry. Among symbols that have a real gene name, coverage is **97.3%**.

**The caveat that matters, and it is not a small one.** HGNC *names* a pseudogene after the parent it
descends from. The curated relation is therefore **not an independent source of truth** about
parentage — it is substantially the same information the symbol heuristic was reading, curated. High
agreement between the two is **expected, and is not corroboration**. What the curated mapping actually
buys is a *versioned, adjudicated* relation with correct handling of renames and family cases — which
is what the reviewer asked for — and nothing more. The **test** itself is untainted: whether a peptide
occurs inside the parent's protein is a pure sequence question, decided independently of any naming.

### The measurement

A processed pseudogene is a retro-copy of a parent gene, so a peptide encoded by a pseudogene ORF can
be an exact substring of the **parent protein**. Of the pseudogene-labelled IEAtlas peptides that are
canonical substrings, **133 are testable**:

| | authoritative | (symbol heuristic) |
|---|---:|---:|
| land in the pseudogene's **own curated parent** | **77 / 133 = 57.9%** | 79 / 132 = 59.8% |
| at gene level, pseudogenes with ≥1 peptide in their own parent | **33 / 35** | 33 / 34 |

**How wrong was the heuristic?** Where both mappings yield a parent, they **agree in 124 / 127 = 97.6%**
of cases; it disagrees 3 times, and in 7 further cases the heuristic yields no parent where a curated
one exists. **Seven peptide-level verdicts flip.** The headline moves from 59.8% to **57.9%** — so the
reviewer's objection was **correct in principle and small in effect**. We report both.

### The null — and two that turned out to be degenerate

The reviewer's objection was that permuting parent labels freely ignores gene-family structure, making
"matched the parent" trivially easy. We ran four nulls (10,000 permutations, seed 0) and report all of
them, including the ones that failed:

| null | construction | result |
|---|---|---|
| **A** — naive / free | parents exchangeable across the whole pool (*the objectionable one*) | 77 obs vs **4.4** null mean; *p* < 1e-4 |
| **B** — HGNC gene-family | shuffle parents only within a curated family | **DEGENERATE** — only 5 permutable items |
| **C** — shared-9-mer component | shuffle only within a sequence-sharing component | **DEGENERATE** — see below |
| **D** — family-decoy swap | replace the true parent with a random **close paralog of that parent** | see below — *the null to read* |

**Nulls B and C are degenerate, and we do not quote their *p*-values.** Null C collapses for a reason
that is itself the answer to the objection: **the 35 parent proteins are pairwise 9-mer-disjoint**, so
every homology stratum is a singleton and the "family-respecting" permutation reduces to the identity.
There is no family structure *among the parents* for a shuffle to exploit. Reporting C's *p* = 1.0 as
if it were evidence would be as dishonest as hiding it.

**Null D is the family-respecting null that actually runs.** It asks the sharp version of the
objection: *could a merely homologous protein have been hit instead of the true parent?* Each parent is
replaced by a random **close paralog of itself** — the hardest decoys available — under two pools:

| decoy pool | true parent hit | random paralog hit | *p* |
|---|---:|---:|---:|
| any shared 9-mer (permissive) | 53 / 101 = **52.5%** | **7.1%** | < 1e-4 |
| **strong paralogs (≥10 shared 9-mers)** | 34 / 65 = **52.3%** | **16.6%** | < 1e-4 |

The hit is **parent-specific, not family-generic** — it survives even against the parent's genuine close
homologs. Without any null at all: **58 / 77 = 75.3%** of parent-hits are compatible with the parent and
with **no** close paralog. The remaining 24.7% are not, and we say so; family structure is real and
quantified rather than denied.

**The objection was right for the reason it gave, and we caught the heuristic in the act.** The old
symbol rule treated two symbols as the same gene if they shared a ≥3-character prefix with a numeric
remainder — so it judged `ZNF720` and `ZNF135` the same gene. `ZNF720P1`'s peptide `KSFSHSSSL` occurs
in *ZNF135*, *ZNF256* and *ZNF483*, and **not** in its true curated parent *KRABD5* (of which `ZNF720`
is merely a previous symbol). The heuristic scored a parent hit by zinc-finger string collision; the
curated relation deletes it. **Aggregate cost of that failure mode: one hit.**

### What it does and does not license

**It does** show that canonical-sequence compatibility is *concentrated in the curated parent* — the
ambiguity sits exactly where descent predicts, and is parent-specific rather than family-generic
(null D). The pseudogene class label is corroborated by sequence structure rather than being arbitrary.

**It does not resolve provenance.** `DEVAFRKF` is encoded by both *RPS3AP12* and *RPS3A*; MS identifies
the sequence, not the locus. A parent hit is not even a *unique* assignment — many parent hits also
match further canonical genes. This is **more** ambiguity, not less.

**It is not a surprise, and we do not dress it up as one.** A processed pseudogene *is* a degenerate
copy of its parent; that its peptides match the parent's protein is the mechanism working as expected.
The value is explanatory — it says *why* a large part of the pseudogene class is source-ambiguous — not
evidentiary.

**The 56 non-parent matches (42.1% of 133) are an unresolved residual.** Out-of-frame retro-copy ORFs
and incidental short-peptide matches are both plausible. **We report this and do not explain it.**

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
