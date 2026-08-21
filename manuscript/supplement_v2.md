# Supplement

Companion to the main manuscript. Three items sit here rather than in the main text, for three different
reasons: **S1** is presented as an **exploratory supplementary analysis**, not a headline result;
**S2** is an
elementary derivation of a phenomenon Woo et al. 2014 established empirically — the closed form
itself is not established prior art (checked directly against Woo et al. 2014) — included because
the main text depends on its exact form; **S3** is a frozen, version-specific worked example from a
separate atlas, included to test whether sequence recurrence survives aggregation into a peptide–HLA
recurrence claim.

---

## S1. Where it can be interrogated, the ambiguity is structured by homology

**This quantifies a known phenomenon.** That a processed pseudogene shares sequence with its parent
gene, and that the resulting peptides are not attributable to one locus, is the classical
shared-peptide / protein-inference problem (Nesvizhskii & Aebersold). We measure it; we do not
discover it.

**Why this is in the Supplement.** Two methodological concerns apply to any symbol-heuristic parent
mapping and a naive permutation null: parentage **derived from the gene symbol** (strip a trailing `P` +
digits) is not authoritative, and a permutation null that treats parents as freely exchangeable ignores
**gene-family structure**. Both are addressed below — with an authoritative annotation and a
family-respecting null — and the result survives both. It stays here anyway, for two reasons the reader should weigh: the testable set is
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
peptide, **36 resolve to a curated parent (58.1%)**. The remaining 26 split into two different reasons,
not one: **25** are clone-style identifiers (`AC005262.1`, `AL158050.1`, …) absent from **every**
registry, and **1** (`GNG5P2`) resolves cleanly to an NCBI GeneID but that GeneID carries no curated
`gene_group` parent relation at all — a real gene symbol with no annotated parent, distinct from the
clone-style cases, which have no symbol-to-GeneID resolution to begin with. 36 + 1 + 25 = 62. Among
symbols that have a real gene name (37 = 36 + 1), coverage is **97.3%**.

**Of the 36 resolved-parent symbols, only 35 carry a testable peptide.** One symbol's curated parent is
`H2BC8`; Swiss-Prot has no reviewed entry under that symbol — not because the protein is missing, but
because UniProt collapses several histone-H2B genes that encode an identical protein into one reviewed
entry filed under a different symbol (`H2BC4`), which gene-ID-level matching cannot see through. That
symbol is therefore excluded from the denominator (the conservative call — it is not scored as a miss),
leaving **35** testable genes, matching the "33/35" gene-level row below.

**The caveat that matters, and it is not a small one.** HGNC *names* a pseudogene after the parent it
descends from. The curated relation is therefore **not an independent source of truth** about
parentage — it is substantially the same information the symbol heuristic was reading, curated. High
agreement between the two is **expected, and is not corroboration**. What the curated mapping actually
buys is a *versioned, adjudicated* relation with correct handling of renames and family cases — which
is the substantive improvement over the symbol heuristic — and nothing more. The **test** itself is untainted: whether a peptide
occurs inside the parent's protein is a pure sequence question, decided independently of any naming.

### The measurement

A processed pseudogene is a retro-copy of a parent gene, so a peptide encoded by a pseudogene ORF can
be an exact substring of the **parent protein**. Of the pseudogene-labelled IEAtlas peptides that are
canonical substrings, **133 are testable**:

| | authoritative | (symbol heuristic) |
|---|---:|---:|
| land in the pseudogene's **own curated parent** | **77 / 133 = 57.9%** | 79 / 132 = 59.8% |
| at gene level, pseudogenes with ≥1 peptide in their own parent | **33 / 35** | 33 / 34 |

**Table S1.** Pseudogene-to-parent hit rate under the authoritative NCBI Gene mapping versus the
symbol-heuristic mapping it replaces.

**How wrong was the heuristic?** Where both mappings yield a parent, they **agree in 124 / 127 = 97.6%**
of cases; it disagrees 3 times, and in 7 further cases the heuristic yields no parent where a curated
one exists. **Seven peptide-level verdicts flip.** The headline moves from 59.8% to **57.9%** — a small,
real correction: the symbol heuristic was directionally right, but not authoritative. We report both.

This head-to-head comparison population is **124 + 3 + 7 = 134 peptides, not 133**: it counts every
peptide whose pseudogene *symbol* resolves to a curated parent, regardless of whether that parent has
its own reviewed SwissProt entry. The **133**-peptide "testable" denominator used everywhere else in
this section additionally requires the parent to have a SwissProt entry, which excludes exactly the one
peptide noted above (`HIST1H2BPS2`, parent `H2BC8`) — so 134 = 133 + 1. That one peptide still has a
well-defined curated-parent *symbol* to compare against the heuristic even though it cannot be scored
for the hit-rate test itself.

### The null — and one that turned out to be degenerate

Permuting parent labels freely ignores gene-family structure, which can make "matched the parent"
trivially easy to satisfy by chance. We ran four nulls (10,000 permutations, seed 0) and report all of
them, including the one that failed:

| null | construction | result |
|---|---|---|
| **A** — naive / free | parents exchangeable across the whole pool (*the naive null*) | 77 obs vs **4.4** null mean; *p* < 1e-4 |
| **B** — HGNC gene-family | shuffle parents only within a curated family | restricted to its **5** permutable items: 5/5 (100.0%) obs vs **3.4** (67.8%) null mean; *p* = 0.20 (**not significant at this n**) |
| **C** — shared-9-mer component | shuffle only within a sequence-sharing component | **DEGENERATE** — see below |
| **D** — family-decoy swap | replace the true parent with a random **close paralog of that parent** | see below — *the null to read* |

**Table S2.** Four family-respecting nulls for the pseudogene→parent hit rate; only Null D is
non-degenerate.

**Null B is not degenerate — it is small, and we report it rather than folding it in with C's true
degeneracy.** Of the 35 parent proteins, only 5 belong to a curated HGNC family with another parent in
the pool at all; restricted to that permutable subset, all 5 (100%) land in their own parent against a
null mean of 3.4 (67.8%). At *n* = 5 this cannot reach significance (*p* = 0.20) and we do not claim it
does — it is reported for completeness, not as independent support. **Null C is the one that is
genuinely degenerate**, and for a reason that is itself the answer to the objection: **the 35 parent
proteins are pairwise 9-mer-disjoint**, so every homology stratum is a singleton and the
"family-respecting" permutation reduces to the identity. There is no family structure *among the
parents* for a shuffle to exploit. Reporting C's *p* = 1.0 as if it were evidence would be as dishonest
as hiding it.

**Null D is the family-respecting null that actually runs.** It asks the sharp version of the
objection: *could a merely homologous protein have been hit instead of the true parent?* Each parent is
replaced by a random **close paralog of itself** — the hardest decoys available — under two pools:

| decoy pool | true parent hit | random paralog hit | *p* |
|---|---:|---:|---:|
| any shared 9-mer (permissive) | 53 / 101 = **52.5%** | **7.1%** | < 1e-4 |
| **strong paralogs (≥10 shared 9-mers)** | 34 / 65 = **52.3%** | **16.6%** | < 1e-4 |

**Table S3.** Null D, the family-decoy swap, under two paralog-decoy pools.

The hit is **parent-specific, not family-generic** — it survives even against the parent's genuine close
homologs. Without any null at all: **58 / 77 = 75.3%** of parent-hits are compatible with the parent and
with **no** close paralog. The remaining 24.7% are not, and we say so; family structure is real and
quantified rather than denied.

**The objection was right for the reason it gave: curated mapping corrected one heuristic false
match.** The old
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

## S2. Set-identification of a class-specific FDR (elementary derivation; the underlying phenomenon is Woo et al. 2014)

**What kind of quantity *q* and *FDP_N* are, stated before the derivation, because "FDR" invites
conflating two different things.** Everything below is an algebraic identity about **realized
proportions in one fixed, already-accepted set of identifications**: it partitions an actual count of
false discoveries into two actual per-class counts, and needs no assumption about repeated draws,
exchangeability, or *f* being non-random. For one realized accepted set there is only one set — the
identity below holds by the arithmetic of a weighted average of two group proportions, not by a
probabilistic argument about FDR as an expectation over a hypothetical ensemble of experiments.
Separately: *q*, as a pipeline reports it, is itself typically a target–decoy **estimate** of that
realized proportion (e.g. a decoy-count-based ratio at a chosen threshold), not a directly observed
ground truth, and decoy-based FDR estimation carries its own assumptions and estimation noise. Θ_N(*q*,
*f*) below is sharp *given* *q* and *f* taken as exact; if *q* is itself estimated with error, that
error propagates into Θ_N and is not additionally bounded by anything here.

Let a target–decoy procedure accept a set of identifications at a reported pooled false-discovery rate
*q*. Partition the accepted set into a non-canonical class *N* and its complement, and let *f* be the
fraction of accepted identifications in class *N*.

Write *FDP_N* for the **realized false-discovery proportion within class *N*** — the actual (unknown)
fraction of accepted class-*N* identifications that are false, in this one accepted set — and *FDP_C*
for that in the complement. We use *FDP*, not *FDR*, throughout this derivation: *FDR* is an
expectation over a hypothetical ensemble of experiments, and *q* itself is typically a target–decoy
**estimate** of a realized proportion, not that proportion directly. The identity below is arithmetic on
realized counts and needs neither name to hold, but calling the unknowns *FDP* keeps it from reading as
a claim about the expectation. The pooled realized proportion is the *f*-weighted mixture:

    p = f · FDP_N + (1 − f) · FDP_C

with the only constraints being that both class-specific proportions are in [0, 1]:

    0 ≤ FDP_N ≤ 1,   0 ≤ FDP_C ≤ 1

Solving for *FDP_N* and imposing those bounds on *FDP_C* gives

    FDP_N = ( p − (1 − f) · FDP_C ) / f

which is decreasing in *FDP_C*. Substituting the two extremes *FDP_C* = 1 and *FDP_C* = 0, and
intersecting with [0, 1]:

    Θ_N(p, f) = [ max(0, (p − (1 − f)) / f),  min(1, p / f) ]

**Sharpness.** Every point in Θ_N is attained by some admissible *FDP_C* ∈ [0, 1], and no point outside
it is: the endpoints correspond exactly to *FDP_C* = 1 and *FDP_C* = 0. So **given *p* and *f* alone,
Θ_N is the sharp identified set for the realized class-*N* false-discovery proportion** — the
data-generating process is not pinned down more tightly by that information. Substituting a *reported* q
for *p* evaluates Θ_N at an estimate of the realized proportion, not at the proportion itself; if *q* has
estimation error, that error is not additionally bounded by anything here (see below).

**What this does *not* say.** Θ_N is sharp *with respect to p and f*, evaluated in practice at *q* as a
stand-in for *p*. It is **not** the claim that no further information could narrow it. Calibrated
class-specific posterior error probabilities, entrapment measurements, or a validated mixture model
would each add information and could tighten the interval. The per-class accepted decoy count *D_N* is
simply the **cheapest sufficient** such object, and one the pipeline already computes: with *T_N* and
*D_N* and the stated threshold, unit and convention, the selected class-specific target–decoy estimate
becomes **reconstructible**. *D_N* does not identify the true class-specific false-discovery
*proportion*, and we do not claim it does.

**Why it bites here.** IEAtlas reports *q* = 0.05 as a PSM-level target–decoy threshold, but publishes
no PSM-level canonical/non-canonical accept counts at all — not even the class-summed, cross-tissue
245,870 epitope total is at the right unit for this (main text, R4). *f* is therefore unknown by any
route we have found, and Θ_N(0.05, *f*) is **unconstrained** — for small *f* the upper endpoint
min(1, *q*/*f*) reaches 1. The interval cannot be evaluated at all from what the resource publishes.
That is the reporting gap, and it is closed by one small table.

**Worked illustration.** At *q* = 0.05:

| *f* (non-canonical share of accepted IDs) | Θ_N(0.05, *f*) |
|---:|---|
| 0.50 | [0.00, 0.10] |
| 0.20 | [0.00, 0.25] |
| 0.05 | [0.00, 1.00] |
| 0.01 | [0.00, 1.00] |

**Table S4.** Worked values of the set-identified class-specific FDR interval Θ_N(0.05, f) at four
plausible class fractions.

The scarcer the class, the less a pooled threshold says about it — which is the entire point of Woo et
al. 2014, restated here only because the resource under audit does not report the quantity that would
resolve it.

---

## S3. Recurrence and scan-level provenance in a frozen ImmunoVerse pilot

### Scope and estimand

This is a **worked external example**, not a pan-atlas prevalence estimate and not an audit of the
live 2026 ImmunoVerse portal. The catalogue and sample map are frozen to the 2025-07-07 preprint
supplements: Table S7 (`ORF_antigen`) and Table S3 (source study, raw file, biological label and
reported HLA genotype). Table S7 contains 17,741 source rows representing **7,770 distinct peptides**.
The current raw-result release post-dates those tables and is used only to
test evidence recoverability in two bounded cancer subsets.

Recurrence is reported at four distinct units: sample/condition label, published source-study
identifier, peptide sequence, and predicted or genotype-compatible peptide–HLA. A source-study label
is not assumed to be an independent patient, cohort or experiment. Table S7's HLA entries are binding
predictions, not direct allele assignments.

### Catalogue-to-study lineage

Table S3 contains **1,771 rows**, **1,679 distinct raw-file names**, **498 biological labels**, and
**47 source-study identifiers**. All **451 / 451** biological labels used by Table S7, and all
**465 / 465** cancer–label pairs, joined exactly to Table S3. Six used labels map to more than one
source study; their study recurrences were bounded rather than resolved by convention. Two labels
carry conflicting reported genotype variants. After formatting normalization, 33,707 Table S7
prediction entries match every reported genotype for their label, 41 match some but not every
reported genotype variant, and none mismatch all reported genotypes.

| Quantity, top 1,000 sequence-recurrent peptides | Result |
|---|---:|
| median sample/condition-label recurrence | **11** |
| median compatible source-study recurrence | **7–7** |
| median conservative label/source-study ratio | **1.667** |
| label recurrence exceeds the study upper bound | **97.2%** |
| median conservative sequence-study/best-predicted-pHLA-study ratio | **1.250** |
| conservative pHLA ratio > 1 | **72.2%** |
| conservative pHLA ratio ≥ 2 | **16.3%** |

**Table S5.** Source-study normalization under the primary strong-or-weak predicted-binder rule. The
conservative pHLA ratio divides the minimum compatible sequence-study recurrence by the maximum
compatible study recurrence of the best predicted pHLA.

The top-1,000 boundary cuts through a 240-peptide tie at seven labels. In the tie-complete set
(*n* = 1,184), the median conservative pHLA ratio remains 1.250 and 72.1% remain above one. In the
732-peptide exact-study-mapping subset, the median sequence recurrence is 6 studies, the best predicted
pHLA recurs in 5, the median ratio is 1.333, and 81.6% exceed one. Under strong binders alone, 913 of
the top 1,000 have a predicted pHLA; the median conservative ratio is 1.333 and 73.8% exceed one.

### Current raw-result recovery: a positive control

The current release publishes MaxQuant, rescoring and consolidated scan outputs separately from the
catalogue. In BLCA, the consolidated table contains 106,026 scan rows. Every historical Table S7
sequence and count was recoverable:

| BLCA recovery check | Result |
|---|---:|
| exact sequence recovery | **16 / 16 peptides** |
| historical total `n_psm` / current exact-sequence rows | **18 / 18** |
| exact biological-label-set recovery | **16 / 16 peptides** |
| exact per-peptide PSM-count agreement | **16 / 16 peptides** |
| at least one current `Identified = +` scan | **15 / 16 peptides** |
| at least one current `Reverse = +` scan | **1 / 16 peptides** |

**Table S6.** BLCA positive-control recovery. Current decision flags are reported as current state,
not assumed to reproduce the historical acceptance policy.

This complete recovery shows why ImmunoVerse is a useful positive implementation example: the
catalogue itself projects away the evidence joins, but the separately released raw-result chain makes
them reconstructible.

### DLBC: recurrence of a sequence is not recurrence of one pHLA

The current DLBC consolidated table contains **51,446 scan rows**. Of 73 historical DLBC Table S7
peptides, **72 / 73** were recovered by exact sequence and exact biological-label set. Historical
`n_psm` summed to **194**, compared with **199** current exact-sequence scan rows; 68 / 73 peptides had
exact per-peptide count agreement, 64 / 73 had at least one current `Identified = +` scan, and 11 / 73
had at least one current `Reverse = +` scan. These differences are decision/version provenance, not
evidence that the historical catalogue was wrong.

All three multi-label sequences were recovered in every catalogue-listed label. Their complete
reported genotype intersections and their common predicted-binder intersections were empty:

| peptide | Table S7 labels | historical `n_psm` | current exact-sequence rows | current identified rows by label | common reported HLA |
|---|---|---:|---:|---|---|
| `AEGPDHHSL` | DOHH2, SUDHL4 | 5 | 5 | DOHH2: 4; SUDHL4: 1 | **none** |
| `VPHTRPVSL` | DOHH2, HBL1, SUDHL4 | 24 | 24 | DOHH2: 6; HBL1: 7; SUDHL4: 7 | **none** |
| `LASPHSPIL` | DOHH2, HBL1, SUDHL4 | 13 | 15 | DOHH2: 0; HBL1: 0; SUDHL4: 0 | **none** |

**Table S7.** Raw-linked recurrent DLBC sequences. `AEGPDHHSL` and `VPHTRPVSL` each have a current
accepted scan in every listed cell line, but the cell lines share no reported HLA allele. Conditional
on those reported genotypes, one fixed peptide–HLA allomorph cannot explain recurrence across all
contexts. `LASPHSPIL` is retained as a decision-version example: exact scans remain, but none carries
the current consolidated acceptance flag.

### Limits

This pilot does not establish independent-patient recurrence, the presenting allele for an individual
scan, source-locus identity, malignant-cell specificity, T-cell recognition, safety, efficacy, or a
catalogue-level FDR. It does establish the narrower estimand correction: a peptide sequence can recur
across contexts in which one common pHLA target is impossible under the published genotypes. The
separately published scan evidence makes this correction auditable; it should be propagated into the
catalogue rather than left for every downstream user to reconstruct.
