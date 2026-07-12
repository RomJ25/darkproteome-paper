# The field's standard for non-canonical peptides, applied to the resources that feed it: one atlas fails it by an order of magnitude

*Working draft. Supersedes `manuscript.md`, whose headline was withdrawn (`REVISION_NOTICE.md`).
Every number is regenerated from the analysis artifacts by `verify_manuscript.py`, which fails the
build on drift. Every quotation from a primary Methods section has been verified in the fetched full
text.*

> ## Framing note — load-bearing, and it should be read first
>
> **We claim no conceptual novelty. Every idea in this paper is already in the literature.**
>
> - That a peptide matching a canonical protein does not identify a non-canonical source, and must be
>   **excluded**: **Bedran et al. 2023** (*Cancer Immunol Res*) — *"Sequences perfectly matching any
>   protein sequence were considered exonic… required to have at least three mismatches with any known
>   protein sequence before being considered noncanonical."* They also measured the residual canonical
>   overlap of four published catalogues (1.4–5%).
> - That shared peptides create source ambiguity and should be dropped: **Kumar et al. 2022**
>   (*Brief. Bioinform.*) — *"Most shared peptides should be dropped if defining… novel-coding
>   regions."*
> - The underlying protein-inference / shared-peptide problem: **Nesvizhskii & Aebersold**, textbook.
> - That a pooled FDR under-controls the minority class, and class-specific FDR is required:
>   **Woo et al. 2014**; Zhang et al.; **pXg** (Choi & Paek 2024).
>
> **What has not been done is to test whether the public resources meet the standard.** Kumar et al.
> state plainly that the problem is *not quantified* and that *no specific databases are criticised by
> name*; Bedran et al., who did quantify, did not include IEAtlas. **This paper is that test.** Its
> contribution is measurements, not ideas.

---

## Abstract

Catalogues of non-canonical (ncORF-derived) HLA-presented peptides are used to nominate
cancer-vaccine targets. A peptide sequence also encoded by a canonical human protein does not
uniquely identify a non-canonical source: tandem MS identifies a *sequence*, and when the same
sequence is encoded by both a canonical protein and an ncORF, the spectrum does not choose between
them. The field's standard is therefore to **exclude** such sequences before calling a peptide
non-canonical (Bedran et al. 2023; Kumar et al. 2022), and the residual overlap of four published
catalogues has been measured at 1.4–5%.

We apply that standard to **IEAtlas**, a public atlas of 245,870 non-canonical HLA epitopes that has
not previously been audited. **At least 56.3%** of its unique cancer-catalogued peptide sequences
(98,193 / 174,465) also occur in reviewed canonical human proteins under an explicitly defined
reference *R* — **11–40× higher than any ncORF catalogue previously audited**, and 300–2,000× higher
than two resources that apply an explicit exclusion rule (CrypticProteinDB, 0.026%; Raja et al.,
0.17%). Neither ORF-class composition nor the false-discovery rate accounts for it.

**The library does.** Measuring the HLA-I-length peptide space of the ncORF libraries themselves,
**34.1%** of nuORFdb v1.2's distinct 9-mers also occur in canonical human proteins, against **1.0–2.4%**
for the GENCODE Ribo-seq ORF sets — a 14–34× difference in *latent* ambiguity. IEAtlas integrates
nuORFdb and applies no exclusion rule. Ouspenskaia et al. searched the *same* nuORFdb and published a
catalogue at 3%, so the library is necessary but not sufficient: the exclusion step does the work.

**The consequence is internal to the resource.** If a catalogued "cancer epitope" is a peptide of an
abundant canonical protein, it should also be presented on normal tissue. It is: canonical-overlapping
cancer epitopes appear in IEAtlas's **own** normal-tissue set at **22.4%**, versus **9.1%** for the
non-overlapping epitopes of the same catalogue (risk ratio 2.4, *z* = 74). **22,003 entries — 12.6% of
every cancer epitope the atlas catalogues — are both canonical-compatible and already observed on
normal tissue by the atlas's own measurement.** No external reference is required to see this.

Applying the field's own exclusion standard would remove 56.3% of the atlas. For scale, class-specific
FDR control cost Ouspenskaia et al. 24% of their nuORF peptides and up to 76% of one ORF class.
Applying an established standard to an ncORF catalogue is expected to be expensive; that is not an
argument against applying it.

**Nothing here shows that any individual entry is canonically derived, or that non-canonical
translation is not real.** It shows that a widely-used catalogue cannot, from what it publishes, tell a
reader which of its entries are source-ambiguous — and that more than half are.

---

## Introduction

Non-canonical open reading frames (ncORFs) — in long non-coding RNAs, pseudogenes, untranslated
regions, and alternative frames of coding genes — yield peptides that are presented on HLA molecules
and recognised by T cells. Because their products are absent from the annotated proteome, they have
been proposed as an unusually attractive class of tumour antigen: potentially tumour-restricted, shared
across patients, and not subject to central tolerance. Several public catalogues now aggregate tens or
hundreds of thousands of such epitopes, and these catalogues are the practical input to target
selection.

Building such a catalogue requires answering one question for every identified peptide: **is this
peptide non-canonical?** The question is harder than it appears, and the field has known why for two
decades.

**Tandem MS identifies a sequence, not a source.** A peptide-spectrum match establishes an amino-acid
sequence. It does not establish which genomic locus produced it. Where two loci encode the same
sequence — a canonical protein and an ncORF — no spectral evidence distinguishes them. This is the
protein-inference problem (Nesvizhskii & Aebersold), and in proteogenomics it is *"exacerbated by the
mapping complexity where many identified peptides map to several loci, both novel and known"*
(Kumar et al. 2022).

**The remedy is established and unambiguous.** Kumar et al. (2022) state that *"most shared peptides
should be dropped if defining… novel-coding regions."* Bedran et al. (2023) implement the rule
directly: sequences perfectly matching any known protein are treated as *exonic* and excluded, with a
stringent requirement of *"at least three mismatches with any known protein sequence"* before a peptide
is called non-canonical. A parallel and additive standard governs statistical confidence: a pooled
target–decoy threshold under-controls a minority class, so class-specific FDR estimation is required
(Woo et al. 2014, who measured 36% novel-class FDR at a 1% combined threshold; Zhang et al.; Choi &
Paek 2024).

**What has not been done is to check whether the public resources meet these standards.** Kumar et al.
review the problem without quantifying it, and note that no specific databases are criticised by name.
Bedran et al. did quantify — reporting residual canonical overlap of 1.4% (Erhard et al. 2020), 3%
(Ouspenskaia et al. 2021), 4% (Chong et al. 2020) and 5% (Laumont et al. 2016) — but their comparison
did not include IEAtlas, the largest such atlas and the one most readily used as an off-the-shelf
source of candidate antigens.

This paper is that check. We apply the field's own metric, against an explicitly stated reference, to
IEAtlas and to two catalogues that apply an explicit exclusion rule; we apply the same measurement to
the ncORF *libraries* those catalogues were searched against; and we ask what the resulting ambiguity
implies for a user selecting tumour-specific targets. **We contribute no new concept, method or
standard. We contribute measurements that had not been made.**

---

## Results

### R1. IEAtlas is an order-of-magnitude outlier on the established metric

**The measurement** (exact substring of a reviewed canonical human protein; unique peptide sequences;
reference *R* pinned by release and hash in Methods):

| catalogue | exclusion rule applied? | canonical-sequence overlap |
|---|---|---:|
| Erhard et al. 2020 † | no | 1.4% |
| Ouspenskaia et al. 2021 † | no | 3% |
| Chong et al. 2020 † | no | 4% |
| Laumont et al. 2016 † | no | 5% |
| CrypticProteinDB | **yes** — *"BLASTP… eliminate all proteins with alignment to canonical proteins"* | **1 / 3,810 = 0.026%** |
| Raja et al. (ovarian) | **yes** — *"peptides mapping to 'protein_coding'… were excluded"* | **5 / 2,979 = 0.17%** |
| **IEAtlas** | **no** — *"only epitopes derived from non-coding regions were retained"* | **98,193 / 174,465 = 56.3%** |

† as reported by Bedran et al. 2023; not re-measured here.

**IEAtlas is 11–40× higher than any ncORF catalogue previously audited**, and 300–2,000× higher than
the two resources that apply an explicit exclusion rule.

#### It is not the ORF-class composition

IEAtlas's ncORF library explicitly includes pseudogenes, the class with the highest canonical
compatibility. Holding the class fixed:

| pseudogene-derived peptides | rule | canonical overlap |
|---|---|---:|
| IEAtlas | retention | **9,874 / 16,323 = 60.5%** |
| Raja et al. | exclusion | **0 / 98 = 0.0%** |

Raja et al. *do* report 98 pseudogene-ORF peptides; none overlaps a canonical protein, because their
exclusion rule removed those that did. Same ORF class, opposite rule, 60.5% versus 0.0%. Moreover the
overlap is **not** a pseudogene artifact: within IEAtlas, **non**-pseudogene ORF types show **56.0%**
overlap (88,827 / 158,688), within a few points of its pseudogenes. Class composition explains almost
none of the 56.3%.

#### It is not the false-discovery rate

A looser FDR admits more *false* identifications, and a false identification is not preferentially a
canonical substring: a composition-matched shuffle places chance canonical overlap near 0.1%. No FDR
can produce a 56% canonical-substring rate. IEAtlas does apply a permissive threshold — a 5% PSM FDR
and, explicitly, *"no protein FDR was set"* — but this cannot be the explanation.

### R2. The magnitude is explained by the library

Applying the same measurement to the **search space** rather than the output — the distinct 9-mers of
each ncORF library, and how many also occur in reviewed canonical human proteins (**full libraries, no
sampling**):

| ncORF library | ORFs | distinct 9-mers | also canonical |
|---|---:|---:|---:|
| **nuORFdb v1.2** — integrated by IEAtlas | 229,251 | 8,448,245 | **34.1%** |
| GENCODE Ribo-seq ORFs (phase 1) | 7,264 | 245,094 | **2.4%** |
| GENCODE Ribo-seq ORFs (phase 2) | 28,359 | 502,528 | **1.0%** |

**ncORF libraries differ by 14–34× in latent canonical ambiguity** — and the GENCODE range (1.0–2.4%)
is precisely where the catalogues Bedran et al. audited sit (1.4–5%). Whole-ORF containment is low
throughout (0.2–0.8%), so this is **extensive partial sharing**, not whole ncORFs nested inside
canonical proteins.

**Independent corroboration.** These figures were obtained twice, by separate implementations over
different k-mer windows: an 8–11mer candidate-universe enumeration gives nuORFdb **34.4%** and GENCODE
Ribo-seq (phase 1) **2.5%**; the 9-mer enumeration reported above gives **34.1%** and **2.4%**. The
agreement is close, and the measurement is not sensitive to the window.

**Two factors are jointly necessary.** The **library creates** the ambiguity: a third of nuORFdb's
peptide space is canonical, and nothing downstream of the search except a peptide-level exclusion can
undo it. The **pipeline decides whether it survives**: Ouspenskaia et al. searched the *same* nuORFdb
and published a catalogue at 3%. IEAtlas is the only resource examined that combines a high-ambiguity
library with no exclusion rule.

*(Why the catalogued 56.3% exceeds the library's 34.1% is a hypothesis we do not test: canonically
encoded peptides derive from abundant proteins and are plausibly over-detected in an immunopeptidome
relative to their share of the search space.)*

### R3. The consequence is internal to the resource

An audit that stops at *"56.3% of these entries are source-ambiguous"* invites the only question that
matters: **so what?**

If a catalogued "cancer epitope" is in fact a peptide of an abundant canonical protein, it should also
be presented on **normal tissue**, because that protein is expressed there too. IEAtlas publishes its
own normal-tissue epitope set (94,375 unique peptides), so the prediction is testable **inside the
resource**, with the non-overlapping epitopes of the same catalogue as an internal control.

| IEAtlas cancer epitopes | also in IEAtlas's **own** normal-tissue set |
|---|---:|
| **canonical-overlapping** (98,193) | **22,003 = 22.4%** |
| non-overlapping — internal control (76,272) | 6,976 = **9.1%** |

**Risk ratio 2.4** (two-proportion *z* = 74; *n* = 174,465). A source-ambiguous "cancer" epitope is more
than twice as likely to appear in the atlas's own normal-tissue set — exactly what is expected if it is
a peptide of a canonical protein that is also expressed in normal tissue.

**And a subset requires no inference at all.** **22,003 entries — 12.6% of every cancer epitope the
atlas catalogues — are both canonical-compatible and already observed on normal tissue by the atlas's
own measurement.** No external reference is needed to see this; it is internal to the resource. A
target-selection pipeline drawing from this catalogue without an exclusion step will encounter them.

This does **not** show that any individual epitope is canonically derived. Presence in the normal set is
evidence about **presentation**, not about **source**. What it shows is that the source-ambiguous
fraction *behaves* as canonically derived peptides would.

### R4. Where it can be interrogated, the ambiguity is structured by homology

**This quantifies a known phenomenon.** That a processed pseudogene shares sequence with its parent
gene, and that the resulting peptides are not attributable to one locus, is the classical shared-peptide
/ protein-inference problem (Nesvizhskii & Aebersold). We measure it; we do not discover it.

A processed pseudogene is a retro-copy of a parent gene, so a peptide encoded by a pseudogene ORF can be
an exact substring of the **parent protein**. Of 213 pseudogene-labelled peptides that are canonical
substrings, 132 are testable (78 carry clone-style symbols with no derivable parent; 3 have no parent in
the reference). Of those, **79 (59.8%)** land in the pseudogene's **own parent gene** (73 by exact
symbol; 6 under a documented rename, e.g. `FAM8A6P → FAM8A1`). At gene level, **33 of 34 pseudogenes
(97.1%)** have at least one peptide in their own parent.

**Null.** A uniform "one gene in ~20,400" null is invalid: pseudogene parent annotations are themselves
homology-derived, and we select peptides already known to match *something* canonical. We instead hold
each peptide's canonical hit set **fixed** — preserving the selection, protein lengths, paralogy and
k-mer structure — and permute only the pseudogene→parent pairing. Observed 79; null mean 4.9 (max 14 in
10,000 permutations); **p < 10⁻⁴**.

**This does not resolve provenance.** `DEVAFRKF` is encoded by both *RPS3AP12* and *RPS3A*; MS
identifies the sequence, not the locus. Even a parent hit is rarely a unique assignment: of the 79, only
52 (65.8%) match the parent *alone*, while 27 (34.2%) match the parent **and further canonical genes** —
more ambiguity, not less. The remaining **53 / 132 (40.2%)** land in a canonical gene that is not the
derived parent; out-of-frame retro-copy ORFs and incidental short-peptide matches are both plausible,
and parent-hits and non-parent hits share a median length of 9 aa, so length does not separate them.
**We report this as an unresolved residual and do not explain it.**

Across classes we report only descriptive heterogeneity in canonical overlap. The retro-copy mechanism
is claimed for **processed pseudogenes only** — the one class where the label is corroborated by
sequence. Low overlap in lncRNA-ORF (0.5%) and altORF (0.2%) classes argues against one specific failure
mode — wholesale pseudogene-like contamination of those classes — and nothing more; it does not validate
those annotations.

### R5. The additive statistical problem, and the remedy the field already demonstrated

**Not novel** (Woo et al. 2014; Choi & Paek 2024). Stated because it is *additive* to R1–R3, and because
IEAtlas reports nothing that would allow it to be assessed.

From a reported pooled FDR *q* and class fraction *f*, the class-specific FDR is only **set-identified**:

    Θ_N(q, f) = [ max(0, (q − (1 − f)) / f),  min(1, q / f) ]

The interval is sharp and cannot be tightened without the per-class accepted decoy count *D_N*. IEAtlas
reports *q* = 0.05 and 245,870 non-canonical epitopes, but **no canonical count**. *f* is therefore
unknown and the interval is **unconstrained**.

**The field has already demonstrated the remedy, and its cost.** Ouspenskaia et al. searched a combined
annotated-ORF/nuORF database and reported that a 1% global FDR gave *"4.6% overall, and as high as 14%
for 3′ dORFs"* among nuORF peptides; group-based filtering *"removed 24% of nuORF peptides overall, and
up to 76% of peptides assigned to 3′ overlap dORFs."* This is the predicted inflation, measured in the
literature by authors who then corrected it. **The 4.6% is a pre-correction diagnostic, not the error
rate of their published catalogue.** They are a positive exemplar.

### R6. Two distinct problems; one reporting remedy; and what it costs

The two defects are **complementary, not the same**, and merging them is a category error:

- **Source ambiguity (R1–R4).** The peptide is *correctly identified as a sequence*; its **source
  locus** is unresolved. Not an FDR problem — no identification is wrong.
- **Class-specific FDR under-control (R5).** A pooled threshold under-controls the minority class, so
  some ncORF **identifications are wrong** — the spectrum was not that peptide.

**A minimal standard addresses both.**

**(a) Per peptide — an exclusivity flag and all compatible source loci.** State whether the sequence is
unique to the nominated ncORF within the searched space, and if not, list every compatible source. This
preserves the peptide *and* the ambiguity, rather than discarding either, and it is strictly more
informative than the exclusion rule it generalises.

**(b) Per class — a class-decoy ledger.** Accepted target and decoy counts for each class, the class
definitions, the thresholding stage, and the formula used. This is exactly what collapses Θ_N to a
point, and Ouspenskaia et al. demonstrate it is achievable.

**(c) Per library — publish the latent canonical ambiguity.** One number (R2), cheap to compute, and it
tells every downstream user how much an exclusion rule matters for the library they are about to search.
**No ncORF library currently publishes it.**

**The cost, stated plainly.** Applying the exclusion standard to IEAtlas removes **98,193 of 174,465
entries (56.3%)**. For scale, class-specific FDR control cost Ouspenskaia et al. 24% of their nuORF
peptides and up to 76% of one ORF class. **Applying an established standard to an ncORF catalogue is
expected to be expensive. That is not an argument against applying it** — it is an argument for applying
it before, rather than after, a candidate enters a clinical pipeline.

---

## Discussion

**This is not an accusation, and should not be read as one.** Every resource examined describes its own
procedures accurately and in public; that is the only reason this audit was possible at all. Two of the
catalogues measured here already apply an exclusion rule. Ouspenskaia et al. already solved the
statistical half and published the cost of doing so. **The field knows how to do this.** What has not
happened is a check that the resources in use actually do it — and the reason is prosaic: nobody had
looked.

**The unexplained becomes explained, and then becomes actionable.** The 56.3% is not a property of
IEAtlas's filtering alone; it is inherited from a library whose peptide space is a third canonical by
sequence. That reframing matters, because it moves the remedy upstream. A downstream exclusion rule
works — CrypticProteinDB and Raja et al. demonstrate it, and Ouspenskaia et al. demonstrate it on the
*same* library — but a library that publishes its own latent ambiguity lets every downstream group know
in advance how much that rule will cost them. **That number does not currently exist for any ncORF
library, and it takes minutes to compute.**

**The stakes are not abstract.** A peptide catalogued as a cancer epitope, which is compatible with an
abundant canonical protein, and which the same atlas has already observed on normal tissue, is not a
promising tumour-restricted target: it is an on-target/off-tumour risk. **22,003 entries in IEAtlas meet
all three conditions.** We are not claiming that any of them *is* canonically derived — MS cannot say —
but a target-selection pipeline that cannot distinguish them is selecting under an ambiguity it has not
been told about.

**Limits, stated plainly.**

- The overlap is **reference-relative**. It is *N*(*R*) for the reference stated in Methods. A broader
  reference can only *raise* it, so 56.3% is a **lower bound**; it cannot be lowered by a different
  reference choice.
- **No individual peptide's provenance is resolved.** Sequence identity is symmetric. Presence in a
  normal-tissue set is evidence about presentation, not source.
- We hold nuORFdb but **not** RPFdb or Translnc, the other two libraries IEAtlas integrates. **34.1% is
  therefore a lower bound** on its library's latent ambiguity.
- The abundance-bias explanation for why the catalogued rate (56.3%) exceeds the library rate (34.1%) is
  a **hypothesis we do not test**.
- The 40.2% of pseudogene peptides landing in a non-parent canonical gene is an **unresolved residual**.
- Class labels for non-pseudogene classes are **source-supplied and uncorroborated**.

**What we would like to be wrong about.** If IEAtlas's pipeline resolves source attribution in a way its
Methods do not describe — for instance through a protein-inference step that assigns shared sequences to
the canonical protein and excludes them — then the 56.3% has an innocent explanation, and this paper is a
correction to a misreading of the Methods rather than a finding about the resource. The search outputs
(`peptides.txt`: `Proteins`, `Leading razor protein`, `Unique (Proteins)`) would settle it in an
afternoon. **We invite the authors to publish them.**

---

## Methods

### Reference

Canonical human proteome: UniProt/Swiss-Prot reviewed human sequences (*R*), release and SHA-256 hash
recorded in `data/SOURCES.md` and printed by `reference_provenance.py`. **All overlap statistics are
`N(R)` and are reference-relative.** IEAtlas searched a February-2022 Swiss-Prot release; our reference
is later, so we make no claim about which canonical proteins were present in *their* search space. A
broader reference can only add matches, so every figure reported here is a **lower bound**.

### The overlap measurement

A catalogued peptide sequence counts as canonical-overlapping if it is an **exact substring** of at least
one protein in *R*. Sequences are compared as unmodified amino-acid strings, uppercased, with inline PTM
annotations stripped. Peptides are **deduplicated to unique sequences** before rates are computed; a rate
is `unique canonical-overlapping / unique scored`. The identical procedure is applied to every catalogue,
so the comparison is internally consistent even where the source studies used different references
themselves.

### The library measurement

For each ncORF library, all distinct **9-mers** (a representative HLA-I ligand length) are enumerated
across the **full library — no sampling** — and intersected with the set of all distinct 9-mers in *R*.
Sampling is invalid here and the bias is large and upward: sampling 4,000 of nuORFdb's 229,251 ORFs gives
43.6%, 20,000 gives 40.7%, and the full library gives 34.1%, because a small sample contains fewer
distinct ncORF-specific k-mers and so over-weights the canonical-shared ones. Whole-ORF containment *is*
sampled (4,000 ORFs, seed 0); that statistic is an O(*n·m*) substring scan and is unbiased under random
sampling, unlike a k-mer union.

### The normal-tissue consequence

IEAtlas's cancer and normal epitope exports are compared as unique sequences. The two-proportion *z* is
computed on the pooled proportion. The **non-overlapping epitopes of the same catalogue** serve as the
internal control, so the comparison is within-resource and requires no external normal-tissue reference.

### The pseudogene→parent test

The putative parent is derived from the **gene symbol** (strip a trailing `P` + digits), never from the
ORF-class label. A hit counts as the parent if the matched gene equals the derived parent, or if the two
share a root of ≥3 characters with a purely numeric remainder (accommodating documented renames:
`PHB → PHB1`, `FAM8A6 → FAM8A1`, `MKRN4 → MKRN1`). Peptides whose symbol yields no parent, or whose
derived parent is absent from *R*, are **untestable and excluded from the denominator** — scoring them as
misses would invent evidence against the class.

**Null model.** Each peptide's canonical hit set is held **fixed** (preserving the selection, protein
lengths, paralogy and k-mer structure) and only the pseudogene→parent pairing is permuted, 10,000 times,
seed 0.

### Class-specific FDR identifiability

Derivation of Θ_N(*q*, *f*) and its sharpness in the Supplement. The result is not new (Woo et al. 2014);
it is stated here because IEAtlas reports *q* but not *f*, so the interval is unconstrained for that
resource.

### Reproducibility

Every headline number in this manuscript is regenerated from the committed artifacts by
`manuscript/verify_manuscript.py`, **which fails the build on drift** and additionally fails if the paper
drops its required prior-art citations or asserts any of the phrasings retracted during review. The
analysis code is guarded by `src/darkproteome/scoring_conformance.py`. All primary Methods quotations
were verified in the fetched full text (EuropePMC / NCBI E-utilities); the fetch scripts and document
hashes are in the repository.

---

## What this paper does not claim

- **Any conceptual novelty whatsoever.** The exclusion standard is Bedran et al. 2023; the shared-peptide
  remedy is Kumar et al. 2022; the inference problem is Nesvizhskii & Aebersold; the class-FDR problem is
  Woo et al. 2014. **This paper contributes measurements, not ideas.**
- That the biology is fake, or that any specific ncORF antigen is not real.
- That any resource "manufactures", "re-labels", or "discards" anything.
- That the canonical source is the correct one for any peptide. **Sequence identity is symmetric; MS
  identifies the sequence, never the locus.**
- That any resource acted improperly. All describe their own procedures accurately.
- That presence in a normal-tissue set establishes a peptide's *source*. It is evidence about
  **presentation**.
- That the 56.3% is reference-independent. It is `N(R)`; a broader reference can only raise it.
- That source ambiguity and class-FDR under-control are the same phenomenon.
- That the abundance-bias explanation for 56.3% > 34.1% is tested. **It is a hypothesis.**
- That we have measured IEAtlas's whole library. We hold nuORFdb but not RPFdb or Translnc, so **34.1%
  is a lower bound.**
