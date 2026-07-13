# Extensive canonical-sequence overlap and unresolved source attribution in a public non-canonical HLA-peptide atlas and its search library

*Working draft. Supersedes `manuscript.md`, whose headline was withdrawn (`REVISION_NOTICE.md`).
Every number is regenerated from the analysis artifacts by `verify_manuscript.py`, which fails the
build on drift. Every quotation from a primary Methods section has been verified in the fetched full
text.*

> ## Prior art, stated once
>
> The underlying principles here are established, and we claim none of them.
>
> - That a peptide matching a canonical protein does not identify a non-canonical source, and should
>   be **excluded** before a peptide is called non-canonical: **Bedran et al. 2023** (*Cancer Immunol
>   Res*) — *"Sequences perfectly matching any protein sequence were considered exonic… required to
>   have at least three mismatches with any known protein sequence before being considered
>   noncanonical."* They also measured the residual canonical overlap of four published catalogues
>   (1.4–5%).
> - That shared peptides create source ambiguity and should be dropped when defining novel coding
>   regions: **Kumar et al. 2022** (*Brief. Bioinform.*) — *"Most shared peptides should be dropped
>   if defining… novel-coding regions."*
> - The underlying protein-inference / shared-peptide problem: **Nesvizhskii & Aebersold**, textbook.
> - That a pooled FDR under-controls a minority class, so class-specific FDR is required: **Woo et al.
>   2014**; Zhang et al.; **pXg** (Choi & Paek 2024).
>
> **Our contribution is empirical:** an audit of whether a major public resource satisfies these
> principles, a measurement of the sequence ambiguity latent in its peptide catalogue and in one of
> its source libraries, a within-resource test of the consequence, and a concrete reporting field
> that would fix it. Kumar et al. state that the problem is *not quantified* and that *no specific
> databases are criticised by name*; Bedran et al., who did quantify, did not include IEAtlas.

---

## Abstract

Catalogues of non-canonical (ncORF-derived) HLA-presented peptides are used to nominate
cancer-vaccine targets. A peptide sequence that is also encoded by a canonical human protein does not
uniquely identify a non-canonical source: tandem MS identifies a *sequence*, and where the same
sequence is encoded by both a canonical protein and an ncORF, the spectrum does not choose between
them. A published and recommended criterion is therefore to **exclude** such sequences before calling
a peptide non-canonical (Bedran et al. 2023; Kumar et al. 2022).

We audit **IEAtlas**, a public atlas of non-canonical HLA epitopes that has not previously been
examined against this criterion. **98,193 of 174,465 unique cancer-catalogued peptide sequences
(56.3%)** exactly match at least one protein in a frozen reviewed human reference *R*.

**The matching canonical proteins were not external to IEAtlas — they were inside its own search
database.** Its Methods state that spectra were searched against *"both our integrated benchmarked
ncORF library **and the canonical human proteome** … (downloaded in February 2022),"* after which
*"only epitopes derived from non-coding regions were retained."* Re-scored against that same
canonical proteome (Swiss-Prot **2022_01**), **97,999 sequences (56.2%)** match a database entry that
was present when the spectrum was assigned, competing with the ncORF the peptide was attributed to.
The finding is therefore **internal to the resource**, not a retrospective comparison against a
reference it never consulted: only 231 sequences (0.24%) are overlaps a February-2022 analyst could
not have made. It is not an artifact of peptide length (**56.3%** length-standardized), nor of
ORF-class composition (**55.8%** with no pseudogene ORFs at all). Under one pipeline, one reference
and one peptide unit, two catalogues that apply an explicit exclusion rule sit at **1 / 3,810
(0.026%)** and **5 / 2,979 (0.17%)**.

Such matches do not establish canonical production. They do prevent unique attribution to the
nominated non-canonical source from peptide sequence and ordinary tandem-MS evidence alone.

**Latent ambiguity differs enormously between search libraries.** Enumerating the HLA-I-length peptide
space of the libraries themselves, **34.1%** of nuORFdb v1.2's distinct 9-mers also occur in canonical
human proteins, against **1.0–2.4%** for the GENCODE Ribo-seq ORF sets — a 14–34× difference in
*latent* ambiguity, measured under our own pipeline for both. **A high-ambiguity library does not
compose into a high-ambiguity one**: Translnc, a second of IEAtlas's three sources, is at **0.6%**, and
the union nuORFdb ∪ Translnc **falls to 20.2%**. With RPFdb v2.0 unobtainable, the complete integrated
library's proportion is **not determined** — it is capped at **47.6%** but has no positive floor.
Ouspenskaia et al. searched the same nuORFdb and published a catalogue at 3%, so a high-ambiguity
library is not sufficient to produce a high-ambiguity catalogue: the exclusion step does the work.

**The consequence is observable inside the resource.** Canonical-overlapping cancer-catalogued
sequences appear in IEAtlas's **own** normal-tissue export at **22.4%**, versus **9.1%** for the
non-overlapping sequences of the same catalogue. The association is not a length artifact
(length-standardized risk ratio **2.42×**) and survives inference that respects the catalogue's
clustering (gene-clustered bootstrap, **95% CI [2.32, 2.52]**). **22,003 unique sequences — 12.6% of
the cancer catalogue — are both canonical-compatible and already present in the atlas's own normal
export.** This is consistent with, but not specific to, greater detectability or expression of
canonical-compatible sequences; it warrants normal-presentation review before such a sequence is
treated as a tumour-restricted target.

Separately and additively, the atlas does not publish the class-resolved target–decoy information
needed to reconstruct a class-specific FDR estimate. We propose retaining shared sequences with an
explicit source-compatibility annotation rather than treating them as uniquely non-canonical.

---

## Introduction

Non-canonical open reading frames (ncORFs) — in long non-coding RNAs, pseudogenes, untranslated
regions, and alternative frames of coding genes — yield peptides that are presented on HLA molecules
and recognised by T cells. Because their products are absent from the annotated proteome, they have
been proposed as an unusually attractive class of tumour antigen: potentially tumour-restricted,
shared across patients, and not subject to central tolerance. Several public catalogues now aggregate
tens or hundreds of thousands of such epitopes, and these catalogues are the practical input to
target selection.

Building such a catalogue requires answering one question for every identified peptide: **is this
peptide non-canonical?** The question is harder than it appears, and the field has known why for two
decades.

**Tandem MS identifies a sequence, not a source.** A peptide-spectrum match establishes an amino-acid
sequence. It does not establish which genomic locus produced it. Where two loci encode the same
sequence — a canonical protein and an ncORF — no spectral evidence distinguishes them. This is the
protein-inference problem (Nesvizhskii & Aebersold), and in proteogenomics it is *"exacerbated by the
mapping complexity where many identified peptides map to several loci, both novel and known"*
(Kumar et al. 2022).

**A criterion is published and recommended.** Kumar et al. (2022) state that *"most shared peptides
should be dropped if defining… novel-coding regions."* Bedran et al. (2023) implement the rule
directly: sequences perfectly matching any known protein are treated as *exonic* and excluded, with a
stringent requirement of *"at least three mismatches with any known protein sequence"* before a
peptide is called non-canonical. We call this a **sequence-exclusivity criterion**. It is published
and recommended; we do not assert it is a universal rule binding every atlas. An atlas may
legitimately retain shared observations — provided it labels them as source-ambiguous rather than as
sequence-unique ncORF products. That distinction is the subject of this paper.

A parallel and additive standard governs statistical confidence: a pooled target–decoy threshold
under-controls a minority class, so class-specific FDR estimation is required (Woo et al. 2014, who
measured 36% novel-class FDR at a 1% combined threshold; Zhang et al.; Choi & Paek 2024).

**What has not been done is to check whether the public resources satisfy these criteria.** Kumar et
al. review the problem without quantifying it, and note that no specific databases are criticised by
name. Bedran et al. did quantify — reporting residual canonical overlap of 1.4% (Erhard et al. 2020),
3% (Ouspenskaia et al. 2021), 4% (Chong et al. 2020) and 5% (Laumont et al. 2016) — but their
comparison did not include IEAtlas, the largest such atlas and the one most readily used as an
off-the-shelf source of candidate antigens.

This paper is that check.

---

## Results

### R1. Most of IEAtlas's cancer-catalogued sequences also occur in canonical proteins

**The measurement.** A catalogued peptide counts as canonical-overlapping if it is an **exact
substring** of at least one protein in a frozen reviewed human reference *R* (release and hash in
Methods). Peptides are deduplicated to **unique sequences**; the unit is a sequence, not an atlas row.

**98,193 of 174,465 unique cancer-catalogued sequences (56.3%)** are canonical-overlapping.

Three checks establish that this is not an artifact of how we measured it.

| robustness check | result |
|---|---:|
| headline, reference *R* (current reviewed human proteome) | **56.3%** (98,193 / 174,465) |
| **era-correct**: Swiss-Prot **2022_01**, the release IEAtlas searched | **56.2%** (97,999 / 174,465) |
| **length-standardized** to the catalogue's own length distribution | **56.3%** |
| **class composition**: if the atlas contained no pseudogene ORFs at all | **55.8%** |

**The era check settles the anachronism question empirically.** Sequence novelty is
reference-relative, so an overlap we score today might reflect a canonical protein that entered the
reference *after* IEAtlas was built — in which case faulting the atlas would be anachronistic. It does
not. Rebuilding the reference at Swiss-Prot 2022_01 (20,376 human proteins) moves the rate by
**0.1 percentage points**; only **231 sequences (0.24% of the overlap set)** are matches that a 2022
analyst could not have made. The overlap was almost entirely visible at build time.

#### The canonical proteins were in IEAtlas's own search database

The era check above understates the point, because it treats Swiss-Prot 2022_01 as a *reference we
chose*. It is not. It is **the canonical half of IEAtlas's own search database.** From its Methods:

> Files were searched against **both** our integrated benchmarked ncORF library **and the canonical
> human proteome** obtained from the UniProt database with Swiss-Prot protein evidence (downloaded in
> February 2022). […] **Only epitopes derived from non-coding regions were retained.**

MaxQuant, given several FASTA inputs, concatenates them into a single target database; every spectrum
is scored against canonical and non-canonical candidates **together**. So for each of the **97,999
(56.2%)** catalogued sequences that occur in a canonical human protein, an entry carrying that
identical sequence was **physically present in the database the spectrum was matched against**,
alongside the ncORF to which the peptide was ultimately attributed.

Nothing here depends on that being one search rather than two. Were the two FASTAs instead searched
separately and compared, the canonical proteome would still have been consulted by the pipeline that
produced the catalogue, and the overlap would still have been visible to it. The claim is only that
**the canonical sequences were inside the procedure, not outside it** — which is what the Methods say.

This changes what kind of claim the paper is making. The source ambiguity is not something an external
auditor reconstructed after the fact against a reference the resource never consulted. **It is a
property of the search itself**, and the search engine's own output records it: MaxQuant's
`peptides.txt` carries a `Proteins` column listing every database entry that contains each peptide.
The information the remedy in §5 asks for already existed, in an intermediate file, and was discarded
at the step described as *"only epitopes derived from non-coding regions were retained."*

The exclusion test is also not a foreign operation. Building the library, IEAtlas reports that *"FASTA
files of peptides were merged, and **peptides entirely contained within other peptides were
removed**."* That is an exact substring-containment test — the same test this paper applies — run
across the ncORF library to remove internal redundancy. It is not described as having been run
**against the canonical proteome**, though that proteome was loaded into the same search. The remedy
is that operation applied once more, in a direction the pipeline had already implemented.

The reading is falsifiable and we state the check: had the canonical half of that database been used
to **exclude** shared sequences — as CrypticProteinDB and Raja et al. describe doing, and as their
rates of 0.026% and 0.17% reflect — this rate would be near zero. It is 56.2%.

**This does not show that any peptide is canonically derived, and it does not show that anyone acted
improperly.** Sequence identity is symmetric, and the canonical entry is not "the right answer" either;
MS identifies the sequence, never the locus. What it shows is that the evidence needed to mark these
sequences *source-ambiguous* was inside the pipeline that produced the catalogue.

**Under one pipeline**, applied identically to every catalogue we could reprocess — same reference,
same exact-substring criterion, same unique-sequence unit:

| catalogue | exclusion rule in its Methods? | canonical-sequence overlap |
|---|---|---:|
| CrypticProteinDB | **yes** — *"BLASTP… eliminate all proteins with alignment to canonical proteins"* | **1 / 3,810 = 0.026%** |
| Raja et al. (ovarian) | **yes** — *"peptides mapping to 'protein_coding'… were excluded"* | **5 / 2,979 = 0.17%** |
| **IEAtlas** | **not described** — *"only epitopes derived from non-coding regions were retained"* | **98,193 / 174,465 = 56.3%** |

For context, Bedran et al. 2023 report residual canonical overlap of 1.4% (Erhard), 3% (Ouspenskaia),
4% (Chong) and 5% (Laumont). **We do not compute a fold-change against those values.** They were
produced with a different reference, normalization, deduplication and peptide unit, and a ratio across
pipelines is arithmetic rather than measurement. They are cited as published context. The controlled
comparison is the table above.

#### It is not the ORF-class composition

IEAtlas's ncORF library explicitly includes pseudogenes, the class with the highest canonical
compatibility. **A peptide may carry more than one ORF label** — 1,801 sequences (1.0%) map to more
than one gene, and 546 carry both a pseudogene and a non-pseudogene ORF — so pseudogene and
non-pseudogene sets are *not* complements. Reported as mutually exclusive strata, which do partition:

| stratum (mutually exclusive) | *n* | canonical overlap |
|---|---:|---:|
| pseudogene-only | 15,777 | **59.4%** |
| non-pseudogene-only | 158,142 | **55.8%** |
| both labels (source-ambiguous within the atlas) | 546 | **93.0%** |
| **total** | **174,465** | **56.3%** |

If the atlas contained no pseudogene ORFs whatever, the rate would be **55.8%**. Class composition
moves the headline by half a percentage point. Separately, Raja et al. report 98 pseudogene-ORF
peptides and **none** overlaps a canonical protein, because their exclusion rule removed those that
did — same ORF class, different rule.

That 1.0% of sequences map to several genes *within the atlas's own annotations* is itself a small,
direct instance of the ambiguity this paper is about.

### R2. Latent canonical ambiguity differs enormously between libraries — and does not compose

Applying the same measurement to the **search space** rather than the output — the distinct 9-mers of
each ncORF library, and how many also occur in reviewed canonical human proteins (**full libraries, no
sampling**):

| ncORF library | ORFs | distinct 9-mers | also canonical |
|---|---:|---:|---:|
| **nuORFdb v1.2** — integrated by IEAtlas | 229,251 | 8,448,245 | **34.1%** |
| GENCODE Ribo-seq ORFs (phase 1) | 7,264 | 245,094 | **2.4%** |
| GENCODE Ribo-seq ORFs (phase 2) | 28,359 | 502,528 | **1.0%** |

Both libraries were measured by us, under one pipeline, so the contrast *is* a controlled one:
**ncORF libraries differ by 14–34× in latent canonical ambiguity.** Whole-ORF containment is low
throughout (0.2–0.8%), so this is **extensive partial sharing**, not whole ncORFs nested inside
canonical proteins.

**Independent corroboration.** These figures were obtained twice, by separate implementations over
different k-mer windows: an 8–11mer candidate-universe enumeration gives nuORFdb **34.4%** and GENCODE
Ribo-seq (phase 1) **2.5%**; the 9-mer enumeration above gives **34.1%** and **2.4%**.

#### The integrated library does not inherit nuORFdb's ambiguity

An earlier draft claimed 34.1% was a **lower bound** on IEAtlas's complete integrated library
(nuORFdb + RPFdb + Translnc). **That was false and is withdrawn.** Adding a library *B* to nuORFdb *A*
changes the rate to |(*A* ∪ *B*) ∩ *C*| / |*A* ∪ *B*|, which is **not monotone**: if *B* contributes
mostly non-canonical *k*-mers, the combined proportion **falls**.

It falls. We obtained Translnc — a second of IEAtlas's three sources, in the version it cites — and
measured the union directly:

| library | distinct 9-mers | also canonical | rate |
|---|---:|---:|---:|
| nuORFdb v1.2 | 8,448,245 | 2,884,119 | **34.1%** |
| Translnc | 6,164,584 | 39,382 | **0.6%** |
| **nuORFdb ∪ Translnc** | 14,364,357 | 2,907,391 | **20.2%** |

Translnc is almost free of latent canonical ambiguity, and the two libraries are nearly disjoint (only
1.7% of the union's 9-mers occur in both). Adding it therefore enlarges the denominator far faster than
the numerator, and the rate **drops by 13.9 pp, from 34.1% to 20.2%**. The withdrawn "lower bound" was
not merely unproven; it is **false in exactly the direction the objection predicted, on IEAtlas's own
second source**. (Re-scored against Swiss-Prot 2022_01: unchanged, 20.2%.)

The gap between the two is not arbitrary, and it generalises. Translnc catalogues peptides from
**lncRNA loci**, which are by construction distinct from protein-coding genes and therefore share
little sequence with them. nuORFdb catalogues ORFs of **coding genes** — uORFs, dORFs, out-of-frame
and in-frame alternative ORFs — whose reading frames overlap or abut the very proteins in the canonical
reference. **A library's latent ambiguity is largely determined by whether its ORFs sit inside coding
genes**, which is a property its builders know and could report at zero cost. This is the concrete
content of recommendation (c) in §5.

**What is and is not now bounded.** RPFdb v2.0 remains genuinely unavailable — the live site serves
only v3.0, and it distributes RibORF genomic *coordinates*, not amino-acid sequences, so
back-translation would produce *our* library rather than IEAtlas's. Carrying its contribution as the
single unknown *m* (novel 9-mers it adds, of which *x* are canonical), the combined rate is
(*h* + *x*) / (*u* + *m*). Maximising over both, the three-source library's 9-mer ambiguity **cannot
exceed 47.6%**, distribution-free and whatever RPFdb contains. That cap is a **hard ceiling, not an
estimate**: it is attained only in the corner where every canonical 9-mer not already in the union is
contributed by RPFdb *and* RPFdb contributes nothing else. There is **no corresponding floor**: as *m*
grows, *h*/(*u* + *m*) → 0. **The combined rate is not determined by our measurement**, and no
measurement of a subset of the libraries could determine it.

**The library is not sufficient on its own.** The evidence that the *exclusion step*, not the library,
is what governs a catalogue's overlap rate is our own same-pipeline reprocessing: CrypticProteinDB and
Raja et al., which describe explicit exclusion rules, sit at **0.026%** and **0.17%** (§1). Ouspenskaia
et al. add the one control neither of those provides — they searched the **same nuORFdb** and their
published catalogue is nonetheless low (3%, Bedran et al.). We attach **no arithmetic** to that 3%: it
is a published cross-pipeline figure, and the unit caveat above applies to it as much as to any other.
The qualitative point is all we need and all we claim — a high-ambiguity library does not force a
high-ambiguity catalogue.

#### A detection effect, tested — and what we refuse to infer

**We do not compare the catalogue's rate to the library's as levels.** An earlier draft opened this
subsection by observing that the catalogued 56.3% "exceeds" nuORFdb's 34.1%, and read the difference as
an excess arising during detection. **That inference is invalid and is withdrawn.** The two numbers are
different objects: 56.3% is over distinct catalogued *peptides* at native lengths, after search, FDR
and deduplication; 34.1% is over distinct *9-mers* of a candidate search space in which nothing has
been detected. Different units, different denominators, different lengths — the same cross-unit error
as the fold-change against other catalogues that we withdraw in §1, and it is not rescued by the union
figure either (56.3% must not be read against 20.2%, or against the 47.6% cap).

The detection hypothesis is nonetheless real and testable *without* that comparison — a peptide of an
abundant, ubiquitous protein should be over-detected in an immunopeptidome relative to its share of the
search space. Earlier drafts asserted this without testing it. Both tests below are internal to a single
space or are ratios of ratios, so neither requires a cross-unit subtraction.

**Prediction 1 — breadth of detection.** A peptide of an abundant, ubiquitous protein should be
detected across more of IEAtlas's 15 cancer types. Canonical-overlapping sequences are seen in a mean
of **1.62** cancer types versus **1.33** for non-overlapping ones, and **28.3%** appear in ≥2 types
versus **16.8%**. Because short peptides both match canonical proteins more readily and recur more
often, we stratified by length: the effect holds in **18 of 18** length strata (8–25 aa). It is not a
length artifact.

**Prediction 2 — enrichment for the abundant housekeeping class.** Ribosomal proteins are the textbook
abundant, ubiquitous class. In the catalogue, canonical-overlapping sequences are **2.51×** enriched
for ribosomal-gene ORFs (2.89% vs 1.15%). This could be pure library composition, so we measured the
same enrichment **in the library**, over (ORF, 9-mer) candidate pairs, where nothing has been detected
yet: there, ribosomal ORFs are slightly **depleted** among canonical-overlapping candidates (**0.91×**).
The catalogue's enrichment is therefore a **2.77× excess over its own search space** — it was not
inherited from the library, and arose during detection.

Note what Prediction 2 does and does not compare. It is a **ratio of ratios** — an *enrichment* measured
within the catalogue, set against the *same enrichment* measured within the library — and a ratio of
ratios is dimensionless, so it survives the unit mismatch that makes a comparison of the two raw rates
meaningless. It says the association between being canonical-overlapping and being ribosomal is
stronger in the catalogue than in the search space it was drawn from. It does **not** license any claim
that the catalogue's overlap rate sits some number of percentage points above the library's — that
quantity is not defined.

**Prediction 3 — the direct test.** Predictions 1 and 2 use *breadth of detection* and *ribosomal
membership* as **proxies** for abundance. We replace them with a measurement: **PaxDb v6.1** (human,
whole-organism integrated, ppm), joined to the canonical proteins the catalogued sequences actually
match (**91.4%** of the reviewed proteome joins; 96,210 of the 98,193 overlapping sequences carry an
abundance).

*Which canonical proteins do the catalogued sequences hit?* **Abundant ones.** At the protein level —
where there is no peptide-clustering problem at all, because the unit *is* the protein — canonical
proteins hit by a catalogued sequence have a median abundance of **0.872 ppm** against **0.086 ppm**
for those never hit, a **10.14×** difference (AUC **0.679**).

*Does abundance predict how broadly a peptide was detected?* **Yes — but weakly, and we state the
weakness rather than the headline.** Binning canonical-overlapping sequences by the abundance of their
matched canonical protein, mean detection breadth rises across every quintile after length
standardization, from **1.40** cancer types (Q1) to **1.78** (Q5) — a Q5−Q1 gap of **0.377**, with a
gene-clustered bootstrap **95% CI [0.335, 0.42]** that excludes zero under **both** cluster definitions
(matched canonical gene, and IEAtlas source gene), and which reproduces on the previous PaxDb release
(**0.38**). **The crude, unstandardized trend saturates** — it climbs through the low quintiles and
flattens across the top — and is monotone *only* after length standardization; the crude gap is
**0.23**. The direction and the interval are solid. **A strong dose–response is not**, and we do not
claim one.

Three controls, each of which could have killed it:

| control | the attack it answers | result |
|---|---|---|
| **peptide length** | short peptides both match canonical proteins more readily *and* recur more | trend holds in **18 / 18** length strata |
| **protein length** | longer proteins are hit more often *by chance*, so "hit" may just mean "long" | hit proteins **are** longer (median 490 vs 360 aa) — but the abundance effect holds in **10 / 10** protein-length deciles |
| **placebo** | does the machinery invent trends? | breaking the peptide→protein link **collapses** the gap to **0.0** (0 of 200 draws reach the observed 0.377) |

**What this licenses, and what it does not.** The abundance explanation is now **measured, not
proxied**: canonical-overlapping sequences preferentially match abundant canonical proteins, and that
abundance predicts detection breadth, surviving peptide length, protein length and a placebo. But the
per-sequence association is **weak** (Spearman ρ = **0.061**): abundance is *one* contributor to which
sequences get detected, **not** the whole explanation, and we do not claim otherwise. Nor does any of
this speak to provenance — it is evidence about **what gets detected**, and MS still identifies the
sequence, never the locus.

### R3. The consequence is observable inside the resource

An audit that stops at *"56.3% of these sequences are source-ambiguous"* invites the only question
that matters: **so what?**

If a catalogued "cancer epitope" is in fact a peptide of an abundant canonical protein, it should also
be presented on **normal tissue**, because that protein is expressed there too. IEAtlas publishes its
own normal-tissue epitope export (94,375 unique peptides), so the prediction is testable **inside the
resource**, using the non-overlapping sequences of the same catalogue as a **within-resource
comparator**.

| IEAtlas cancer-catalogued sequences | also in IEAtlas's **own** normal-tissue export |
|---|---:|
| **canonical-overlapping** (98,193) | **22,003 = 22.4%** |
| non-overlapping — within-resource comparator (76,272) | 6,976 = **9.1%** |

**The inference, done correctly.** These 174,465 sequences are *not* independent Bernoulli
observations — they are clustered by source gene, gene family, dataset, tissue, donor, HLA allele and
pipeline — and peptide length confounds both arms. Accordingly:

- **Per length, unpooled**, the risk ratio ranges **1.55–3.08** and the effect holds at **all 18**
  peptide lengths (8–25 aa).
- **Length-standardized** (direct standardization to the catalogue's own length distribution) the risk
  ratio is **2.42×** — essentially the crude 2.45×, so length is not driving it.
- **Gene-clustered bootstrap** (resampling 22,765 source-gene clusters with replacement, *B* = 2,000):
  **95% CI [2.32, 2.52]**.
- **Within every ORF-class stratum**: 2.39× (pseudogene-only), 2.40× (non-pseudogene-only).

Earlier drafts attached a two-proportion *z* = 74 to this contrast. **That statistic was invalid** — it
treated a heavily structured catalogue as 174,465 independent experiments — and it is withdrawn. Being
directionally right does not rescue it. The clustered interval above replaces it.

**A subset requires no inference at all.** **22,003 unique sequences — 12.6% of the cancer catalogue —
are both canonical-compatible and already present in the atlas's own normal-tissue export.** No
external reference is needed to see this; it is internal to the resource.

**What this means, bounded.** This is **consistent with, but not specific to**, greater detectability
or expression of canonical-compatible sequences. It does not show that any individual sequence is
canonically derived: presence in a normal-tissue export is evidence about **presentation**, not about
**source**. Its practical consequence is that such a sequence **warrants normal-presentation review
before it is treated as a tumour-restricted target** — clinical risk depends additionally on allele
matching, tissue context, abundance and TCR avidity, none of which we assess.

### R4. The additive statistical problem, and the remedy the field already demonstrated

**Not novel** (Woo et al. 2014; Choi & Paek 2024). Stated because it is *additive* to R1–R3, and
because IEAtlas reports nothing that would allow it to be assessed.

From a reported pooled FDR *q* and class fraction *f*, the class-specific FDR is only
**set-identified**:

    Θ_N(q, f) = [ max(0, (q − (1 − f)) / f),  min(1, q / f) ]

**From *q* and *f* alone, this identified set is sharp.** IEAtlas reports *q* = 0.05 and 245,870
non-canonical epitopes, but **no canonical count**. *f* is therefore unknown and the interval is
**unconstrained**.

Reporting the per-class accepted target and decoy counts (*T_N*, *D_N*), the threshold, the unit and
the convention would make the selected class-specific target–decoy estimate **reconstructible**. This
is not the only information that could tighten the interval — calibrated class-specific posterior error
probabilities or entrapment measurements could also do so — and *D_N* does not identify the true
class-specific false-discovery proportion. It is simply the cheapest sufficient object, and one the
pipeline already computes.

**The field has demonstrated the remedy, and its cost.** Ouspenskaia et al. searched a combined
annotated-ORF/nuORF database and reported that a 1% global FDR gave *"4.6% overall, and as high as 14%
for 3′ dORFs"* among nuORF peptides; group-based filtering *"removed 24% of nuORF peptides overall, and
up to 76% of peptides assigned to 3′ overlap dORFs."* **The 4.6% is a pre-correction diagnostic, not
the error rate of their published catalogue.** They are a positive exemplar.

### R5. Two distinct problems; one reporting remedy; and what it costs

The two defects are **complementary, not the same**, and merging them is a category error:

- **Source ambiguity (R1–R3).** The peptide is *correctly identified as a sequence*; its **source
  locus** is unresolved. This is not an FDR problem — no identification is wrong.
- **Class-specific FDR under-control (R4).** A pooled threshold under-controls the minority class, so
  some ncORF **identifications are wrong** — the spectrum was not that peptide.

An earlier draft argued that FDR could not *explain* the 56.3%, on the grounds that a
composition-matched shuffle places chance canonical overlap near 0.1%. **That argument was wrong and is
withdrawn.** A false target PSM is not an arbitrary shuffled string; it is an accepted, incorrect
candidate drawn from the actual search database, and its class composition is not described by a
shuffle null. The correct statement needs no such argument: source ambiguity is present *even when the
sequence is correctly identified*. FDR concerns whether the spectrum was assigned to the right
sequence; canonical overlap concerns whether that correctly-identified sequence determines a source.
They are different objects.

**A minimal standard addresses both.**

**(a) Per peptide — an exclusivity flag and all compatible source loci.** State whether the sequence is
unique to the nominated ncORF within the searched space, and if not, list every compatible source. This
preserves the peptide *and* the ambiguity, rather than discarding either, and it is strictly more
informative than the exclusion rule it generalises.

*And it is cheap.* A remedy that demanded a long list of loci per peptide would be glib, so we measured
the list. Across the 97,999 catalogued sequences that match the canonical half of IEAtlas's own search
database, the **median number of compatible canonical genes is 1**; **93.1% are compatible with exactly
one** canonical gene and **98.0% with at most two** (maximum 22). The label is, for the large majority
of ambiguous sequences, **a single gene symbol** — one column, not a redesign. And as §1 showed, a
concatenated search already computes it: MaxQuant's `peptides.txt` lists every database entry
containing each peptide. The remedy asks resources to **publish a column their own pipelines already
produce**.

**(b) Per class — a class-decoy ledger.** Accepted target and decoy counts per class, the class
definitions, the thresholding stage, and the formula used — enough to reconstruct the class-specific
estimate. Ouspenskaia et al. demonstrate it is achievable.

**(c) Per library — publish the latent canonical ambiguity.** One number (R2), cheap to compute, and it
tells every downstream user how much an exclusion rule matters for the library they are about to
search. **No ncORF library currently publishes it.**

**The cost, stated plainly.** Under a sequence-exclusive definition, **98,193 of 174,465 unique
cancer-catalogued sequences (56.3%)** would be ineligible for designation as uniquely non-canonical.
That does **not** require deleting them from the atlas — the remedy in (a) is to retain them and label
them source-ambiguous with their compatible loci. For scale, class-specific FDR control cost
Ouspenskaia et al. 24% of their nuORF peptides and up to 76% of one ORF class. **Applying an
established criterion to an ncORF catalogue is expected to be expensive. That is not an argument
against applying it** — it is an argument for applying it before, rather than after, a candidate enters
a clinical pipeline.

---

## Discussion

**This is not an accusation, and should not be read as one.** Every resource examined describes its own
procedures accurately and in public; that is the only reason this audit was possible at all. Two of the
catalogues measured here already apply an exclusion rule. Ouspenskaia et al. already solved the
statistical half and published the cost of doing so. **The field knows how to do this.** What has not
happened is a check that the resources in use actually do it — and the reason is prosaic: nobody had
looked.

**Where the ambiguity comes from, and where the remedy belongs.** A third of nuORFdb's peptide space is
canonical by sequence. We cannot say what the corresponding figure is for IEAtlas's full integrated
library, and we do not claim the library *quantitatively explains* 56.3%. But a library that publishes
its own latent ambiguity lets every downstream group know in advance how much an exclusion rule will
cost them, and the detection effect measured in R2 shows that what a catalogue reports is not simply a
readout of what its library contains. **That number does not currently exist for any ncORF library, and
it takes minutes to compute.**

**The stakes are not abstract.** A sequence catalogued as a cancer epitope, which is compatible with an
abundant canonical protein, and which the same atlas has already observed on normal tissue, is not a
promising tumour-restricted target without further review. **22,003 unique sequences in IEAtlas meet
all three conditions.** We are not claiming that any of them *is* canonically derived — MS cannot say —
but a target-selection pipeline that cannot distinguish them is selecting under an ambiguity it has not
been told about.

**Limits, stated plainly.**

- The overlap is **reference-relative** — it is *N*(*R*). For a fixed query set, overlap is monotone
  under **nested expansion of the same *R***; 56.3% is a lower bound with respect to supersets of *this*
  reference, **not** with respect to every possible reference definition. A narrower or differently
  defined reference can lower it. The era-correct check (56.2% against Swiss-Prot 2022_01) is the one
  that matters for judging the resource at build time.
- **No individual peptide's provenance is resolved.** Sequence identity is symmetric.
- We hold nuORFdb but **not** RPFdb or Translnc. The latent ambiguity of IEAtlas's complete integrated
  library is **unknown**, and 34.1% is **not** a lower bound on it.
- The 1.4–5% values for four other catalogues are **published figures from a different pipeline**, cited
  as context. We report **no fold-change** against them.
- Class labels for non-pseudogene classes are **source-supplied and uncorroborated**.
- The pseudogene→parent homology analysis (Supplement) now uses a **curated, versioned** annotation
  (NCBI Gene `gene_group`) rather than symbol-stripping, and the parent hit survives a
  **family-respecting** null even when the decoys are the parent's own **strong paralogs** (52.3%
  observed vs 16.6%, *p* < 1e-4) — the two objections that put it there. It stays in the Supplement
  anyway: **133 testable peptides**, and the curated mapping is not
  independent of the symbol it replaces, because HGNC names a pseudogene *after* its parent. It
  explains *why* part of the pseudogene class is source-ambiguous; it does not measure the headline.
- The detection-bias result (R2) is now a **direct measurement** (PaxDb v6.1), not a proxy, and it
  survives peptide length, protein length and a placebo. But the **per-sequence association is weak**
  (Spearman ρ = 0.061). Abundance is **one** contributor to what gets detected, not the explanation.

**What we would like to be wrong about.** If IEAtlas's pipeline resolves source attribution in a way its
Methods do not describe — for instance through a protein-inference step that assigns shared sequences to
the canonical protein — then the 56.3% has an innocent explanation, and this paper is a correction to a
misreading of the Methods rather than a finding about the resource. The search outputs (`peptides.txt`:
`Proteins`, `Leading razor protein`, `Unique (Proteins)`) would settle it in an afternoon. **We invite
the authors to publish them.**

---

## Methods

### Reference

Canonical human proteome: UniProt/Swiss-Prot reviewed human sequences (*R*), release and SHA-256 hash
recorded in `data/SOURCES.md` and printed by `reference_provenance.py`. **All overlap statistics are
`N(R)` and are reference-relative.** Because IEAtlas searched a February-2022 Swiss-Prot release, every
headline is additionally re-derived against **Swiss-Prot 2022_01** (566,996 entries; human subset by
`OX=9606`, 20,376 proteins), which is the estimand relevant to whether the overlap was detectable when
the resource was built. Both are reported.

### The overlap measurement

A catalogued peptide sequence counts as canonical-overlapping if it is an **exact substring** of at
least one protein in *R*. Sequences are compared as unmodified amino-acid strings, uppercased, with
inline PTM annotations stripped. Peptides are **deduplicated to unique sequences** before rates are
computed; a rate is `unique canonical-overlapping / unique scored`. The unit is a peptide **sequence**,
never an atlas row. The identical procedure is applied to every catalogue we reprocess, so that
comparison is internally consistent. Rates published by other groups under other pipelines are cited as
context and **never** combined into a ratio with ours.

**Length standardization.** Because exact-substring probability is strongly length-dependent, rates are
also reported after **direct standardization** to a common length distribution (IEAtlas's own), with
per-length rates given unpooled. A standardized rate is not reported for a catalogue with fewer than 20
overlap events, where it would be falsely precise.

### The library measurement

For each ncORF library, all distinct **9-mers** (a representative HLA-I ligand length) are enumerated
across the **full library — no sampling** — and intersected with the set of all distinct 9-mers in *R*.
Sampling is invalid here and the bias is large and upward: sampling 4,000 of nuORFdb's 229,251 ORFs
gives 43.6%, 20,000 gives 40.7%, and the full library gives 34.1%, because a small sample contains
fewer distinct ncORF-specific k-mers and so over-weights the canonical-shared ones. Whole-ORF
containment *is* sampled (4,000 ORFs, seed 0); that statistic is an O(*n·m*) substring scan and is
unbiased under random sampling, unlike a k-mer union.

**Union caveat, and the union itself.** For libraries *A*, *B* and canonical set *C*, the combined rate
is |(*A* ∪ *B*) ∩ *C*| / |*A* ∪ *B*|, which is **not** monotone in the addition of *B*. No bound on the
combined library's ambiguity is claimed from nuORFdb alone. The union of the two sources we hold is
measured directly (nuORFdb ∪ Translnc = 20.2%, *below* nuORFdb's 34.1%). RPFdb v2.0 is not measured
and is **not approximated**: the live site serves only v3.0, and it distributes RibORF genomic
*coordinates* rather than amino-acid sequences, so back-translation would require many free parameters
and would yield *our* library rather than IEAtlas's. Its contribution is carried as a single free
unknown *m*, giving a distribution-free cap of **47.6%** on the three-source library and **no positive
floor**.

**A comparison we do not make.** The library rate (distinct 9-mers of an undetected candidate space)
and the catalogue rate (distinct peptides at native lengths, after search, FDR and deduplication) are
**different objects with different units and denominators**. We report no difference, ratio or excess
between them. The detection test is stated only as a within-catalogue contrast (breadth) and a **ratio
of ratios** (ribosomal enrichment in the catalogue against the same enrichment in the library), both of
which are immune to the unit mismatch.

### The detection-bias test

IEAtlas records the cancer type of each observation (15 types). Breadth of detection (number of distinct
cancer types per unique sequence) is compared between canonical-overlapping and non-overlapping
sequences, **stratified by peptide length**. The ribosomal-enrichment test compares the share of
ribosomal-gene ORFs (`RPL*`/`RPS*`/`MRPL*`/`MRPS*`) between the two groups *in the catalogue*, against
the same share computed over **(ORF, 9-mer) candidate pairs in nuORFdb**, which is the composition of
the search space prior to any detection. The reported effect is the **excess** of the former over the
latter.

### The normal-tissue consequence

IEAtlas's cancer and normal exports are compared as unique sequences. The non-overlapping sequences of
the same catalogue are a **within-resource comparator**, not a control: they do not control abundance,
detectability, HLA coverage or study composition. Inference uses a **gene-clustered bootstrap** —
source-gene clusters resampled with replacement, *B* = 2,000, seed 20260713 — of the
**length-standardized** risk ratio, reported as a percentile interval. A two-proportion *z*-test is
**not** valid here (the observations are clustered) and is not used.

### Class strata

ORF-class strata are **mutually exclusive** (`pseudogene-only`, `non-pseudogene-only`, `both`), because
a peptide may carry several ORF labels; a pseudogene / non-pseudogene split is *not* a partition and
double-counts the 546 sequences carrying both.

### Class-specific FDR identifiability

Derivation of Θ_N(*q*, *f*) and its sharpness **given *q* and *f*** in the Supplement. The result is not
new (Woo et al. 2014); it is stated here because IEAtlas reports *q* but not *f*.

### Reproducibility

Every headline number in this manuscript is regenerated from the committed artifacts by
`manuscript/verify_manuscript.py`, **which fails the build on drift** and additionally fails if the
paper drops its required prior-art citations or asserts any of the phrasings retracted during review.
The analysis code is guarded by `src/darkproteome/scoring_conformance.py`. All primary Methods
quotations were verified in the fetched full text (EuropePMC / NCBI E-utilities); the fetch scripts and
document hashes are in the repository.

---

## What this paper does not claim

- That the underlying principles are ours. The exclusion criterion is Bedran et al. 2023; the
  shared-peptide remedy is Kumar et al. 2022; the inference problem is Nesvizhskii & Aebersold; the
  class-FDR problem is Woo et al. 2014. **Our contribution is empirical.**
- That the biology is fake, or that any specific ncORF antigen is not real.
- That any resource "manufactures", "re-labels", or "discards" anything.
- That the canonical source is the correct one for any peptide. **Sequence identity is symmetric; MS
  identifies the sequence, never the locus.**
- That any resource acted improperly. All describe their own procedures accurately. We say that
  IEAtlas's **published Methods do not describe** a peptide-level canonical-exclusion rule — not that
  its pipeline applies none.
- That the sequence-exclusivity criterion is a universal field-wide rule. It is published and
  recommended. An atlas may retain shared sequences if it **labels** them source-ambiguous.
- That presence in a normal-tissue export establishes a peptide's *source*, or that the 22,003
  sequences are established on-target/off-tumour risks. They **warrant normal-presentation review**.
- That the 56.3% is reference-independent, or a lower bound against *any* reference. It is `N(R)`,
  monotone only under nested expansion of *R*.
- That IEAtlas is "11–40×" any other catalogue. **We report no fold-change against externally
  published rates.**
- That we have measured IEAtlas's whole library. We hold nuORFdb and Translnc but **not RPFdb v2.0**,
  which is unobtainable. The combined proportion is **not determined** — capped at 47.6%, with no
  positive floor, and certainly **not bounded below by 34.1%**.
- That the catalogue's overlap rate can be compared to a library's **as levels**. They are different
  units (peptides vs 9-mers) over different denominators (detected output vs candidate space). We
  report no difference, ratio or excess between them — **including against the 20.2% union and the
  47.6% cap.**
- That detection bias is *proven*, or that abundance is the *explanation*. It is now **directly
  measured** (PaxDb v6.1) rather than proxied, and survives length, protein length and a placebo — but
  the per-sequence association is **weak** (ρ = 0.061) and the crude trend **saturates**. Abundance is
  **one** contributor to what gets detected.
- That the pseudogene→parent analysis is headline-grade. It now uses a **curated, versioned** parent
  annotation (NCBI Gene `gene_group`) and survives a family-respecting null, but it rests on **133
  testable peptides**, and the curated mapping is **not independent** of the gene symbol it replaces.
  It stays in the Supplement as a **mechanistic vignette**, not a measurement.
