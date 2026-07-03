# Canonical self in cryptic cancer-epitope catalogues and the class-decoy ledger needed to verify non-canonical antigen claims

**Rom Jan**¹  
¹ Independent Researcher. ORCID: 0009-0004-1015-1627. Correspondence: romjacobjan2005@gmail.com.

*All analyses use public data and are fully scripted and reproducible.*

---

## Abstract
Non-canonical open reading frames (ncORFs) are increasingly reported as sources of tumor-specific HLA antigens, but their clinical reuse depends on whether each published claim stays checkable after publication. We audited **306,844 published non-canonical tumor-antigen claims** (180,650 unique peptides) from two flagship end-to-end cohorts and the two largest public cryptic-peptide atlases, scoring only evidence in the reported record. The largest measured defect is sequence-level: **56.3% of catalogued cryptic "cancer" epitopes in IEAtlas are exact substrings of the canonical human proteome** — canonical *self* by sequence, which cannot serve as evidence for the non-canonical, tumor-specific antigen claim the catalogue nominates it for, whatever ORF it is attributed to. The effect is class-resolved, not field-wide: primary-cohort altORF/lncRNA peptides are clean (≈0.2%), pseudogene-derived claims inherit canonical parent sequences (37.1%), and a clean atlas (CrypticProteinDB, 0.0%) shows the rate is avoidable; of the genuinely non-canonical remainder, 9.1% appear on normal tissue. Separately, of all 306,844 claims **none** carries a machine-readable per-claim evidence package — a **reusability** gap that even known-real MAGE/SSX antigens share, and that downstream pipelines inherit when they treat a catalogue row as a validated target — locating the problem in reporting, not biology. The gap is structural: from reported summary statistics the non-canonical class false-discovery rate is only *set*-identified. The single field that closes it — the per-class accepted decoy count **`D_N`**, a "class-decoy ledger" — makes future claims re-verifiable.

## Significance
A clinically important class of cancer antigens is being nominated faster than it can be checked from the published record. We show that a majority of catalogued cryptic "cancer" epitopes in the largest public atlas are canonical by sequence, that the published record rarely carries the per-claim evidence to re-verify a claim, and that the gap has a provable statistical core — then give the minimal reporting standard (a class-decoy ledger) that closes it.

---

## Introduction

The human genome is translated far more broadly than its ≈20,000 canonical protein-coding genes.
Ribosome profiling and mass-spectrometry immunopeptidomics have shown that non-canonical open
reading frames (ncORFs) — in long non-coding RNAs, pseudogenes, 5′/3′ untranslated regions, and
alternative reading frames of coding genes — give rise to short peptides that are presented on
human leukocyte antigen (HLA) class I molecules and recognized by T cells [13].
Because these "dark-proteome" or "cryptic" peptides can be tumor-restricted, they have become a
prominent candidate source of tumor-specific antigens for cancer vaccines and T-cell therapies,
and the published catalogue of such claims now numbers in the hundreds of thousands [14,15].

A genuine tumor antigen must clear a four-step chain: its source ORF must be genuinely translated;
its peptide must be genuinely HLA-presented; the peptide must be tumor-specific (absent from normal
presentation); and T cells must genuinely respond. As the field moves toward the clinic, a
prerequisite — and, to our knowledge, unmeasured — question is whether the claims already in the
literature carry, *per claim*, the evidence needed for an independent party to re-verify them along
this chain. This is a question about the reusability of the published record, not about the
underlying biology.

One link in that chain rests on a statistical control that is known to be fragile for this class of
peptide. Non-canonical peptides are identified by searching tandem mass spectra against vastly
inflated databases (multi-frame translations of non-coding transcripts, pseudogenes, ncORF
catalogues), with the false-discovery rate (FDR) controlled by target-decoy competition [6]. It is well
established that a single *global* FDR pooled over canonical and non-canonical candidates
underestimates the true FDR of the non-canonical subclass: because that subclass is typically a
small fraction (<1–5%) of accepted identifications, a global threshold is excessively permissive for
it [4,5,10,11]. This has been quantified for over a decade — a nominal 1% combined FDR was shown to
correspond to a 36% FDR for the novel-peptide class while the canonical class sat at 0.03% [19] — and
recent immunopeptidomics pipelines adopt separate class-specific FDR for exactly this reason [20]. The
field's recommended remedies — class-specific (group) FDR and multistage search — concern *how to
compute* the subclass FDR when one holds the raw spectra and the per-class decoy counts [5,10,11]. What has not been examined is the
complementary, reproducibility-facing question: given only the statistics that papers actually
*report*, what can a reader independently conclude about a claim's class-specific FDR — and is that
enough to re-verify the claim at all?

Here we make three contributions. **First**, we conduct what is, to our knowledge, the first
systematic, leakage-controlled, public-data audit of the published non-canonical tumor-antigen
record — 306,844 claims (180,650 unique peptides) from two flagship end-to-end cohorts and the two
largest public cryptic-peptide atlases — and report two findings. (i) A direct measurement: **56.3% of
catalogued cryptic "cancer" epitopes in IEAtlas are exact substrings of the canonical proteome** (≈800× a
composition-matched null), strongly class-dependent (pseudogene-ORF 37.1%, altORF 0.2%; a second atlas
that applies a canonical filter is clean) and, to our knowledge, unquantified. (ii) A re-verifiability result: scored — exactly as reported — on
the four evidence axes, none of the claims carries reusable per-claim evidence on all four at once and
only seven are auditable on all four; critically, known-real canonical cancer-testis antigens
(MAGE/SSX) fail the identical bar, establishing this as a reporting and reusability gap rather than a
statement about the biology. **Second**, we show the gap is structural, not incidental, by
characterizing what the class-specific FDR can be from the reported statistics: it is only
*set*-identified, to the sharp interval [0, min(1, α/(T_N/T))]; the common upper "budget" bound is
provably as tight as summary statistics allow; the quantity that governs the class decoy
allocation at a threshold is a score-weighted *effective ρ* that can depart from the global
candidate-space ratio (in a regime we characterize);
and, propagated to the claim level, a calibrated per-claim truth probability is **not identifiable**
— the honest output is a partial order over claims. **Third**, we give the minimal fix: a per-claim
reporting standard — chiefly, report the per-class decoy count `D_N` on the claim's unit — that
provably collapses the identifiability gap, together with a proposal for a neutral, leakage-
controlled benchmark for the field.

Throughout, our claims concern the reporting and the statistics, never the biology: we do not assert
that any specific antigen is false, and the failure of the canonical controls under the same audit
is our central guard against that misreading. All inputs are public and all analyses are scripted
and openly reproducible (see Data and code availability).

## Results I — The audit: sequence non-novelty, specificity, and re-verifiability

### A four-axis, as-reported audit of 306,844 published claims, stratified by resource type
We assembled 306,844 published non-canonical tumor-antigen claims (180,650 unique peptides) from
two flagship end-to-end immunopeptidomics cohorts — ovarian [1] and hepatocellular
carcinoma [2] — together with the two largest public cryptic-peptide resources,
IEAtlas (293,222 cancer epitopes) and CrypticProteinDB (9,101). These are different kinds of object,
and we keep the strata separate throughout: a cohort paper *could in principle* report a complete
per-claim evidence package, whereas an atlas/catalogue record was never built to carry one. Each claim
was scored, exactly as the source reports it, on four independent evidence axes: (i) *source-ORF
translation* (is the ORF plausibly translated — Ribo-seq periodicity or the protein-existence consensus
bar, applied to the ORF and not to the single ligand); (ii) *HLA presentation* (an eluted 8–12-mer with
an assigned allele); (iii) *tumor specificity* (checked against a broad normal panel); and (iv)
*immunogenicity* (an autologous or HLA-matched human T-cell response). A claim is a *strict survivor*
only if it strict-passes all four; anything the source leaves underspecified is labelled *unverifiable*
rather than failed. The auditor never recomputes FDR from spectra — it scores the evidence papers
actually publish (Methods). We report two distinct findings: a direct, positive *measurement* of
sequence non-novelty, and a measurement of *re-verifiability* — whether the reported record carries a
reusable per-claim evidence package.

### More than half of catalogued cryptic "cancer" epitopes are canonical by sequence
The first finding is a positive, *measured* rate, and — to our knowledge — the first time it has been
quantified for these widely reused resources against the field's own sequence-novelty bar.
Cross-referencing IEAtlas's own cancer epitope list against the canonical human proteome,
**56.3% (98,193 / 174,465)** of catalogued cryptic "cancer" epitopes are *exact substrings of canonical
proteins* and are therefore canonical *self* by sequence (an identical sequence yields an identical
spectrum), whatever ORF they are attributed to — **≈800× a composition-matched shuffle null** (Fig 1),
so the overlap is structural, not short-peptide coincidence. Nor is it an artifact of peptide length:
restricting to the canonical HLA-I window (8–11mers) the rate is **53.9% (67,852 / 125,999)**, and it
rises monotonically with length (46.2% at 8mer to 58.1% at 11mer, 62.6% beyond) — as expected, since a
longer exact canonical substring is stronger, not weaker, evidence of canonical origin. Mechanistically these are predominantly **in-frame alternative ORFs of coding genes**: 88.8% of the canonical-self epitopes are exact substrings of the canonical protein of *their own annotated gene* (e.g. the IEAtlas ORF `RBM47_223aa` → RBM47), whereas the out-of-frame remainder that is *not* canonical-self maps to its own gene's canonical protein in **0.0%** of cases — an internal negative control confirming the substring test does not match spuriously (`src/darkproteome/ieatlas_frame_audit.py`). We use "canonical self" strictly as a
*sequence / MS-origin* criterion: an exact canonical-proteome substring is not sequence-novel, and the
spectrum alone cannot establish a non-canonical origin. This is not a claim that the source ORF is
untranslated — but the sequence alone cannot support the non-canonical, tumor-specific claim built on
it, because the T-cell receptor and the mass spectrometer both read the sequence, not the locus. That non-canonical sets can harbour
canonical sequences, and that pseudogene peptides are frequently identical to their parent protein, is
recognised, and canonical-subtraction or "uniquely-mapping" requirements are established best practice
[3]; a related large-scale non-canonical proteome/immunopeptidome study resolves ambiguous ORF
assignment by priority-mapping (start-codon confidence, Kozak strength, transcript expression) but does
not report a canonical-sequence overlap rate [24]; what has not been reported is *how large the rate
is* in these catalogues.

It is strongly
**class-resolved**, not uniform (Fig 1): in the primary cohorts, altORF and lncRNA-ORF peptides are
essentially never canonical substrings (Raja altORF 0.2%, 5/2,592 — a composition-matched null gives
≈0.8×, i.e. no enrichment; these are genuinely non-canonical). This is not in tension with the 88.8%
in-frame figure above: Raja's altORF class is defined as *out-of-frame* ORFs within a coding gene body,
whereas IEAtlas's canonical-self burden is concentrated specifically in *in-frame* ORFs of coding genes
— the same frame distinction, not a contradiction between "altORF" classes, so the contamination has a
precise mechanistic locus rather than being a property of altORFs generally. Pseudogene-ORF peptides,
by contrast, carry an intrinsic risk of this failure — realized at 37.1% (43/116) in the HCC cohort, each an exact
substring of its parent canonical gene (e.g. RPS3AP12→RPS3A), a consequence of pseudogene biology, not
negligence, that renders the peptide MS-unfalsifiable as pseudogene-derived rather than parent-derived.
The risk is avoidable by upstream filtering, not absent in careful studies: Raja's own pseudogene-ORF
claims are 0/98, because that pipeline already excludes parent-gene-matching peptides before
publication; HCC's did not, so the same biological hazard reached the published claim set. The
canonical-self burden is thus concentrated in the aggregator atlas and in the pseudogene class; the primary-cohort
altORF/lncRNA peptides are clean. Critically, it is **resource-specific and avoidable, not intrinsic to
cryptic-peptide catalogues**: at a comparable peptide-length profile, CrypticProteinDB carries **0.0%
(1 / 3,774)** canonical-self peptides versus IEAtlas's 56.3% — i.e. IEAtlas omits a canonical-overlap
filter that another public atlas already applies, making CrypticProteinDB a clean negative control for
the substring test.

Note the nulls differ by stratum: the ≈800-fold enrichment is an atlas-level
statement, whereas the per-cohort altORF enrichment is ≈1×. These four resources are not independent
corroborations, and should not be read as such: IEAtlas aggregates public immunopeptidome deposits that
include both cohorts, and 45 of the 48 cohort canonical-self peptides (Raja's 5 + HCC's 43) are already
present in IEAtlas's own list. More broadly, across all ORF classes (not just canonical-self peptides),
823 peptides overlap in total between {cohorts, CrypticProteinDB} and IEAtlas (743 of all cohort
peptides of any class; 80 of CrypticProteinDB's). The 56.3% and 37.1% figures are thus substantially
overlapping measurements of the same underlying peptides, not separate confirmations across resources —
each remains individually correct over its own stated denominator, but the four-resource framing should
be read as ORF-class-resolved slices of a largely shared corpus, not four independent measurements.

### A measurable tumour-specificity floor
Among the genuinely non-canonical remainder (after removing the canonical-self peptides), **9.1%
(6,976 / 76,272)** appear in IEAtlas's *own* normal-tissue epitope list, and nearly two-thirds of
those map to a critical normal organ (brain/heart/liver/lung/kidney) — a direct on-target/off-tumour
signal that any pipeline reusing these atlases as a tumour-antigen source inherits. For the
canonical-self pseudogene peptides specifically, an independent normal reference (the HLA Ligand
Atlas [12] benign HLA-I ligandome) directly observes 16 of 43 on normal tissue — a hard,
finite-sampling-conservative floor on non-specificity. Separately, normal-tissue expression (GTEx v8)
of the canonical parent gene itself is measured, not assumed, for all 43 — each is broadly expressed
(median breadth: all 54 tissues) — a specificity *concern* at the RNA-expression level, not a
presentation-level confirmation. Together, the non-novelty floor and this parent-gene expression
signal flag all 37% (43/116) of the HCC pseudogene "antigen" claims with at least the RNA-level
concern; 16 of those 43 additionally carry the stronger, presentation-level confirmation of direct
observation on the benign ligandome itself.

### No published claim carries a reusable four-axis evidence package
Turning from sequence to re-verifiability: of 306,844 claims, **none** strict-pass all four axes (95%
Wilson upper bound ≈ 0.0013%), and only **seven** carry enough reported evidence to be even *decidable*
on all four at once (Fig 2). This bar measures the *machine-readable reusability of the reported
record* — not whether an antigen is biologically real — and reusability is exactly the property a
downstream pipeline depends on when it treats a catalogue row as a validated target. The gap holds
across both strata, for different reasons: the atlas records were never built to carry a per-claim
T-cell assay or allele call, while the end-to-end cohorts complete the chain in their figures and prose
but not in reusable, per-claim form. Critically, the bar is not unfairly harsh — a panel of known-real
canonical cancer-testis antigens (MAGE/SSX), scored through the identical pipeline, **also fails** (0
strict survivors, N=144; 95% Wilson upper bound ≈2.6%). That a positive control fails is expected, and is the central point: it locates the
problem in the *form of the reporting* rather than the truth of the biology, and makes the remedy a
reporting standard rather than a wave of corrections. The published path from "a peptide was observed"
to "a validated tumor antigen" is rarely available in reusable, per-claim, machine-readable form.
The bar is satisfiable by construction, not vacuous: a claim reported with the class-decoy ledger
(Table 1) — a per-class decoy count, an assigned allele, quantified translation evidence, and a
queried normal-ligandome search space — would strict-pass the three reporting-addressable axes
(translation, presentation, specificity). The immunogenicity axis additionally requires the T-cell
assay itself to be reported in reusable form, which no reporting standard can manufacture on its
own; a claim carrying both the ledger and a reusable assay result would strict-pass all four. The
bar is therefore satisfiable, not
vacuous — the gap is unmet reporting, not an impossible standard.

### Where the evidence package breaks down
The 0% survival is forced jointly by two axes, not by knob-twisting on one (Fig 2). Source-ORF
translation rigour is reported, per claim, for ≈0% of the corpus. The HLA-presentation axis hits an
allele-assignment cliff: 99.96% of unique peptides (all but ≈80) carry no per-peptide HLA allele in
reusable form. The immunogenicity axis is empty of machine-readable human positives — across both
cohorts the maximal *reusable* human-validated set tops out at ≈40 peptides, each requiring
humanised-mouse data or manual figure digitisation to recover. An escalating-leniency analysis
confirms the result is structural: strict survivors remain 0 even after fully relaxing the source-ORF
and allele requirements, and reach only 41 when an entire axis is dropped (Methods).

## Results II — The class-specific FDR is only set-identifiable; the per-claim truth is a partial order

### The subsidy, made precise
At a fixed acceptance threshold and a fixed statistical unit (PSM-, peptide-, or unique-peptide-level —
mixing them invalidates what follows), let `T_N`, `T_C` be the accepted target identifications in the
non-canonical and canonical classes, `T = T_N + T_C`, and class fraction `f = T_N/T`. Let `θ_N`, `θ_C`
be the (unobserved) class false-discovery proportions and `q` the global target-level false-discovery
proportion on the same unit. By the law of total proportion,
`q = f·θ_N + (1-f)·θ_C`. A reported global FDR controlled at α gives only `q ≤ α`; because the canonical
class dominates the accepted set (`1-f ≈ 1`), it absorbs most of the false budget and leaves `θ_N` nearly
free. This is the formal content of the established observation that a pooled FDR underestimates the rare
class's FDR [4,5,10,11].

### The class FDR is set-identified — and the budget bound is as tight as the data allow
Solving the identity for `θ_N` and letting the unobserved `θ_C` range over [0,1] gives the **sharp
identified set** `Θ_N(q,f) = [ max(0,(q-(1-f))/f), min(1, q/f) ]` — the same mixture-bound construction
as the classical ecological-inference "method of bounds" [23], applied here to a setting that field
lacks: the class rate is *directly observable* via the accepted decoy count (below), so the
identifiability gap is a reporting omission rather than the permanent impossibility ecological
inference confronts. When only `q ≤ α` is reported it widens
to **`Θ_N = [0, min(1, α/f)]`**. The upper endpoint is what we term the *budget bound* (all false budget assigned to the
non-canonical class). We prove it **sharp**: for any `θ*` in the set, the implied canonical rate
`θ_C* = (q - f·θ*)/(1-f)` lies in [0,1], and a law in which non-canonical targets are false at rate `θ*` and
canonical targets at rate `θ_C*` reproduces the same reported `(q,f)` (Methods, Prop. 1). Two consequences:
**(a)** no estimator can tighten the budget bound from the reported summary statistics alone — it is the best
possible; **(b)** the lower endpoint is 0 — a global FDR statement forces no false discovery into the
non-canonical class, so summary statistics can never show a class FDR is *high*, only bound how high it
could be.

For the ovarian cohort [1] (α = 3%): the headline noncoding-cryptic subset is <1% of accepted
identifications (`f < α`), so its class FDR set is **[0, 1]** — entirely unconstrained by the reported 3%.
The broader cryptic class (`f = 4.4–7.1%`) gives `Θ_N = [0, 42–68%]`, an upper bound 14–23× the reported
global figure (Fig 3a). This bounds what the *published statistics* support, not the realized class FDR —
which direct re-searches with class-specific decoys have shown can in fact be large (e.g. 36% at a nominal
1% [19,20]); that the realized value can be high yet remain invisible to the reported global figure is
precisely the identifiability gap, not a claim that this cohort's realized class FDR is high.

### A diagnostic: the class split follows an *effective* ρ, not the global candidate ratio
One might predict the class split of false (decoy) discoveries from search-space sizes. Under a common
null with per-spectrum top-hit competition, the expected accepted-decoy ratio is a *threshold-weighted*
average of per-spectrum class fractions, **`ρ_eff(τ) = Σ_s w_s π_sN / Σ_s w_s π_sC`** (`w_s = 1-F(τ)^{m_s}`,
`π_sg = m_sg/m_s`; derivation in Methods) — **not** the global ratio `ρ_global = Σ m_sN / Σ m_sC`. The two
can diverge arbitrarily (e.g. one spectrum `(m_N,m_C)=(900,100)` plus a hundred `(1,9)` spectra give
`ρ_global = 1` but `ρ_eff = 0.21` at a per-candidate tail probability `t = 1-F(τ) = 0.01` — the same
construction gives `ρ_eff ≈ 0.96` at `t = 10⁻⁴` and `ρ_eff ≈ 0.12` as `t → 1`, so the toy value is
threshold-specific, not construction-intrinsic), and a direct simulation confirms the accepted `D_N/D_C`
tracks `ρ_eff`, not `ρ_global`, as the threshold relaxes (Fig 3b — a separate, larger
heterogeneous-population simulation, not the two-spectrum toy example above; Methods). We treat `ρ_eff`
as a diagnostic, not a measured quantity: the common-null assumption it requires is least secure exactly
where `ρ_eff` departs from `ρ_global` (non-tryptic cryptics plausibly have `F_N ≠ F_C`), and in the
stringent-threshold regime where immunopeptidomics operates `ρ_eff → ρ_global` — but the reported record
does not establish which per-candidate tail probability `t` a given cryptic search actually corresponds
to, so whether a cohort sits in that stringent regime or the divergent one is itself an empirical
question, resolved only by the class-labelled decoy counts `D_N(τ), D_C(τ)`, absent from the reported
record.

### A single statistic collapses the uncertainty
Both results point to one missing datum. The class FDR becomes **consistently estimable** — by the ordinary
target-decoy estimate `θ̂_N = (D_N+1)/T_N` [6,7] — exactly when the per-class accepted decoy count `D_N` is
reported on the claim's unit and convention (equivalently, when the per-class decoy *split* is recoverable,
since the reported global FDR and class fraction fix the total); absent it, the class FDR is **not
identifiable** from the published record (Methods, Prop. 2). `D_N` is the single number that turns an
unfalsifiable interval into a checkable estimate. We stress the scope of that claim: `D_N` closes the
*reconstructibility* gap — it makes the class estimate *computable* — but whether the estimate is
well-*calibrated* is a separate question, since the target-decoy equal-chance assumption can be fragile for the
inflated, often non-tryptic cryptic search space; only an entrapment-based check settles it. The reporting
standard below therefore asks for `D_N` together with a matched entrapment result.

### Per-claim truth is an upper-bound ranking, not a probability
Propagated to the object a reader ultimately wants — the probability a single claim is true — the
set-identification has a sharp consequence. Lacking reusable immunogenicity evidence (Results I), we
bound *presented-ligand* truth `Y = M·G·S` (correctly-identified × genuinely non-canonical ×
tumor-specific), not full antigen truth. The three channels are coupled (search-space inflation ties
M–G, pseudogene-parent identity ties G–S, detectability ties M–S), so the per-claim probability is not a
product of marginals; the sharp Fréchet bounds (Methods) give an interval. The M-marginal lower bound is
`1 - min(1,α/f)`, which reaches 0 only in the limit `f ≤ α` (the headline noncoding-cryptic subset); more
generally it is strictly positive (e.g. `0.32–0.58` for the broader cryptic class). On the present data,
though, every claim's *joint* probability still has lower bound 0 — not because of the M-marginal, but
because the non-novelty and specificity floors (Results I) each independently reach a lower bound of 0
for any single claim, and that alone is enough to collapse the Fréchet joint lower bound. So every such
claim's probability is **not point-identified** — claims are comparable only by their sharp *upper* bound
(a partial order). On our data a random Raja altORF claim has `P ≤ 0.998` versus an HCC pseudogene claim
`P ≤ 0.629` (the 37.1% non-novelty floor), and an exact-canonical-substring claim collapses to `P = 0`
(Fig 4): an upper-bound ranking, not a calibrated probability. Crucially, our three audit instruments *are* the marginals of this
object — the non-novelty floor is `1-g⁺`, the normal-presentation floor is `1-s⁺`, and the class-FDR
budget bound is the M-marginal — so the empirical audit (Results I) and the identifiability theory are one
construction, not two.

### Acting on an identified set: a triage rule
The partial order has an operational reading. With each claim's truth probability bounded to a sharp identified
interval `P(Y=1) ∈ [L, U]` (Methods), a user selecting targets at a benefit/cost break-even `p* = C/(B+C)` can act
by interval dominance: **reject** if `U < p*` (even the best compatible probability is too low), **accept** if
`L ≥ p*` (even the worst clears the bar), and otherwise treat the claim as **not identified** — collecting the one
missing statistic (the per-class decoy count `D_N` for the M-factor; a normal-ligandome search for the S-factor)
rather than guessing a probability. Ranking the not-identified claims requires a stated risk posture (optimistic,
by `U`; conservative, by `L`; minimax-regret weighs both): the reported record yields the *set*, and the user
chooses the rule [18]. On the present data this rarely licenses acceptance — class-FDR lower bounds are zero — so
its value is a valid rejection/triage rule and a precise statement of which single measurement would change the decision.

## Results III — A reporting standard, a neutral benchmark, and its feasibility

### A minimal reporting standard makes claims re-verifiable
The identifiability results convert directly into a short, actionable reporting standard — a minimum-information **class-decoy ledger** that travels with each non-canonical antigen claim and supplies exactly the
statistics that collapse the identified sets to checkable estimates (Table 1; reference implementation `class_decoy_ledger.py`). The single most
important field is the **per-class accepted decoy count `D_N`**, on the same unit (PSM / peptide) and
FDR convention as the claim; with it, the class FDR becomes the ordinary target-decoy estimate rather
than an unfalsifiable interval. The remaining fields close the other axes: the class target counts and
fraction (`T_N, T_C, f`); the FDR convention and unit, stated explicitly; the per-peptide HLA allele;
machine-readable source-ORF translation evidence (periodicity or protein-level statistics) and the
peptide's canonical-substring status; and, for tumor specificity, the normal-ligandome search space
actually queried with its class-specific decoy counts. None of these requires new experiments — they
are quantities the original analysis already computed and discarded at reporting time.

### A neutral, leakage-controlled benchmark
A reporting standard fixes new claims but cannot, by itself, adjudicate the existing record or grow
the missing ground truth. We therefore propose a neutral, openly governed benchmark for non-canonical
tumor antigens, on the model of community benchmarks in adjacent fields where the party controlling
the ground truth is distinct from those being evaluated (organizer ≠ competitor ≠ assessor). Two tiers
follow the evidence: a **presentation tier**, buildable now from public data, scoring claims on the
re-verifiable axes under one uniform class-specific-FDR pipeline; and a small **immunogenicity tier**,
grown prospectively with time-separated or embargoed ground truth to prevent the circularity that
would otherwise let a method be tuned on its own test set. The audit of Results I is the natural seed —
it defines the corpus, the axes, and the failure modes the benchmark must measure.

### The empirical check is feasible for one researcher
The standard asks authors for `D_N`; for the already-published record, the same statistic is
recoverable independently. The class-labelled decoy count `D_N(τ)` can be obtained by re-searching
public raw spectra under a class-preserving target-decoy design, and an entrapment search
[8,9] independently tests whether the target-decoy control is actually
calibrated for the cryptic class. A power analysis shows this is within reach of a single researcher on
public data: for an entrapment search at multiplicity `k` over a class of `n` accepted non-canonical
discoveries with true FDR `θ`, the estimate has standard error `√(θ/(kn))`, so for the ovarian cohort's
headline class (`n ≈ 311`) a single matched entrapment (`k = 1`) already resolves the class FDR to a
95% half-width of ≈9 percentage points — enough to distinguish, e.g., a 40% from a 70% class FDR. This precision is governed by the absolute number of accepted non-canonical discoveries `n`, not by the class's vanishing fraction of the catalogue, so the rarity of the class is not itself an obstacle to estimating its FDR. The
decisive comparison is then whether the measured ratio `D_N(τ)/D_C(τ)` tracks `ρ_eff(τ)` or the global
candidate ratio across thresholds — the test that converts the identifiability bound into a number.
This re-analysis is the natural next study; it is not required for the present results, which stand on
the published record and the identifiability theory.

---

**Table 1. The class-decoy ledger: a minimal per-claim reporting standard.** Each field is a quantity the original analysis already computes; reporting it makes the corresponding identified set checkable. All FDR fields are on one stated statistical unit (PSM / peptide / unique-peptide) and threshold.

| Field to report (per claim) | What it makes re-verifiable |
|---|---|
| Statistical unit + FDR convention (D/T, (D+1)/T, …) | Makes every FDR figure interpretable; without it the class-FDR identity cannot be applied at all |
| Class target counts `T_N`, `T_C` (⇒ fraction `f`) | Yields the sharp budget bound `[0, min(1, α/f)]` |
| **Per-class accepted decoy count `D_N`** (same unit + threshold) | **Collapses the class FDR to the target-decoy estimate `(D_N+1)/T_N` — the single highest-value field (computability)** |
| *(paired)* a matched entrapment class-FDP check | Calibration: confirms the `(D_N+1)/T_N` estimate is unbiased — target-decoy alone gives computability, not calibration |
| Per-peptide HLA allele | The presentation axis (resolves the allele-assignment cliff) |
| Source-ORF translation evidence (Ribo-seq periodicity % or protein-level FDR + ≥2 unique ≥9-aa peptides), machine-readable | The translation axis |
| Peptide canonical-substring status vs a named proteome version | The non-novelty floor (the `G` factor) |
| Normal-ligandome search space queried + its class-specific decoy counts | The specificity axis (the `S` factor) — distinguishes "absent" from "not searched" |
| *(recommended)* per-spectrum candidate counts `m_sN`, `m_sC` | Recovers the effective ρ; diagnoses the global-ρ shortcut |

## Discussion

We set out to measure whether the published dark-proteome tumor-antigen record can be independently
re-verified. We found, first, a concrete measured defect — **56.3% of catalogued cryptic "cancer"
epitopes are canonical by sequence** (≈800-fold over a composition-matched null, strongly
class-dependent, and to our knowledge unquantified by any atlas) — and, more broadly, that the record as
reported essentially cannot be re-verified: of 306,844 claims, none carry reusable per-claim evidence
across all four evidence axes, and known-real canonical antigens fail the same bar. The contribution of this work is to show that this is not merely uneven reporting but has a
structural core. The class-specific false-discovery rate underwriting the presentation axis is only
*set*-identified from the statistics papers report; the common budget bound is provably the tightest
such bound; the quantity that governs the class allocation at a threshold can be a
threshold-weighted *effective ρ* rather than the global candidate ratio (in a regime we characterize); and, propagated to the claim
level, a calibrated per-claim truth probability is not identifiable at all — the honest referee output
is a partial order. A single reported statistic, the per-class decoy count `D_N`, collapses the central
uncertainty.

That the canonical controls fail the same audit is the result's most important feature, not a caveat
(Results I) — the remedy is accordingly a standard, not a wave of corrections. Our recommendation mirrors and extends the field's own
reform trajectory — the HUPO/HPP [17] and GENCODE Ribo-seq ORF [16] reporting standards, entrapment-based
assessment of FDR control [8], the established recognition that class-undiscriminated FDR
underestimates rare-class error [4,5,19], and recent calls for verifiable FDR reporting more broadly [21].
Concurrent work proposes general verifiable-FDR diagnostics — the scope, calibration, and stability of any
reported FDR, with decoy-inclusive outputs for independent checking [21]; our contribution is the
complementary, *class-specific* core — a proof that the non-canonical class FDR is unidentifiable from
summary statistics without the per-class decoy count, and the single field (`D_N`) that closes it. Where
prior work asks *how to compute* a class-specific FDR from raw data [5,10,11,22], we ask what a reader can
conclude from the published *summary statistics*, and answer it with a partial-identification analysis [18]
and a concrete class-decoy ledger that makes the needed quantities travel with the claim. The set-identification is not
specific to immunopeptidomics: the same mixture identity bounds any rare relevant subclass's FDR from a pooled
target-decoy threshold — proteogenomics, rare post-translational modifications, any setting where a global
multiple-testing threshold is applied to a small relevant subset — so the contribution is the publication-record
estimand (what a reader can identify from the summaries papers report) and the class-decoy ledger that closes it, not the
elementary mixture algebra itself.

Several limitations bound the scope. The full four-axis chain is measured on the two cohorts that
report it end to end; the atlas scale-up tests the generality of the reporting gap but not of the
immunogenicity axis, whose reusable human ground truth is genuinely thin (≈40 peptides, largely
figure-locked). The identifiability results are statements about what the published record *supports*;
they neither assert nor deny that any specific antigen is real. Establishing the realized class FDR for
a given study — and thereby testing empirically whether `D_N(τ)/D_C(τ)` tracks the effective or the
global ρ — requires re-analysis of raw spectra, which we show is feasible for one researcher on public
data but leave to a dedicated study. Finally, the composition-matched null underlying the non-novelty
enrichment is stochastic and is reported as a bound.

As non-canonical antigens move toward the clinic, the cost of a claim that cannot be re-verified rises. The fix
we propose is cheap — a short class-decoy ledger of statistics the original analysis already computed — and it
converts a class of claims that today cannot be independently checked into one that can.

## Methods

### Corpus assembly and statistical unit
We assembled 306,844 claims (180,650 unique peptides): two end-to-end immunopeptidomics cohorts —
ovarian [1] and HCC [2] — plus IEAtlas (293,222 cancer epitopes) and
CrypticProteinDB (9,101; 8,879 cryptic immunopeptides + 222 allele-resolved epitopes), with 144
canonical cancer-testis antigen claims (MAGE/SSX; 37 unique peptides) as a positive control. IEAtlas
reanalyzes public immunopeptidome deposits and therefore overlaps the other resources at the peptide
level (823 peptides shared between {cohorts, CrypticProteinDB} and IEAtlas: 743 cohort, 80
CrypticProteinDB — see Results I for the reading this supports). Claims are deduplicated and stored in a 19-column contract
(Supplement). Headline rates use resource-specific
denominators, stated at each use: the 56.3% canonical-self rate is over the 174,465 unique IEAtlas
cancer peptides entering the substring test (the 293,222 catalogue rows after deduplication and the
8–12-mer presentation filter), the 9.1% normal-tissue rate is over the 76,272 genuinely non-canonical
remainder, and 306,844 / 180,650 count claims / unique peptides across all resources. **Unit
discipline:** every
class-FDR quantity below is computed at one fixed statistical unit (PSM-, peptide-, or
unique-peptide-level) and one acceptance threshold; mixing units (e.g. PSM-level decoys against
peptide-level claims) invalidates the identity in §"Identified set."

### Four-axis, as-reported scoring
Each claim is scored on four axes with the field's standards: source-ORF translation (Ribo-seq ≥70%
periodicity *or* the protein-existence consensus bar [3] — applied to the ORF, never to
a single 8–12-mer ligand); HLA presentation (eluted 8–12-mer with an assigned allele); tumor
specificity (broad normal panel); immunogenicity (autologous/HLA-matched human T-cell assay).
Underspecified fields are labelled *unverifiable*, not failed. The auditor independently recomputes
the consensus-bar verdict from raw fields. Robustness is assessed with an escalating-leniency ladder
(Supplement): strict survivors remain 0 through relaxation of the source-ORF and allele requirements
and reach 41 only when an entire axis is dropped.

### Non-novelty floor
A peptide that is an exact substring of the SwissProt (reviewed) human proteome is canonical-self by
sequence (identical sequence ⇒ identical spectrum). We compute exact and I/L-collapsed substring rates
per ORF class, with a composition-matched shuffle null (seeded) and Wilson 95% intervals; pseudogene
parent genes are mapped from the `GN=` field. The shuffle-null enrichment multiplier is stochastic and
reported as a bound (≈800×, i.e. >several-hundred-fold).

### Specificity floor
For canonical-substring peptides the correct normal reference is canonical-space: the HLA Ligand Atlas
benign HLA-I ligandome. We test exact membership and map hits to tissue criticality. Absence from a
finite atlas is treated as non-informative (a conservative floor), never as positive evidence of
specificity.

### Class-specific FDR — identified set and identifiability
Fix the unit and threshold, and **condition on the observed accepted counts** `T_N, T_C` (`T=T_N+T_C`,
`f=T_N/T`); let `θ_N=E[V_N/T_N]`, `θ_C=E[V_C/T_C]` be the (conditional) class false-discovery proportions
and `q=E[(V_N+V_C)/T]` the global one. Because `f` is then fixed, **`q = f·θ_N + (1-f)·θ_C`** holds exactly.
We take a reported "global FDR controlled at α" to mean the target-decoy estimate bounds the expected
global FDP, `q ≤ α` (the standard, calibration-dependent reading).

**Proposition 1 (sharp identified set).** The set of `θ_N` consistent with a reported `(q,f)` and
`θ_C∈[0,1]` is `Θ_N(q,f)=[max(0,(q-(1-f))/f), min(1,q/f)]`, and every value is attained.
*Proof.* `θ_N=(q-(1-f)θ_C)/f` is affine-decreasing in `θ_C`; over `θ_C∈[0,1]` it spans
`[(q-(1-f))/f, q/f]`, intersected with `[0,1]`. Attainment: given `θ*∈Θ_N`, set
`θ_C*=(q-f·θ*)/(1-f)`, which lies in `[0,1]` by construction; a law in which each accepted
non-canonical (resp. canonical) target is independently false with probability `θ*` (resp. `θ_C*`) has
expected/conditional class FDPs `(θ*,θ_C*)` and global `q`. ∎ When only `q≤α` is reported, the union
over `q∈[0,α]` is `[0, min(1,α/f)]`; the upper endpoint (the "budget bound") is, by Prop 1, the tightest
upper bound that is a function of `(α,f)` alone.

**Proposition 2 (reconstructibility).** Given `(q,f)`, `θ_N` is consistently estimable from the record iff
the record determines the per-class accepted decoy *split* — equivalently the count `D_N`, since (under
the stated unit, threshold, and FDR convention) `(q,f)` fix the total decoy count. *Proof.* Sufficiency:
under concatenated target-decoy competition the
class-restricted estimator `θ̂_N=(D_N+1)/T_N` consistently estimates `θ_N` (conservatively, under the
standard equal-chance assumption).
Necessity: if the record fixes only `(q,f)` with `f>0` and `Θ_N` nondegenerate, the laws "all false
discoveries canonical" (`θ_N=0`) and "all non-canonical" (`θ_N=min(1,q/f)`) share the reported `(q,f)`;
e.g. `T_N=100, T_C=9900, q=0.01` admits both `θ_N=0` and `θ_N=1`. ∎

### The effective ρ
Under a common null in which a spectrum's candidate scores are i.i.d. with CDF `F` and per-spectrum
top-hit competition, the winning candidate for spectrum `s` lies in class `g` with probability
`m_sg/m_s` — independent of the winning score, because the index of the argmax of exchangeable
variables is independent of its value — and the winner exceeds `τ` with probability `1-F(τ)^{m_s}`.
Hence the expected accepted decoy count in class `g` is `μ_g(τ)=Σ_s (m_sg/m_s)(1-F(τ)^{m_s})`, and
`ρ_eff(τ)=μ_N/μ_C=Σ_s w_s π_sN / Σ_s w_s π_sC`, `w_s=1-F(τ)^{m_s}`, `π_sg=m_sg/m_s`. In the rare-tail
regime `1-F(τ)^{m_s}≈m_s(1-F(τ))`, `w_s∝m_s` and `ρ_eff→ρ_global=Σ m_sN/Σ m_sC`; otherwise they differ.
`ρ_eff(τ)` is recoverable from class-labelled decoy counts `D_N(τ),D_C(τ)` or per-spectrum candidate
counts plus the null law — neither in the reported record. This assumes a single null `F` shared across
classes; if class-specific score distributions differ (`F_N ≠ F_C`, plausible for non-tryptic cryptic
candidates) the winner-class probability `m_sg/m_s` no longer holds and `ρ_eff` must be measured rather
than derived.

### Per-claim Fréchet bounds
With `Y=M·G·S` and interval marginals only, the sharp Fréchet bounds on the joint are
`P(Y=1)∈[max(0,m⁻+g⁻+s⁻-2), min(m⁺,g⁺,s⁺)]`. The M-marginal lower bound is `m⁻=1-min(1,α/f)`, which
reaches 0 only when `f≤α` (the headline noncoding-cryptic subset); more generally `m⁻>0` (e.g.
`0.32–0.58` for the broader cryptic class). In the present worked examples the joint lower bound is 0
regardless of `m⁻`, because `g⁻=s⁻=0` (any single claim could be the non-novel or normal-presented case
within its own floor) already forces `max(0,·)=0`; claims are then comparable only by upper bound — a
partial order. Worked example: Raja altORF `g⁺=1-5/2592=0.998 ⇒ P≤0.998`; HCC pseudogene `g⁺=1-43/116=0.629`,
`s⁺≤1-16/116 ⇒ P≤0.629`; strict upper-bound ordering (not set-dominance; both sets reach 0), gap 0.37; an exact canonical-substring HCC claim has `G=0 ⇒ P=0`.

### Simulation
`scripts/verify_effective_rho.py` (numpy, seed 0). At per-candidate tail probability `t=1-F(τ)`, a spectrum's
decoy top-hit is accepted with probability `1-(1-t)^{m_s}` and is non-canonical with probability
`m_sN/m_s`; `D_N(t),D_C(t)` accumulate over Monte-Carlo replicates. We confirm (i) the two-spectrum
construction (`ρ_global=1`, `ρ_eff=0.213` at `t=0.01`; `4.69` swapped) exactly, and (ii) in a heterogeneous regime
empirical `D_N/D_C` tracks `ρ_eff(t)` to within Monte-Carlo error while both drift from `ρ_global` (both
checks validate the derivation under the common-null model, not its fit to real spectra).

### Entrapment power
For an entrapment search at multiplicity `k` over `n` accepted non-canonical discoveries with true FDR
`θ`, accepted entrapment hits `E~Poisson(knθ)` and `θ̂=E/(kn)`, so `SE(θ̂)=√(θ/(kn))` and a 95%
half-width `≈1.96√(θ/(kn))`. At `n=311, θ=0.7, k=1`: `≈9.3` percentage points.

### Software, data, reproducibility
All inputs are public (PRIDE PXD055609; the two cohort supplements; IEAtlas; CrypticProteinDB; the HLA Ligand
Atlas; UniProt/SwissProt) — accessions in Supplement. All analyses are scripted (`src/darkproteome/`,
`scripts/verify_effective_rho.py`); code and derived tables are included in this repository.
A reference implementation of the reporting standard (Table 1), `class_decoy_ledger.py`,
emits the per-class decoy ledger from standard target-decoy outputs (a mokapot/Percolator table or a
search-engine pepXML) and interoperates with general verifiable-FDR tooling [21]; a worked example on the
PXD055609 deposit is included. FDR conventions and the statistical unit are stated for every reported quantity.

---

## Data and code availability
All inputs are public: PRIDE PXD055609; the two cohort supplements; IEAtlas; CrypticProteinDB; the HLA Ligand Atlas; UniProt/SwissProt (accessions in Supplement). All analyses are scripted (`src/darkproteome/`, `scripts/verify_effective_rho.py`); code and derived tables are included in this repository. A reference implementation of the Table 1 reporting standard (`src/darkproteome/class_decoy_ledger.py`) emits the per-class decoy ledger from standard target-decoy outputs (mokapot/Percolator or pepXML) and interoperates with general verifiable-FDR tooling [21]; a worked example on the PXD055609 deposit is provided (`examples/`).

## Author contributions / Competing interests
R.J. designed and performed all analyses and wrote the manuscript. The author declares no competing interests.

## Figure legends

**Figure 1. Non-novelty and inherited specificity floors.** (a) Exact canonical-substring rate by ORF class and atlas: altORF/lncRNA and Raja pseudogene ≈0%, HCC pseudogene-ORF 37.1%, and the IEAtlas catalogue 56.3% versus a clean CrypticProteinDB at 0.0% — the clean atlas applies the canonical-overlap filter IEAtlas omits. (b) Of IEAtlas cancer "cryptic" epitopes, 56.3% are canonical-self by sequence; of the genuine remainder, 9.1% appear on normal tissue. Separately, of the HCC cohort's own canonical-self pseudogene peptides (panel a), 16/43 are directly observed in the benign HLA Ligand Atlas.

**Figure 2. The re-verifiability audit.** Of 306,844 published claims (180,650 unique peptides), none strict-pass all four evidence axes; the known-real canonical positive control (MAGE/SSX) collapses identically. Counts surviving each successive axis are shown on a linear scale so each step's true attrition is visually comparable; the two terminal zero-survivor counts are plotted at their true value (x = 0). The first step's ratio (180,650→82,453) is pooled across the full corpus and is a distinct measurement from the IEAtlas-specific 56.3% canonical-self headline reported in Results I — different denominator and scope, not a discrepancy (Methods). The specificity step here (82,453→75,471, 8.5%) is pooled across all sources; it is a distinct measurement from the IEAtlas-specific 9.1% normal-tissue rate reported in Results I, which uses IEAtlas's own genuinely-non-canonical remainder (76,272) as its denominator.

**Figure 3. The identifiability theory.** (a) The class-specific FDR is only set-identified; the sharp upper bound is `min(1, α/f)`, giving ≤ 42–68% for the broad cryptic class (`f` = 4.4–7.1%) at the reported α = 3%. (b) A diagnostic — the class decoy split follows a threshold-weighted *effective ρ*, not the global candidate ratio: simulated `D_N/D_C` tracks `ρ_eff(τ)` and drifts from `ρ_global` as the threshold relaxes; the shaded band marks this simulation's stringent end, where `ρ_eff ≈ ρ_global` — consistent with, but not a numeric calibration of, the stringent-threshold regime Results II argues immunopeptidomics operates in.

**Figure 4. Per-claim truth is a partial order.** Claims admit only a sharp upper bound on P(true), not a point value: a Raja altORF claim P ≤ 0.998, an HCC pseudogene claim P ≤ 0.629, and an exact-canonical-substring claim P = 0. Both non-degenerate intervals share the lower bound 0, so this is an upper-bound ranking, not set-dominance. Intervals are drawn open-capped at the upper bound (a bound, not a confirmed value); the filled diamond marks the point-identified P = 0 case.

## References
*DOIs / PMCIDs verified for this manuscript; author lists abbreviated. In-text citations use numbered
references keyed to this list.*

1. Raja R, Mangalaparthi KK, Madugundu AK, et al.; Curtis M. Immunogenic cryptic peptides dominate the antigenic landscape of ovarian cancer. *Sci Adv* 2025;11:eads7405. doi:10.1126/sciadv.ads7405. (PMC11837991).
2. Camarena ME, Theunissen P, Ruiz M, et al.; Albà MM. Microproteins encoded by noncanonical ORFs are a major source of tumor-specific antigens in a liver cancer patient meta-cohort. *Sci Adv* 2024;10:eadn3628. doi:10.1126/sciadv.adn3628.
3. Prensner JR, et al. What can Ribo-seq, immunopeptidomics, and proteomics tell us about the noncanonical proteome? *Mol Cell Proteomics* 2023;22:100631. doi:10.1016/j.mcpro.2023.100631. (PMC10506109).
4. Zhang B, Bassani-Sternberg M. Current perspectives on mass spectrometry-based immunopeptidomics: the computational angle to tumor antigen discovery. *J Immunother Cancer* 2023;11:e007073. doi:10.1136/jitc-2023-007073. (PMC10619091). — explicitly states that class-undiscriminated target-decoy competition underestimates the non-canonical-peptide FDR.
5. Fu Y, Qian X. Transferred subgroup false discovery rate for rare post-translational modifications detected by mass spectrometry. *Mol Cell Proteomics* 2014;13:1359–1368. (PMC4014291).
6. Elias JE, Gygi SP. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nat Methods* 2007;4:207–214. doi:10.1038/nmeth1019.
7. Keich U, Kertesz-Farkas A, Noble WS. Improved false discovery rate estimation procedure for shotgun proteomics. *J Proteome Res* 2015;14(8):3148–3161. doi:10.1021/acs.jproteome.5b00081. (PMC4533616). — rigorous target-decoy FDR estimation; grounding for the conservative (D+1)/T estimate.
8. Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U. Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nat Methods* 2025. doi:10.1038/s41592-025-02719-x. (PMC12240826).
9. FDRBench: entrapment-based FDP estimation (Noble Lab). https://github.com/Noble-Lab/FDRBench (Apache-2.0).
10. Fondrie WE, Noble WS. mokapot: fast and flexible semisupervised learning for peptide detection. *J Proteome Res* 2021;20:1966–1971. — grouped/group-specific FDR.
11. Kong AT, Leprevost FV, Avtonomov DM, et al. MSFragger: ultrafast and comprehensive peptide identification. *Nat Methods* 2017;14:513–520. — no-enzyme search; FragPipe group FDR.
12. Marcu A, et al. HLA Ligand Atlas: a benign reference of HLA-presented peptides to improve T-cell-based cancer immunotherapy. *J Immunother Cancer* 2021;9:e002071. doi:10.1136/jitc-2020-002071. (PMC8054196).
13. Laumont CM, et al.; Perreault C. Global proteogenomic analysis of human MHC class I-associated peptides derived from non-canonical reading frames. *Nat Commun* 2016;7:10238. doi:10.1038/ncomms10238.
14. Cai Y, Lv D, Li D, et al.; Xu J. IEAtlas: an atlas of HLA-presented immune epitopes derived from non-coding regions. *Nucleic Acids Res* 2023;51(D1):D409–D417. doi:10.1093/nar/gkac776. (PMC9825419).
15. Othoum G, Maher CA. CrypticProteinDB: an integrated database of proteome and immunopeptidome derived non-canonical cancer proteins. *NAR Cancer* 2023;5:zcad024.
16. Mudge JM, et al. Standardized annotation of translated open reading frames. *Nat Biotechnol* 2022;40:994–999. doi:10.1038/s41587-022-01369-0.
17. Deutsch EW, et al. Human Proteome Project Mass Spectrometry Data Interpretation Guidelines 3.0. *J Proteome Res* 2019;18(12):4108–4116. doi:10.1021/acs.jproteome.9b00542. (PMC6986310).
18. Manski CF. *Partial Identification of Probability Distributions.* Springer; 2003. — methodological grounding for the sharp identified-set and Fréchet-bound analysis.
19. Woo S, Cha SW, Na S, Guest C, Liu T, Smith RD, Rodland KD, Payne SH, Bafna V. Proteogenomic strategies for identification of aberrant cancer peptides using large-scale next-generation sequencing data. *Proteomics* 2014;14(23–24):2719–2730. doi:10.1002/pmic.201400206. (PMID 25263569; PMC4256132). — original quantified demonstration that a 1% combined target-decoy FDR corresponds to ≈36% FDR for the novel-peptide class (0.03% canonical); recommends a two-stage/separate FDR.
20. Choi S, Paek E. pXg: Comprehensive Identification of Noncanonical MHC-I-Associated Peptides From De Novo Peptide Sequencing Using RNA-Seq Reads. *Mol Cell Proteomics* 2024;23(4):100743. doi:10.1016/j.mcpro.2024.100743. (PMC10979277; PMID 38403075). — modern immunopeptidomics pipeline using separate class-specific target-decoy FDR for canonical vs non-canonical peptides.
21. Chion M, Godmer A, Douché T, Matondo M, Giai Gianetto Q. Verifiable False Discovery Rate Reporting in Proteomics via Scope, Calibration, and Stability Diagnostics. *bioRxiv* 2026. doi:10.64898/2026.04.16.718468. (R package diagFDR). — concurrent proposal for general verifiable-FDR reporting (scope/calibration/stability; decoy-inclusive outputs).
22. Lin A, Plubell DL, Keich U, Noble WS. Improving power while controlling the false discovery rate when only a subset of peptides are relevant. *J Proteome Res* 2021;20(8):4153–4164. doi:10.1021/acs.jproteome.1c00483. (PMC8489664). — subset-FDR estimation and the "neighbor peptide" failure mode.
23. Duncan OD, Davis B. An alternative to ecological correlation. *Am Sociol Rev* 1953;18(6):665–666. — the classical ecological-inference "method of bounds"; Proposition 1's mixture-bound identity is this construction's instance for target-decoy class-specific FDR.
24. Ruiz Cuevas MV, Hardy M-P, Hollý J, et al.; Perreault C, Yewdell JW. Most non-canonical proteins uniquely populate the proteome or immunopeptidome. *Cell Rep* 2021;34(10):108815. doi:10.1016/j.celrep.2021.108815. (PMC8040094). — resolves ambiguous ORF assignment by priority-mapping (start-codon confidence, Kozak motif, transcript expression), not canonical-sequence overlap; reports no comparable contamination rate.
