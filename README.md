# darkproteome

**An audit of a public non-canonical HLA-peptide atlas, and of the search library that feeds it.**

Code, derived data tables, and manuscript source accompanying *"Extensive canonical-sequence overlap
and unresolved source attribution in a public non-canonical HLA-peptide atlas and its search library"*
(Rom Jan).

The thesis, in one line: **before you trust a catalogue of candidate antigens, check whether its
entries can even be attributed to the source it names them for.**

---

> ### ⚠️ The paper in `manuscript/manuscript.md` is SUPERSEDED. Do not cite it.
>
> Its headline — *"0 of 306,844 claims pass all four evidence axes"* — **was true by construction**
> and has been withdrawn: the scorer could only record a pass via statistics that no source publishes,
> so no claim *could* pass. A second defect scored *unassayed* claims as immunogenicity **failures**.
> Both are the same error — absence of evidence recorded as evidence of absence — and both are
> precisely what this project exists to criticise. See `manuscript/REVISION_NOTICE.md`.
>
> ### The current paper is `manuscript/manuscript_v2.md` (+ `supplement_v2.md`).
>
> **The principles it applies are not ours, and it says so in its first paragraph.** That a peptide
> matching a canonical protein does not identify a non-canonical source, and should be **excluded**, is
> **Bedran et al. 2023** (*Cancer Immunol Res*), who also measured four catalogues at 1.4–5%. That
> *"most shared peptides should be dropped"* is **Kumar et al. 2022** (*Brief. Bioinform.*). The
> underlying protein-inference problem is **Nesvizhskii & Aebersold**. That a pooled FDR under-controls
> a minority class is **Woo et al. 2014**.
>
> **What had not been done was to test whether the public resources satisfy them.** This repository is
> that test. **The contribution is empirical.**

## The result

Tandem MS identifies a **sequence**, never a **locus**. Where the same sequence is encoded by both a
canonical protein and a non-canonical ORF, the spectrum does not choose between them.

**1 — Most of IEAtlas's cancer-catalogued sequences also occur in canonical proteins.**
**98,193 of 174,465 unique sequences (56.3%)** exactly match at least one protein in a frozen reviewed
human reference. Three ways of trying to make that go away all fail:

| robustness check | result |
|---|---:|
| headline (current reviewed reference) | **56.3%** |
| **era-correct** — Swiss-Prot **2022_01**, the release IEAtlas actually searched | **56.2%** |
| **length-standardized** | **56.3%** |
| if the atlas contained **no pseudogene ORFs at all** | **55.8%** |

Only **231 sequences (0.24%)** are overlaps a February-2022 analyst could not have seen, so the audit is
not anachronistic. Under one pipeline — same reference, same criterion, same peptide unit — two
catalogues that apply an explicit exclusion rule sit at **1 / 3,810** (CrypticProteinDB) and
**5 / 2,979** (Raja et al.).

**A canonical match does not establish canonical production.** It means the sequence cannot be *uniquely
attributed* to its nominated non-canonical source from sequence and ordinary tandem-MS evidence alone.
Sequence identity is symmetric.

**2 — The search library already carries extensive canonical overlap.** **34.1%** of nuORFdb v1.2's
distinct 9-mers are already canonical, against **1.0–2.4%** for the GENCODE Ribo-seq ORF sets (both
measured here, under one pipeline). But the library is **not sufficient** — Ouspenskaia et al. searched
the *same* nuORFdb and published a catalogue at 3%, so the exclusion step does the work — and it is not
the whole story: the catalogued rate *exceeds* the library's, and that excess arises **during
detection**. Canonical-overlapping sequences are detected across more cancer types (holding in 18 of 18
peptide-length strata), and ribosomal-gene ORFs are **2.51×** enriched among them in the catalogue while
being slightly **depleted** (0.91×) in the library they were drawn from.

**3 — The consequence is observable inside the resource.** Canonical-overlapping sequences appear in
IEAtlas's **own normal-tissue export** at **22.4%**, versus **9.1%** for the non-overlapping sequences of
the same catalogue — **length-standardized risk ratio 2.42×**, gene-clustered bootstrap **95% CI
[2.32, 2.52]**. **22,003 unique sequences — 12.6% of the cancer catalogue — are both canonical-compatible
and already present in the atlas's own normal export.** No external reference is needed to see this.

This is *consistent with, but not specific to*, greater detectability or expression of
canonical-compatible sequences. It does not show any sequence is canonically derived; it means such a
sequence **warrants normal-presentation review** before being treated as a tumour-restricted target.

## The proposed fix: label, don't delete

For each peptide, state whether the sequence is unique to the nominated ncORF within the searched space,
and if not, **list every compatible source locus**. That preserves the peptide *and* the ambiguity, and
is strictly more informative than the exclusion rule it generalises. Separately, publishing the per-class
accepted target/decoy counts makes a class-specific FDR estimate reconstructible — and every ncORF
library could publish its own latent canonical ambiguity: one number, minutes to compute, and **no
library currently does**.

## Reproduce

```bash
python3 manuscript/verify_manuscript.py         # every headline number, regenerated from artifacts
python3 src/darkproteome/scoring_conformance.py

python3 scripts/era_correct_reference.py        # 56.2% vs the reference IEAtlas actually searched
python3 scripts/consequence_robust.py           # RR 2.42x, gene-clustered 95% CI [2.32, 2.52]
python3 scripts/abundance_bias.py               # the detection effect, with the control that could kill it
python3 scripts/cross_catalogue.py              # same-pipeline comparison, length-standardized
python3 scripts/library_ambiguity.py            # nuORFdb 34.1% vs GENCODE 1.0-2.4%
python3 scripts/ouspenskaia_verify.py           # the published remedy, verified at source
```

`verify_manuscript.py` **fails the build on drift**. It also fails if the paper asserts any of **23
retracted phrasings or stale values**, drops its required prior-art citations, or if the ORF-class strata
fail to partition. The large public inputs (atlas exports, proteome FASTAs) are **not redistributed**; on
a clean checkout the checks that need them are reported as *skipped* rather than crashing. See
`data/SOURCES.md` and `data/external/README.md` to fetch them.

**Two scripts deliberately refuse to run** — `scripts/consequence.py` and `scripts/rule_predicts_rate.py`.
They produced results that were **retracted during review**: an invalid significance test that treated
clustered observations as independent, and an FDR argument built on the wrong null object (a false PSM is
drawn from the search database, not from a shuffle). They are kept for the record and print a retraction
notice instead of output, so their numbers cannot be quoted.

## Layout

```
manuscript/manuscript_v2.md      the paper      (manuscript.md is SUPERSEDED — see REVISION_NOTICE.md)
manuscript/supplement_v2.md      S1 pseudogene homology · S2 the class-FDR derivation
manuscript/verify_manuscript.py  regenerates every headline number; fails on drift
manuscript/figures_v2/           Figures 1–3
scripts/                         one script per result — each reproduces its own numbers
data/                            derived tables + the JSON artifacts the paper is verified against
src/darkproteome/                the analysis package (stdlib-only core)
```

## What this work does not claim

- Not that the biology is fake, or that any ncORF antigen is not real.
- Not that the canonical source is the correct one for any peptide. **MS identifies the sequence, never
  the locus.**
- Not that any resource acted improperly. All describe their own procedures accurately, which is the
  only reason this audit was possible. IEAtlas's published **Methods do not describe** a peptide-level
  canonical-exclusion rule — which is different from saying its pipeline applies none.
- Not that the exclusion criterion is a universal, field-wide rule. It is *published and recommended*. An
  atlas may legitimately retain shared sequences **if it labels them source-ambiguous**.
- **No fold-change** against the 1.4–5% rates published for four other catalogues. Those came from a
  different pipeline, reference and peptide unit; a ratio across pipelines is arithmetic, not
  measurement.
- Not that IEAtlas's *combined* library (nuORFdb + RPFdb + Translnc) has been measured. It has not — and
  34.1% is **not** a lower bound on it, because union rates are not monotone.

## License

MIT (`LICENSE`). Public data only; sources and versions in `data/SOURCES.md`.
