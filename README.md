# darkproteome

**An audit-first assessment of the dark-proteome tumor-antigen literature.**

Code, derived data tables, and manuscript source accompanying the preprint *"Canonical self in
cryptic cancer-epitope catalogues and the class-decoy ledger needed to verify non-canonical
antigen claims"* (Rom Jan).

Most "dark proteome" cancer-antigen claims, non-canonical ORFs (ncORFs) / microproteins
presented on HLA, are nominated on evidence that isn't independently re-verifiable from what
gets published. This project audits the published record against the field's own evidence bar,
characterizes the resulting gap as a statistical identifiability problem, and proposes the
minimal reporting fix.

The thesis, in one line: **before you trust a catalogue of candidate antigens, audit whether the
claims in it can even be re-verified from what was reported.**

---

> ### ⚠️ The paper in `manuscript/manuscript.md` is SUPERSEDED. Do not cite it.
>
> Its headline — *"0 of 306,844 claims pass all four evidence axes"* — **was true by construction**
> and has been withdrawn: the scorer could only record a pass via statistics that no source publishes,
> so no claim *could* pass. A second defect scored *unassayed* claims as immunogenicity **failures**.
> Both are the same error — absence of evidence recorded as evidence of absence — and both are
> precisely what this project exists to criticise. See `manuscript/REVISION_NOTICE.md`.
>
> ### The current paper is `manuscript/manuscript_v2.md`.
>
> **It claims no conceptual novelty, and says so in its first paragraph.**
>
> The standard it applies is **not ours**. That a peptide matching a canonical protein does not
> identify a non-canonical source, and must be **excluded**, is **Bedran et al. 2023**
> (*Cancer Immunol Res*), who also measured four catalogues at 1.4–5%. That *"most shared peptides
> should be dropped"* is **Kumar et al. 2022** (*Brief. Bioinform.*). The underlying
> protein-inference problem is **Nesvizhskii & Aebersold**. That a pooled FDR under-controls a
> minority class is **Woo et al. 2014**.
>
> **What had not been done was to test whether the public resources meet the standard.** This
> repository is that test. Its contribution is **three measurements, and no new ideas.**

## Headline result

**1 — IEAtlas is an order-of-magnitude outlier.** At least **56.3%** of its unique cancer-catalogued
peptide sequences (98,193 / 174,465) are exact substrings of reviewed canonical human proteins under
an explicitly defined reference *R* — against **1.4–5%** for every ncORF catalogue previously audited,
and **0.026%** / **0.17%** for two resources that apply an explicit exclusion rule (CrypticProteinDB;
Raja et al.). Neither ORF-class composition nor the false-discovery rate accounts for it.

**2 — The library explains it.** **34.1%** of nuORFdb v1.2's distinct 9-mers already occur in the
canonical proteome, against **1.0–2.4%** for the GENCODE Ribo-seq ORF sets. IEAtlas integrates
nuORFdb and applies no exclusion rule. *Ouspenskaia et al. searched the same nuORFdb and published a
catalogue at 3%* — so the library is necessary but not sufficient; the exclusion step does the work.
(Corroborated by an independent 8–11mer enumeration: 34.4% / 2.5%.)

**3 — The consequence is internal to the resource.** Canonical-overlapping "cancer" epitopes appear
in IEAtlas's **own normal-tissue set** at **22.4%**, versus **9.1%** for the non-overlapping epitopes
of the same catalogue (risk ratio 2.4; *z* = 74). **22,003 entries — 12.6% of every cancer epitope the
atlas catalogues — are both canonical-compatible and already observed on normal tissue by the atlas's
own measurement.** No external reference is needed to see this.

**What this does not claim.** That the biology is fake; that any resource manufactures, re-labels or
discards anything; that any individual peptide is canonically derived (**MS identifies the sequence,
never the locus**); or that any resource acted improperly — all describe their own procedures
accurately, which is the only reason this audit was possible.

**Reproduce:**
```bash
python3 scripts/rule_predicts_rate.py     # 1 — the outlier, with both confounds controlled
python3 scripts/library_ambiguity.py      # 2 — latent canonical ambiguity of the ncORF libraries
python3 scripts/consequence.py            # 3 — normal-tissue presentation, internal control
python3 scripts/ouspenskaia_verify.py     # the published remedy, verified at source
python3 manuscript/verify_manuscript.py   # regenerates every headline number; fails on drift
```

## What's here

- **`manuscript/`** — the paper itself. `manuscript.md` is the canonical source; `tex/` builds
  the typeset PDF (`tectonic manuscript/tex/manuscript.tex --outdir manuscript/tex`);
  `manuscript.html` is a browser-viewable render (`python3 manuscript/md_to_html.py` to
  regenerate); `figures/` holds all 4 figures as PDF + PNG, built by `make_figures.py`.
- **`src/darkproteome/`** — the analysis code, stdlib-first. Each entry point prints the
  headline numbers it's responsible for (see "Run the audit" below).
- **`tests/`** — regression tests for the class-decoy ledger and the ecological-inference
  diagnostic.
- **`examples/`** — a worked example of the class-decoy ledger tool run against the real
  PXD055609 deposit (Raja et al., ovarian immunopeptidome).
- **`data/`** — the derived/scored claim tables the audit runs on, plus every data source's
  license and provenance (`data/SOURCES.md`, `data/external/README.md`). Large public inputs
  (SwissProt, IEAtlas, CrypticProteinDB, the HLA Ligand Atlas, raw PRIDE deposits) are not
  included here for size reasons; `data/external/README.md` documents exactly how to
  re-download each one.
- **`scripts/verify_effective_rho.py`** — reproduces the Methods "Simulation" subsection cited
  throughout the manuscript.

## Six evidence dimensions (encoded in `src/darkproteome/evidence_dimensions.py`)

A claim is scored on **six evidence dimensions**, and **never as a single pass/fail survivor
count**. A joint pass/fail across dimensions the reported record cannot even decide measures the
scorer, not the claims.

Each dimension is scored on a **cumulative reporting ladder** — rung *k* counts claims that clear
rungs 1..*k*:

| Rung | Question |
|---|---|
| `asserted` | Is the experiment or analysis named for this claim? |
| `claim_linked` | Does an individual result travel *with* the claim (not just a study-level method)? |
| `quantitative` | Is that result a value a prespecified criterion can be applied to? |
| `modality_appropriate` | Does the measurement answer the endpoint actually being claimed? |
| **`adjudicable`** | **All of the above → the criterion can be applied independently.** |

| Dimension | Question |
|---|---|
| `source_translation` | Is the nominated ORF translated? |
| `hla_elution` | Was the peptide observed on HLA by MS? |
| `allele_restriction` | *Which allele* presented it? |
| `normal_presentation` | Is it absent from **normal HLA presentation**? |
| `human_tcell_assay` | Do human T cells respond? |
| `class_fdr_reconstructible` | Can the class-specific identification error be recomputed? |

Two rules govern every dimension, and violating either manufactures results:

- **A zero is not a finding until the record could have said otherwise.** A *structural* zero (no
  claim reports the field the criterion needs) and an *empirical* zero (claims report it, none
  clear it) look identical in a survivor count and mean opposite things. `adjudicable` separates
  them.
- **Absence of evidence is never evidence of absence.** An MS-observed peptide with no T-cell
  assay is not a negative immunogenicity result — it is *unassayed*. A ligand outside the HLA-I
  length window is not a contradiction — no source here states the MHC class, and 13–25mers are
  ordinary HLA-II ligands.

`adjudicable` is deliberately **not** source-attribution resolution (whether the record leaves
only the nominated source compatible). A record can carry excellent claim-linked translation
evidence while a canonical source remains perfectly compatible.

The auditor never recomputes FDR from raw spectra; it scores **as reported**, and always reports
**stratified** (atlas records / end-to-end cohorts / T-cell-tested), because a pooled denominator
dominated by atlas records hides the cohorts.

**Is the bar rigged to fail?** No, and this is enforced: `data/sample_claims.csv` carries
hand-built claims that report everything the criteria ask for, and `scoring_conformance.py`
asserts they come out adjudicable-and-supporting on *every* dimension — one witness per
independent route to a pass. If a criterion ever becomes unsatisfiable, the guard fails the build.

## Run the audit

Core auditor uses only the Python standard library. Three of the four commands below need no
external data; `tier1_nonnovelty.py` is the exception -- it reads the SwissProt FASTA (and the
IEAtlas normal-tissue table) from `data/external/`, so populate that first (see "Reproducing the
derived data tables" below) or it will exit with a `FileNotFoundError`.

```bash
python3 src/darkproteome/audit.py data/sample_claims.csv          # smoke test on synthetic rows
python3 src/darkproteome/probe_report.py data/claim_catalog_real.csv   # cohort headline + MAGE/SSX control
python3 src/darkproteome/tier1_nonnovelty.py                       # class-resolved non-novelty floor (needs data/external/)
python3 src/darkproteome/robustness.py data/claim_catalog_real.csv # leniency ladder + reusable-positive ceiling
```

Self-check that the environment is intact (with `data/external/` populated): `tier1_nonnovelty.py`
must print `5/2979`, `213/369 = 57.7%`, `0 mismatches`, `54.4%`.

The class-decoy ledger tool (`src/darkproteome/class_decoy_ledger.py`) is demonstrated end to
end on a real deposited run in `examples/`; see `examples/README.md`.

## Run the tests

Stdlib only, no test framework needed; each file is also its own runner:

```bash
python3 tests/test_class_decoy_ledger.py
python3 tests/test_eco_diagnostic.py
```

### Reproducing the derived data tables (needs the external inputs; see below)

All of these read from `data/external/`; set `$DARKPROTEOME_DATA` if you keep that data
elsewhere (see `data/external/README.md`).

```bash
python3 src/darkproteome/ingest_cohorts.py       # -> data/claim_catalog_real.csv (needs openpyxl)
python3 src/darkproteome/ingest_atlases.py       # -> data/claim_catalog_scaled.csv (307,318 claims; gitignored, large)
python3 src/darkproteome/reference_model.py      # -> data/claim_catalog_scored.csv + Fig. 1 survivorship funnel
python3 src/darkproteome/ieatlas_frame_audit.py  # the 56.3% canonical-self headline (Results I, Fig. 2)
python3 src/darkproteome/deepen_specificity.py           # falsification-tested normal-tissue overlap (Results I)
python3 src/darkproteome/pseudogene_specificity.py       # -> data/pseudogene_specificity.csv, the 16/43 HLA Ligand Atlas floor
python3 src/darkproteome/gtex_specificity.py             # -> data/gtex_pseudogene_specificity.csv, measured GTEx floor (43/43 parents expressed)
python3 src/darkproteome/gtex_class_specificity.py       # -> data/gtex_class_specificity.csv, class-resolved (altORF/lncRNA-ORF) GTEx floor
python3 src/darkproteome/lncrna_ensg_specificity.py      # -> data/lncrna_ensg_specificity.csv, lncRNA-ORF floor at ENSG resolution
```

## Data sources (all public)

PRIDE PXD055609 (raw immunopeptidomics) · two flagship cohort supplements (ovarian, HCC) ·
IEAtlas · CrypticProteinDB · GENCODE Ribo-seq ORFs · nuORFdb · HLA Ligand Atlas · GTEx v8 ·
UniProt/SwissProt. Full accessions, versions, and licenses in `data/SOURCES.md`.

## License

MIT (see `LICENSE`). All inputs are public data; see `data/SOURCES.md` for source-specific
licensing terms.
