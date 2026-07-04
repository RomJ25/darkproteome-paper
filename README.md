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

## Four independent evidence axes (encoded in `src/darkproteome/axes.py`)

A claim is scored on **four axes**, not one blended bar, because a single HLA-I antigen is
one short 8–12mer ligand, so a protein-existence standard is the wrong test for it. A claim is
a **strict survivor** only if it strict-passes all four:

| Axis | Question | Standard |
|---|---|---|
| **source_orf** | Is the source ORF plausibly translated? | Ribo-seq ≥70% periodicity **or** the Prensner protein bar (≤0.1% FDR, ≥2 unique ≥9-aa peptides) |
| **hla_presentation** | Is the peptide actually presented? | eluted-ligand MS, plausible 8–12mer, allele assigned (class-specific FDR confirmed downstream) |
| **tumor_specificity** | Absent from normal presentation/expression? | broad normal panel (HLA Ligand Atlas / GTEx-Recount3); "no public normal evidence", never "safe" |
| **immunogenicity** | Do T cells respond? | autologous / HLA-matched T-cell assay (in-vivo = weaker) |

The auditor never recomputes FDR from raw spectra; it scores each axis **as reported**, and
labels anything underspecified `unverifiable`. Headline = **strict survivor fraction** (passes
all four) with a Wilson 95% CI, plus per-axis survival so you see exactly where claims die.

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
must print `5/2979`, `43/116 = 37.1%`, `0 mismatches`, `54.4%`.

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
python3 src/darkproteome/ingest_atlases.py       # -> data/claim_catalog_scaled.csv (306,844 claims; gitignored, large)
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

## Headline result

**0 of 306,844** published non-canonical tumor-antigen claims carry per-claim reusable evidence
on all four axes at once; a known-real canonical control (MAGE/SSX) fails the identical bar,
locating the gap in reporting, not the underlying biology. **56.3%** of catalogued cryptic
"cancer" epitopes in the largest public atlas (IEAtlas) are exact substrings of the canonical
human proteome: canonical *self* by sequence. The class-specific false-discovery rate
underwriting these claims is only *set*-identified from what papers report; the fix is a
one-line reporting standard, a **class-decoy ledger**, that makes it re-verifiable.

## License

MIT (see `LICENSE`). All inputs are public data; see `data/SOURCES.md` for source-specific
licensing terms.
