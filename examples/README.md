# Worked example — the class-decoy ledger on a real deposited run

This shows the companion tool [`../src/darkproteome/class_decoy_ledger.py`](../src/darkproteome/class_decoy_ledger.py)
emitting the manuscript's proposed **class-decoy ledger** — the minimal reporting block that makes a
non-canonical FDR claim re-verifiable — directly from a **published** target-decoy deposit, in **one
command, with no re-search**. It is the concrete feasibility proof behind the reporting standard: the one
field that closes the identifiability gap (`D_N`, the per-class accepted decoy count) is a *free output* of
data the field already deposits.

## Command

```bash
python3 src/darkproteome/class_decoy_ledger.py \
    --pepxml T5-1546_Ovarian_MSFragger_MHC_ClassI_101821.pepXML \
    --alpha 0.03 \
    --canonical-fasta data/external/swissprot_human.fasta \
    --out examples/PXD055609_T5_ledger
```

**Input:** one MSFragger `pepXML` from **PXD055609** (Raja et al., ovarian immunopeptidome — the cohort the
manuscript audits end to end). The `pepXML` *intermediate* is the right input precisely because it **retains
decoys** (`rev_`) and class-labelable accessions; the FDR-filtered *final* reports strip both. mokapot /
Percolator PSM tables (q-values + decoy labels) are the other clean input (`--psms --q-col`).

## Output (`PXD055609_T5_ledger.tsv` / `.json`)

| class | T (targets) | D (decoys) | class-FDR `(D+1)/T` | f |
|---|---|---|---|---|
| canonical | 1133 | 8 | 0.79% | 21.2% |
| noncanonical (broad) | 3822 | 150 | 3.95% | 71.7% |
| variant | 379 | 2 | 0.79% | 7.1% |
| **global** | **5334** | **160** | **3.02%** | — |

Plus the non-novelty axis: **93.0%** (3315/3563) of accepted "noncanonical" target peptides are exact
substrings of canonical SwissProt. The `.json` additionally records `analysis_id`, software, unit, α, the FDR
convention, and the input SHA-256 for provenance.

## What this is — and is NOT (read before citing)

- **It demonstrates the tool and the feasibility claim.** `D_N` per class is recovered from a single pass over
  data already in PRIDE — no co-author, no raw-spectra re-search. That is the whole point of the reporting
  standard: the number that makes the claim checkable is essentially free.
- **The "noncanonical" row here is the BROAD accession class** — the entire 3-frame RNA-seq translation
  (`1546T*`), which is **93% canonical-coding** (the 93.0% self figure). It is **NOT** the genuinely-noncanonical
  (cryptic / ncORF) class the manuscript's bound is about. Splitting the decoys to the *genuinely*-noncanonical
  class requires the class map from a database rebuild (the natural next study); the deposit alone cannot do it.
  **So 3.95% is not "Raja's cryptic FDR"** — do not read it as a finding about that cohort's biology.
- **Computability, not calibration.** `(D+1)/T` is what the target-decoy record *implies*; whether it is
  well-calibrated for an inflated, often non-tryptic cryptic search space needs a *matched entrapment check*
  (manuscript Table 1, paired row). The ledger reports the computable estimate and says so.
- **This example is one deposited run** — the tool itself runs per-file. A separate pooled check across
  all 5 deposited PXD055609 samples (T1–T5, not just T5) confirms the broad noncanonical class-FDR is
  stable cohort-wide (3.99% pooled vs. 3.95% T5-only); the caveat above (broad, not
  genuinely-noncanonical, class) still applies at cohort scale too.

In short: the ledger is the reporting *vehicle*; the genuinely-noncanonical *number* is the open measurement.
The tool makes the vehicle frictionless today.

## PSM-replication-depth stratification (`--stratify-multiplicity`)

`--stratify-multiplicity` reports
class-FDR broken out by how many accepted PSMs support each peptide (n=1/n=2/n≥3). Pooled across all 5
PXD055609 samples, singleton (n=1) IDs carry a materially higher class-FDR than multi-PSM peptides **in every
class** (canonical 1.68%→0.29% n=1-to-n≥3; noncanonical 6.25%→0.66%) — a general target-decoy hygiene signal,
not cryptic-specific — while the noncanonical class stays elevated over canonical by a consistent ~2.3–3.7× at
every stratum, showing that gap isn't an artifact of cryptic peptides being disproportionately low-multiplicity
(they aren't — canonical and noncanonical have statistically indistinguishable PSM-multiplicity distributions).
`--unit unique-peptide` dedupes to one row per (class, peptide, decoy-status) as documented.

## Guarding against the ecological-inference mistake (`eco_diagnostic.py`)

On real data, trying to recover a class's true FDR by regressing pooled per-sample FDR on the class fraction
(Goodman-style ecological regression) fails here in a specifically *deceptive* way: the pooled FDR is
clamped near the target α by the procedure under audit, so the fit reports a tiny, confident-looking
standard error while being off by 2.5×+ from the truth. `src/darkproteome/eco_diagnostic.py` makes
this failure visible instead of silent — given several samples, it prints the directly-measured
per-class FDR (the trustworthy number), a constancy check across samples, and what the naive
regression would have said, flagging explicitly when the regression's SE is deceptively small rather
than honestly wide. Worked example: `PXD055609_eco_diagnostic.txt` (all 5 T1–T5 samples; reproduces
the "θ_C_hat=0.0297 vs measured 0.0112" deceptive-failure numbers exactly).

```bash
python3 src/darkproteome/eco_diagnostic.py data/external/pxd055609_pepxml/*.pepXML --alpha 0.03
```

## Interoperating with diagFDR

This is a companion to diagFDR (the general verifiable-FDR R package on CRAN [21]), not a rival — three bridges:
- **It reads diagFDR's input shape.** `--psms` consumes a PSM table with id / decoy-label / q-value / accession
  columns — the same per-PSM contract diagFDR expects.
- **It emits that contract.** `--emit-diagfdr out.diagfdr.tsv` writes the per-PSM table
  (`id, is_decoy, q, pep, run, score`), so a deposit available only as a **pepXML** — which diagFDR cannot read —
  can be passed *through* this tool into diagFDR. This file carries every PSM (target and decoy, accepted and
  rejected, matching the threshold sweep diagFDR's own diagnostics need) but no class label — it lets diagFDR
  diagnose the whole analysis's scope/calibration/stability, not a class-specific split.
- **It emits the class-stratified contract diagFDR needs for that split.** `--emit-diagfdr-by-class DIR` writes
  one `<class>.tsv` file per class (same six-column contract) plus `class_manifest.json`, reusing the exact
  per-PSM classification already computed for the main ledger — no new heuristic, and every class file's
  alpha-accepted subset reproduces the main ledger's own `T_class`/`D_class` exactly (verified, both on this
  fixture and live on PXD055609 T5). Load the resulting files as named universes —
  `dfdr_run_all(xs = list(canonical = ..., noncanonical = ...))` — for a per-class diagFDR diagnostic; diagFDR's
  own `dfdr_run_all()` requires a *named* list and already runs every diagnostic independently per name, so no
  diagFDR-side change is needed. One caveat: `dfdr_run_all()` also auto-generates a scope-disagreement (Jaccard)
  plot whenever more than one universe is supplied — for mutually exclusive class partitions (not diagFDR's own
  run-wise-vs-global use case) that plot is expected to show near-zero overlap by construction, not a real
  disagreement; read it accordingly.
