# Data sources

Every URL / schema / count below was checked against the live source. This is the
executable manifest the ingestion layer targets.
**The T-cell validation labels are thinner and messier than the headline counts suggest — see
the "label reality" notes below for each cohort.**

---

## A. Ground-truth cohorts (the T-cell labels)

### HCC — Camarena, Albà et al., Sci Adv 2024, `10.1126/sciadv.adn3628`
- Full text (open): https://pmc.ncbi.nlm.nih.gov/articles/PMC11235171/
- Supplement (ALL tables S1–S26 in one xlsx, CC BY 4.0):
  - Figshare DOI: https://doi.org/10.6084/m9.figshare.24448723.v5
  - Direct: https://ndownloader.figshare.com/files/44916451  (`...SupplementaryTables_S1_S26.xlsx`, 76 MB)
  - Figshare API (bot-friendly): https://api.figshare.com/v2/articles/24448723
  - science.org / bioRxiv HTML return 403 to scripts; use Figshare API + PMC.
- Tables that matter: **S16** candidates (13 synthesized peptides), **S17** HLA-A*02:01
  binding assay, **S18** mouse ELISPOT immunogenicity, **S19** per-sample tumor-specificity
  + normal-adjacent FPKM, **S23** immunopeptidomics evidence for 1196 tumor-specific genes
  (gene-level), **S26** their own MHCquant MS (36 unique peptides, has `q-value`), **S2** HLA typing.
- **Label reality (important):** the T-cell validation is **mouse HHD-DR1, not human PBMC**,
  and tiny: 4 peptides ELISPOT-tested → **2 reactive** (`WMSLDWELYV`, `GLFHIYHKI`),
  **2 non-reactive** (`HLWHSATSL`, `FLTLQVHGA`). Plus **4 HLA non-binders** (binding-level
  negatives) in S17. No per-peptide HLA allele column anywhere (whole validation set is A*02:01 by design).
  Ribo-seq is a binary RibORF call (≥5 footprints & score ≥0.5), not a periodicity %.
- **Where HCC is actually rich:** the *audit* corpus — 1196 tumor-specific genes with
  immunopeptidomics evidence + source provenance (S23), normal-adjacent expression per
  tumor (S19), and a real FDR column on their own MS (S26). Use HCC for the AUDIT and the
  PRESENTATION tier, not as a T-cell benchmark.

### Ovarian — Raja et al., Sci Adv 2025, `10.1126/sciadv.ads7405`
- Full text (open, CC BY-NC): https://pmc.ncbi.nlm.nih.gov/articles/PMC11837991/
- Supplement (works, no auth):
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11837991/supplementaryFiles
  → zip → `sciadv.ads7405_tables_s1_to_s4.xlsx` (sheets S1–S4) + `sciadv.ads7405_sm.pdf`.
- Raw: PRIDE **PXD055609** (MS) · SRA **PRJNA1160863** (RNA).
- **Deposited pepXML (used by `class_decoy_ledger.py` / `psm_multiplicity_probe.py`):** 5 MSFragger
  pepXML intermediates, one per patient sample T1–T5, at
  `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD055609/` — retain decoys (`rev_` prefix)
  and class-labelable accessions (`sp|`=canonical, `1546T*`/etc.=RNA-seq-derived, `ENSP*-Mut`=variant).
  Kept locally at `data/external/pxd055609_pepxml/` (gitignored, ~190 MB; see
  `data/external/README.md`).
- Tables: **S2** all 311 cryptic peptides (Sample, Peptide, Gene Symbol, Biotype,
  Transcript, coords, Gen_Location/Region — no HLA/FDR/T-cell). **S3** the 38 tested
  candidates (Sample, Peptide-id, Sequence, Gene Symbol, %Rank EL, %Rank BA,
  Predicted HLA allele). **S1** per-patient 6-allele HLA typing. **Data S1** = 26-peptide
  MS spectral validation (NOT the reactive count — the two are easy to conflate).
- **Label reality (critical):** funnel **311 → 38 tested → ~70% reactive (~26) / ~12 non-reactive**
  is right on 311→38, but the **per-peptide reactive/non-reactive call is NOT in any table.**
  It lives only in **Fig 6B** (non-reactive = black dots) + Fig 6C / figs S5–S7 bar graphs.
  The 38 sequences + HLA + locus are machine-readable (join S3↔S2 by sequence); the
  positive/negative *label* requires **figure digitization**, cross-checked by EL/BA rank.
  Authoritative per-peptide values may need a request to the corresponding author of
  the Raja et al. 2025 paper above. Treat the ~26/~12 split as approximate until digitized.

**Net labels available:** ~28 human immunogenic positives (ovarian, figure-locked) + ~12 human
negatives + 2 mouse positives / 2 mouse + 4 binding negatives (HCC). This is the TESLA regime
(37 positives) — small but *normal* for this field.

---

## B. Background / integration resources (presentation tier + normal filtering)

| Resource | Access | License | Best use |
|---|---|---|---|
| **GENCODE Ribo-seq ORFs** | https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/ (phase1 7,264 / phase2 10,127–28,359); code MIT github.com/jorruior/gencode-riboseqORFs | EMBL open (cite Mudge 2022) | **ncORF catalog + Ribo-seq evidence (primary).** Carries HLA peptide-support tiers. |
| **nuORFdb v1.2** | https://proteomics.broadinstitute.org/nuORFdb/ (BED + 229k-entry FASTA + xlsx) | unstated (cite Ouspenskaia 2022) | ncORF catalog w/ ORF-class deflines (uORF/dORF/lncRNA-ORF/pseudogene). |
| **CrypticProteinDB** | https://www.maherlab.com/crypticproteindb-download (5 CSVs) | article CC BY; files unstated | MS-evidenced cryptic proteins + **epitopes with HLA + score** (14 cancers). |
| **IEAtlas** | http://bio-bigdata.hrbmu.edu.cn/IEAtlas/download.jsp (TSV; needs browser UA) | CC BY 4.0 | **Observed non-canonical HLA epitopes, cancer + normal** (245k; normal-tissue file = direct subtraction). |
| **HLA Ligand Atlas** | https://hla-ligand-atlas.org/rel/hla_2020.12.zip | **CC-BY 4.0 (cleanest)** | **Benign normal-tissue HLA ligandome** = self-peptide negatives / specificity filter. Match by exact peptide seq. |
| **Recount3** | Bioconductor `recount3`; AWS open data | Artistic-2.0 / open | GTEx + huge normal RNA background; quantify any ncORF locus (bigWig). |
| **GTEx v8 gene-median TPM** | https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz (~7 MB) | GTEx open-access (cite GTEx Consortium, *Science* 2020) | **Measured normal-tissue expression** (56,200 genes × 54 tissues). The specificity floor: pseudogene PARENT expression (`gtex_specificity.py` → 43/43) + lncRNA-ORF / altORF source expression (`gtex_class_specificity.py`, `lncrna_ensg_specificity.py`). On disk: `data/external/gtex/`. |
| **GENCODE v26 lncRNA GTF** | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.long_noncoding_RNAs.gtf.gz (~2.6 MB) | EMBL open (cite Frankish 2019) | **ENST↔ENSG↔gene_name for lncRNA genes** (v26 == GTEx v8). Maps lncRNA-ORF antigens to their lncRNA gene's ENSG → GTEx (coverage 32%→84%); lncRNA-only ⇒ no coding-neighbour contamination. On disk: `data/external/gencode/`. |
| **PRIDE / MassIVE** | https://www.ebi.ac.uk/pride/ws/archive/v2/ ; massive.ucsd.edu | CC0 default | Raw HLA-elution to re-search vs a custom ncORF DB (must reprocess; deposits used canonical DBs). |
| caAtlas | zhang-lab.org/caatlas | **CC BY-NC-ND** (watch for commercial) | Canonical/PTM only — comparator, not dark-proteome. |

No single resource gives catalog + Ribo-seq + raw HLA + normal background. Stack:
GENCODE+nuORFdb (catalog/translation) · PRIDE-reprocess + CrypticProteinDB + IEAtlas (presentation)
· HLA Ligand Atlas + GTEx-median (+ Recount3) + IEAtlas-normal (specificity).
