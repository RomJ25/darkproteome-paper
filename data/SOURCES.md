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
| **UniProt/Swiss-Prot reviewed human proteome** (`swissprot_human.fasta`, the canonical reference *R* behind the manuscript's headline 56.3% figure) | `https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta` (live REST query, unversioned) | UniProt open | **The canonical reference** for every overlap statistic in `manuscript_v2.md`. 20,431 entries, 11,418,104 residues as of the last fetch; since the query has no pinned release string, `reference_provenance.py`'s printed SHA-256 hash is the reproducibility anchor, not a dated release. The era-correct comparator (Swiss-Prot 2022_01, 566,996 entries / 20,376 human) is separately pinned as `data/external/uniprot_sprot.fasta.gz`. |
| **GENCODE Ribo-seq ORFs** | https://ftp.ebi.ac.uk/pub/databases/gencode/riboseq_orfs/ (phase1 7,264 / phase2 10,127–28,359); code MIT github.com/jorruior/gencode-riboseqORFs | EMBL open (cite Mudge 2022) | **ncORF catalog + Ribo-seq evidence (primary).** Carries HLA peptide-support tiers. |
| **nuORFdb v1.2** | https://proteomics.broadinstitute.org/nuORFdb/ (BED + 229k-entry FASTA + xlsx) | unstated (cite Ouspenskaia 2022) | ncORF catalog w/ ORF-class deflines (uORF/dORF/lncRNA-ORF/pseudogene). |
| **Translnc** (`lncRNA_peptide_AA_seq.fasta`, IEAtlas's second integrated ncORF source; cite Lv/Chang et al. 2022, *Nucleic Acids Research* 50(D1):D413, doi:10.1093/nar/gkab928) | http://bio-bigdata.hrbmu.edu.cn/TransLnc/download.jsp (same host/group as IEAtlas) | not separately stated on the download page | Latent-canonical-ambiguity measurement for R2's headline library-union rate. **Multi-species by design** (its own paper title: "in multiple species") — 583,840 headers mix human (ALL-CAPS gene-symbol convention) and mouse (mixed-case, e.g. `Tug1-...`) entries. Species-filtered to human-only via `looks_human_translnc()` (duplicated in `abundance_bias.py`/`abundance_direct.py`/`library_union.py`), which also excludes two confirmed non-human GenBank WGS-scaffold prefixes (`AABR07*` = rat, `CAAA01*` = mouse — verified against their NCBI nuccore records) that the naive all-caps heuristic alone would have wrongly kept; a third WGS-style prefix, `AUXG01*`, is genuinely human and is correctly kept. |
| **CrypticProteinDB** | https://www.maherlab.com/crypticproteindb-download (5 CSVs) | article CC BY; files unstated | MS-evidenced cryptic proteins + **epitopes with HLA + score** (14 cancers). |
| **IEAtlas** | http://bio-bigdata.hrbmu.edu.cn/IEAtlas/download.jsp (TSV; needs browser UA) | CC BY 4.0 | **Observed non-canonical HLA epitopes, cancer + normal** (245k; normal-tissue file = direct subtraction). |
| **HLA Ligand Atlas** | https://hla-ligand-atlas.org/rel/hla_2020.12.zip | **CC-BY 4.0 (cleanest)** | **Benign normal-tissue HLA ligandome** = self-peptide negatives / specificity filter. Match by exact peptide seq. |
| **Recount3** | Bioconductor `recount3`; AWS open data | Artistic-2.0 / open | GTEx + huge normal RNA background; quantify any ncORF locus (bigWig). |
| **GTEx v8 gene-median TPM** | https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz (~7 MB) | GTEx open-access (cite GTEx Consortium, *Science* 2020) | **Measured normal-tissue expression** (56,200 genes × 54 tissues). The specificity floor: pseudogene PARENT expression (`gtex_specificity.py` → 43/43) + lncRNA-ORF / altORF source expression (`gtex_class_specificity.py`, `lncrna_ensg_specificity.py`). On disk: `data/external/gtex/`. |
| **GENCODE v26 lncRNA GTF** | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.long_noncoding_RNAs.gtf.gz (~2.6 MB) | EMBL open (cite Frankish 2019) | **ENST↔ENSG↔gene_name for lncRNA genes** (v26 == GTEx v8). Maps lncRNA-ORF antigens to their lncRNA gene's ENSG → GTEx (coverage 32%→84%); lncRNA-only ⇒ no coding-neighbour contamination. On disk: `data/external/gencode/`. |
| **PRIDE / MassIVE** | https://www.ebi.ac.uk/pride/ws/archive/v2/ ; massive.ucsd.edu | CC0 default | Raw HLA-elution to re-search vs a custom ncORF DB (must reprocess; deposits used canonical DBs). |
| caAtlas | zhang-lab.org/caatlas | **CC BY-NC-ND** (watch for commercial) | Canonical/PTM only — comparator, not dark-proteome. |
| **COD-dipp** (Bedran et al. 2023, *Cancer Immunol Res* 11(6):747–62, `10.1158/2326-6066.CIR-22-0621`) | github.com/immuno-informatics/COD-dipp; data at Figshare `10.6084/m9.figshare.16538097`; 26-study MS reanalysis (772 samples, 11 cancers) | CC BY 4.0 | Independent ncMAP atlas + tumor-selectivity funnel — see note below, don't skip it. |
| **ImmunoVerse** (Li, Guzmán-Bringas, et al.; corr. Yarmarkovich; NYU), "A pan-cancer atlas of therapeutic T cell targets," bioRxiv `10.1101/2025.01.22.634237`, PMC12265682 | Supp. Table S7 direct: `https://www.biorxiv.org/content/biorxiv/early/2025/07/07/2025.01.22.634237/DC7/embed/media-7.xlsx?download=true`; code (no peptide-level data) at github.com/frankligy/pan_cancer_intracellular_antigen_atlas | bioRxiv default (preprint, not yet peer-reviewed as of this writing — recheck publication status before citing) | 17,741-row nuORF-derived peptide catalogue ("ORF_antigen" sheet) with per-peptide ORF-class label (Pseudogene/lncRNA/uORF/dORF/Out-of-Frame/Other). Used as a third, independently-computed cross-study data point (this project's own analysis code, not shipped in this release) (0.1% canonical-substring overlap overall, 0.5% pseudogene-ORF — traced to inheriting Ouspenskaia et al. 2022's already-curated ORF coordinates, not a canonical-exclusion step of its own). |
| **Deutsch, Prensner, et al. 2026** (TransCODE/GENCODE/PeptideAtlas/HUPO-HPP/HUPO-HIPP consortium), "Expanding the human proteome with microproteins and peptideins," *Nature* 654:813–825, `10.1038/s41586-026-10459-x`, PMC13275300 | Full text (paywalled on nature.com) via Cambridge repository `https://www.repository.cam.ac.uk/items/33098e00-9718-4bab-af50-6fac63593db4`; Supp. Info `https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-026-10459-x/MediaObjects/41586_2026_10459_MOESM1_ESM.pdf`; companion preprints bioRxiv `10.1101/2024.09.09.612016` (methodology detail) and `10.1101/2025.02.19.639069` (community assessment); peptide/ORF-level bulk data via PeptideAtlas build 568, `https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptides?atlas_build_id=568&redundancy_constraint=4&biosequence_name_constraint=CONTRIB_GENCODE%25&apply_action=QUERY&output_mode=tsv` (no login) | CC BY-NC-ND 4.0 (published article, per Cambridge repository deposit); PeptideAtlas bulk export open | 2026-05-06, the field's most current classification-standardization effort for ncORF/microprotein evidence, ~half the paper devoted to HLA immunopeptidomics. Explicit own canonical-exclusion criterion ("peptides ... mapping to ncORFs and not canonical proteins, with at most 10 distinct mappings" — bioRxiv 612016 Methods). Used as a fourth, independently-computed cross-study data point (analysis code not shipped in this release): 1.7% canonical-substring overlap (53/3,061) after approximating the paper's own filters — squarely in the "confirmed mechanism, low residual" bucket alongside Ouspenskaia/Chong/CrypticProteinDB/Raja, not the IEAtlas/HCC outlier bucket. **IEAtlas, nuORFdb, OpenProt, sORFs.org, CrypticProteinDB and ImmunoVerse are not mentioned anywhere in the paper, its Supplement, or either companion preprint** (verified by direct full-text grep) — this project's IEAtlas/HCC finding is not redundant with this consortium's own current work. States a verbatim "7-point research agenda" of open questions (bioRxiv 612016, Box 1); question 2, "Should HLA immunopeptidomics be used as evidence that a ncORF encodes a protein-coding gene?", is directly relevant to this project's core finding. Prensner JR (this project's own planned domain-read gate) is a co-author, and ncORF/microprotein annotation is a sustained, current line of his research (also a 2025 *Trends in Genetics* review "Microproteins in cancer" and co-authorship on the companion Wacholder et al. *Nat Commun* 2026 community-assessment paper). No dedicated public channel exists for submitting evidence against the stated research agenda — only generic feedback forms (PeptideAtlas `feedback.php`, GENCODE contact page). |
| **SysteMHC Atlas v2.0** (Huang, Gan, Cui, Lan, Liu, Caron, Shao), *Nucleic Acids Research* 2024, 52(D1):D1062–D1071, `10.1093/nar/gkad1068` | **Trusted host only: `https://systemhc.sjtu.edu.cn`** — verified against the paper's own Data Availability statement on `academic.oup.com` (an unaffected host). **`systemhcatlas.org` (a different, similarly-named domain that surfaced during research) 302-redirects to an unrelated non-academic domain (`survey-smiles.com`) and must never be used for anything** — reconnaissance, collection, or casual browsing. The two hosts' `robots.txt` also differ materially: `systemhcatlas.org` has `Disallow: /` (hard no); `systemhc.sjtu.edu.cn` has no `robots.txt` at all. General search: `https://systemhc.sjtu.edu.cn/explore` (GET params `protein_id`, `search_hit`, `top_allele`, `PTM`, `mhc_class`); dedicated `/Non-UniProt` page is broken (freezes the browser tab on every load, confirmed twice — do not attempt it again) | CC BY-NC (per the publishing article's copyright notice; the site itself states no separate data license) | Its own paper defines a **"non-UniProt peptides"** category (78,959 total, 4,442 binders, across 2,447 samples/7,190 MS files/~303 allotypes) as peptides matched only via each original study's own customized search database, not UniProt — explicitly **mixed canonical + non-canonical with no described cross-check** ("non-UniProt canonical peptides may represent tumor-specific mutated peptides, whereas non-UniProt non-canonical peptides may represent tumor-specific non-mutated peptides from non-coding regions... or other non-canonical genomic sources"), and states finer classification awaits **"a future version of the atlas"** (both quoted verbatim in `systemhc_check.py`'s docstring) — the same self-flagged gap already found in IEAtlas/HCC. No bulk peptide export exists (`/download` offers only ~300 allele-specific spectral libraries); the general `/explore` search has no `Protein_Database` column, but its `protein_id` field's **format** is a usable proxy: `sp\|ACCESSION\|NAME_ORGANISM` = UniProt-matched, `ucNNNxxx.V` (a UCSC Known Gene transcript accession) = non-UniProt-matched (confirmed by tracing which underlying dataset each format traces to, via `/datasets`). Used as a fifth, independently-computed cross-study data point (analysis code not shipped in this release): **98.9% canonical-substring overlap (442/447)** on a small, pre-registered, bounded sample (8 common HLA-I alleles + 1 allele-blind pull, one query each, ~0.57% coverage of the 78,959 total) — even higher than IEAtlas/HCC, but concentrated in just 2 underlying datasets, one of which (`SYSMHC00027` = Sarkizova et al. 2020, PMID 31844290, a **healthy mono-allelic reference panel**, not a tumor cohort) plausibly reflects a benign accession-namespace gap (canonical proteins simply best-matched to a UCSC transcript instead of a UniProt one in that study's reprocessing) rather than the dark-proteome mechanism this project audits elsewhere — reported plainly with that caveat, not treated as a clean sixth "confirmed instance." Full collection log, the classification regex, and known dead ends (non-deterministic query ordering, no dataset-ID filter in the search form) are documented in the script's own docstring. **This healthy-panel speculation was then tested directly** in a follow-up analysis (same day, same raw data, zero new web queries — the TSVs already carry the SysteMHC dataset ID per row): SYSMHC00027 (healthy panel) 98.7% (313/317, Wilson 95% CI 96.8–99.5%) vs SYSMHC00057 (melanoma tumor cohort, PMID 34391888) 99.2% (130/131, Wilson 95% CI 95.8–99.9%) — comparable, CIs overlapping almost completely, the tumor cohort's rate numerically not lower. **This undermines rather than supports the healthy-panel-confound explanation**: whatever produces the high overlap is at least as present in genuine tumor-derived data, not confined to the non-tumor reference panel — the more interesting, less comforting outcome. Still a small sample (n=131 for the tumor cohort, drawn incidentally from an allele-stratified pull, not a dataset-targeted design), not treated as a confirmed sixth instance of the IEAtlas/HCC pattern. |

No single resource gives catalog + Ribo-seq + raw HLA + normal background. Stack:
GENCODE+nuORFdb (catalog/translation) · PRIDE-reprocess + CrypticProteinDB + IEAtlas (presentation)
· HLA Ligand Atlas + GTEx-median (+ Recount3) + IEAtlas-normal (specificity).

**COD-dipp is worth reading closely, not just cataloguing (checked):**
- **The preprint's headline numbers did not survive peer review — don't cite them.** The 2022
  bioRxiv abstract advertised "140,966 immune-visible genomic regions" and a "7.8×"
  immunogenicity claim; the *published* 2023 version reports **8,601 ncMAPs (1.7% of 516,382
  peptides)** and drops the 7.8× framing (the real HLA-supertype finding is "A03 shows 5%
  noncanonical presentation vs. ~1% average" — a narrower, different claim). Always cite the
  published numbers.
- **Independent corroboration of our headline pattern.** A 3-step tumor-selectivity funnel
  (exclude ncMAPs MS-detected in a normal panel → GTEx v8 parent-gene TPM <1 across 29 tissues,
  testis excluded → zero Human Protein Atlas protein detection) collapses their 8,601 ncMAPs to
  **17 stringently cancer-selective candidates** (Table 1) — a ~500-fold funnel, run
  independently by a different group on different data, landing on the same "the big number
  doesn't survive a real specificity check" shape as our own audit.
- **A real caveat in general, but checked live against our own floor — doesn't apply here.**
  Their Table 1 lists two histone-gene ncMAPs (`HIST1H4L`, `HIST1H2BB`) with low GTEx TPM
  (0.04, 0.26) that *fail* to qualify as cancer-selective once checked against the Human Protein
  Atlas — protein detected in 43–44/56 healthy tissues despite the low transcript signal. That's
  documented mRNA/protein decoupling, and our `gtex_specificity.py` / `gtex_class_specificity.py`
  floor is GTEx-TPM-only. **Checked: it doesn't create a gap in our results.** The
  43/43 pseudogene-parent floor is saturated at the ceiling already — every one of the 43
  peptides clears >=1 TPM in all 54 GTEx tissues, most by orders of magnitude (one, `REIQTAVRL`,
  matches the H2BC/HIST1H2B family itself — `H2BC3`≡`HIST1H2BB` — at max 200.2 TPM across all
  54 tissues via the existing max-over-matched-genes logic), so there is no currently-"passes as
  specific" claim in that floor for a protein check to flip. The altORF and lncRNA-ORF class
  results (`gtex_class_specificity.py`) already stop at "source transcribed in normal tissue,
  not proof of normal presentation" — they never claim RNA absence proves protein absence, so
  they don't inherit this gap either. No pipeline change needed; this stays a documented
  external caveat, not an open action item.
- No T-cell/immunogenicity validation anywhere in the paper — their own Discussion states
  "identified ncMAPs require further validation... immunogenicity prediction is still in its
  infancy." Same presented-≠-validated gap our audit already flags.


---

## E. Manuscript reference-list verification notes

The bibliography itself lists entries plainly, without dagger marks or inline asides; the
verification methodology and citation-context notes live here instead.

- **Verified directly against fetched full-text XML** held in this repository
  (`data/external/fulltext/`): Bedran et al. 2023 (*Cancer Immunol Res*), Cai et al. 2023 (IEAtlas,
  *Nucleic Acids Research*), Chong et al. 2020 (*Nature Communications*), Ouspenskaia et al. 2022
  (*Nature Biotechnology*), Raja et al. 2025 (*Science Advances*).
- The remaining references were independently confirmed against publisher/PubMed records during
  this revision, or are corroborated by this project's own prior citation-verification records
  (see `scripts/CITATION_ATTRIBUTION_AUDIT_2026-07-22.md`).
- **Erhard et al. 2018** (*Nature Methods*) is cited only via Bedran et al. 2023's re-analysis of
  its published non-canonical-peptide canonical-overlap rate, not for a figure taken directly
  from it.
- **Nesvizhskii & Aebersold 2005** (*Molecular & Cellular Proteomics*) is cited as the standard
  reference for shared-peptide ambiguity in shotgun proteomics.
