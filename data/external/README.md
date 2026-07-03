# data/external — large published inputs (NOT in git)

These files are **published, large (~358 MB), and gitignored**. They belong in this
directory (or wherever `$DARKPROTEOME_DATA` points); `src/darkproteome/paths.py`
resolves them.

**This README is the only file in this directory tracked by git** — so the
recovery instructions survive even though the data does not.

## What must be here

| path | what | size | re-download |
|---|---|---|---|
| `swissprot_human.fasta` | UniProt SwissProt human reviewed proteome (canonical reference for the non-novelty substring test) | 14 MB | UniProt: `https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta` |
| `atlases/IEAtlas_Epitopes_In_Cancer_Tissues.txt` + `..._In_Normal_Tissues.txt` (+ other IEAtlas/CPDB files) | IEAtlas non-canonical epitopes (cancer+normal) + CrypticProteinDB | ~69 MB | IEAtlas `http://bio-bigdata.hrbmu.edu.cn/IEAtlas/download.jsp` (needs browser UA); CPDB `https://www.maherlab.com/crypticproteindb-download` |
| `hla_ligand_atlas/hla_2020.12/HLA_aggregated.tsv` | HLA Ligand Atlas benign normal-tissue HLA-I ligandome (canonical-space normal reference) | 66 MB | `https://hla-ligand-atlas.org/rel/hla_2020.12.zip` (CC-BY 4.0) |
| `eupmc_extracted/tables/sciadv.ads7405_tables_s1_to_s4.xlsx` | Raja 2025 ovarian supplement (S1–S4) | <1 MB | `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11837991/supplementaryFiles` → unzip |
| `supp_tables.xlsx` | HCC 2024 (Camarena/Albà) supplement (S1–S26) | 76 MB | `https://ndownloader.figshare.com/files/44916451` (Figshare 24448723, CC-BY 4.0) |
| `gtex/GTEx_v8_gene_median_tpm.gct.gz` | GTEx v8 gene-median TPM (56,200 genes × 54 normal tissues) — measured normal-expression specificity floor | 7 MB | `https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz` |
| `gencode/gencode.v26.long_noncoding_RNAs.gtf.gz` | GENCODE v26 lncRNA GTF (ENST↔ENSG↔name; v26 == GTEx v8) — maps lncRNA-ORF antigens to their lncRNA gene's ENSG | 2.6 MB | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.long_noncoding_RNAs.gtf.gz` |
| `pxd055609_pepxml/T{1..5}-*.pepXML` | Raja 2025 ovarian cohort's 5 deposited MSFragger pepXML intermediates (decoys retained, class-labelable accessions) — input to `class_decoy_ledger.py` and `psm_multiplicity_probe.py` | ~190 MB | PRIDE FTP: `ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD055609/` (5 files, `T1-6472_*_082020.pepXML` through `T5-1546_*_101821.pepXML`; CC0) |

Full provenance, licenses, and table-level notes: **`../SOURCES.md`**.

## If this directory is empty/lost
1. Re-download each file above into this exact relative layout, **or**
2. set `DARKPROTEOME_DATA` to wherever they live.

Self-check after restoring: `python3 src/darkproteome/tier1_nonnovelty.py` must
print `5/2979`, `43/116 = 37.1%`, `0 mismatches`, `54.4%`.
