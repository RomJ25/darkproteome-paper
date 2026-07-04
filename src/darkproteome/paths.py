"""Centralized data-path resolution.

The large external inputs (SwissProt, IEAtlas, CrypticProteinDB, HLA Ligand Atlas, the
two cohort supplements) are PUBLISHED + LARGE, so they are not committed to git -- they
live in a repo-local, gitignored `data/external/`.

Resolution order (first that exists wins):
  1. $DARKPROTEOME_DATA            -- explicit override
  2. <repo>/data/external          -- the default

If `data/external/` is ever lost, re-download every file from `data/SOURCES.md`
(URLs + licenses) and `data/external/README.md`.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(_HERE))          # .../darkproteome


def _resolve_base():
    env = os.environ.get("DARKPROTEOME_DATA")
    if env and os.path.isdir(env):
        return env
    return os.path.join(REPO, "data", "external")


BASE = _resolve_base()

SPROT          = os.path.join(BASE, "swissprot_human.fasta")
ATLAS_DIR      = os.path.join(BASE, "atlases")
IEATLAS_CANCER = os.path.join(ATLAS_DIR, "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
IEATLAS_NORMAL = os.path.join(ATLAS_DIR, "IEAtlas_Epitopes_In_Normal_Tissues.txt")
HLALA          = os.path.join(BASE, "hla_ligand_atlas", "hla_2020.12", "HLA_aggregated.tsv")
RAJA_SI        = os.path.join(BASE, "eupmc_extracted", "tables", "sciadv.ads7405_tables_s1_to_s4.xlsx")
HCC_SI         = os.path.join(BASE, "supp_tables.xlsx")
GTEX_MEDIAN    = os.path.join(BASE, "gtex", "GTEx_v8_gene_median_tpm.gct.gz")  # genes x 54 tissues
GENCODE_LNC    = os.path.join(BASE, "gencode", "gencode.v26.long_noncoding_RNAs.gtf.gz")  # ENST<->ENSG (v26 = GTEx v8)
CPDB_IMMUNO    = os.path.join(ATLAS_DIR, "immunopeptides_cryptic.csv")
CPDB_EPITOPES  = os.path.join(ATLAS_DIR, "epitopes_cryptic.csv")


def require(*required):
    """Raise a clear, actionable error if any required input file is missing."""
    missing = [p for p in required if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            "Missing required data input(s):\n  " + "\n  ".join(missing) +
            f"\n\nResolved data base = {BASE}\n"
            "Populate data/external/ (see data/external/README.md + data/SOURCES.md) "
            "or set $DARKPROTEOME_DATA to where the inputs live.")
