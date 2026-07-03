"""Scale the audit corpus from 2 cohorts to the major public cryptic-antigen atlases.

Ingests IEAtlas (cancer non-canonical epitopes) + CrypticProteinDB (cryptic immunopeptides
and allele-resolved cryptic epitopes), mapped as-reported into the Claim Catalog, and
appends them to the cohort corpus -> data/claim_catalog_scaled.csv.

Purpose: show the reporting/reusability gap generalises from 2 flagship papers (see
probe_report.py) to hundreds of thousands of published claims, not a 2-paper quirk. Also
computes the IEAtlas cancer-vs-normal presentation overlap as a bonus stat.

    python3 src/darkproteome/ingest_atlases.py
"""
import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import COLUMNS, NOT_REPORTED  # noqa: E402

import paths  # noqa: E402  -- centralized data paths

NR = NOT_REPORTED
ATLAS = paths.ATLAS_DIR
REAL = os.path.join(paths.REPO, "data", "claim_catalog_real.csv")
OUT = os.path.join(paths.REPO, "data", "claim_catalog_scaled.csv")

IEATLAS_DOI = "10.1093/nar/gkac776"
CPDB_DOI = "10.1093/narcan/zcad024"
FIELDS = COLUMNS + ["_canonical", "_source"]

csv.field_size_limit(10_000_000)


def base(**kw):
    r = {c: NR for c in COLUMNS}
    r["meets_consensus_bar_as_reported"] = "insufficient-info"
    r["source_provenance"] = "primary"
    r["extraction_confidence"] = "med"   # atlas bulk-extraction
    r["_canonical"] = "no"
    r.update(kw)
    return r


def atlas_class(t):
    t = (t or "").strip().lower()
    if "utr" in t and "3" in t:
        return "3'UTR-ORF"
    if "utr" in t and "5" in t:
        return "uORF"
    if "frame" in t:
        return "altORF"
    if "pseudogene" in t:
        return "pseudogene-ORF"
    if "intron" in t:
        return "intronic-ORF"
    if t in ("ncrna", "noncoding_rna", "lncrna", "lincrna", "antisense"):
        return "lncRNA-ORF"
    return "other"


def ingest_ieatlas(rows_out):
    path = f"{ATLAS}/IEAtlas_Epitopes_In_Cancer_Tissues.txt"
    n = 0
    cancer_seqs = set()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)  # header Sequence/Length/ORF_ID/tissue
        for r in rd:
            if len(r) < 4 or not r[0]:
                continue
            seq = r[0].strip().upper()
            if not seq.isalpha():
                continue
            cancer_seqs.add(seq)
            try:
                ln = int(float(r[1]))
            except (ValueError, TypeError):
                ln = len(seq)
            rows_out.append(base(
                peptide_sequence=seq,
                orf_id_or_locus=r[2] or NR,
                orf_class="other",                  # IEAtlas = non-coding-derived; class not per-row
                antigen_type="TSA",
                cancer_type=r[3] or NR,
                evidence_types="immunopeptidomics",  # MS-observed HLA-presented
                min_peptide_len=str(ln),
                tumor_specificity_basis=NR,          # IEAtlas does not report specificity per epitope
                validation_level="MS-presented",
                underlying_dataset_accession="IEAtlas",
                citation_doi_pmid=IEATLAS_DOI,
                citation_location="Epitopes In Cancer Tissues",
                _source="IEAtlas",
            ))
            n += 1
    return n, cancer_seqs


def normal_overlap(cancer_seqs):
    path = f"{ATLAS}/IEAtlas_Epitopes_In_Normal_Tissues.txt"
    normal = set()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if r and r[0]:
                normal.add(r[0].strip().upper())
    both = cancer_seqs & normal
    return len(normal), len(both)


def ingest_cpdb_immuno(rows_out):
    path = f"{ATLAS}/immunopeptides_cryptic.csv"
    n = 0
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            seq = (r.get("Peptide") or "").strip().upper()
            if not seq or not seq.isalpha():
                continue
            rows_out.append(base(
                peptide_sequence=seq,
                orf_id_or_locus=(r.get("gene.symbol") or NR),
                orf_class=atlas_class(r.get("type")),
                antigen_type="TSA",
                cancer_type=NR,                      # not in this table
                evidence_types="immunopeptidomics",
                min_peptide_len=(r.get("Peptide length") or str(len(seq))),
                tumor_specificity_basis=NR,
                validation_level="MS-presented",
                underlying_dataset_accession="CrypticProteinDB",
                citation_doi_pmid=CPDB_DOI,
                citation_location="immunopeptides_cryptic.csv",
                _source="CrypticProteinDB-immuno",
            ))
            n += 1
    return n


def ingest_cpdb_epitopes(rows_out):
    """epitopes_cryptic.csv has an off-by-one header; parse by position (17 fields)."""
    path = f"{ATLAS}/epitopes_cryptic.csv"
    n = 0
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh)
        next(rd, None)  # malformed header
        for r in rd:
            if len(r) < 17 or not r[0]:
                continue
            seq = r[0].strip().upper()
            if not seq.isalpha():
                continue
            allele = (r[2] or "").strip()
            gtex_low = (r[10] or "").strip().lower() == "yes"   # Expression <1 TPM in GTEx
            try:
                ln = int(float(r[6]))
            except (ValueError, TypeError):
                ln = len(seq)
            rows_out.append(base(
                peptide_sequence=seq,
                orf_id_or_locus=(r[8] or NR),
                orf_class=atlas_class(r[16]),
                antigen_type="TSA",
                cancer_type=(r[1] or NR),
                hla_allele=allele or NR,
                evidence_types="immunopeptidomics, MHC binding prediction",
                min_peptide_len=str(ln),
                tumor_specificity_basis="broad-normal-panel" if gtex_low else "not stated",
                validation_level="MS-presented",
                underlying_dataset_accession="CrypticProteinDB",
                citation_doi_pmid=CPDB_DOI,
                citation_location="epitopes_cryptic.csv",
                _source="CrypticProteinDB-epitopes",
            ))
            n += 1
    return n


def main():
    rows = []
    # carry forward the cohort corpus, tagging its source
    with open(REAL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            r.setdefault("_source", "cohort:" + ("HCC" if "Hepato" in r.get("cancer_type", "")
                                                 else "ovarian"))
            rows.append(r)
    n_cohort = len(rows)

    n_ie, cancer_seqs = ingest_ieatlas(rows)
    n_cp1 = ingest_cpdb_immuno(rows)
    n_cp2 = ingest_cpdb_epitopes(rows)
    n_norm, n_both = normal_overlap(cancer_seqs)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, NR) for k in FIELDS})

    uniq = len({r["peptide_sequence"] for r in rows if r.get("_canonical") != "yes"})
    print(f"cohorts carried:            {n_cohort}")
    print(f"+ IEAtlas cancer epitopes:  {n_ie}")
    print(f"+ CrypticProteinDB immuno:  {n_cp1}")
    print(f"+ CrypticProteinDB epitope: {n_cp2}")
    print(f"= TOTAL claims:             {len(rows)}  (unique non-canonical peptides: {uniq})")
    print(f"by source: {dict(Counter(r.get('_source') for r in rows))}")
    print(f"wrote -> {OUT}")
    print("\nBONUS — IEAtlas presentation specificity:")
    print(f"  IEAtlas normal-tissue epitopes: {n_norm}")
    print(f"  cancer epitopes ALSO presented in normal tissue: {n_both} "
          f"({100*n_both/max(1,len(cancer_seqs)):.1f}% of {len(cancer_seqs)} unique cancer epitopes)")


if __name__ == "__main__":
    main()
