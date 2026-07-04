"""Ingest the two landmark cohort supplements into the Claim Catalog schema.

Reads the REAL published supplementary tables (HCC Camarena/Albà 2024 + ovarian Raja 2025)
and emits one row per published dark-proteome tumor-antigen *claim*, mapped faithfully to
the 19-column contract. Everything the paper does not machine-report becomes `not reported`
(a finding, not a blank).

Depends on openpyxl (ingestion only; the core auditor stays stdlib). Inputs are the
downloaded supplements (see data/SOURCES.md for provenance/URLs).

    python3 src/darkproteome/ingest_cohorts.py [HCC.xlsx] [OVARIAN.xlsx] > /dev/null
"""

import csv
import os
import sys

import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import COLUMNS, NOT_REPORTED, validate_row  # noqa: E402

import paths  # noqa: E402  -- centralized data paths

NR = NOT_REPORTED
HCC_DEFAULT = paths.HCC_SI
OVA_DEFAULT = paths.RAJA_SI
OUT = os.path.join(paths.REPO, "data", "claim_catalog_real.csv")

HCC_DOI = "10.1126/sciadv.adn3628"
OVA_DOI = "10.1126/sciadv.ads7405"

# Verified from HCC S18 (mouse HHD-DR1 IFN-γ ELISPOT) by direct extraction.
HCC_REACTIVE = {"WMSLDWELYV", "GLFHIYHKI"}          # in-vivo (humanized mouse) reactive
HCC_TCELL_NEG = {"HLWHSATSL", "FLTLQVHGA"}          # assayed, NOT reactive -> validated_negative


def blank_row():
    r = {c: NR for c in COLUMNS}
    r["meets_consensus_bar_as_reported"] = "insufficient-info"
    r["source_provenance"] = "primary"
    r["extraction_confidence"] = "high"
    return r


def map_class(biotype, gen_location=None):
    g = (biotype or "").strip().lower()
    loc = (gen_location or "").strip().lower()
    if "pseudogene" in g:
        return "pseudogene-ORF", False
    if g in ("lncrna", "lincrna", "antisense", "bidirectional_promoter_lncrna",
             "sense_intronic", "processed_transcript", "macro_lncrna",
             "3prime_overlapping_ncrna", "sense_overlapping"):
        return "lncRNA-ORF", False
    if g in ("protein_coding", "protein-coding"):
        if loc == "genebody":          # ovarian: cryptic peptide off-frame in a coding gene body
            return "altORF", False
        return "other", True           # HCC: a canonical cancer-testis antigen (MAGE/SSX...) = control
    return "other", False


def peptides_from_cell(cell):
    if not cell:
        return []
    out = []
    for tok in str(cell).replace(";", ",").split(","):
        s = tok.strip().upper()
        if s and s.isalpha() and len(s) >= 7:
            out.append(s)
    return out


def build_gene_type_map(wb):
    """gene_name -> gene_type, harvested from HCC S19 (richest gene_type coverage)."""
    m = {}
    if "S19" not in wb.sheetnames:
        return m
    ws = wb["S19"]
    rows = ws.iter_rows(values_only=True)
    next(rows, None)  # header
    for r in rows:
        if not r:
            continue
        name, gtype = r[2], r[3]
        if name and gtype and name not in m:
            m[str(name)] = str(gtype)
    return m


def ingest_hcc(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    gtmap = build_gene_type_map(wb)
    rows = []

    # --- S23: gene-level immunopeptidomics presentation claims (the bulk corpus) ---
    ws = wb["S23"]
    it = ws.iter_rows(values_only=True)
    next(it, None)  # header
    for r in it:
        if not r or not r[0]:
            continue
        gene_id, _tx, gene_name, gene_type, n, _pct, source, pepcell, _hi = (list(r) + [None] * 9)[:9]
        peps = peptides_from_cell(pepcell)
        cls, canonical = map_class(gene_type)
        row = blank_row()
        row["peptide_sequence"] = peps[0] if peps else NR
        row["orf_id_or_locus"] = str(gene_name or gene_id or NR)
        row["orf_class"] = cls
        row["antigen_type"] = "TAA" if canonical else "TSA"
        row["cancer_type"] = "Hepatocellular Carcinoma"
        row["evidence_types"] = "immunopeptidomics"
        row["n_unique_peptides"] = str(len(peps)) if peps else NR
        row["min_peptide_len"] = str(min(len(p) for p in peps)) if peps else NR
        row["tumor_specificity_basis"] = "broad-normal-panel"
        row["validation_level"] = "MS-presented"
        row["citation_doi_pmid"] = HCC_DOI
        row["citation_location"] = "Table S23"
        row["_canonical"] = "yes" if canonical else "no"
        rows.append(row)

    # --- S16/S17/S18: the 13 experimentally-tested peptides (immunogenicity seed) ---
    ws = wb["S16"]
    it = ws.iter_rows(values_only=True)
    next(it, None)  # header
    for r in it:
        if not r or not r[0]:
            continue
        pep, _gid, gene_name, _orf, _pos, _eid, _rank, _aff, binding, tumorspec, _np, riborf = (list(r) + [None] * 12)[:12]
        pep = str(pep).strip().upper()
        cls, canonical = map_class(gtmap.get(str(gene_name), ""))
        if pep in HCC_REACTIVE:
            val = "in-vivo"
        elif pep in HCC_TCELL_NEG:
            val = "validated_negative"
        else:
            val = "MS-presented"
        ev = "immunopeptidomics, MHC binding assay, ELISPOT"
        if str(riborf).strip().upper() == "YES":
            ev = "Ribo-seq, " + ev
        row = blank_row()
        row["peptide_sequence"] = pep
        row["orf_id_or_locus"] = str(gene_name or NR)
        row["orf_class"] = cls
        row["antigen_type"] = "TSA"
        row["cancer_type"] = "Hepatocellular Carcinoma"
        row["hla_allele"] = "HLA-A*02:01"
        row["evidence_types"] = ev
        row["min_peptide_len"] = str(len(pep))
        row["tumor_specificity_basis"] = "broad-normal-panel" if str(tumorspec).strip().lower() == "yes" else "not stated"
        row["validation_level"] = val
        row["citation_doi_pmid"] = HCC_DOI
        row["citation_location"] = "Tables S16/S17/S18"
        row["_canonical"] = "no"
        rows.append(row)

    wb.close()
    return rows


def _find_header(ws, needles):
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        vals = [str(c).strip().lower() if c is not None else "" for c in r]
        if all(any(n in v for v in vals) for n in needles):
            return i, vals
    return None, None


def ingest_ovarian(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = []
    pep_biotype = {}

    # --- S2: all identified cryptic peptides (presentation claims) ---
    ws = wb["Supplementary Table S2"]
    all_rows = list(ws.iter_rows(values_only=True))
    hdr = all_rows[0]
    cols = {str(c).strip(): j for j, c in enumerate(hdr) if c}
    ci_pep = cols.get("Peptide")
    ci_bio = cols.get("Biotype")
    ci_loc = cols.get("Gen_Location")
    ci_gene = cols.get("Gene Symbol")
    s2_peps = []
    for r in all_rows[1:]:
        if not r or ci_pep is None or not r[ci_pep]:
            continue
        pep = str(r[ci_pep]).strip().upper()
        biotype = r[ci_bio] if ci_bio is not None else None
        genloc = r[ci_loc] if ci_loc is not None else None
        gene = r[ci_gene] if ci_gene is not None else None
        pep_biotype[pep] = (biotype, genloc)
        s2_peps.append((pep, biotype, genloc, gene))

    # --- S3: the 38 prioritized + immunogenicity-tested candidates ---
    ws3 = wb["Supplementary Table S3"]
    hi, hvals = _find_header(ws3, ["sequence", "predicted hla", "gene symbol"])
    s3_set = set()
    if hi is not None:
        idx = {name: j for j, name in enumerate(hvals)}
        def col(part):
            for name, j in idx.items():
                if part in name:
                    return j
            return None
        j_seq, j_gene, j_hla = col("sequence"), col("gene symbol"), col("predicted hla")
        for r in list(ws3.iter_rows(values_only=True))[hi + 1:]:
            if not r or j_seq is None or not r[j_seq]:
                continue
            pep = str(r[j_seq]).strip().upper()
            if not pep.isalpha():
                continue
            s3_set.add(pep)
            biotype, genloc = pep_biotype.get(pep, (None, None))
            cls, _ = map_class(biotype, genloc)
            row = blank_row()
            row["peptide_sequence"] = pep
            row["orf_id_or_locus"] = str(r[j_gene]) if j_gene is not None and r[j_gene] else NR
            row["orf_class"] = cls
            row["antigen_type"] = "TSA"
            row["cancer_type"] = "Metastatic Ovarian Cancer"
            row["hla_allele"] = str(r[j_hla]) if j_hla is not None and r[j_hla] else NR
            row["evidence_types"] = "immunopeptidomics, T-cell assay (result figure-locked)"
            row["min_peptide_len"] = str(len(pep))
            row["tumor_specificity_basis"] = "broad-normal-panel"   # prioritized via HLA Ligand Atlas + GTEx
            row["validation_level"] = NR  # T-cell tested, but per-peptide reactive/non-reactive only in Fig 6
            row["citation_doi_pmid"] = OVA_DOI
            row["citation_location"] = "Table S3 + Fig 6"
            row["_canonical"] = "no"
            rows.append(row)

    # S2 claims, excluding the 38 already added at higher resolution from S3
    for pep, biotype, genloc, gene in s2_peps:
        if pep in s3_set:
            continue
        cls, _ = map_class(biotype, genloc)
        row = blank_row()
        row["peptide_sequence"] = pep
        row["orf_id_or_locus"] = str(gene) if gene else NR
        row["orf_class"] = cls
        row["antigen_type"] = "TSA"
        row["cancer_type"] = "Metastatic Ovarian Cancer"
        row["evidence_types"] = "immunopeptidomics"
        row["min_peptide_len"] = str(len(pep))
        row["tumor_specificity_basis"] = NR   # full identified set, not per-peptide specificity-filtered
        row["validation_level"] = "MS-presented"
        row["citation_doi_pmid"] = OVA_DOI
        row["citation_location"] = "Table S2"
        row["_canonical"] = "no"
        rows.append(row)

    wb.close()
    return rows


def main():
    hcc_path = sys.argv[1] if len(sys.argv) > 1 else HCC_DEFAULT
    ova_path = sys.argv[2] if len(sys.argv) > 2 else OVA_DEFAULT
    paths.require(hcc_path, ova_path)
    rows = ingest_hcc(hcc_path) + ingest_ovarian(ova_path)

    bad = 0
    for r in rows:
        probs = validate_row({k: v for k, v in r.items() if not k.startswith("_")})
        if probs:
            bad += 1
            if bad <= 5:
                print(f"  !! invalid: {r.get('peptide_sequence')}: {probs}", file=sys.stderr)

    fieldnames = COLUMNS + ["_canonical"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, NR) for k in fieldnames})

    from collections import Counter
    print(f"wrote {len(rows)} claims to {OUT}  (schema-invalid: {bad})")
    print("by cohort:", dict(Counter(r["cancer_type"] for r in rows)))
    print("by orf_class:", dict(Counter(r["orf_class"] for r in rows)))
    print("by validation_level:", dict(Counter(r["validation_level"] for r in rows)))
    print("canonical control rows:", sum(1 for r in rows if r["_canonical"] == "yes"))


if __name__ == "__main__":
    main()
