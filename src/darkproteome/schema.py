"""The Claim Catalog data contract.

One row per distinct published claim that a non-canonical ORF / microprotein is an
HLA-presented tumor antigen. Any claim catalog (e.g. the cohort-ingested
`claim_catalog_real.csv`) must conform to this exactly, so it drops straight into the auditor.

`NOT_REPORTED` is a first-class value, not a blank: "the paper did not state this" is a
finding, and it is what pushes a claim to `insufficient-info` in the audit.
"""

NOT_REPORTED = "not reported"

# The 19 columns, in order.
COLUMNS = [
    "peptide_sequence",
    "orf_id_or_locus",
    "orf_class",
    "antigen_type",
    "cancer_type",
    "hla_allele",
    "evidence_types",
    "reported_fdr",
    "n_unique_peptides",
    "min_peptide_len",
    "periodicity_pct",
    "tumor_specificity_basis",
    "validation_level",
    "meets_consensus_bar_as_reported",
    "underlying_dataset_accession",
    "source_provenance",
    "citation_doi_pmid",
    "citation_location",
    "extraction_confidence",
]

# Controlled vocabularies for the categorical fields. Anything outside these is flagged by
# the validator (so the corpus stays clean and groupable).
CONTROLLED_VOCAB = {
    "orf_class": {
        "uORF", "lncRNA-ORF", "altORF", "3'UTR-ORF", "intronic-ORF",
        "pseudogene-ORF", "dORF", "other", NOT_REPORTED,
    },
    "antigen_type": {"TSA", "TAA", "normal-or-other", NOT_REPORTED},
    "tumor_specificity_basis": {
        "broad-normal-panel", "matched-normal-only", "not stated", NOT_REPORTED,
    },
    "validation_level": {
        "nominated", "MS-presented", "T-cell-validated", "in-vivo",
        "validated_negative", NOT_REPORTED,
    },
    "meets_consensus_bar_as_reported": {"yes", "no", "insufficient-info"},
    "source_provenance": {"primary", "review"},
    "extraction_confidence": {"high", "med", "low"},
}


def validate_row(row):
    """Return a list of human-readable problems with a claim row (empty == valid)."""
    problems = []
    for col in COLUMNS:
        if col not in row:
            problems.append(f"missing column: {col}")
    for col, allowed in CONTROLLED_VOCAB.items():
        val = (row.get(col) or "").strip()
        if val and val not in allowed:
            problems.append(f"{col}={val!r} not in controlled vocabulary")
    if not (row.get("citation_doi_pmid") or "").strip():
        problems.append("citation_doi_pmid is required (no claim without a primary source)")
    return problems
