"""The Claim Catalog data contract.

One row per distinct published claim that a non-canonical ORF / microprotein is an
HLA-presented tumor antigen. Any claim catalog (e.g. the cohort-ingested
`claim_catalog_real.csv`) must conform to this exactly, so it drops straight into the auditor.

`NOT_REPORTED` is a first-class value, not a blank: "the paper did not state this" is a
finding, and it is what pushes a claim to `insufficient-info` in the audit.
"""

NOT_REPORTED = "not reported"

# The columns, in order.
#
# Two distinctions are load-bearing and must never collapse back into one field:
#
#   ligand_len  vs  source_pep_len
#       The HLA ligand (8-12aa) and the tryptic peptide supporting the SOURCE PROTEIN (>=9aa) are
#       different measurements. Most HLA ligands are >= 9aa, so a shared column would satisfy the
#       protein-existence rule as a coincidence of the number rather than as evidence.
#
#   tumor_specificity_modality  vs  _scope  vs  _result
#       A normal-ligandome search, a GTEx RNA threshold and a cohort inclusion criterion are not
#       equivalent evidence and must not earn the same verdict. Transcript abundance does not
#       determine HLA presentation, so an RNA threshold cannot answer a presentation claim.
COLUMNS = [
    "peptide_sequence",
    "orf_id_or_locus",
    "orf_class",
    "antigen_type",
    "cancer_type",
    "hla_allele",
    "evidence_types",
    "reported_fdr",             # protein/ORF-level FDR, as reported
    "psm_qvalue",               # PSM-level q-value — a DIFFERENT unit; never read as an FDR
    # D_N: accepted DECOY PSMs in this claim's class, at the accepted threshold. Without it a
    # class-specific FDR is not reconstructible. It is populated on ZERO rows of the corpus, and
    # that is the point — HCC's S26 even ships a `target_decoy` column whose 43 rows are ALL
    # targets. The column exists here so the absence is explicit and machine-checkable rather
    # than merely narrated, and so any future source that DOES report it drops straight in.
    "accepted_decoys_in_class",
    "n_unique_peptides",
    "ligand_len",               # the HLA ligand itself (8–12aa)
    "source_pep_len",           # tryptic peptide supporting the SOURCE protein (≥9aa)
    "periodicity_pct",
    "tumor_specificity_modality",
    "tumor_specificity_scope",
    # WHAT the check actually returned. The old code recorded only the "yes, tumour-specific"
    # answers and threw the "no" away -- but an explicit detection in normal tissue is the ONE
    # true empirical negative for tumour absence, and discarding it made the record look merely
    # silent where it had in fact spoken.
    "tumor_specificity_result",
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
    # (what was measured) x (which normals). Only a BROAD normal LIGANDOME is matched to a
    # presentation claim: RNA abundance does not determine presentation, and adjacent tissue is
    # one organ from one donor.
    "tumor_specificity_modality": {
        "normal-ligandome-broad", "normal-ligandome-matched",
        "normal-rna-broad", "normal-rna-matched", "not stated", NOT_REPORTED,
    },
    # WHETHER it was reported per claim, or is merely how the authors filtered the table.
    "tumor_specificity_scope": {
        "per-claim-reported", "cohort-inclusion-criterion", "not stated", NOT_REPORTED,
    },
    # WHAT it returned. `detected-in-normal` is the only true empirical negative in the corpus.
    "tumor_specificity_result": {
        "absent-from-normal", "detected-in-normal", "not stated", NOT_REPORTED,
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
