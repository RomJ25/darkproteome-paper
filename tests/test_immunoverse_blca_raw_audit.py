import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "darkproteome"))

from immunoverse_blca_raw_audit import audit_records, normalize_raw_name  # noqa: E402


def test_raw_name_normalization_removes_only_known_terminal_extensions():
    assert normalize_raw_name("sample.raw") == "sample"
    assert normalize_raw_name("sample.mzML") == "sample"
    assert normalize_raw_name("sample.d") == "sample"
    assert normalize_raw_name("sample.raw.extra") == "sample.raw.extra"


def test_raw_crosswalk_preserves_label_psm_and_engine_statuses():
    s7 = [{
        "pep": "PEPTIDEAA",
        "samples": "BIO1",
        "n_psm": 2,
        "presented_by_each_sample_hla": repr({
            "BIO1": [("HLA-A*02:01", 0.2, 50.0, "SB")]
        }),
    }]
    s3 = [{"sample": "RAW1.raw", "biology": "BIO1", "HLA": "A*02:01"}]
    scans = [
        {
            "Raw file": "RAW1", "Scan number": 1, "Sequence": "PEPTIDEAA",
            "Identified": "+", "Reverse": "", "Identified_vanilla": True,
            "Identified_rescore": False, "Tesorai": False, "final_identity": "nuORF",
        },
        {
            "Raw file": "RAW1", "Scan number": 2, "Sequence": "PEPTIDEAA",
            "Identified": "", "Reverse": "+", "Identified_vanilla": False,
            "Identified_rescore": True, "Tesorai": True, "tesorai_sequence": "PEPTIDEAA",
            "final_identity": "nuORF",
        },
        {"Raw file": "RAW1", "Scan number": 3, "Sequence": "OTHER"},
    ]
    result = audit_records(s7, s3, scans)
    assert result["peptides_with_exact_sequence_recovery"] == 1
    assert result["peptides_with_exact_biological_label_set_match"] == 1
    assert result["peptides_with_exact_psm_count_match"] == 1
    assert result["peptides_with_current_identified_scan"] == 1
    assert result["peptides_with_reverse_flagged_scan"] == 1
    assert result["peptides_all_matched_rows_final_identity_nuorf"] == 1


def test_recurrent_sequence_with_disjoint_genotypes_has_no_common_phla():
    s7 = [{
        "pep": "PEPTIDEAA",
        "samples": "BIO1,BIO2",
        "n_psm": 2,
        "presented_by_each_sample_hla": repr({
            "BIO1": [("HLA-A*02:01", 0.2, 50.0, "SB")],
            "BIO2": [("HLA-B*07:02", 0.3, 70.0, "SB")],
        }),
    }]
    s3 = [
        {"sample": "RAW1.raw", "biology": "BIO1", "HLA": "A*02:01"},
        {"sample": "RAW2.raw", "biology": "BIO2", "HLA": "B*07:02"},
    ]
    scans = [
        {
            "Raw file": "RAW1", "Scan number": 1, "Sequence": "PEPTIDEAA",
            "Identified": "+",
        },
        {
            "Raw file": "RAW2", "Scan number": 2, "Sequence": "PEPTIDEAA",
            "Identified": "+",
        },
    ]
    result = audit_records(s7, s3, scans)
    record = result["peptide_records"][0]
    assert record["common_reported_hla_alleles_across_expected_labels"] == []
    assert record["common_predicted_binder_alleles_across_expected_labels"] == []
    assert result[
        "recurrent_peptides_with_complete_genotypes_and_no_common_reported_hla_allele"
    ] == 1
    assert result[
        "recurrent_peptides_with_current_identified_scan_in_every_label_and_no_common_reported_hla"
    ] == 1
    assert result["recurrent_peptides_with_no_common_predicted_binder_allele"] == 1


def test_missing_raw_mapping_does_not_silently_claim_label_match():
    s7 = [{"pep": "PEPTIDEAA", "samples": "BIO1", "n_psm": 1}]
    scans = [{"Raw file": "UNKNOWN", "Scan number": 1, "Sequence": "PEPTIDEAA"}]
    result = audit_records(s7, [], scans)
    record = result["peptide_records"][0]
    assert record["observed_s3_biological_labels"] == []
    assert record["exact_biological_label_set_match"] is False
