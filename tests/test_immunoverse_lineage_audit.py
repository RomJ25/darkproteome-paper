import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "darkproteome"))

from immunoverse_lineage_audit import (  # noqa: E402
    analyze_tables,
    minimum_compatible_study_count,
    normalize_hla_allele,
    study_bounds,
)


def _s7(peptide, labels, mapping, cancer="TEST"):
    return {
        "cancer": cancer,
        "pep": peptide,
        "samples": ",".join(labels),
        "presented_by_each_sample_hla": repr(mapping),
    }


def _s3(label, study, raw_file, hla="A*02:01"):
    return {"biology": label, "study": study, "sample": raw_file, "HLA": hla}


def test_hla_normalization_changes_format_not_allele_identity():
    assert normalize_hla_allele("HLA-A*02:01") == "A02:01"
    assert normalize_hla_allele("C14:02") == "C14:02"
    assert normalize_hla_allele(" B*13:02 ") == "B13:02"
    assert normalize_hla_allele(None) is None


def test_minimum_compatible_study_count_is_exact():
    assert minimum_compatible_study_count([{"A"}, {"A", "B"}, {"B", "C"}]) == 2
    assert minimum_compatible_study_count([{"A", "B"}, {"B", "C"}]) == 1
    assert minimum_compatible_study_count([]) == 0


def test_study_bounds_preserve_multi_study_ambiguity():
    lineage = {
        "L1": {"studies": {"A"}},
        "L2": {"studies": {"A", "B"}},
        "L3": {"studies": {"C"}},
    }
    bounds = study_bounds(["L1", "L2", "L3"], lineage)
    assert bounds["lower"] == 2
    assert bounds["upper"] == 3
    assert bounds["exact"] is False


def test_full_join_computes_conservative_study_phla_ratio_and_hla_match():
    s3 = {
        "TEST": [
            _s3("L1", "S1", "R1.raw"),
            _s3("L2", "S1", "R2.raw"),
            _s3("L3", "S2", "R3.raw"),
            _s3("L3", "S3", "R4.raw"),
        ]
    }
    s7 = [
        _s7(
            "PEPTIDEAA",
            ["L1", "L2", "L3"],
            {
                "L1": [("HLA-A*02:01", 0.2, 50.0, "SB")],
                "L2": [("HLA-A*02:01", 0.3, 60.0, "SB")],
                "L3": [(None, None, None, None)],
            },
        )
    ]
    result = analyze_tables(s7, s3, top_n=1)
    record = result["sensitivities"]["strong_or_weak_predicted_binder"]["top_examples"][0]
    assert record["sequence_recurrence_sample_condition_labels"] == 3
    assert record["sequence_recurrence_source_studies_lower"] == 2
    assert record["sequence_recurrence_source_studies_upper"] == 3
    assert record["best_predicted_phla_source_studies_upper"] == 1
    assert record["sequence_study_to_predicted_phla_ratio_lower_bound"] == 2.0
    assert result["hla_genotype_consistency"][
        "prediction_entries_matching_every_reported_s3_genotype"
    ] == 2


def test_incomplete_s7_to_s3_join_fails_closed():
    s3 = {"TEST": [_s3("L1", "S1", "R1.raw")]}
    s7 = [_s7("PEPTIDEAA", ["L2"], {"L2": []})]
    try:
        analyze_tables(s7, s3, top_n=1)
    except ValueError as error:
        assert "lineage join incomplete" in str(error)
    else:
        raise AssertionError("incomplete lineage join must fail closed")
