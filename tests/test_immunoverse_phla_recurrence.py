import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "darkproteome"))

from immunoverse_phla_recurrence import (  # noqa: E402
    aggregate_rows,
    analyze_rows,
    parse_predicted_hla,
    parse_sample_conditions,
    select_rank_sets,
)


def _row(peptide, samples, mapping):
    return {
        "pep": peptide,
        "samples": samples,
        "presented_by_each_sample_hla": repr(mapping),
    }


def test_parsers_preserve_sample_order_and_validate_hla_shape():
    assert parse_sample_conditions("A, B,A") == ("A", "B")
    parsed = parse_predicted_hla("{'A': [('HLA-A*02:01', 0.2, 50.0, 'SB')], 'B': []}")
    assert parsed["A"][0][0] == "HLA-A*02:01"
    assert parsed["B"] == []


def test_duplicate_source_rows_union_samples_and_deduplicate_hla_entries():
    rows = [
        _row(
            "PEPTIDEAA",
            "S1,S2",
            {
                "S1": [("HLA-A*02:01", 0.2, 50.0, "SB")],
                "S2": [("HLA-B*07:02", 1.2, 500.0, "WB")],
            },
        ),
        _row(
            "peptideaa",
            "S2,S3",
            {
                "S2": [
                    ("HLA-B*07:02", 1.2, 500.0, "WB"),
                    ("HLA-B*07:02", 1.2, 500.0, "WB"),
                ],
                "S3": [("HLA-A*02:01", 0.4, 80.0, "SB")],
            },
        ),
    ]
    aggregated = aggregate_rows(rows)
    record = aggregated["records_by_sensitivity"]["strong_or_weak_predicted_binder"][0]
    assert record["sequence_recurrence_sample_conditions"] == 3
    assert record["max_predicted_phla_recurrence_sample_conditions"] == 2
    assert record["top_predicted_hla_alleles"] == ["HLA-A*02:01"]
    assert math.isclose(record["sequence_to_predicted_phla_ratio"], 1.5)


def test_strong_only_sensitivity_does_not_credit_weak_binders():
    rows = [
        _row(
            "PEPTIDEAA",
            "S1,S2",
            {
                "S1": [("HLA-A*02:01", 0.2, 50.0, "SB")],
                "S2": [("HLA-A*02:01", 1.2, 500.0, "WB")],
            },
        )
    ]
    result = analyze_rows(rows, top_n=1)
    any_binder = result["sensitivities"]["strong_or_weak_predicted_binder"]["all_peptides"]
    strong_only = result["sensitivities"]["strong_predicted_binder_only"]["all_peptides"]
    assert any_binder["median_max_predicted_phla_fraction"] == 1.0
    assert strong_only["median_max_predicted_phla_fraction"] == 0.5
    assert strong_only["median_sequence_to_predicted_phla_ratio"] == 2.0


def test_none_tuple_is_an_explicit_no_prediction_sentinel():
    rows = [
        _row(
            "PEPTIDEAA",
            "S1,S2",
            {
                "S1": [(None, None, None, None)],
                "S2": [("HLA-A*02:01", 0.2, 50.0, "SB")],
            },
        )
    ]
    aggregated = aggregate_rows(rows)
    record = aggregated["records_by_sensitivity"]["strong_or_weak_predicted_binder"][0]
    assert record["max_predicted_phla_recurrence_sample_conditions"] == 1
    assert record["sample_conditions_without_predicted_hla"] == 1
    assert aggregated["audit"]["no_predicted_hla_sentinel_entries"] == 1
    assert aggregated["binder_class_entry_counts"] == {
        "NO_PREDICTED_HLA_SENTINEL": 1,
        "SB": 1,
    }


def test_top_n_records_tie_cutoff_instead_of_hiding_boundary_ties():
    records = [
        {"peptide": "A", "sequence_recurrence_sample_conditions": 5},
        {"peptide": "B", "sequence_recurrence_sample_conditions": 4},
        {"peptide": "C", "sequence_recurrence_sample_conditions": 4},
        {"peptide": "D", "sequence_recurrence_sample_conditions": 1},
    ]
    selected, tie_complete, metadata = select_rank_sets(records, 2)
    assert [r["peptide"] for r in selected] == ["A", "B"]
    assert [r["peptide"] for r in tie_complete] == ["A", "B", "C"]
    assert metadata["cutoff_sequence_recurrence_sample_conditions"] == 4
    assert metadata["peptides_tied_at_cutoff"] == 2
