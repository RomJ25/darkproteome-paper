"""Trace historical ImmunoVerse catalogue peptides back to current scan-level output.

The default bounded audit crosswalks the 16 BLCA ncORF peptides in the 2025-07-07 Supplementary
Table S7 to the authors' current consolidated scan table.  ``--cancer`` supports another cancer with
the same schema.  The audit checks exact sequence recovery, raw-file-to-biological-label agreement,
PSM-count agreement and current engine decision flags.  Because the consolidated tables post-date
S7, status differences are versioned observations, not evidence that a historical identification
was wrong.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping

from immunoverse_lineage_audit import (
    DEFAULT_S3,
    DEFAULT_S7,
    normalize_hla_allele,
    parse_hla_genotype,
)
from immunoverse_phla_recurrence import (
    ANALYSIS_AUTHOR,
    REPO,
    normalize_peptide,
    parse_predicted_hla,
    sha256_file,
)


CANCER = "BLCA"
DEFAULT_CONSOLIDATED = (
    REPO / "data" / "external" / "immunoverse" / "raw_ms_result" / "consolidated" /
    CANCER / "msmsScans_all_add_tesorai.txt"
)
DEFAULT_OUTPUT = REPO / "data" / "derived_immunoverse_blca_raw_audit.json"
CONSOLIDATED_BASE_URL = (
    "https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/raw_MS_result/consolidated"
)
SCAN_COLUMNS = [
    "Raw file", "Scan number", "Identified", "Reverse", "Sequence", "Proteins", "Score", "PEP",
    "Identified_vanilla", "qval_vanilla", "Identified_rescore", "qval_rescore",
    "tesorai_sequence", "tesorai_proteins", "tesorai_score", "Tesorai", "final_identity",
]


def _clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return _clean(value).lower() in {"true", "1", "+", "yes"}


def _number(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def normalize_raw_name(value: object) -> str:
    name = _clean(value)
    lowered = name.lower()
    for suffix in (".raw", ".mzml", ".d"):
        if lowered.endswith(suffix):
            return name[:-len(suffix)]
    return name


def _split_labels(value: object) -> set[str]:
    return {part.strip() for part in _clean(value).split(",") if part.strip()}


def audit_records(
    s7_rows: Iterable[Mapping[str, object]],
    s3_rows: Iterable[Mapping[str, object]],
    scan_rows: Iterable[Mapping[str, object]],
) -> dict[str, object]:
    expected = defaultdict(lambda: {
        "labels": set(),
        "n_psm": 0,
        "source_rows": 0,
        "predicted_alleles_by_label": defaultdict(set),
    })
    for row in s7_rows:
        peptide = normalize_peptide(row.get("pep"))
        if peptide is None:
            continue
        expected[peptide]["labels"].update(_split_labels(row.get("samples")))
        expected[peptide]["n_psm"] += int(row.get("n_psm", 0))
        expected[peptide]["source_rows"] += 1
        try:
            prediction_mapping = parse_predicted_hla(
                row.get("presented_by_each_sample_hla", "{}")
            )
        except (SyntaxError, ValueError):
            prediction_mapping = {}
        for label, entries in prediction_mapping.items():
            for entry in entries:
                allele = normalize_hla_allele(entry[0])
                binder_class = _clean(entry[3]).upper()
                if allele and binder_class in {"SB", "WB"}:
                    expected[peptide]["predicted_alleles_by_label"][label].add(allele)

    raw_to_labels: dict[str, set[str]] = defaultdict(set)
    label_to_reported_hla: dict[str, set[str]] = defaultdict(set)
    for row in s3_rows:
        raw_name = normalize_raw_name(row.get("sample"))
        label = _clean(row.get("biology"))
        if raw_name and label:
            raw_to_labels[raw_name].add(label)
        if label:
            label_to_reported_hla[label].update(parse_hla_genotype(row.get("HLA")))

    matched = defaultdict(list)
    scan_source_rows = 0
    for row in scan_rows:
        scan_source_rows += 1
        peptide = normalize_peptide(row.get("Sequence"))
        if peptide in expected:
            matched[peptide].append(row)

    peptide_records: list[dict[str, object]] = []
    for peptide in sorted(expected):
        rows = matched.get(peptide, [])
        observed_raw_files = sorted({normalize_raw_name(row.get("Raw file")) for row in rows})
        observed_labels = sorted(set().union(*(
            raw_to_labels.get(raw_file, set()) for raw_file in observed_raw_files
        )) if observed_raw_files else set())
        expected_labels = sorted(expected[peptide]["labels"])
        reported_hla_by_label = {
            label: sorted(label_to_reported_hla.get(label, set())) for label in expected_labels
        }
        predicted_hla_by_label = {
            label: sorted(expected[peptide]["predicted_alleles_by_label"].get(label, set()))
            for label in expected_labels
        }
        common_reported_hla = sorted(set.intersection(*(
            set(reported_hla_by_label[label]) for label in expected_labels
        ))) if expected_labels else []
        common_predicted_hla = sorted(set.intersection(*(
            set(predicted_hla_by_label[label]) for label in expected_labels
        ))) if expected_labels else []
        labels_without_reported_hla = sorted(
            label for label in expected_labels if not reported_hla_by_label[label]
        )
        final_identities = sorted({_clean(row.get("final_identity")) for row in rows if _clean(row.get("final_identity"))})
        exact_scan_rows_by_label = Counter()
        identified_scan_rows_by_label = Counter()
        for row in rows:
            row_labels = raw_to_labels.get(normalize_raw_name(row.get("Raw file")), set())
            for label in row_labels:
                exact_scan_rows_by_label[label] += 1
                if _clean(row.get("Identified")) == "+":
                    identified_scan_rows_by_label[label] += 1
        evidence = []
        for row in rows:
            evidence.append({
                "raw_file": normalize_raw_name(row.get("Raw file")),
                "scan_number": int(row["Scan number"]) if _number(row.get("Scan number")) is not None else None,
                "identified_current_consolidated": _clean(row.get("Identified")) == "+",
                "reverse_flag": _clean(row.get("Reverse")) == "+",
                "identified_vanilla": _truthy(row.get("Identified_vanilla")),
                "qval_vanilla": _number(row.get("qval_vanilla")),
                "identified_rescore": _truthy(row.get("Identified_rescore")),
                "qval_rescore": _number(row.get("qval_rescore")),
                "identified_tesorai": _truthy(row.get("Tesorai")),
                "tesorai_sequence_exact": normalize_peptide(row.get("tesorai_sequence")) == peptide,
                "score": _number(row.get("Score")),
                "pep": _number(row.get("PEP")),
                "final_identity": _clean(row.get("final_identity")),
                "proteins": _clean(row.get("Proteins")),
            })
        peptide_records.append({
            "peptide": peptide,
            "s7_source_row_count": expected[peptide]["source_rows"],
            "s7_expected_sample_condition_labels": expected_labels,
            "observed_raw_files": observed_raw_files,
            "observed_s3_biological_labels": observed_labels,
            "exact_biological_label_set_match": observed_labels == expected_labels,
            "exact_sequence_scan_rows_by_label": {
                label: exact_scan_rows_by_label[label] for label in expected_labels
            },
            "current_identified_scan_rows_by_label": {
                label: identified_scan_rows_by_label[label] for label in expected_labels
            },
            "every_expected_label_has_exact_sequence_scan": all(
                exact_scan_rows_by_label[label] > 0 for label in expected_labels
            ),
            "every_expected_label_has_current_identified_scan": all(
                identified_scan_rows_by_label[label] > 0 for label in expected_labels
            ),
            "reported_hla_alleles_by_label": reported_hla_by_label,
            "expected_labels_without_reported_hla": labels_without_reported_hla,
            "common_reported_hla_alleles_across_expected_labels": common_reported_hla,
            "predicted_binder_alleles_by_label": predicted_hla_by_label,
            "common_predicted_binder_alleles_across_expected_labels": common_predicted_hla,
            "s7_n_psm": expected[peptide]["n_psm"],
            "consolidated_exact_sequence_scan_rows": len(rows),
            "exact_psm_count_match": len(rows) == expected[peptide]["n_psm"],
            "currently_identified_scan_rows": sum(
                _clean(row.get("Identified")) == "+" for row in rows
            ),
            "has_current_identified_scan": any(
                _clean(row.get("Identified")) == "+" for row in rows
            ),
            "has_reverse_flagged_scan": any(_clean(row.get("Reverse")) == "+" for row in rows),
            "final_identities": final_identities,
            "evidence": evidence,
        })

    recurrent = [
        row for row in peptide_records if len(row["s7_expected_sample_condition_labels"]) > 1
    ]
    return {
        "scan_source_rows": scan_source_rows,
        "s7_peptide_count": len(peptide_records),
        "s7_total_n_psm": sum(int(row["s7_n_psm"]) for row in peptide_records),
        "consolidated_exact_sequence_scan_rows": sum(
            int(row["consolidated_exact_sequence_scan_rows"]) for row in peptide_records
        ),
        "peptides_with_exact_sequence_recovery": sum(
            int(row["consolidated_exact_sequence_scan_rows"]) > 0 for row in peptide_records
        ),
        "peptides_with_exact_biological_label_set_match": sum(
            bool(row["exact_biological_label_set_match"]) for row in peptide_records
        ),
        "peptides_with_exact_psm_count_match": sum(
            bool(row["exact_psm_count_match"]) for row in peptide_records
        ),
        "peptides_with_current_identified_scan": sum(
            bool(row["has_current_identified_scan"]) for row in peptide_records
        ),
        "peptides_with_reverse_flagged_scan": sum(
            bool(row["has_reverse_flagged_scan"]) for row in peptide_records
        ),
        "peptides_all_matched_rows_final_identity_nuorf": sum(
            bool(row["evidence"]) and row["final_identities"] == ["nuORF"]
            for row in peptide_records
        ),
        "recurrent_multi_label_peptides": len(recurrent),
        "recurrent_peptides_with_exact_biological_label_set_match": sum(
            bool(row["exact_biological_label_set_match"]) for row in recurrent
        ),
        "recurrent_peptides_with_exact_sequence_scan_in_every_expected_label": sum(
            bool(row["every_expected_label_has_exact_sequence_scan"]) for row in recurrent
        ),
        "recurrent_peptides_with_current_identified_scan_in_every_expected_label": sum(
            bool(row["every_expected_label_has_current_identified_scan"]) for row in recurrent
        ),
        "recurrent_peptides_with_complete_reported_hla_genotypes": sum(
            not row["expected_labels_without_reported_hla"] for row in recurrent
        ),
        "recurrent_peptides_with_complete_genotypes_and_no_common_reported_hla_allele": sum(
            not row["expected_labels_without_reported_hla"]
            and not row["common_reported_hla_alleles_across_expected_labels"]
            for row in recurrent
        ),
        "recurrent_peptides_with_current_identified_scan_in_every_label_and_no_common_reported_hla": sum(
            bool(row["every_expected_label_has_current_identified_scan"])
            and not row["expected_labels_without_reported_hla"]
            and not row["common_reported_hla_alleles_across_expected_labels"]
            for row in recurrent
        ),
        "recurrent_peptides_with_no_common_predicted_binder_allele": sum(
            not row["common_predicted_binder_alleles_across_expected_labels"] for row in recurrent
        ),
        "peptide_records": peptide_records,
    }


def analyze_files(
    s7_path: Path, s3_path: Path, consolidated_path: Path, *, cancer: str = CANCER
) -> dict[str, object]:
    cancer = cancer.strip().upper()
    for path in (s7_path, s3_path, consolidated_path):
        if not path.exists():
            raise FileNotFoundError(f"missing required input: {path}")
    import pandas as pd

    s7 = pd.read_excel(
        s7_path,
        sheet_name="ORF_antigen",
        usecols=["cancer", "pep", "samples", "n_psm", "presented_by_each_sample_hla"],
    )
    s7 = s7[s7["cancer"].astype(str).map(lambda value: cancer in _split_labels(value))]
    s3 = pd.read_excel(
        s3_path, sheet_name=cancer, usecols=["study", "sample", "biology", "HLA"]
    )
    targets = {normalize_peptide(value) for value in s7["pep"]}
    scan_rows: list[dict[str, object]] = []
    scan_source_rows = 0
    for chunk in pd.read_csv(
        consolidated_path, sep="\t", usecols=SCAN_COLUMNS, chunksize=100_000, low_memory=False
    ):
        scan_source_rows += len(chunk)
        subset = chunk[chunk["Sequence"].astype(str).str.upper().isin(targets)]
        scan_rows.extend(subset.to_dict(orient="records"))
    result = audit_records(
        s7.to_dict(orient="records"), s3.to_dict(orient="records"), scan_rows
    )
    result["scan_source_rows"] = scan_source_rows
    return {
        "schema_version": 1,
        "analysis_author": ANALYSIS_AUTHOR,
        "generated_by": "src/darkproteome/immunoverse_blca_raw_audit.py",
        "sources": {
            "table_s7": {"path": os.path.relpath(s7_path, REPO), "sha256": sha256_file(s7_path)},
            "table_s3": {"path": os.path.relpath(s3_path, REPO), "sha256": sha256_file(s3_path)},
            "consolidated_scan_table": {
                "path": os.path.relpath(consolidated_path, REPO),
                "url": (
                    f"{CONSOLIDATED_BASE_URL}/{cancer}/msmsScans_all_add_tesorai.txt"
                ),
                "sha256": sha256_file(consolidated_path),
                "size_bytes": consolidated_path.stat().st_size,
                "server_last_modified": "11:50 (directory listing)",
            },
        },
        "scope": {
            "cancer": cancer,
            "historical_catalogue_version": "bioRxiv supplements dated 2025-07-07",
            "scan_table_version": "current consolidated release, last modified 2026-01-23",
        },
        "limitations": [
            "The consolidated scan table post-dates the historical S7 catalogue.",
            "Current Identified/Reverse/engine flags are not assumed to reproduce historical policy.",
            "A single-cancer subset cannot establish pan-cancer behavior.",
            "Exact sequence recovery confirms evidence linkage, not biological source-locus identity.",
        ],
        "result": result,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-s7", type=Path, default=DEFAULT_S7)
    parser.add_argument("--table-s3", type=Path, default=DEFAULT_S3)
    parser.add_argument("--cancer", default=CANCER)
    parser.add_argument("--consolidated", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    cancer = args.cancer.strip().upper()
    consolidated = args.consolidated or (
        REPO / "data" / "external" / "immunoverse" / "raw_ms_result" / "consolidated" /
        cancer / "msmsScans_all_add_tesorai.txt"
    )
    output = args.output or (
        REPO / "data" / f"derived_immunoverse_{cancer.lower()}_raw_audit.json"
    )
    artifact = analyze_files(args.table_s7, args.table_s3, consolidated, cancer=cancer)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result = artifact["result"]
    print(f"IMMUNOVERSE {cancer} RAW-LINEAGE AUDIT")
    print(
        f"  exact sequence recovery: {result['peptides_with_exact_sequence_recovery']}/"
        f"{result['s7_peptide_count']} peptides"
    )
    print(
        f"  S7 PSM total vs exact scan rows: {result['s7_total_n_psm']}/"
        f"{result['consolidated_exact_sequence_scan_rows']}"
    )
    print(
        f"  exact biological-label sets: {result['peptides_with_exact_biological_label_set_match']}/"
        f"{result['s7_peptide_count']}"
    )
    print(
        f"  current Identified support: {result['peptides_with_current_identified_scan']}/"
        f"{result['s7_peptide_count']}; reverse-flagged: "
        f"{result['peptides_with_reverse_flagged_scan']}/{result['s7_peptide_count']}"
    )
    print(f"  wrote {os.path.relpath(output, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
