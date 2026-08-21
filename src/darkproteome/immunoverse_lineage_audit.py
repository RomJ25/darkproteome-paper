"""Normalize ImmunoVerse peptide recurrence to published source-study metadata.

Supplementary Table S7 collapses detections to biological sample/condition labels.  Supplementary
Table S3 maps those labels to source-study identifiers, raw files and reported HLA genotypes.  This
audit joins the same-version tables and asks two conservative questions:

1. How much sequence recurrence remains after labels from the same source study are collapsed?
2. Does one predicted peptide-HLA pairing span that study-level recurrence?

Six published biological labels map to more than one source study.  Because S7 does not identify
which study supplied a peptide occurrence for those labels, the audit reports sharp compatible-study
bounds.  It does not call source studies independent patients, infer malignant-cell presentation, or
turn binding predictions into direct allele assignments.

Run from the repository root:

    python3 src/darkproteome/immunoverse_lineage_audit.py
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping

from immunoverse_phla_recurrence import (
    ANALYSIS_AUTHOR,
    REPO,
    SHEET,
    SENSITIVITIES,
    _linear_percentile,
    _rounded,
    normalize_peptide,
    parse_predicted_hla,
    parse_sample_conditions,
    sha256_file,
)


DEFAULT_S7 = REPO / "data" / "external" / "immunoverse" / "media-7_table_s7.xlsx"
DEFAULT_S3 = REPO / "data" / "external" / "immunoverse" / "media-3_table_s3.xlsx"
DEFAULT_OUTPUT = REPO / "data" / "derived_immunoverse_lineage_audit.json"
S3_URL = (
    "https://www.biorxiv.org/content/biorxiv/early/2025/07/07/"
    "2025.01.22.634237/DC3/embed/media-3.xlsx?download=true"
)
S7_URL = (
    "https://www.biorxiv.org/content/biorxiv/early/2025/07/07/"
    "2025.01.22.634237/DC7/embed/media-7.xlsx?download=true"
)


def _clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def normalize_hla_allele(value: object) -> str | None:
    """Normalize formatting only: ``HLA-A*02:01`` and ``A02:01`` become ``A02:01``."""
    allele = _clean(value).upper().replace(" ", "")
    if not allele:
        return None
    if allele.startswith("HLA-"):
        allele = allele[4:]
    allele = allele.replace("*", "")
    return allele or None


def parse_hla_genotype(value: object) -> set[str]:
    alleles: set[str] = set()
    for raw in _clean(value).split(","):
        allele = normalize_hla_allele(raw)
        if allele:
            alleles.add(allele)
    return alleles


def minimum_compatible_study_count(study_sets: Iterable[set[str] | frozenset[str]]) -> int:
    """Minimum studies needed to intersect every compatible-study set (exact hitting set)."""
    normalized = tuple(sorted({frozenset(values) for values in study_sets}, key=lambda s: sorted(s)))
    if any(not values for values in normalized):
        raise ValueError("cannot bound recurrence with an empty compatible-study set")
    if not normalized:
        return 0

    forced = frozenset(next(iter(values)) for values in normalized if len(values) == 1)
    remaining = tuple(values for values in normalized if not values.intersection(forced))

    @lru_cache(maxsize=None)
    def solve(uncovered: tuple[frozenset[str], ...]) -> int:
        if not uncovered:
            return 0
        branch = min(uncovered, key=lambda values: (len(values), sorted(values)))
        return 1 + min(
            solve(tuple(values for values in uncovered if study not in values))
            for study in sorted(branch)
        )

    return len(forced) + solve(remaining)


def study_bounds(labels: Iterable[str], lineage: Mapping[str, Mapping[str, set[str]]]) -> dict[str, object]:
    label_list = sorted(set(labels))
    compatible = [set(lineage[label]["studies"]) for label in label_list]
    if not compatible:
        return {"lower": 0, "upper": 0, "exact": True, "studies": []}
    union = set().union(*compatible)
    lower = minimum_compatible_study_count(compatible)
    return {
        "lower": lower,
        "upper": len(union),
        "exact": all(len(values) == 1 for values in compatible),
        "studies": sorted(union),
    }


def build_lineage(rows_by_cancer: Mapping[str, Iterable[Mapping[str, object]]]) -> dict[str, object]:
    lineage = defaultdict(lambda: {
        "studies": set(),
        "raw_files": set(),
        "cancers": set(),
        "hla_alleles": set(),
        "hla_genotypes": set(),
    })
    cancer_label_pairs: set[tuple[str, str]] = set()
    source_rows = 0
    for cancer, rows in rows_by_cancer.items():
        for row in rows:
            source_rows += 1
            label = _clean(row.get("biology"))
            study = _clean(row.get("study"))
            raw_file = _clean(row.get("sample"))
            if not label or not study:
                raise ValueError(f"S3 row lacks biology/study in cancer sheet {cancer!r}")
            record = lineage[label]
            record["studies"].add(study)
            record["cancers"].add(cancer)
            if raw_file:
                record["raw_files"].add(raw_file)
            genotype = frozenset(parse_hla_genotype(row.get("HLA")))
            if genotype:
                record["hla_alleles"].update(genotype)
                record["hla_genotypes"].add(genotype)
            cancer_label_pairs.add((cancer, label))
    return {
        "lineage": dict(lineage),
        "source_rows": source_rows,
        "cancer_label_pairs": cancer_label_pairs,
    }


def aggregate_s7(rows: Iterable[Mapping[str, object]]) -> dict[str, object]:
    peptide_labels: dict[str, set[str]] = defaultdict(set)
    allele_labels = {
        sensitivity: defaultdict(lambda: defaultdict(set)) for sensitivity in SENSITIVITIES
    }
    label_cancer_pairs: set[tuple[str, str]] = set()
    prediction_entries: set[tuple[str, str, str, str]] = set()
    audit = Counter()

    for row in rows:
        audit["source_rows"] += 1
        peptide = normalize_peptide(row.get("pep"))
        if peptide is None:
            audit["invalid_peptide_rows"] += 1
            continue
        labels = set(parse_sample_conditions(row.get("samples", "")))
        cancers = set(parse_sample_conditions(row.get("cancer", "")))
        if not labels:
            audit["rows_without_sample_condition"] += 1
            continue
        peptide_labels[peptide].update(labels)
        label_cancer_pairs.update((cancer, label) for cancer in cancers for label in labels)

        try:
            mapping = parse_predicted_hla(row.get("presented_by_each_sample_hla", "{}"))
        except (SyntaxError, ValueError):
            audit["predicted_hla_parse_errors"] += 1
            continue
        audit["sample_conditions_missing_mapping_key"] += len(labels - set(mapping))
        audit["mapping_keys_absent_from_samples_field"] += len(set(mapping) - labels)
        for label in labels:
            for entry in mapping.get(label, []):
                allele = normalize_hla_allele(entry[0])
                binder_class = _clean(entry[3]).upper()
                if allele is None and not binder_class:
                    audit["no_predicted_hla_sentinel_entries"] += 1
                    continue
                if allele is None:
                    audit["prediction_entries_without_allele"] += 1
                    continue
                prediction_entries.add((peptide, label, allele, binder_class))
                for sensitivity, allowed in SENSITIVITIES.items():
                    if binder_class in allowed:
                        allele_labels[sensitivity][peptide][allele].add(label)

    return {
        "peptide_labels": peptide_labels,
        "allele_labels": allele_labels,
        "label_cancer_pairs": label_cancer_pairs,
        "prediction_entries": prediction_entries,
        "audit": dict(sorted(audit.items())),
    }


def lineage_records(
    peptide_labels: Mapping[str, set[str]],
    allele_labels: Mapping[str, Mapping[str, set[str]]],
    lineage: Mapping[str, Mapping[str, set[str]]],
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for peptide, labels in peptide_labels.items():
        sequence = study_bounds(labels, lineage)
        allele_bounds = {
            allele: study_bounds(predicted_labels, lineage)
            for allele, predicted_labels in allele_labels.get(peptide, {}).items()
        }
        best_upper = max((int(bounds["upper"]) for bounds in allele_bounds.values()), default=0)
        best_lower = max((int(bounds["lower"]) for bounds in allele_bounds.values()), default=0)
        top_upper_alleles = sorted(
            allele for allele, bounds in allele_bounds.items() if int(bounds["upper"]) == best_upper
        ) if best_upper else []
        label_count = len(labels)
        sequence_lower = int(sequence["lower"])
        sequence_upper = int(sequence["upper"])
        records.append({
            "peptide": peptide,
            "sequence_recurrence_sample_condition_labels": label_count,
            "sequence_recurrence_source_studies_lower": sequence_lower,
            "sequence_recurrence_source_studies_upper": sequence_upper,
            "study_mapping_exact": bool(sequence["exact"]),
            "compatible_source_studies": sequence["studies"],
            "best_predicted_phla_source_studies_lower": best_lower,
            "best_predicted_phla_source_studies_upper": best_upper,
            "top_predicted_hla_alleles_by_study_upper": top_upper_alleles,
            "sequence_label_to_study_ratio_lower_bound": label_count / sequence_upper,
            "sequence_label_to_study_ratio_upper_bound": label_count / sequence_lower,
            "sequence_study_to_predicted_phla_ratio_lower_bound": (
                sequence_lower / best_upper if best_upper else None
            ),
            "sequence_study_to_predicted_phla_ratio_upper_bound": (
                sequence_upper / best_lower if best_lower else None
            ),
        })
    return records


def _median(records: list[dict[str, object]], field: str) -> float | None:
    values = [float(row[field]) for row in records if row[field] is not None]
    return _rounded(statistics.median(values) if values else None)


def summarize(records: list[dict[str, object]]) -> dict[str, object]:
    with_phla = [
        row for row in records
        if row["sequence_study_to_predicted_phla_ratio_lower_bound"] is not None
    ]
    exact = [row for row in records if row["study_mapping_exact"]]
    exact_with_phla = [
        row for row in exact
        if row["sequence_study_to_predicted_phla_ratio_lower_bound"] is not None
    ]
    conservative_ratios = [
        float(row["sequence_study_to_predicted_phla_ratio_lower_bound"]) for row in with_phla
    ]
    exact_ratios = [
        float(row["sequence_study_to_predicted_phla_ratio_lower_bound"])
        for row in exact_with_phla
    ]
    return {
        "peptide_count": len(records),
        "peptides_with_any_predicted_phla": len(with_phla),
        "peptides_with_exact_study_mapping": len(exact),
        "median_sequence_recurrence_sample_condition_labels": _median(
            records, "sequence_recurrence_sample_condition_labels"
        ),
        "median_sequence_recurrence_source_studies_lower": _median(
            records, "sequence_recurrence_source_studies_lower"
        ),
        "median_sequence_recurrence_source_studies_upper": _median(
            records, "sequence_recurrence_source_studies_upper"
        ),
        "median_sequence_label_to_study_ratio_lower_bound": _median(
            records, "sequence_label_to_study_ratio_lower_bound"
        ),
        "fraction_with_label_recurrence_exceeding_study_upper_bound": _rounded(
            sum(
                int(row["sequence_recurrence_sample_condition_labels"])
                > int(row["sequence_recurrence_source_studies_upper"])
                for row in records
            ) / len(records) if records else None
        ),
        "median_conservative_sequence_study_to_best_predicted_phla_ratio": _rounded(
            statistics.median(conservative_ratios) if conservative_ratios else None
        ),
        "p90_conservative_sequence_study_to_best_predicted_phla_ratio": _rounded(
            _linear_percentile(conservative_ratios, 0.9)
        ),
        "fraction_conservatively_above_one": _rounded(
            sum(value > 1 for value in conservative_ratios) / len(conservative_ratios)
            if conservative_ratios else None
        ),
        "fraction_conservatively_at_least_two": _rounded(
            sum(value >= 2 for value in conservative_ratios) / len(conservative_ratios)
            if conservative_ratios else None
        ),
        "exact_mapping_subset": {
            "peptide_count": len(exact),
            "peptides_with_any_predicted_phla": len(exact_with_phla),
            "median_sequence_source_studies": _median(
                exact, "sequence_recurrence_source_studies_lower"
            ),
            "median_best_predicted_phla_source_studies": _median(
                exact_with_phla, "best_predicted_phla_source_studies_upper"
            ),
            "median_sequence_study_to_best_predicted_phla_ratio": _rounded(
                statistics.median(exact_ratios) if exact_ratios else None
            ),
            "fraction_above_one": _rounded(
                sum(value > 1 for value in exact_ratios) / len(exact_ratios)
                if exact_ratios else None
            ),
            "fraction_at_least_two": _rounded(
                sum(value >= 2 for value in exact_ratios) / len(exact_ratios)
                if exact_ratios else None
            ),
        },
        "bounds_note": (
            "Conservative pHLA ratios divide the minimum compatible sequence-study recurrence "
            "by the maximum compatible recurrence of the best predicted pHLA."
        ),
    }


def select_rank_sets(
    records: list[dict[str, object]], top_n: int
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    ranked = sorted(
        records,
        key=lambda row: (
            -int(row["sequence_recurrence_sample_condition_labels"]), str(row["peptide"])
        ),
    )
    selected = ranked[: min(top_n, len(ranked))]
    cutoff = int(selected[-1]["sequence_recurrence_sample_condition_labels"]) if selected else 0
    tie_complete = [
        row for row in ranked
        if int(row["sequence_recurrence_sample_condition_labels"]) >= cutoff
    ] if selected else []
    return selected, tie_complete, {
        "requested_top_n": top_n,
        "selected_count": len(selected),
        "cutoff_sample_condition_recurrence": cutoff,
        "peptides_above_cutoff": sum(
            int(row["sequence_recurrence_sample_condition_labels"]) > cutoff for row in ranked
        ),
        "peptides_tied_at_cutoff": sum(
            int(row["sequence_recurrence_sample_condition_labels"]) == cutoff for row in ranked
        ),
        "tie_complete_count": len(tie_complete),
        "tie_break_rule": "descending sample/condition recurrence, then peptide ascending",
    }


def analyze_tables(
    s7_rows: Iterable[Mapping[str, object]],
    s3_rows_by_cancer: Mapping[str, Iterable[Mapping[str, object]]],
    *,
    top_n: int = 1000,
) -> dict[str, object]:
    built = build_lineage(s3_rows_by_cancer)
    lineage = built["lineage"]
    s7 = aggregate_s7(s7_rows)
    used_labels = set().union(*s7["peptide_labels"].values()) if s7["peptide_labels"] else set()
    missing_labels = sorted(used_labels - set(lineage))
    missing_pairs = sorted(s7["label_cancer_pairs"] - built["cancer_label_pairs"])
    if missing_labels or missing_pairs:
        raise ValueError(
            f"S7-to-S3 lineage join incomplete: {len(missing_labels)} labels and "
            f"{len(missing_pairs)} cancer-label pairs missing"
        )

    hla_audit = Counter()
    mismatch_examples: list[dict[str, str]] = []
    for peptide, label, allele, binder_class in sorted(s7["prediction_entries"]):
        genotypes = lineage[label]["hla_genotypes"]
        if not genotypes:
            hla_audit["prediction_entries_without_reported_s3_genotype"] += 1
        elif all(allele in genotype for genotype in genotypes):
            hla_audit["prediction_entries_matching_every_reported_s3_genotype"] += 1
        elif any(allele in genotype for genotype in genotypes):
            hla_audit[
                "prediction_entries_matching_some_but_not_every_reported_s3_genotype"
            ] += 1
        else:
            hla_audit["prediction_entries_not_in_any_reported_s3_genotype"] += 1
            if len(mismatch_examples) < 20:
                mismatch_examples.append({
                    "peptide": peptide,
                    "sample_condition": label,
                    "predicted_allele": allele,
                    "reported_genotypes": " | ".join(
                        ",".join(sorted(genotype)) for genotype in sorted(
                            genotypes, key=lambda values: sorted(values)
                        )
                    ),
                    "binder_class": binder_class,
                })

    sensitivities: dict[str, object] = {}
    for sensitivity in SENSITIVITIES:
        records = lineage_records(
            s7["peptide_labels"], s7["allele_labels"][sensitivity], lineage
        )
        top, tie_complete, ranking = select_rank_sets(records, top_n)
        sensitivities[sensitivity] = {
            "all_peptides": summarize(records),
            "top_n_by_sample_condition_recurrence": summarize(top),
            "tie_complete_at_top_n_cutoff": summarize(tie_complete),
            "ranking": ranking,
            "top_examples": top[:20],
        }

    all_studies = set().union(*(record["studies"] for record in lineage.values()))
    all_raw_files = set().union(*(record["raw_files"] for record in lineage.values()))
    return {
        "lineage_inventory": {
            "s3_source_rows": built["source_rows"],
            "s3_biological_labels": len(lineage),
            "s3_source_study_identifiers": len(all_studies),
            "s3_distinct_raw_file_names": len(all_raw_files),
            "s7_used_biological_labels": len(used_labels),
            "s7_labels_with_exact_s3_match": len(used_labels) - len(missing_labels),
            "s7_cancer_label_pairs": len(s7["label_cancer_pairs"]),
            "s7_cancer_label_pairs_with_exact_s3_match": (
                len(s7["label_cancer_pairs"]) - len(missing_pairs)
            ),
            "used_labels_mapping_to_multiple_source_studies": sum(
                len(lineage[label]["studies"]) > 1 for label in used_labels
            ),
            "multi_study_labels": {
                label: sorted(lineage[label]["studies"])
                for label in sorted(used_labels)
                if len(lineage[label]["studies"]) > 1
            },
            "used_labels_with_multiple_reported_hla_genotypes": sum(
                len(lineage[label]["hla_genotypes"]) > 1 for label in used_labels
            ),
            "multiple_reported_hla_genotype_labels": {
                label: [
                    sorted(genotype) for genotype in sorted(
                        lineage[label]["hla_genotypes"], key=lambda values: sorted(values)
                    )
                ]
                for label in sorted(used_labels)
                if len(lineage[label]["hla_genotypes"]) > 1
            },
        },
        "s7_audit": s7["audit"],
        "hla_genotype_consistency": {
            **dict(sorted(hla_audit.items())),
            "mismatch_examples": mismatch_examples,
            "normalization": "uppercase; remove HLA- prefix and asterisk only",
        },
        "sensitivities": sensitivities,
    }


def analyze_workbooks(s7_path: Path, s3_path: Path, *, top_n: int = 1000) -> dict[str, object]:
    for path in (s7_path, s3_path):
        if not path.exists():
            raise FileNotFoundError(f"missing required source workbook: {path}")
    import pandas as pd

    s7_frame = pd.read_excel(
        s7_path,
        sheet_name=SHEET,
        usecols=["cancer", "pep", "samples", "presented_by_each_sample_hla"],
    )
    s3_book = pd.ExcelFile(s3_path)
    s3_rows_by_cancer = {
        sheet: pd.read_excel(
            s3_path, sheet_name=sheet, usecols=["study", "sample", "biology", "HLA"]
        ).to_dict(orient="records")
        for sheet in s3_book.sheet_names
    }
    result = analyze_tables(
        s7_frame.to_dict(orient="records"), s3_rows_by_cancer, top_n=top_n
    )
    return {
        "schema_version": 1,
        "analysis_author": ANALYSIS_AUTHOR,
        "generated_by": "src/darkproteome/immunoverse_lineage_audit.py",
        "sources": {
            "table_s7": {
                "path": os.path.relpath(s7_path, REPO),
                "sheet": SHEET,
                "url": S7_URL,
                "sha256": sha256_file(s7_path),
                "size_bytes": s7_path.stat().st_size,
                "source_row_count": len(s7_frame),
            },
            "table_s3": {
                "path": os.path.relpath(s3_path, REPO),
                "sheets": s3_book.sheet_names,
                "url": S3_URL,
                "sha256": sha256_file(s3_path),
                "size_bytes": s3_path.stat().st_size,
            },
        },
        "estimand": {
            "lineage_unit": "published source-study identifier in Supplementary Table S3",
            "ambiguous_label_policy": "sharp lower/upper bounds; no forced study assignment",
            "phla_layer": "published SB/WB predictions, not direct allele-specific evidence",
        },
        "limitations": [
            "Source-study identifiers are not assumed to be independent patients or cohorts.",
            "Table S7 does not retain the raw-file evidence supporting each peptide occurrence.",
            "Six biological labels map to multiple source studies; results retain that ambiguity.",
            "HLA annotations are binding predictions, not direct allele-specific presentation.",
            "The analysis does not identify malignant versus nonmalignant presenting cells.",
        ],
        "result": result,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table-s7", type=Path, default=DEFAULT_S7)
    parser.add_argument("--table-s3", type=Path, default=DEFAULT_S3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-n", type=int, default=1000)
    args = parser.parse_args(argv)
    if args.top_n <= 0:
        parser.error("--top-n must be positive")
    artifact = analyze_workbooks(args.table_s7, args.table_s3, top_n=args.top_n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    result = artifact["result"]
    summary = result["sensitivities"]["strong_or_weak_predicted_binder"][
        "top_n_by_sample_condition_recurrence"
    ]
    print("IMMUNOVERSE SOURCE-STUDY LINEAGE AUDIT")
    print(
        f"  S7 labels joined to S3: {result['lineage_inventory']['s7_labels_with_exact_s3_match']:,}/"
        f"{result['lineage_inventory']['s7_used_biological_labels']:,}"
    )
    print(
        f"  top {args.top_n:,}: median labels {summary['median_sequence_recurrence_sample_condition_labels']:.1f}; "
        f"median source studies {summary['median_sequence_recurrence_source_studies_lower']:.1f}-"
        f"{summary['median_sequence_recurrence_source_studies_upper']:.1f}"
    )
    print(
        "  conservative median sequence-study / best predicted-pHLA-study ratio: "
        f"{summary['median_conservative_sequence_study_to_best_predicted_phla_ratio']:.3f}"
    )
    print(
        "  interpretation: source-study lineage bounds, not patient reach or direct pHLA evidence"
    )
    print(f"  wrote {os.path.relpath(args.output, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
