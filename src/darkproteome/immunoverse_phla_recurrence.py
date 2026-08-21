"""Audit sequence recurrence against predicted peptide-HLA recurrence in ImmunoVerse.

This module deliberately measures a *prediction-layer feasibility diagnostic*.  ImmunoVerse
Supplementary Table S7 reports peptide detections by sample/condition and, separately, candidate
HLA binders encoded as tuples containing an allele, rank, affinity and ``SB``/``WB`` label.  Those
annotations are not allele-specific immunoprecipitation or direct experimental HLA assignments.

The analysis therefore asks a bounded question:

    How often does recurrence of a peptide string exceed recurrence of its most widespread
    *predicted* peptide-HLA pairing across the published sample/condition labels?

It does not estimate independent-patient reach, observed pHLA recurrence, immunogenicity, or
malignant-cell presentation.  Duplicate peptide rows (for example, multiple source annotations)
are collapsed by taking the union of their sample/condition and predicted-HLA annotations.

Run from the repository root:

    python3 src/darkproteome/immunoverse_phla_recurrence.py
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


REPO = Path(__file__).resolve().parents[2]
DEFAULT_XLSX = REPO / "data" / "external" / "immunoverse" / "media-7_table_s7.xlsx"
DEFAULT_OUTPUT = REPO / "data" / "derived_immunoverse_phla_recurrence.json"
SHEET = "ORF_antigen"
SENSITIVITIES = {
    "strong_or_weak_predicted_binder": frozenset({"SB", "WB"}),
    "strong_predicted_binder_only": frozenset({"SB"}),
}
ANALYSIS_AUTHOR = "Rom Jan"


def normalize_peptide(value: object) -> str | None:
    peptide = str(value).strip().upper()
    if not peptide or not peptide.isalpha():
        return None
    return peptide


def parse_sample_conditions(value: object) -> tuple[str, ...]:
    """Parse the comma-separated ImmunoVerse sample/condition field, preserving first order."""
    seen: set[str] = set()
    parsed: list[str] = []
    for raw in str(value).split(","):
        sample = raw.strip()
        if sample and sample not in seen:
            seen.add(sample)
            parsed.append(sample)
    return tuple(parsed)


def parse_predicted_hla(value: object) -> dict[str, list[tuple[object, ...]]]:
    """Safely parse the published Python-literal mapping from sample to HLA tuples."""
    parsed = ast.literal_eval(str(value))
    if not isinstance(parsed, dict):
        raise ValueError("predicted-HLA field is not a dictionary")
    out: dict[str, list[tuple[object, ...]]] = {}
    for raw_sample, raw_entries in parsed.items():
        sample = str(raw_sample).strip()
        if not sample:
            raise ValueError("predicted-HLA mapping contains an empty sample key")
        if raw_entries is None:
            raw_entries = []
        if not isinstance(raw_entries, (list, tuple)):
            raise ValueError(f"predicted-HLA entries for {sample!r} are not a list")
        entries: list[tuple[object, ...]] = []
        for entry in raw_entries:
            if not isinstance(entry, (list, tuple)) or len(entry) < 4:
                raise ValueError(f"malformed predicted-HLA tuple for {sample!r}: {entry!r}")
            entries.append(tuple(entry))
        out[sample] = entries
    return out


def _linear_percentile(values: list[float], probability: float) -> float | None:
    """Excel/Pandas-style linearly interpolated percentile for a finite numeric list."""
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _rounded(value: float | None) -> float | None:
    return None if value is None else round(float(value), 12)


def aggregate_rows(rows: Iterable[Mapping[str, object]]) -> dict[str, object]:
    """Collapse source-duplicated rows into peptide-level sample and predicted-HLA sets."""
    peptide_samples: dict[str, set[str]] = defaultdict(set)
    peptide_allele_samples = {
        sensitivity: defaultdict(lambda: defaultdict(set)) for sensitivity in SENSITIVITIES
    }
    peptide_sample_alleles = {
        sensitivity: defaultdict(lambda: defaultdict(set)) for sensitivity in SENSITIVITIES
    }
    audit = Counter()
    binder_classes = Counter()

    for row in rows:
        audit["source_rows"] += 1
        peptide = normalize_peptide(row.get("pep"))
        if peptide is None:
            audit["invalid_peptide_rows"] += 1
            continue
        samples = set(parse_sample_conditions(row.get("samples", "")))
        if not samples:
            audit["rows_without_sample_condition"] += 1
            continue
        peptide_samples[peptide].update(samples)

        try:
            mapping = parse_predicted_hla(row.get("presented_by_each_sample_hla", "{}"))
        except (SyntaxError, ValueError):
            audit["predicted_hla_parse_errors"] += 1
            continue

        audit["sample_conditions_missing_mapping_key"] += len(samples - set(mapping))
        audit["mapping_keys_absent_from_samples_field"] += len(set(mapping) - samples)

        for sample in samples:
            entries = mapping.get(sample, [])
            if not entries:
                audit["sample_condition_entries_with_empty_hla_list"] += 1
            for entry in entries:
                raw_allele = entry[0]
                raw_binder_class = entry[3]
                if raw_allele is None and raw_binder_class is None:
                    audit["no_predicted_hla_sentinel_entries"] += 1
                    binder_classes["NO_PREDICTED_HLA_SENTINEL"] += 1
                    continue
                allele = "" if raw_allele is None else str(raw_allele).strip()
                binder_class = (
                    "" if raw_binder_class is None else str(raw_binder_class).strip().upper()
                )
                binder_classes[binder_class] += 1
                if not allele:
                    audit["empty_allele_entries"] += 1
                    continue
                if binder_class not in {"SB", "WB"}:
                    audit["unexpected_binder_class_entries"] += 1
                for sensitivity, allowed in SENSITIVITIES.items():
                    if binder_class not in allowed:
                        continue
                    peptide_allele_samples[sensitivity][peptide][allele].add(sample)
                    peptide_sample_alleles[sensitivity][peptide][sample].add(allele)

    records_by_sensitivity: dict[str, list[dict[str, object]]] = {}
    for sensitivity in SENSITIVITIES:
        records: list[dict[str, object]] = []
        for peptide, samples in peptide_samples.items():
            allele_counts = {
                allele: len(allele_samples)
                for allele, allele_samples in peptide_allele_samples[sensitivity][peptide].items()
            }
            max_recurrence = max(allele_counts.values(), default=0)
            top_alleles = sorted(
                allele for allele, count in allele_counts.items() if count == max_recurrence
            ) if max_recurrence else []
            sample_count = len(samples)
            annotated_samples = sum(
                bool(peptide_sample_alleles[sensitivity][peptide].get(sample)) for sample in samples
            )
            singleton_samples = sum(
                len(peptide_sample_alleles[sensitivity][peptide].get(sample, set())) == 1
                for sample in samples
            )
            records.append({
                "peptide": peptide,
                "sequence_recurrence_sample_conditions": sample_count,
                "max_predicted_phla_recurrence_sample_conditions": max_recurrence,
                "top_predicted_hla_alleles": top_alleles,
                "candidate_predicted_hla_allele_count": len(allele_counts),
                "max_predicted_phla_fraction": max_recurrence / sample_count,
                "sequence_to_predicted_phla_ratio": (
                    sample_count / max_recurrence if max_recurrence else None
                ),
                "predicted_hla_annotated_sample_fraction": annotated_samples / sample_count,
                "single_predicted_hla_sample_fraction": singleton_samples / sample_count,
                "sample_conditions_without_predicted_hla": sample_count - annotated_samples,
            })
        records_by_sensitivity[sensitivity] = records

    return {
        "peptide_samples": peptide_samples,
        "records_by_sensitivity": records_by_sensitivity,
        "audit": dict(sorted(audit.items())),
        "binder_class_entry_counts": dict(sorted(binder_classes.items())),
    }


def summarize(records: list[dict[str, object]]) -> dict[str, object]:
    with_hla = [r for r in records if r["max_predicted_phla_recurrence_sample_conditions"] > 0]
    fractions = [float(r["max_predicted_phla_fraction"]) for r in with_hla]
    ratios = [float(r["sequence_to_predicted_phla_ratio"]) for r in with_hla]
    recurrences = [int(r["sequence_recurrence_sample_conditions"]) for r in records]
    return {
        "peptide_count": len(records),
        "peptides_with_any_predicted_hla": len(with_hla),
        "peptides_without_any_predicted_hla": len(records) - len(with_hla),
        "peptides_with_at_least_one_unannotated_sample_condition": sum(
            int(r["sample_conditions_without_predicted_hla"]) > 0 for r in records
        ),
        "median_sequence_recurrence_sample_conditions": _rounded(
            statistics.median(recurrences) if recurrences else None
        ),
        "maximum_sequence_recurrence_sample_conditions": max(recurrences, default=None),
        "fraction_with_one_predicted_hla_covering_every_sample_condition": _rounded(
            sum(f == 1 for f in fractions) / len(fractions) if fractions else None
        ),
        "fraction_with_max_predicted_phla_fraction_below_one": _rounded(
            sum(f < 1 for f in fractions) / len(fractions) if fractions else None
        ),
        "fraction_with_max_predicted_phla_fraction_at_most_half": _rounded(
            sum(f <= 0.5 for f in fractions) / len(fractions) if fractions else None
        ),
        "median_max_predicted_phla_fraction": _rounded(
            statistics.median(fractions) if fractions else None
        ),
        "median_sequence_to_predicted_phla_ratio": _rounded(
            statistics.median(ratios) if ratios else None
        ),
        "p90_sequence_to_predicted_phla_ratio": _rounded(_linear_percentile(ratios, 0.9)),
        "fraction_denominator_note": (
            "Fraction and ratio summaries exclude peptides with no predicted HLA under this "
            "sensitivity; their count is reported separately."
        ),
    }


def select_rank_sets(
    records: list[dict[str, object]], top_n: int
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, int]]:
    """Return deterministic top-N and the tie-complete set at the Nth recurrence threshold."""
    ranked = sorted(
        records,
        key=lambda row: (-int(row["sequence_recurrence_sample_conditions"]), str(row["peptide"])),
    )
    selected = ranked[: min(top_n, len(ranked))]
    cutoff = int(selected[-1]["sequence_recurrence_sample_conditions"]) if selected else 0
    tie_complete = [
        row for row in ranked if int(row["sequence_recurrence_sample_conditions"]) >= cutoff
    ] if selected else []
    above = sum(
        int(row["sequence_recurrence_sample_conditions"]) > cutoff for row in ranked
    )
    tied = sum(
        int(row["sequence_recurrence_sample_conditions"]) == cutoff for row in ranked
    )
    return selected, tie_complete, {
        "requested_top_n": top_n,
        "selected_count": len(selected),
        "cutoff_sequence_recurrence_sample_conditions": cutoff,
        "peptides_above_cutoff": above,
        "peptides_tied_at_cutoff": tied,
        "tie_complete_count": len(tie_complete),
        "tie_break_rule": "descending sequence recurrence, then peptide string ascending",
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze_rows(
    rows: Iterable[Mapping[str, object]], *, top_n: int = 1000
) -> dict[str, object]:
    aggregated = aggregate_rows(rows)
    peptide_samples = aggregated["peptide_samples"]
    records_by_sensitivity = aggregated["records_by_sensitivity"]
    sample_conditions = set().union(*peptide_samples.values()) if peptide_samples else set()
    result: dict[str, object] = {
        "unique_peptide_sequences": len(peptide_samples),
        "unique_sample_condition_labels": len(sample_conditions),
        "recurrent_peptides_n_ge_2_sample_conditions": sum(
            len(samples) >= 2 for samples in peptide_samples.values()
        ),
        "audit": aggregated["audit"],
        "binder_class_entry_counts": aggregated["binder_class_entry_counts"],
        "sensitivities": {},
    }
    for sensitivity, records in records_by_sensitivity.items():
        recurrent = [
            row for row in records if int(row["sequence_recurrence_sample_conditions"]) >= 2
        ]
        top, tie_complete, rank_metadata = select_rank_sets(records, top_n)
        result["sensitivities"][sensitivity] = {
            "all_peptides": summarize(records),
            "recurrent_peptides_n_ge_2": summarize(recurrent),
            "top_n_by_sequence_recurrence": summarize(top),
            "tie_complete_at_top_n_cutoff": summarize(tie_complete),
            "ranking": rank_metadata,
            "top_examples": top[:20],
        }
    return result


def analyze_workbook(xlsx: Path, *, top_n: int = 1000) -> dict[str, object]:
    if not xlsx.exists():
        raise FileNotFoundError(f"missing {xlsx}; see data/SOURCES.md for the source")
    import pandas as pd

    frame = pd.read_excel(
        xlsx,
        sheet_name=SHEET,
        usecols=["pep", "samples", "presented_by_each_sample_hla"],
    )
    result = analyze_rows(frame.to_dict(orient="records"), top_n=top_n)
    return {
        "schema_version": 1,
        "generated_by": "src/darkproteome/immunoverse_phla_recurrence.py",
        "analysis_author": ANALYSIS_AUTHOR,
        "source": {
            "path": os.path.relpath(xlsx, REPO),
            "sheet": SHEET,
            "sha256": sha256_file(xlsx),
            "size_bytes": xlsx.stat().st_size,
            "source_row_count": len(frame),
        },
        "estimand": {
            "sequence_recurrence_unit": "distinct published sample/condition labels per peptide",
            "phla_layer": "published predicted HLA binders (SB/WB), not direct allele-specific evidence",
            "primary_comparison": (
                "sequence recurrence divided by recurrence of the most widespread predicted "
                "peptide-HLA pairing"
            ),
        },
        "limitations": [
            "Sample/condition labels are not normalized to independent patients or studies.",
            "HLA annotations are binding predictions, not direct allele-specific immunoprecipitation.",
            "The analysis does not identify the malignant or nonmalignant presenting cell.",
            "The analysis does not estimate T-cell recognition, safety, or therapeutic efficacy.",
            "Duplicate peptide rows are unioned across source annotations before recurrence is counted.",
        ],
        "result": result,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-n", type=int, default=1000)
    args = parser.parse_args(argv)
    if args.top_n <= 0:
        parser.error("--top-n must be positive")

    artifact = analyze_workbook(args.xlsx, top_n=args.top_n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    result = artifact["result"]
    primary = result["sensitivities"]["strong_or_weak_predicted_binder"][
        "top_n_by_sequence_recurrence"
    ]
    strong = result["sensitivities"]["strong_predicted_binder_only"][
        "top_n_by_sequence_recurrence"
    ]
    print("IMMUNOVERSE pHLA RECURRENCE FEASIBILITY DIAGNOSTIC")
    print(f"  source rows: {artifact['source']['source_row_count']:,}")
    print(f"  unique peptides: {result['unique_peptide_sequences']:,}")
    print(f"  recurrent in >=2 sample/conditions: {result['recurrent_peptides_n_ge_2_sample_conditions']:,}")
    print(
        f"  top {args.top_n:,}, strong+weak prediction: median best-pHLA fraction "
        f"{primary['median_max_predicted_phla_fraction']:.3f}, median sequence/pHLA ratio "
        f"{primary['median_sequence_to_predicted_phla_ratio']:.3f}"
    )
    print(
        f"  top {args.top_n:,}, strong-only prediction: median best-pHLA fraction "
        f"{strong['median_max_predicted_phla_fraction']:.3f}, median sequence/pHLA ratio "
        f"{strong['median_sequence_to_predicted_phla_ratio']:.3f}"
    )
    print("  interpretation: prediction-layer sample/condition diagnostic, not observed pHLA reach")
    print(f"  wrote {os.path.relpath(args.output, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
