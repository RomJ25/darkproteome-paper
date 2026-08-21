"""Crosswalk IEAtlas cancer sequences to PXD038782 benign processed peptide tables.

This is deliberately a staged external recurrence analysis, not a raw-spectrum re-search.
PXD038782 is a partial PRIDE submission whose processed tables came from a canonical-only
PEAKS DB search and do not expose donor HLA alleles or USIs. The artifact therefore uses
``processed-table recurrence`` language and records those limitations explicitly.

Run: ``python3 src/darkproteome/pxd038782_benign_crosswalk.py``
"""
import csv
import io
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

PROJECT = "PXD038782"
API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{PROJECT}"
FILES_API = API + "/files"
OUT = os.path.join(paths.REPO, "data", "derived_pxd038782_benign_crosswalk.json")
SCORED = os.path.join(paths.REPO, "data", "claim_catalog_scored.csv")
MOD = re.compile(r"\([^)]*\)")
USER_AGENT = "darkproteome-provenance-audit/1.0"


def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as response:
        return json.load(response)


def peptide_set(path):
    values = set()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if row:
                peptide = row[0].strip().upper()
                if peptide.isalpha():
                    values.add(peptide)
    return values


def https_location(file_record):
    for loc in file_record["publicFileLocations"]:
        value = loc.get("value", "")
        if "pride/data/archive/" in value:
            suffix = value.split("pride/data/archive/", 1)[1]
            return "https://ftp.pride.ebi.ac.uk/pride/data/archive/" + urllib.parse.quote(
                suffix, safe="/"
            )
    raise ValueError(f"no public archive location for {file_record['fileName']}")


def file_labels(name):
    parts = name.split("_")
    benign_i = parts.index("benign")
    tissue = parts[benign_i + 2]
    if "W6-32" in name:
        hla_class = "I"
    elif "Tue39L243" in name:
        hla_class = "II"
    else:
        hla_class = "unresolved"
    return tissue, hla_class


def main():
    paths.require(paths.IEATLAS_CANCER, paths.IEATLAS_NORMAL)
    if not os.path.exists(SCORED):
        sys.exit(f"missing {SCORED}")

    project = fetch_json(API)
    files = fetch_json(FILES_API)
    benign = sorted(
        (f for f in files
         if f["fileCategory"]["value"] == "SEARCH" and "_benign_" in f["fileName"]),
        key=lambda f: f["fileName"],
    )
    cancer = peptide_set(paths.IEATLAS_CANCER)
    normal = peptide_set(paths.IEATLAS_NORMAL)
    exact = {
        r["peptide"].strip().upper()
        for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))
        if r["canonical_self"] == "1"
    }

    evidence = defaultdict(lambda: {"tissues": set(), "classes": set(), "files": set()})
    valid_rows = 0
    found_by = set()
    file_provenance = []
    for i, record in enumerate(benign, 1):
        name = record["fileName"]
        tissue, hla_class = file_labels(name)
        url = https_location(record)
        print(f"[{i:02d}/{len(benign)}] {name}")
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req) as response:
            text = io.TextIOWrapper(response, encoding="utf-8-sig", errors="replace", newline="")
            for row in csv.DictReader(text):
                peptide = MOD.sub("", (row.get("Peptide") or "")).strip().upper()
                if not peptide.isalpha():
                    continue
                valid_rows += 1
                if row.get("Found By"):
                    found_by.add(row["Found By"].strip())
                if peptide in cancer:
                    evidence[peptide]["tissues"].add(tissue)
                    evidence[peptide]["classes"].add(hla_class)
                    evidence[peptide]["files"].add(name)
        file_provenance.append({
            "file_name": name,
            "size_bytes": record["fileSizeBytes"],
            "sha1_from_pride": record["checksum"],
            "tissue_label": tissue,
            "preparation_class": hla_class,
        })

    matches = set(evidence)
    class_i = {p for p, e in evidence.items() if "I" in e["classes"]}
    class_ii = {p for p, e in evidence.items() if "II" in e["classes"]}
    breadth = [len(e["tissues"]) for e in evidence.values()]
    artifact = {
        "project": PROJECT,
        "project_title": project["title"],
        "submission_type": project["submissionType"],
        "submission_date": project["submissionDate"],
        "publication_date": project["publicationDate"],
        "license": project["license"],
        "source_apis": [API, FILES_API],
        "analysis_unit": "unique unmodified peptide sequence in deposited benign processed tables",
        "n_benign_search_tables": len(benign),
        "n_valid_processed_rows": valid_rows,
        "found_by_values": sorted(found_by),
        "n_ieatlas_cancer_sequences_recurrent": len(matches),
        "n_recurrent_exact_canonical_compatible": len(matches & exact),
        "n_recurrent_already_ieatlas_normal": len(matches & normal),
        "n_recurrent_absent_ieatlas_normal": len(matches - normal),
        "preparation_class_counts": {
            "class_I": len(class_i),
            "class_II": len(class_ii),
            "both": len(class_i & class_ii),
        },
        "tissue_breadth": {
            "n_unique_labels": len({t for e in evidence.values() for t in e["tissues"]}),
            "n_at_least_2": sum(x >= 2 for x in breadth),
            "n_at_least_5": sum(x >= 5 for x in breadth),
            "maximum": max(breadth) if breadth else 0,
        },
        "limitations": [
            "processed-table crosswalk, not a raw-spectrum re-search or USI-level validation",
            "PXD038782 is a PARTIAL submission",
            "deposited results derive from a canonical-proteome PEAKS DB search, not joint canonical/ncORF competition",
            "preparation antibody supports a class-level label, but donor HLA alleles are unavailable here",
            "recurrence establishes sequence-list agreement, not source-locus provenance",
        ],
        "files": file_provenance,
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(artifact, fh, indent=2)
        fh.write("\n")
    print(json.dumps({k: v for k, v in artifact.items() if k != "files"}, indent=2))
    print(f"wrote {os.path.relpath(OUT, paths.REPO)}")


if __name__ == "__main__":
    main()
