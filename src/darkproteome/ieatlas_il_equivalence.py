"""Measure exact and mass-spectrometry-equivalent I/L canonical compatibility.

Routine tandem mass spectrometry does not distinguish leucine from isoleucine by
precursor or fragment mass. The manuscript's literal exact-string result remains
the primary estimand; this script adds the nested sensitivity analysis obtained by
collapsing both residues to ``J`` in the peptide and canonical reference.

Run: ``python3 src/darkproteome/ieatlas_il_equivalence.py``
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

SCORED = os.path.join(paths.REPO, "data", "claim_catalog_scored.csv")
OUT = os.path.join(paths.REPO, "data", "derived_ieatlas_il_equivalence.json")


def collapse_il(sequence):
    return sequence.replace("I", "J").replace("L", "J")


def il_hits(queries):
    """Return queries matching the reference after I/L collapse in one FASTA pass."""
    collapsed_to_original = {}
    for peptide in queries:
        collapsed_to_original.setdefault(collapse_il(peptide), []).append(peptide)
    lengths = sorted({len(p) for p in queries})
    hits, buf = set(), []

    def scan(sequence):
        sequence = collapse_il(sequence)
        for length in lengths:
            for start in range(len(sequence) - length + 1):
                originals = collapsed_to_original.get(sequence[start:start + length])
                if originals:
                    hits.update(originals)

    with open(paths.SPROT, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if buf:
                    scan("".join(buf))
                buf = []
            else:
                buf.append(line.strip().upper())
        if buf:
            scan("".join(buf))
    return hits


def main():
    paths.require(paths.SPROT, paths.IEATLAS_CANCER)
    if not os.path.exists(SCORED):
        sys.exit(f"missing {SCORED}")

    selfmap = {
        r["peptide"].strip().upper(): int(r["canonical_self"])
        for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))
    }
    peptides = set()
    with open(paths.IEATLAS_CANCER, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if row:
                peptide = row[0].strip().upper()
                if peptide.isalpha():
                    peptides.add(peptide)
    exact = {peptide for peptide in peptides if selfmap.get(peptide) == 1}
    il_compatible = exact | il_hits(set(peptides) - exact)

    n = len(peptides)
    added = il_compatible - exact
    artifact = {
        "reference": "Swiss-Prot reviewed human reference R used by claim_catalog_scored.csv",
        "rule": "replace both I and L by J before substring matching",
        "n_total": n,
        "n_exact_compatible": len(exact),
        "pct_exact_compatible": round(100 * len(exact) / n, 1),
        "n_il_equivalent_compatible": len(il_compatible),
        "pct_il_equivalent_compatible": round(100 * len(il_compatible) / n, 1),
        "n_additional_il_equivalent": len(added),
        "delta_percentage_points": round(100 * len(added) / n, 1),
        "interpretation": (
            "Nested sensitivity analysis for sequence compatibility under routine MS measurement; "
            "it does not identify the source locus."
        ),
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(artifact, fh, indent=2)
        fh.write("\n")
    print(json.dumps(artifact, indent=2))
    print(f"wrote {os.path.relpath(OUT, paths.REPO)}")


if __name__ == "__main__":
    main()
