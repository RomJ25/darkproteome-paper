"""Headline report: audit the REAL two-cohort corpus, four axes, as-reported.

Answers the project's core question on real data:
  "What fraction of published dark-proteome tumor-antigen claims even report enough to be
   auditable on all four axes, and of those, what fraction strict-pass?"

    python3 src/darkproteome/probe_report.py data/claim_catalog_real.csv
"""
import csv
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence_dimensions as ed  # noqa: E402  -- THE scorer


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def report(label, rows):
    """The reporting ladder per evidence dimension.

    There is NO "strict survivor" line any more, and reviving one would be a mistake: a joint
    pass/fail count across dimensions the record cannot even decide measures the scorer, not the
    claims. The old version printed `fail` for every MS-presented peptide -- i.e. it scored
    "nobody ran a T-cell assay" as an immunogenicity FAILURE.
    """
    n = len(rows)
    if n == 0:
        print(f"\n## {label}: 0 claims")
        return
    m = ed.matrix(rows)
    print(f"\n## {label}  (N={n})")
    for k in ed.DIMENSIONS:
        d = m[k]
        adj, lo, hi = d["adjudicable"], *wilson(d["adjudicable"], n)
        print(f"  {k:<27} asserted {d['asserted']:>5} | claim-linked {d['claim_linked']:>5} | "
              f"quantitative {d['quantitative']:>5} | ADJUDICABLE {adj:>5} "
              f"({100*adj/n:4.1f}%, CI [{100*lo:.1f}-{100*hi:.1f}]) | supports {d[ed.SUPPORTS]:>4}")


def main(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    noncanon = [r for r in rows if r.get("_canonical") != "yes"]
    canon = [r for r in rows if r.get("_canonical") == "yes"]
    hcc = [r for r in noncanon if "Hepato" in r["cancer_type"]]
    ova = [r for r in noncanon if "Ovarian" in r["cancer_type"]]

    print("=" * 78)
    print("Four-axis audit of the real dark-proteome tumor-antigen corpus")
    print("  sources: HCC Camarena/Alba 2024 (10.1126/sciadv.adn3628) + ovarian Raja 2025")
    print("=" * 78)

    report("ALL claims", rows)
    report("DARK-PROTEOME (noncanonical) claims", noncanon)
    report("  - HCC subset", hcc)
    report("  - Ovarian subset", ova)
    report("CANONICAL CTA control (MAGE/SSX etc., known-real)", canon)

    # immunogenicity reality
    iv = sum(1 for r in noncanon if r["validation_level"] == "in-vivo")
    neg = sum(1 for r in noncanon if r["validation_level"] == "validated_negative")
    tcv = sum(1 for r in noncanon if r["validation_level"] == "T-cell-validated")
    figlock = sum(1 for r in noncanon if "figure-locked" in (r["evidence_types"] or ""))
    print("\n" + "=" * 78)
    print("IMMUNOGENICITY REALITY (noncanonical):")
    print(f"  human T-cell-validated, machine-readable: {tcv}")
    print(f"  humanized-mouse reactive (in-vivo, weak): {iv}")
    print(f"  assayed NON-reactive (validated_negative): {neg}")
    print(f"  T-cell tested but result FIGURE-LOCKED (unextractable per-peptide): {figlock}")
    print("  (the ovarian paper reports ~70% of 38 reactive in Fig 6 — none recoverable as data)")
    print("=" * 78)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/claim_catalog_real.csv")
