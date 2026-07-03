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
from axes import AXES, score_all, is_strict_survivor  # noqa: E402


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def fully_auditable(v):
    return all(v[a] != "unverifiable" for a in AXES)


def report(label, rows):
    n = len(rows)
    if n == 0:
        print(f"\n## {label}: 0 claims")
        return
    verdicts = [score_all(r) for r in rows]
    print(f"\n## {label}  (N={n})")
    for a in AXES:
        c = Counter(v[a] for v in verdicts)
        sp, wp, fl, un = c["strict-pass"], c["weak-pass"], c["fail"], c["unverifiable"]
        print(f"  {a:<18} strict {sp:>4} ({100*sp/n:4.1f}%) | weak {wp:>4} | fail {fl:>4} | unverif {un:>4} ({100*un/n:4.1f}%)")
    fa = sum(fully_auditable(v) for v in verdicts)
    flo, fhi = wilson(fa, n)
    surv = sum(is_strict_survivor(v) for v in verdicts)
    slo, shi = wilson(surv, n)
    print(f"  --> fully auditable (all 4 axes decidable): {fa}/{n} = {100*fa/n:.1f}%  "
          f"(95% CI [{100*flo:.1f}%, {100*fhi:.1f}%])")
    print(f"  --> STRICT SURVIVORS (strict-pass all 4):    {surv}/{n} = {100*surv/n:.2f}%  "
          f"(95% CI [{100*slo:.1f}%, {100*shi:.1f}%])")


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
