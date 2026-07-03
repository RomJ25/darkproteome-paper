"""Run the four-axis evidence audit over a Claim Catalog CSV.

    python3 src/darkproteome/audit.py data/sample_claims.csv

Scores every claim on four INDEPENDENT axes (source ORF validity, HLA presentation,
tumor specificity, immunogenicity) as reported, and reports the headline:

    strict survivor fraction = claims that strict-pass ALL FOUR axes / total claims

with a Wilson 95% confidence interval, plus per-axis survival so you can see exactly
where claims die.
"""

import csv
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from axes import AXES, score_all, is_strict_survivor  # noqa: E402
from schema import validate_row  # noqa: E402

_CODE = {"strict-pass": "P", "weak-pass": "w", "fail": "F", "unverifiable": "?"}


def load_catalog(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def deduplicate(rows):
    """Distinct-evidence count: rows sharing a raw dataset accession collapse to one."""
    seen, distinct = set(), 0
    for r in rows:
        acc = (r.get("underlying_dataset_accession") or "").strip().lower()
        if acc and acc not in ("not reported", "na", "n/a", ""):
            if acc in seen:
                continue
            seen.add(acc)
        distinct += 1
    return distinct


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def run(path):
    rows = load_catalog(path)
    n = len(rows)
    invalid = [(i, p) for i, r in enumerate(rows, 1) if (p := validate_row(r))]

    verdicts = [score_all(r) for r in rows]
    survivors = [is_strict_survivor(v) for v in verdicts]
    axis_pass = {a: sum(1 for v in verdicts if v[a] == "strict-pass") for a in AXES}
    n_surv = sum(survivors)
    lo, hi = wilson_ci(n_surv, n)

    print(f"\n=== darkproteome four-axis audit: {os.path.basename(path)} ===\n")
    print(f"raw claims:               {n}")
    print(f"deduplicated (by raw MS): {deduplicate(rows)}")
    print(f"schema-invalid rows:      {len(invalid)}")

    print("\nper-axis strict-pass survival:")
    for a in AXES:
        print(f"  {a:<18} {axis_pass[a]:>3}/{n}  ({100*axis_pass[a]/n:.0f}%)")

    print("\n>>> STRICT SURVIVOR FRACTION (passes all four axes):")
    print(f"    {n_surv}/{n} = {100*n_surv/n:.0f}%   "
          f"(Wilson 95% CI [{100*lo:.0f}%, {100*hi:.0f}%])")

    print("\nper-claim verdicts  (P=strict-pass  w=weak  F=fail  ?=unverifiable):")
    print(f"  {'claim':<26} {'orf':>4} {'hla':>4} {'tum':>4} {'imm':>4}  survivor")
    for r, v, s in list(zip(rows, verdicts, survivors))[:20]:
        name = (r.get("peptide_sequence") or "?")[:26]
        cells = "".join(f"{_CODE[v[a]]:>5}" for a in AXES)
        print(f"  {name:<26}{cells}   {'YES' if s else ''}")
    if n > 20:
        print(f"  ... ({n - 20} more rows)")

    if invalid:
        print(f"\n!! {len(invalid)} schema-invalid row(s):")
        for i, probs in invalid[:10]:
            print(f"   row {i}: {'; '.join(probs)}")
    print()
    return {"n": n, "survivors": n_surv, "axis_pass": axis_pass}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    run(sys.argv[1])
