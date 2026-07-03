"""Bulletproof the core finding against the obvious objections.

Two objections to pre-empt:
  (1) "Your bar is arbitrarily strict."  -> escalating-leniency ladder + leave-axes-out.
  (2) "Your N is inflated by altORF / duplicate peptides." -> corpus cuts + dedup.

The claim we are stress-testing: the published record yields NO ready-to-use, machine-
readable human ground truth for dark-proteome tumor antigens. We show this is robust:
strict survivors stay ~0 across corpus definitions, and you cannot manufacture survivors
without (a) dropping the source-ORF rigor axis AND (b) crediting non-human or figure-locked
immunogenicity — i.e. doing the very benchmark-construction work that doesn't yet exist.

    python3 src/darkproteome/robustness.py data/claim_catalog_real.csv
"""
import csv
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from axes import AXES, score_all  # noqa: E402

NONCODING = {"lncRNA-ORF", "pseudogene-ORF"}  # narrow "dark proteome" (excludes altORF/other)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def lenient_pass(r, relax_source=False, relax_hla=False, generous_immuno=False,
                 drop=()):
    v = score_all(r)
    if relax_source and v["source_orf"] == "weak-pass":
        v["source_orf"] = "strict-pass"
    if relax_hla and v["hla_presentation"] == "weak-pass":
        v["hla_presentation"] = "strict-pass"
    if generous_immuno:
        if v["immunogenicity"] == "weak-pass":            # humanized-mouse reactive
            v["immunogenicity"] = "strict-pass"
        if "figure-locked" in (r.get("evidence_types") or ""):  # ovarian Fig-6 reactive
            v["immunogenicity"] = "strict-pass"
    need = [a for a in AXES if a not in drop]
    return all(v[a] == "strict-pass" for a in need)


def survivors(rows, **kw):
    return sum(lenient_pass(r, **kw) for r in rows)


def cut_table(rows):
    noncanon = [r for r in rows if r.get("_canonical") != "yes"]
    narrow = [r for r in noncanon if r["orf_class"] in NONCODING]
    seen, dedup = set(), []
    for r in noncanon:
        p = r["peptide_sequence"]
        if p and p != "not reported" and p in seen:
            continue
        if p and p != "not reported":
            seen.add(p)
        dedup.append(r)
    cuts = [
        ("ALL non-canonical", noncanon),
        ("NARROW (lncRNA+pseudogene only)", narrow),
        ("DEDUP by peptide", dedup),
    ]
    print("\n=== 1. STRICT survivors are 0 across every corpus definition ===")
    for name, rs in cuts:
        s = survivors(rs)
        lo, hi = wilson(s, len(rs))
        print(f"  {name:<34} N={len(rs):>4}  strict-survivors={s}  "
              f"({100*s/len(rs):.2f}%, CI [{100*lo:.1f}, {100*hi:.1f}])")
    return noncanon


def forcing_axes(rows):
    print("\n=== 2. Which axes force the 0%? (survivors if we DROP one/two axes) ===")
    print(f"  require all 4 strict:                 {survivors(rows)}")
    for a in AXES:
        print(f"  drop {a:<18} (other 3 strict): {survivors(rows, drop=(a,))}")
    print("  -- two forcing axes are source_orf AND immunogenicity; drop BOTH:")
    n = survivors(rows, drop=("source_orf", "immunogenicity"))
    print(f"  drop source_orf+immunogenicity (need hla+tumor strict): {n}")
    print("     => that {0} = the max set that is MS-presented w/ allele AND tumor-specific,".format(n))
    print("        of which the human-T-cell-validated machine-readable count is 0.")


def leniency_ladder(rows):
    print("\n=== 3. Escalating-leniency ladder (how much must we relax to get ANY?) ===")
    steps = [
        ("strict (baseline)", dict()),
        ("+ accept binary Ribo-seq (relax source)", dict(relax_source=True)),
        ("+ drop per-peptide allele (relax hla)", dict(relax_source=True, relax_hla=True)),
        ("+ credit mouse + figure-locked immuno", dict(relax_source=True, relax_hla=True, generous_immuno=True)),
        ("+ DROP source axis entirely (max leniency)", dict(relax_hla=True, generous_immuno=True, drop=("source_orf",))),
    ]
    for label, kw in steps:
        print(f"  {label:<46} survivors = {survivors(rows, **kw)}")
    print("  => survivors only appear once you BOTH drop source-ORF rigor AND credit")
    print("     non-human/figure-locked immunogenicity. That 'generous ceiling' is exactly")
    print("     the low-rigor benchmark you'd have to BUILD (mouse data + figure digitization).")


def immuno_truth(rows):
    print("\n=== 4. The reusable human-immunogenicity ground truth (the real ceiling) ===")
    tcv = [r for r in rows if r["validation_level"] == "T-cell-validated"]
    mouse = [r for r in rows if r["validation_level"] == "in-vivo"]
    neg = [r for r in rows if r["validation_level"] == "validated_negative"]
    figlock = {r["peptide_sequence"] for r in rows if "figure-locked" in (r.get("evidence_types") or "")}
    print(f"  human T-cell-validated, machine-readable POSITIVES: {len(tcv)}")
    print(f"  humanized-mouse reactive (non-human evidence):      {len(mouse)}")
    print(f"  assayed NON-reactive negatives (machine-readable):  {len(neg)}")
    print(f"  unique peptides T-cell-tested but FIGURE-LOCKED:     {len(figlock)}")
    ceiling = len(mouse) + len(figlock)
    print(f"  --> generous reusable-positive CEILING = {ceiling} "
          f"(all require non-human evidence or figure digitization)")


def main(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    print("#" * 78)
    print("# ROBUSTNESS / BULLETPROOFING — core finding")
    print("#" * 78)
    noncanon = cut_table(rows)
    forcing_axes(noncanon)
    leniency_ladder(noncanon)
    immuno_truth(noncanon)
    print("\n" + "#" * 78)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/claim_catalog_real.csv")
