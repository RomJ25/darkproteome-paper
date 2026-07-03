#!/usr/bin/env python3
"""
psm_multiplicity_probe.py — test whether the class-decoy ledger's class-FDR understates the
cryptic (noncanonical) class's true error because cryptic peptides are disproportionately
single-PSM (low replicate depth), and single-PSM IDs are the least reliable stratum.

Tests three falsifiable predictions about a stratification gap in class_decoy_ledger.py:
it reports one class-FDR per class, blind to PSM-replicate-depth within that class.
  (i)   noncanonical (cryptic) peptides have a markedly lower median PSM count than canonical
  (ii)  unique-peptide-level noncanonical class-FDR exceeds PSM-level noncanonical class-FDR
  (iii) the excess concentrates in the n=1 (single-PSM) stratum

Reuses class_decoy_ledger.py's pepXML loader + q-value sweep + classify() so results are
directly comparable to the shipped ledger tool's own numbers (verified against the T5 baseline
in examples/PXD055609_T5_ledger.json: canonical 0.79%, noncanonical 3.95%, 93.0% canonical-self).

Usage:
    python3 psm_multiplicity_probe.py PEPXML [PEPXML ...] [--alpha 0.03]
"""
import sys, os, argparse, statistics
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from class_decoy_ledger import load_pepxml, classify, assign_qvalues_sweep, accept_by_q, DEFAULT_DECOY_PREFIX


def load_accepted(pepxml_path, alpha):
    psms = load_pepxml(pepxml_path, DEFAULT_DECOY_PREFIX)  # hyperscore: higher = better
    for p in psms:
        p["class"] = classify(p["accession"], None)
        p["sample"] = os.path.basename(pepxml_path).split("-")[0]
    assign_qvalues_sweep(psms, ascending=False)
    return accept_by_q(psms, alpha)


def fdr(T, D):
    return (D + 1) / max(T, 1)


def bucket(n):
    if n == 1:
        return "n=1"
    if n == 2:
        return "n=2"
    return "n>=3"


def analyze(accepted):
    """Returns (psm_level, peptide_level, mult_target, mult_decoy) keyed by class."""
    psm_level = defaultdict(lambda: {"T": 0, "D": 0})
    pep_groups = defaultdict(int)  # (class, peptide, is_decoy) -> PSM count
    for p in accepted:
        psm_level[p["class"]]["D" if p["is_decoy"] else "T"] += 1
        pep_groups[(p["class"], p["peptide"], p["is_decoy"])] += 1

    peptide_level = defaultdict(lambda: {"T": 0, "D": 0})
    mult_target = defaultdict(list)
    mult_decoy = defaultdict(list)
    for (cls, pep, is_dec), n_psm in pep_groups.items():
        peptide_level[cls]["D" if is_dec else "T"] += 1
        (mult_decoy if is_dec else mult_target)[cls].append(n_psm)
    return psm_level, peptide_level, mult_target, mult_decoy


def stratified_counts(mult_target_cls, mult_decoy_cls):
    buckets = defaultdict(lambda: {"T": 0, "D": 0})
    for n in mult_target_cls:
        buckets[bucket(n)]["T"] += 1
    for n in mult_decoy_cls:
        buckets[bucket(n)]["D"] += 1
    return buckets


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pepxml", nargs="+")
    ap.add_argument("--alpha", type=float, default=0.03)
    a = ap.parse_args()

    all_accepted = []
    print(f"# PSM-multiplicity probe  (alpha={a.alpha})\n")
    for path in a.pepxml:
        accepted = load_accepted(path, a.alpha)
        all_accepted.extend(accepted)
        print(f"  {os.path.basename(path)}: {len(accepted)} accepted PSMs")
    print()

    psm_level, peptide_level, mult_target, mult_decoy = analyze(all_accepted)

    print("=== (i) median PSM-multiplicity: canonical vs noncanonical (target peptides only) ===")
    for cls in ("canonical", "noncanonical", "variant"):
        xs = mult_target.get(cls, [])
        if xs:
            med = statistics.median(xs)
            mean = statistics.mean(xs)
            frac_singleton = sum(1 for x in xs if x == 1) / len(xs)
            print(f"  {cls:14s} n_peptides={len(xs):6d}  median_PSM/peptide={med:.1f}  "
                  f"mean={mean:.2f}  %singleton(n=1)={frac_singleton:.1%}")

    print("\n=== (ii) PSM-level vs unique-peptide-level class-FDR ===")
    for cls in ("canonical", "noncanonical", "variant"):
        psm = psm_level.get(cls, {"T": 0, "D": 0})
        pep = peptide_level.get(cls, {"T": 0, "D": 0})
        print(f"  {cls:14s} PSM-level:    T={psm['T']:6d} D={psm['D']:5d}  FDR=(D+1)/T={fdr(psm['T'], psm['D']):.2%}")
        print(f"  {cls:14s} peptide-level: T={pep['T']:6d} D={pep['D']:5d}  FDR=(D+1)/T={fdr(pep['T'], pep['D']):.2%}")
        delta = fdr(pep["T"], pep["D"]) - fdr(psm["T"], psm["D"])
        print(f"  {cls:14s} peptide-level minus PSM-level FDR = {delta:+.2%}")

    print("\n=== (iii) class-FDR stratified by PSM-multiplicity bucket, ALL classes (fair comparison) ===")
    for cls in ("canonical", "noncanonical", "variant"):
        print(f"  -- {cls} --")
        buckets = stratified_counts(mult_target.get(cls, []), mult_decoy.get(cls, []))
        for b in ("n=1", "n=2", "n>=3"):
            c = buckets.get(b, {"T": 0, "D": 0})
            print(f"    {b:6s}  T={c['T']:6d}  D={c['D']:5d}  class-FDR=(D+1)/T={fdr(c['T'], c['D']):.2%}")

    print("\n=== sanity check against known T5-only ledger baseline (if T5 alone was passed) ===")
    if len(a.pepxml) == 1 and "T5" in a.pepxml[0]:
        c = psm_level.get("noncanonical", {"T": 0, "D": 0})
        print(f"  expect T=3822 D=150 (3.95%) -> got T={c['T']} D={c['D']} ({fdr(c['T'], c['D']):.2%})")


if __name__ == "__main__":
    main()
