#!/usr/bin/env python3
"""Regression test for eco_diagnostic.py's pure math (stdlib only; no pytest needed).

The pepXML pipeline itself is validated by examples/PXD055609_eco_diagnostic.txt (a
committed worked-example run against the real T1-T5 deposit). This file tests
the two pure functions (ols, constancy_flag) directly with hand-computed fixtures, the same
pattern test_class_decoy_ledger.py uses for its --psms synthetic-fixture path. Run:
python3 tests/test_eco_diagnostic.py
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src", "darkproteome"))
from eco_diagnostic import ols, constancy_flag  # noqa: E402


def main():
    fails = []

    # ols: a perfect line y = 1 + 2x through (0,1),(1,3),(2,5) -> a=1, b=2, zero residual.
    a, b, se_at = ols([0, 1, 2], [1, 3, 5])
    if abs(a - 1.0) > 1e-9 or abs(b - 2.0) > 1e-9:
        fails.append(f"  ols perfect-line: expected a=1 b=2, got a={a} b={b}")
    if se_at(0.0) > 1e-9 or se_at(1.0) > 1e-9:
        fails.append(f"  ols perfect-line: expected zero SE (zero residual), got se(0)={se_at(0.0)} se(1)={se_at(1.0)}")

    # ols: noisy case with a known residual, sanity-check SE grows away from xbar.
    a2, b2, se_at2 = ols([0, 1, 2, 3], [1.0, 1.9, 3.2, 2.8])
    if not (se_at2(0.0) > 0 and se_at2(10.0) > se_at2(1.5)):
        fails.append("  ols noisy-case: expected SE > 0 and SE to grow further from xbar (extrapolation)")

    # constancy_flag: a clearly-varying series (4x fold) must flag LIKELY VIOLATED.
    fold, cv, verdict = constancy_flag([1.0, 2.0, 4.0])
    if abs(fold - 4.0) > 1e-9 or verdict != "LIKELY VIOLATED":
        fails.append(f"  constancy_flag varying: expected fold=4.0 LIKELY VIOLATED, got fold={fold} verdict={verdict!r}")

    # constancy_flag: a near-constant series must NOT flag.
    fold2, cv2, verdict2 = constancy_flag([1.0, 1.01, 0.99])
    if verdict2 != "plausible, not clearly violated":
        fails.append(f"  constancy_flag near-constant: expected 'plausible...', got {verdict2!r}")

    # the real PXD055609 T1-T5 numbers reproduce the documented failure mode.
    theta_C = [0.0178, 0.0077, 0.0129, 0.0104, 0.0071]
    f_vals = [0.759, 0.688, 0.745, 0.768, 0.788]
    q_vals = [0.0300, 0.0300, 0.0300, 0.0300, 0.0300]
    fold_C, _, verdict_C = constancy_flag(theta_C)
    if verdict_C != "LIKELY VIOLATED":
        fails.append(f"  real-data theta_C: expected LIKELY VIOLATED (known 2.5x spread), got {verdict_C!r}")
    a3, b3, se_at3 = ols(f_vals, q_vals)
    if abs(a3 - 0.0297) > 5e-4 or abs((a3 + b3) - 0.0301) > 5e-4:
        fails.append(f"  real-data ols: expected theta_C_hat~0.0297 theta_N_hat~0.0301, got {a3:.4f}/{a3+b3:.4f}")

    # constancy_flag: a sample with 0 accepted targets for a class makes its rate NaN
    # (per_sample_rates). Python's min/max would silently drop NaN unless it was the first
    # element (order-dependent), and statistics.pstdev would crash outright with ValueError --
    # neither is the right behavior for a tool whose whole purpose is "make failure visible,
    # not silent." Guards against both: a conservative, explicit, position-independent verdict.
    nan = float("nan")
    for label, vals in (("NaN first", [nan, 0.03, 0.03]), ("NaN middle", [0.03, nan, 0.03])):
        fold_n, cv_n, verdict_n = constancy_flag(vals)
        if fold_n != float("inf") or "zero accepted targets" not in verdict_n:
            fails.append(f"  constancy_flag {label}: expected inf fold + an explicit zero-target "
                         f"verdict, got fold={fold_n} verdict={verdict_n!r}")

    # ols: at exactly n=2 (this tool's CLI-enforced minimum), OLS has 0 residual degrees of
    # freedom, so se_at() is ~0 by construction regardless of the data -- verified with two
    # points that are NOT close together, so a nonzero SE would be the naive expectation.
    _, _, se_at_n2 = ols([0.1, 0.9], [0.02, 0.09])
    if se_at_n2(1.0) > 1e-9:
        fails.append(f"  ols n=2: expected ~0 SE (0 residual DOF), got se(1.0)={se_at_n2(1.0)}")

    if fails:
        print("FAIL — eco_diagnostic:\n" + "\n".join(fails)); sys.exit(1)
    print("PASS — ols (perfect-line exact, noisy-case SE grows with extrapolation distance, "
          "n=2 SE~0 by construction), constancy_flag (4x-fold -> LIKELY VIOLATED, near-constant "
          "-> not violated, NaN -> explicit zero-target verdict regardless of position), "
          "and the real PXD055609 T1-T5 numbers reproduce the documented deceptive-failure mode "
          "(theta_C 2.5x fold VIOLATED; OLS fit ~0.0297/0.0301).")


if __name__ == "__main__":
    main()
