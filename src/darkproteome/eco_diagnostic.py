#!/usr/bin/env python3
"""
eco_diagnostic.py — the shipped guard against the ecological-inference mistake.

On real PXD055609 data (all 5 samples, T1-T5),
trying to recover the non-canonical class's true FDR by regressing a pooled per-sample
FDR `q_i` on the non-canonical class fraction `f_i` (Goodman-style ecological regression)
FAILS here: the regressand `q_i` is clamped near the target alpha by the FDR procedure
itself, which zeroes the covariance the estimator needs, so the fitted line collapses to
"both classes have FDR ~= alpha" — a confident, precise-looking, WRONG answer. This tool
makes that mistake visible and hard to make silently: given several samples/runs, it
prints (a) the DIRECTLY MEASURED per-class FDR per sample (the trustworthy number,
already free from any pepXML/PSM table via the decoy count), (b) whether those measured
rates look constant across samples (the assumption ecological regression needs and this
project's own data violates), and (c) what naive ecological regression WOULD have told
you, so the gap between (a) and (c) is explicit rather than silently trusted.

This does not replace class_decoy_ledger.py — it operates on the SAME per-class T/D counts,
just across multiple samples at once, to make the constancy question directly visible.

Stdlib only (no numpy). MIT.

Example
-------
python3 eco_diagnostic.py data/external/pxd055609_pepxml/*.pepXML --alpha 0.03
"""
import argparse, math, os, sys, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from class_decoy_ledger import load_pepxml, classify, assign_qvalues_sweep, accept_by_q, DEFAULT_DECOY_PREFIX

FOLD_FLAG_THRESHOLD = 1.5  # max/min ratio above which we call constancy "LIKELY VIOLATED"


def per_sample_rates(path, alpha):
    """One (pepXML) sample -> {f, q, theta_A (canonical), theta_B (rest)}."""
    psms = load_pepxml(path, DEFAULT_DECOY_PREFIX)  # hyperscore: higher = better (this
    for p in psms:                                  # tool is pepXML/MSFragger-only, no CLI)
        p["class"] = "canonical" if classify(p["accession"], None) == "canonical" else "noncanonical"
    assign_qvalues_sweep(psms, ascending=False)
    accepted = accept_by_q(psms, alpha)
    c = {"canonical": {"T": 0, "D": 0}, "noncanonical": {"T": 0, "D": 0}}
    for p in accepted:
        c[p["class"]]["D" if p["is_decoy"] else "T"] += 1
    T = c["canonical"]["T"] + c["noncanonical"]["T"]
    D = c["canonical"]["D"] + c["noncanonical"]["D"]
    return {
        "sample": os.path.basename(path).split("-")[0],
        "f": c["noncanonical"]["T"] / T if T else float("nan"),
        "q": D / T if T else float("nan"),
        "theta_C": c["canonical"]["D"] / c["canonical"]["T"] if c["canonical"]["T"] else float("nan"),
        "theta_N": c["noncanonical"]["D"] / c["noncanonical"]["T"] if c["noncanonical"]["T"] else float("nan"),
    }


def constancy_flag(values):
    """A sample with 0 accepted targets for a class makes its rate NaN (per_sample_rates).
    Python's min/max don't reject NaN -- they silently drop it unless it's the first element
    (order-dependent), and statistics.pstdev raises outright on NaN input. Neither is the
    right behavior here: a zero-target sample is itself strong evidence AGAINST constancy,
    not a data point to quietly ignore or crash on -- so it forces the conservative verdict,
    with an explicit reason, instead."""
    n_nan = sum(1 for v in values if math.isnan(v))
    if n_nan:
        return (float("inf"), float("nan"),
                f"LIKELY VIOLATED ({n_nan}/{len(values)} sample(s) had zero accepted targets "
                "for this class -- rate is undefined, not constant)")
    lo, hi = min(values), max(values)
    fold = (hi / lo) if lo > 0 else float("inf")
    cv = (statistics.pstdev(values) / statistics.mean(values)) if statistics.mean(values) else float("nan")
    verdict = "LIKELY VIOLATED" if fold >= FOLD_FLAG_THRESHOLD else "plausible, not clearly violated"
    return fold, cv, verdict


def ols(xs, ys):
    """Simple linear regression y ~ a + b*x. Returns (a, b, se_at_x fn)."""
    n = len(xs)
    xbar, ybar = statistics.mean(xs), statistics.mean(ys)
    Sxx = sum((x - xbar) ** 2 for x in xs)
    Sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    b = Sxy / Sxx if Sxx else float("nan")
    a = ybar - b * xbar
    resid2 = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    s2 = resid2 / max(n - 2, 1)

    def se_at(x):
        return (s2 * (1.0 / n + (x - xbar) ** 2 / Sxx)) ** 0.5 if Sxx else float("nan")

    return a, b, se_at


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pepxml", nargs="+", help="2+ MSFragger pepXML files, one per sample/run")
    ap.add_argument("--alpha", type=float, default=0.03)
    a = ap.parse_args()

    if len(a.pepxml) < 2:
        sys.exit("Need >=2 samples to test constancy across units (that's the whole diagnostic).")

    rows = [per_sample_rates(p, a.alpha) for p in a.pepxml]

    print(f"# eco_diagnostic  (alpha={a.alpha}; canonical vs rest)\n")
    print(f"{'sample':8s} {'f_i':>7s} {'q_i':>8s} {'theta_C':>9s} {'theta_N':>9s}")
    for r in rows:
        print(f"{r['sample']:8s} {r['f']:7.3f} {r['q']:8.4f} {r['theta_C']:9.4f} {r['theta_N']:9.4f}")

    print("\n=== (a) MEASURED per-class FDR — the trustworthy number, direct from decoy counts ===")
    theta_C_vals = [r["theta_C"] for r in rows]
    theta_N_vals = [r["theta_N"] for r in rows]
    print(f"  canonical:    mean={statistics.mean(theta_C_vals):.4f}  range=[{min(theta_C_vals):.4f}, {max(theta_C_vals):.4f}]")
    print(f"  noncanonical: mean={statistics.mean(theta_N_vals):.4f}  range=[{min(theta_N_vals):.4f}, {max(theta_N_vals):.4f}]")

    print("\n=== (b) constancy check — does ecological regression's key assumption hold? ===")
    for label, vals in (("theta_C (canonical)", theta_C_vals), ("theta_N (noncanonical)", theta_N_vals)):
        fold, cv, verdict = constancy_flag(vals)
        print(f"  {label:24s} fold(max/min)={fold:.2f}x  CV={cv:.1%}  constancy: {verdict}")

    print("\n=== (c) what naive ecological regression (q_i ~ a + b*f_i) WOULD tell you ===")
    f_vals = [r["f"] for r in rows]
    q_vals = [r["q"] for r in rows]
    q_fold, q_cv, _ = constancy_flag(q_vals) if min(q_vals) > 0 else (float("inf"), float("nan"), "")
    a_hat, b_hat, se_at = ols(f_vals, q_vals)
    theta_C_hat, theta_N_hat = a_hat, a_hat + b_hat
    print(f"  q_i itself: fold(max/min)={q_fold:.2f}x  CV={q_cv:.2%}  "
          f"{'(clamped near alpha -> regression has no honest signal)' if q_fold < 1.05 else ''}")
    print(f"  fitted:  theta_C_hat={theta_C_hat:.4f}  (SE ~ {se_at(0.0):.4f})")
    print(f"           theta_N_hat={theta_N_hat:.4f}  (SE ~ {se_at(1.0):.4f})")
    gap_C = abs(theta_C_hat - statistics.mean(theta_C_vals))
    gap_N = abs(theta_N_hat - statistics.mean(theta_N_vals))
    print(f"  gap vs. measured mean: canonical {gap_C:.4f}, noncanonical {gap_N:.4f}")

    print("\n=== VERDICT ===")
    any_violated = any(constancy_flag(v)[0] >= FOLD_FLAG_THRESHOLD for v in (theta_C_vals, theta_N_vals))
    if len(rows) == 2:
        # At exactly 2 samples (this tool's own enforced minimum), OLS has 0 residual
        # degrees of freedom -- any line fits 2 points exactly, so se_at() is ~0 by
        # construction, REGARDLESS of the data. That's not "the fit looks confident," it's
        # "confidence is undefined here" -- don't let it silently satisfy se_wide==False and
        # feed the DECEPTIVE-FAILURE branch below, which would misdescribe a math artifact
        # as a property of the data.
        print("  n=2 samples: OLS has 0 residual degrees of freedom, so its SE is ~0 by")
        print("  construction (any line fits 2 points exactly) -- NOT evidence of a good fit.")
        print("  Add a 3rd sample before trusting this tool's own confidence read at all.")
        if any_violated:
            print("  Measured rates are already inconstant across these 2 samples (see (b)) --")
            print("  use the MEASURED per-class rates directly; the fit above can't be evaluated.")
        else:
            print("  Measured rates look consistent across these 2 samples, but 2 points can't")
            print("  establish constancy either -- treat as inconclusive, not confirmed.")
        return
    se_wide = se_at(1.0) > 0.05
    deceptive = any_violated and not se_wide and max(gap_C, gap_N) > 0.01
    if deceptive:
        print("  DECEPTIVE FAILURE: constancy is violated, but the regression's own SE is tiny")
        print(f"  (~{se_at(1.0):.4f}) — it will LOOK highly confident while being off by "
              f"{max(gap_C, gap_N):.4f} from the measured rate.")
        print("  This is the clamped-regressand failure mode: q_i barely varies across samples, so the")
        print("  fit has almost no honest signal but reports a tiny SE anyway. DO NOT trust the SE as a")
        print("  sign of accuracy here.")
    elif any_violated or se_wide:
        print("  Constancy is violated and/or the regression's extrapolation SE is honestly wide.")
        print("  DO NOT use the ecological-regression fit above as a point estimate.")
    else:
        print("  Constancy looks plausible here — but this is still an as-reported summary-stat")
        print("  reconstruction, not a direct measurement.")
    if any_violated or se_wide or deceptive:
        print("  Use the MEASURED per-class rates in (a) directly, or report D_N per class")
        print("  (class_decoy_ledger.py) rather than reconstructing it from pooled aggregates.")


if __name__ == "__main__":
    main()
