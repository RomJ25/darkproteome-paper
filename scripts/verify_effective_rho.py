"""
Monte Carlo verification of the manuscript's effective-rho identity (Methods):
  rho_eff(tau) = sum_s w_s(tau) pi_sN / sum_s w_s(tau) pi_sC,
  w_s = 1 - F(tau)^{m_s},  pi_sg = m_sg/m_s  -- NOT the global ratio rho_global.

Model (no special functions): at per-candidate tail prob t = 1-F(tau) (small t = stringent
threshold), a spectrum with m_s iid-null candidates yields an accepted decoy top-hit with
prob w_s = 1-(1-t)^{m_s}; given acceptance the winner's class is N w.p. pi_sN (argmax-class
is independent of the max value under a common iid null). Monte-Carlo D_N(t), D_C(t) and
compare empirical D_N/D_C to rho_eff(t) and to rho_global = sum m_sN / sum m_sC.

RESULT (seed 0): (A) the two-spectrum worked example in the manuscript reproduces exactly
(rho_eff=0.2132 / 4.6902); (B) empirical D_N/D_C tracks rho_eff(t) to 3 decimals across
thresholds and both drift off rho_global. Confirms the identity used in Results II.
"""
import numpy as np
rng = np.random.default_rng(0)

def rho_eff(t, mN, mC):
    m = mN + mC
    w = 1.0 - (1.0 - t)**m
    return (w*(mN/m)).sum() / (w*(mC/m)).sum()

def simulate(t, mN, mC, reps=400):
    m = mN+mC; piN = mN/m; DN=DC=0
    for _ in range(reps):
        accept = rng.random(len(m)) < (1.0-(1.0-t)**m)
        isN    = rng.random(len(m)) < piN
        DN += int(np.count_nonzero(accept &  isN))
        DC += int(np.count_nonzero(accept & ~isN))
    return DN, DC

if __name__ == "__main__":
    print("=== (A) two-spectrum counterexample (manuscript's worked example) ===")
    mN = np.array([900]+[1]*100); mC = np.array([100]+[9]*100)
    print("  rho_global               =", mN.sum()/mC.sum())
    print("  rho_eff(t=0.01)          =", round(rho_eff(0.01,mN,mC),4), " (manuscript: ~0.213)")
    print("  rho_eff(t=0.01) swapped  =", round(rho_eff(0.01,mC,mN),4), " (manuscript: ~4.69)")

    print("\n=== (B) heterogeneous regime: does emp D_N/D_C track rho_eff(t)? ===")
    S_low,S_high = 4000,200
    mN = np.concatenate([rng.integers(1,3,S_low),   rng.integers(400,1200,S_high)])
    mC = np.concatenate([rng.integers(5,20,S_low),  rng.integers(20,120,S_high)])
    rg = mN.sum()/mC.sum(); print(f"  rho_global = {rg:.3f}")
    print(f"  {'t=1-F':>10} {'rho_eff':>9} {'emp DN/DC':>10} {'vs global':>10}")
    for t in [1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1]:
        re = rho_eff(t,mN,mC); DN,DC = simulate(t,mN,mC)
        print(f"  {t:>10.0e} {re:>9.3f} {DN/max(DC,1):>10.3f} {re/rg:>9.2f}x")
