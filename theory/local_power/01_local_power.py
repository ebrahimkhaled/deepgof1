"""01_local_power.py -- W1, part 1: evaluate the limit formulas of Theorem L (THEOREM.md).

Under pi_n = pi_0 + h / sqrt(n) the deployed DeepGOF-1 (B = 199, k = 10) has limit power

    beta(h) = E[ Bin(k-1; B, 1 - F0(S~(delta(h) + Z))) ],   Z ~ N(0, Sigma),
    delta(h) = D^-1/2 (a(h) - G I^-1 b(h)),   Sigma = I - D^-1/2 G I^-1 G' D^-1/2.

All population integrals are computed on an exact product grid in probability space (every cell gets an
exact sub-grid). The rivals are evaluated on the SAME drift: Hosmer-Lemeshow (a Gaussian quadratic form
over the deciles of eta_0) and Stukel (noncentral chi-square, 2 df).

For each departure the script finds tau_50, the amplitude at which each test's limit power is .5, and
writes the limit powers of the three tests at the DeepGOF-1 tau_50 to FORMULA.json for 02_simulate.py.

Design: x1, x2 ~ N(0, 1) independent; eta_0 = -0.5 + x1 + 0.6 x2.

Run:  python -u 01_local_power.py
"""
import os, json
import numpy as np
from scipy.special import ndtri, expit
from scipy.stats import chi2, binom, ncx2
from lp_core import score, check_bench, gfun, K, B, KREJ, ALPHA, BETA0, DEPS

HERE = os.path.dirname(os.path.abspath(__file__))
print("shipped bench reproduced: %.1e" % check_bench())

# ---- 1. population quantities -----------------------------------------------------------------
G_ = 1200
u = (np.arange(G_) + 0.5) / G_
xg = ndtri(u)
bins = np.minimum((K * u).astype(int), K - 1)
X1, X2 = np.meshgrid(xg, xg, indexing="ij"); X1, X2 = X1.ravel(), X2.ravel()
CELL = (bins[:, None] * K + bins[None, :]).ravel()
Xt = np.column_stack([np.ones_like(X1), X1, X2])
eta0 = Xt @ BETA0
pi0 = expit(eta0); v = pi0 * (1 - pi0)
Npt = len(X1)

def cellsum(idx, w, nc):
    return np.bincount(idx, weights=w, minlength=nc) / Npt

Dv = cellsum(CELL, v, 36)
Gm = np.column_stack([cellsum(CELL, v * Xt[:, j], 36) for j in range(3)])
Iinf = (Xt * v[:, None]).T @ Xt / Npt
Iinv = np.linalg.inv(Iinf)
Dm12 = 1 / np.sqrt(Dv)
Sigma = np.eye(36) - (Dm12[:, None] * Gm) @ Iinv @ (Gm.T * Dm12[None, :])
ev, evec = np.linalg.eigh(Sigma)
print("Sigma: eigenvalues in [%.2e, %.4f]; rank %d (35 expected: only the intercept direction is null)"
      % (ev.min(), ev.max(), int((ev > 1e-8).sum())))
Sroot = evec * np.sqrt(np.clip(ev, 0, None))

def drift_parts(name):
    h = v * gfun(name, X1, X2)                         # tau = 1
    return h, Xt.T @ h / Npt

def delta_unit(name):
    h, b = drift_parts(name)
    return Dm12 * (cellsum(CELL, h, 36) - Gm @ Iinv @ b)

# Hosmer-Lemeshow: population deciles of eta_0 ~ N(-0.5, 1.36)
HLG = np.searchsorted(-0.5 + np.sqrt(1.36) * ndtri(np.arange(1, 10) / 10), eta0)
Dhl = cellsum(HLG, v, 10)
Ghl = np.column_stack([cellsum(HLG, v * Xt[:, j], 10) for j in range(3)])
Phl = cellsum(HLG, np.ones(Npt), 10); pbar = cellsum(HLG, pi0, 10) / Phl
whl = Phl * pbar * (1 - pbar)
e2, v2 = np.linalg.eigh(np.diag(Dhl) - Ghl @ Iinv @ Ghl.T)
Rhl = v2 * np.sqrt(np.clip(e2, 0, None))

# Stukel: the two added covariates at eta_0
Zs = np.column_stack([0.5 * eta0 ** 2 * (eta0 >= 0), -0.5 * eta0 ** 2 * (eta0 < 0)])
Ezx = (Zs * v[:, None]).T @ Xt / Npt
Oinv = np.linalg.inv((Zs * v[:, None]).T @ Zs / Npt - Ezx @ Iinv @ Ezx.T)

# ---- 2. limit power of each test ------------------------------------------------------------------
rngF = np.random.default_rng(90210)
NF = 300000
F0 = np.sort(score(rngF.normal(size=(NF, 36)) @ Sroot.T))
ZF = rngF.normal(size=(NF, 36)) @ Sroot.T
ZH = rngF.normal(size=(NF, 10)) @ Rhl.T

def power_dg(d):
    W_ = 1 - np.searchsorted(F0, score(ZF + d[None, :]), side="right") / NF
    return float(binom.cdf(KREJ - 1, B, W_).mean())

def power_hl(tau, name):
    h, b = drift_parts(name)
    U = ZH + tau * (cellsum(HLG, h, 10) - Ghl @ Iinv @ b)[None, :]
    return float(((U ** 2 / whl).sum(1) > chi2.ppf(1 - ALPHA, 8)).mean())

def power_st(tau, name):
    h, b = drift_parts(name)
    s = tau * (Zs.T @ h / Npt - Ezx @ Iinv @ b)
    return float(ncx2.sf(chi2.ppf(1 - ALPHA, 2), 2, float(s @ Oinv @ s))) if tau > 0 else ALPHA

def tau_for(fn, target=0.5):
    hi = 1.0
    while fn(hi) < target:
        hi *= 2
        if hi > 4096: return float("inf")
    lo = 0.0
    for _ in range(25):
        mid = (lo + hi) / 2
        if fn(mid) < target: lo = mid
        else: hi = mid
    return hi

size = dict(DG=power_dg(np.zeros(36)), HL=power_hl(0, "quadratic"),
            Stukel=power_st(0, "quadratic"))
print("\nlimit size at tau = 0 (DeepGOF-1 must equal k/(B+1) = .05):  " +
      "  ".join("%s %.4f" % kv for kv in size.items()))

FORM = {"size": size, "dep": {}}
print("\n%-12s %10s | %8s %8s %8s | limit power of each test at the DeepGOF-1 tau_50"
      % ("departure", "|delta|/tau", "t50 DG", "t50 HL", "t50 St"))
for name in DEPS:
    d1 = delta_unit(name)
    f = dict(delta_norm=float(np.linalg.norm(d1)),
             tau50=dict(DG=tau_for(lambda t: power_dg(t * d1)),
                        HL=tau_for(lambda t: power_hl(t, name)),
                        Stukel=tau_for(lambda t: power_st(t, name))))
    t = f["tau50"]["DG"]
    f["power_at_tauDG"] = dict(DG=power_dg(t * d1), HL=power_hl(t, name),
                               Stukel=power_st(t, name))
    f["ARE_vs_DG"] = {k: (f["tau50"][k] / t) ** 2 for k in ("HL", "Stukel")}
    FORM["dep"][name] = f
    pw = f["power_at_tauDG"]
    print("%-12s %10.4f | %8.2f %8.2f %8.2f | DG %.3f HL %.3f St %.3f"
          % (name, f["delta_norm"], t, f["tau50"]["HL"], f["tau50"]["Stukel"],
             pw["DG"], pw["HL"], pw["Stukel"]))

print("\nasymptotic relative efficiency, (tau_50 rival / tau_50 DeepGOF-1)^2 -- above 1 means DeepGOF-1 needs fewer patients:")
for name in DEPS:
    a = FORM["dep"][name]["ARE_vs_DG"]
    print("  %-12s HL %7.2f   Stukel %7.2f" % (name, a["HL"], a["Stukel"]))
json.dump(FORM, open(os.path.join(HERE, "FORMULA.json"), "w"), indent=2)
print("\nwrote theory/local_power/FORMULA.json")
