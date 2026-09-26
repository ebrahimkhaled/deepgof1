"""02_np_gap.py -- W2: how close is DeepGOF-1 to the most powerful map-based test for its own prior?

THEOREM (THEOREM_NP.md). Let P0, P1 be the null and prior-averaged alternative laws of the map, f =
dP1/dP0, phi* the Neyman-Pearson test at level alpha, and phi_S = 1{S > c} any score test at exact level
alpha. For every b,
        E1 phi* - E1 phi_S  <=  E0 | f - exp(S - b) |,
and phi* beats every map-based level-alpha test.

WHAT THIS SCRIPT MEASURES, on held-out maps from the shipped prior (01_family_maps.R):
  A. SURROGATE WORLD (exact): P0 = N(mu0, Sigma0) fitted to the bootstrap-null maps; P1 = the equal-
     weight mixture of N(m_f, Sigma0) over the families, m_f estimated on half of each family's maps.
     Here f is known in closed form, so the NP power, the DeepGOF-1 power and the theorem's bound are
     all computed exactly (Monte Carlo, SE reported).
  B. REAL MAPS: DeepGOF-1 against the surrogate-LR test on the other half of the real held-out maps,
     both thresholded at the 95% point of the real null maps.
Both for the 100 TRAINING families and for 100 FRESH families from the same prior.

Run:  python -u 02_np_gap.py
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.special import logsumexp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "local_power"))
from lp_core import score

HERE = os.path.dirname(os.path.abspath(__file__))
D = pd.read_csv(os.path.join(HERE, "FAMILY_MAPS.csv"))
cols = ["c%d" % i for i in range(1, 37)]
rng = np.random.default_rng(314159)
A = 0.05
null = D[D.y == 0][cols].to_numpy()
rng.shuffle(null)
null_fit, null_test = null[: len(null) // 2], null[len(null) // 2:]
mu0 = null_fit.mean(0); S0 = np.cov(null_fit.T)
S0i = np.linalg.inv(S0); L0 = np.linalg.cholesky(S0)
print("null maps: %d (fit %d / test %d); Sigma0 eigenvalues %.3f .. %.3f"
      % (len(null), len(null_fit), len(null_test), *np.linalg.eigvalsh(S0)[[0, -1]]))

def lr_stat(M, means):
    """log f(m) for the equal-weight Gaussian mixture with common covariance Sigma0."""
    X = M - mu0
    lin = X @ S0i @ means.T                                    # n x F
    quad = 0.5 * np.einsum("fi,ij,fj->f", means, S0i, means)
    return logsumexp(lin - quad[None, :], axis=1) - np.log(len(means))

def rate(x):
    r = float(np.mean(x)); return r, float(np.sqrt(max(r * (1 - r), 1e-6) / len(x)))

out = {}
for setname in ("train", "fresh"):
    Dm = D[D.set == setname]
    fams = sorted(Dm.fam.unique())
    est, test, means = [], [], []
    for f in fams:
        M = Dm[Dm.fam == f][cols].to_numpy()
        h = len(M) // 2
        est.append(M[:h]); test.append(M[h:]); means.append(M[:h].mean(0) - mu0)
    means = np.array(means)
    # ---- A. surrogate world -----------------------------------------------------------------
    NS = 200000
    Z0 = mu0 + rng.normal(size=(NS, 36)) @ L0.T
    fid = rng.integers(0, len(fams), NS)
    Z1 = mu0 + means[fid] + rng.normal(size=(NS, 36)) @ L0.T
    lf0, lf1 = lr_stat(Z0, means), lr_stat(Z1, means)
    s0, s1 = score(Z0), score(Z1)
    k_np, c_dg = np.quantile(lf0, 1 - A), np.quantile(s0, 1 - A)
    p_np, se_np = rate(lf1 > k_np); p_dg, se_dg = rate(s1 > c_dg)
    # the theorem's bound: min over b of E0 | f - exp(S - b) |
    f0 = np.exp(lf0 - lf0.max()) * np.exp(lf0.max())            # f on null draws
    bs = np.linspace(-10, 10, 401)
    bound = [float(np.mean(np.abs(f0 - np.exp(s0 - b)))) for b in bs]
    bmin = float(min(bound)); b_at = float(bs[int(np.argmin(bound))])
    # the gap measured directly, with a paired SE
    gap_draws = (lf1 > k_np).astype(float) - (s1 > c_dg).astype(float)
    gap, gse = float(gap_draws.mean()), float(gap_draws.std() / np.sqrt(NS))
    # ---- B. real maps ------------------------------------------------------------------------
    R1 = np.vstack(test); lr_null = lr_stat(null_test, means); s_null = score(null_test)
    k_r, c_r = np.quantile(lr_null, 1 - A), np.quantile(s_null, 1 - A)
    pr_lr, se_lr = rate(lr_stat(R1, means) > k_r); pr_dg, se_rdg = rate(score(R1) > c_r)
    per_fam = []
    for f, T in zip(fams, test):
        per_fam.append((float(np.mean(score(T) > c_r)), float(np.mean(lr_stat(T, means) > k_r))))
    per_fam = np.array(per_fam)
    out[setname] = dict(families=len(fams), surrogate=dict(np_power=p_np, np_se=se_np, dg_power=p_dg,
                        dg_se=se_dg, gap=gap, gap_se=gse, share_of_np=p_dg / p_np,
                        theorem_bound=bmin, bound_b=b_at),
                        real=dict(lr_power=pr_lr, lr_se=se_lr, dg_power=pr_dg, dg_se=se_rdg,
                                  dg_beats_lr_in_families=int((per_fam[:, 0] > per_fam[:, 1]).sum()),
                                  maps=len(R1)))
    o = out[setname]
    print("\n=== %s families (%d) ===" % (setname, len(fams)))
    print("  SURROGATE (exact NP available): NP power %.3f (%.3f)   DeepGOF-1 %.3f (%.3f)"
          % (p_np, se_np, p_dg, se_dg))
    print("     measured gap %.3f (SE %.3f) -> DeepGOF-1 attains %.1f%% of the NP-optimal power"
          % (gap, gse, 100 * p_dg / p_np))
    print("     theorem bound min_b E0|f - e^(S-b)| = %.3f (at b = %.2f)" % (bmin, b_at))
    print("  REAL MAPS: surrogate-LR test %.3f (%.3f)   DeepGOF-1 %.3f (%.3f); DeepGOF-1 ahead in %d of %d families"
          % (pr_lr, se_lr, pr_dg, se_rdg, o["real"]["dg_beats_lr_in_families"], len(fams)))
# ---- the FAIR out-of-sample comparison: the LR test built from the TRAINING families' means, scored
# on the FRESH families' held-out maps. DeepGOF-1 never saw the fresh families either.
tr = D[D.set == "train"]; fr = D[D.set == "fresh"]
means_tr = np.array([tr[tr.fam == f][cols].to_numpy()[: (tr.fam == f).sum() // 2].mean(0) - mu0
                     for f in sorted(tr.fam.unique())])
Rf = np.vstack([fr[fr.fam == f][cols].to_numpy()[(fr.fam == f).sum() // 2:] for f in sorted(fr.fam.unique())])
k_fair = np.quantile(lr_stat(null_test, means_tr), 1 - A)
c_fair = np.quantile(score(null_test), 1 - A)
pf_lr, sf_lr = rate(lr_stat(Rf, means_tr) > k_fair); pf_dg, sf_dg = rate(score(Rf) > c_fair)
out["fair_fresh"] = dict(lr_trainmeans_power=pf_lr, lr_se=sf_lr, dg_power=pf_dg, dg_se=sf_dg)
print("\nFAIR (both built on training families, scored on FRESH families' real maps):")
print("  Gaussian-mixture LR from training-family means %.3f (%.3f)   DeepGOF-1 %.3f (%.3f)" % (pf_lr, sf_lr, pf_dg, sf_dg))
json.dump(out, open(os.path.join(HERE, "NP_GAP.json"), "w"), indent=2)
print("\nwrote theory/optimality/NP_GAP.json")
