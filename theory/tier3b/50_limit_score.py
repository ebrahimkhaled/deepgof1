"""50_limit_score.py -- limiting power of the score test for x1^2, x2^2, x1*x2 under the local alternatives
of Theorem 2 (theory/local_power/THEOREM.md, 01_local_power.py), for the five departures of the paper's
Table 2, and the sample-size ratio relative to DeepGOF-1 at limit power .5.

Design and departures exactly as theory/local_power (lp_core.py): x1, x2 ~ N(0,1), eta_0 = -0.5 + x1 + 0.6 x2,
pi_n = pi_0 + tau v g(x) / sqrt(n). Population integrals on the same 1200 x 1200 product grid in probability
space. With two covariates the pair DeepGOF-1 selects is (x1, x2), so the score columns are fixed:
Z = (x1^2, x2^2, x1 x2). The projected score n^-1/2 sum (y - pi_hat) Z is asymptotically
N(tau s, Omega) with s = E[Z h1] - E[v Z x~'] I^-1 E[x~ h1] (h1 = v g, tau = 1) and
Omega = E[v Z Z'] - E[v Z x~'] I^-1 E[v x~ Z'], so the statistic is noncentral chi-square, 3 df,
noncentrality tau^2 s' Omega^-1 s. The bootstrap-calibrated version has the same limit.
tau_50(DG) is read from theory/local_power/FORMULA.json (the numbers behind the paper's Table 2);
n ratio = (tau_50(score) / tau_50(DG))^2: above 1 means DeepGOF-1 needs fewer observations.
Run: python -u 50_limit_score.py
"""
import os, sys, json, csv
import numpy as np
from scipy.special import ndtri, expit
from scipy.stats import chi2, ncx2

LP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/local_power"
sys.path.insert(0, LP)
from lp_core import gfun, BETA0, DEPS, ALPHA   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
G_ = 1200
u = (np.arange(G_) + 0.5) / G_
xg = ndtri(u)
X1, X2 = np.meshgrid(xg, xg, indexing="ij"); X1, X2 = X1.ravel(), X2.ravel()
Xt = np.column_stack([np.ones_like(X1), X1, X2]); N = len(X1)
pi0 = expit(Xt @ BETA0); v = pi0 * (1 - pi0)
I = (Xt * v[:, None]).T @ Xt / N; Iinv = np.linalg.inv(I)
Z = np.column_stack([X1 ** 2, X2 ** 2, X1 * X2])
Ezx = (Z * v[:, None]).T @ Xt / N
Omega = (Z * v[:, None]).T @ Z / N - Ezx @ Iinv @ Ezx.T
Oinv = np.linalg.inv(Omega)
crit = chi2.ppf(1 - ALPHA, 3)
FORM = json.load(open(os.path.join(LP, "FORMULA.json")))

def ncp1(name):
    h1 = v * gfun(name, X1, X2)
    s = Z.T @ h1 / N - Ezx @ Iinv @ (Xt.T @ h1 / N)
    return float(s @ Oinv @ s)

def power(tau, lam1):
    return float(ncx2.sf(crit, 3, tau ** 2 * lam1)) if tau > 0 else ALPHA

# noncentrality giving power .5 for a 3-df noncentral chi-square
lo, hi = 0.0, 50.0
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if ncx2.sf(crit, 3, mid) < 0.5 else (lo, mid)
LAM50 = hi
rows = []
print("noncentrality for limit power .5 at 3 df: %.4f" % LAM50)
print("%-12s %10s %10s %10s %12s %14s | n ratio vs DG: score  HL  Stukel  leCessie" % ("departure", "ncp/tau^2", "t50 score", "t50 DG", "n ratio", "score pw@tDG"))
for name in DEPS:
    l1 = ncp1(name); t50 = float(np.sqrt(LAM50 / l1)) if l1 > 0 else float("inf")
    f = FORM["dep"][name]; tdg = f["tau50"]["DG"]
    r = dict(departure=name, ncp_per_tau2=l1, tau50_score=t50, tau50_DG=tdg,
             n_ratio_score=(t50 / tdg) ** 2, power_score_at_tau50_DG=power(tdg, l1),
             n_ratio_HL=f["ARE_vs_DG"]["HL"], n_ratio_Stukel=f["ARE_vs_DG"]["Stukel"],
             n_ratio_leCessie=f["ARE_vs_DG"]["leCessie"],
             tau50_HL=f["tau50"]["HL"], tau50_Stukel=f["tau50"]["Stukel"], tau50_leCessie=f["tau50"]["leCessie"])
    rows.append(r)
    print("%-12s %10.5f %10.2f %10.2f %12.3f %14.3f | %6.2f %6.2f %6.2f %6.2f" % (name, l1, t50, tdg, r["n_ratio_score"],
          r["power_score_at_tau50_DG"], r["n_ratio_score"], r["n_ratio_HL"], r["n_ratio_Stukel"], r["n_ratio_leCessie"]))
# check: the Stukel column of FORMULA.json is reproduced by the same machinery (2 df, eta_0 columns)
eta0 = Xt @ BETA0
Zs = np.column_stack([0.5 * eta0 ** 2 * (eta0 >= 0), -0.5 * eta0 ** 2 * (eta0 < 0)])
Es = (Zs * v[:, None]).T @ Xt / N
Os = np.linalg.inv((Zs * v[:, None]).T @ Zs / N - Es @ Iinv @ Es.T)
for name in DEPS:
    h1 = v * gfun(name, X1, X2); t = FORM["dep"][name]["tau50"]["DG"]
    s = t * (Zs.T @ h1 / N - Es @ Iinv @ (Xt.T @ h1 / N))
    p = float(ncx2.sf(chi2.ppf(1 - ALPHA, 2), 2, float(s @ Os @ s)))
    print("check Stukel power at tau50(DG), %-11s: recomputed %.6f vs FORMULA.json %.6f" % (name, p, FORM["dep"][name]["power_at_tauDG"]["Stukel"]))
    rows[DEPS.index(name)]["check_stukel_recomputed"] = p
    rows[DEPS.index(name)]["check_stukel_formula_json"] = FORM["dep"][name]["power_at_tauDG"]["Stukel"]
os.makedirs(os.path.join(HERE, "tables"), exist_ok=True)
with open(os.path.join(HERE, "tables", "LIMIT_SCORE_TABLE2.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote tables/LIMIT_SCORE_TABLE2.csv")
