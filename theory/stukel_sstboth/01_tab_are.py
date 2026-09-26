"""01_tab_are.py -- the Stukel column of Table tab:are in the SstBoth form.

Same population grid, design and departures as theory/local_power/01_local_power.py (G = 1200 product grid in
probability space; x1, x2 ~ N(0, 1); eta_0 = -0.5 + x1 + 0.6 x2); DeepGOF-1's tau_50 is read from FORMULA.json
(unchanged). The joint form is recomputed as a check: it must reproduce FORMULA.json's Stukel tau_50.

Run:  python -u 01_tab_are.py     -> TAB_ARE_SSTBOTH.json
"""
import os, sys, json
import numpy as np
from scipy.special import ndtri, expit
from scipy.stats import chi2, ncx2

HERE = os.path.dirname(os.path.abspath(__file__))
LP = os.path.join(os.path.dirname(HERE), "local_power")
sys.path.insert(0, HERE)
from sst_core import SstBoth, tau_for, ALPHA

BETA0 = np.array([-0.5, 1.0, 0.6])
DEPS = ["quadratic", "interaction", "cubic", "bump", "threshold"]

def gfun(name, x1, x2):                                   # verbatim from lp_core.gfun
    g = {"quadratic": x1 ** 2 - 1, "interaction": x1 * x2, "cubic": x1 ** 3 - 3 * x1,
         "bump": 3 * np.exp(-((x1 - 1) ** 2 + (x2 + 0.5) ** 2) / 0.5), "threshold": 2.0 * (x2 > 1)}[name]
    return np.clip(g, -6, 6)

FORM = json.load(open(os.path.join(LP, "FORMULA.json")))
G_ = 1200
u = (np.arange(G_) + 0.5) / G_
xg = ndtri(u)
X1, X2 = np.meshgrid(xg, xg, indexing="ij"); X1, X2 = X1.ravel(), X2.ravel()
Xt = np.column_stack([np.ones_like(X1), X1, X2])
eta0 = Xt @ BETA0; pi0 = expit(eta0); v = pi0 * (1 - pi0); Npt = len(X1)
Iinv = np.linalg.inv((Xt * v[:, None]).T @ Xt / Npt)
Zs = np.column_stack([0.5 * eta0 ** 2 * (eta0 >= 0), -0.5 * eta0 ** 2 * (eta0 < 0)])
Ezx = (Zs * v[:, None]).T @ Xt / Npt
Omega = (Zs * v[:, None]).T @ Zs / Npt - Ezx @ Iinv @ Ezx.T
Oinv = np.linalg.inv(Omega)
S = SstBoth(Omega)
print("Omega = %s; correlation of the two marginal scores rho = %.4f; weights %s" % (Omega.round(5).tolist(), S.rho, S.lam.round(4)))
size = S.size()
print("limit null size of SstBoth referred to chi2(2): %.4f (joint form: .0500)" % size)

out = dict(design="x1, x2 ~ N(0,1); eta0 = -0.5 + x1 + 0.6 x2; grid G = 1200 as 01_local_power.py",
           Omega=Omega.tolist(), rho=S.rho, weights=S.lam.tolist(), size_sstboth=size, dep={})
for name in DEPS:
    h = v * gfun(name, X1, X2)
    s = Zs.T @ h / Npt - Ezx @ Iinv @ (Xt.T @ h / Npt)
    ncp1 = float(s @ Oinv @ s)
    pj = lambda t: float(ncx2.sf(chi2.ppf(1 - ALPHA, 2), 2, t * t * ncp1)) if t > 0 else ALPHA
    pm = lambda t: S.power(t, s)
    tdg = FORM["dep"][name]["tau50"]["DG"]
    tj, tm = tau_for(pj), tau_for(pm)
    rec = dict(s=s.tolist(), marginal_z_drift=(S.dm12 * s).tolist(), ncp_joint=ncp1,
               tau50_DG=tdg, tau50_joint=tj, tau50_joint_formula_json=FORM["dep"][name]["tau50"]["Stukel"],
               tau50_sstboth=tm, ratio_joint=(tj / tdg) ** 2, ratio_sstboth=(tm / tdg) ** 2,
               power_at_tauDG_joint=pj(tdg), power_at_tauDG_sstboth=pm(tdg))
    out["dep"][name] = rec
    print("%-12s s=(%+.4f, %+.4f)  tau50 DG %.2f | joint %.2f (FORMULA %.2f) ratio %.3f | SstBoth %.2f ratio %.3f | power at tauDG joint %.3f SstBoth %.3f"
          % (name, s[0], s[1], tdg, tj, rec["tau50_joint_formula_json"], rec["ratio_joint"], tm, rec["ratio_sstboth"],
             rec["power_at_tauDG_joint"], rec["power_at_tauDG_sstboth"]))
json.dump(out, open(os.path.join(HERE, "TAB_ARE_SSTBOTH.json"), "w"), indent=1)
print("wrote TAB_ARE_SSTBOTH.json")
