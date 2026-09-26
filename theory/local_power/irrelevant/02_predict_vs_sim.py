"""02_predict_vs_sim.py -- the local-power formulas at k = 8 (the niche design, c_irr = 0.15) against the
niche simulation's rows (theory/niche_paper/NICHE_ROWS.csv, READ ONLY; its proj_superseded column is the
benchmark's first projection weight, which this script predicts; Stukel_joint is the joint score form, which
the local formula of power_stukel describes).

Mapping: the study generates eta = eta_0 + C g(x1, x2) at sample size n, so the local amplitude is
tau = C sqrt(n) and h = tau v g (Theorem L, 'lin'). The finite-n drift the study actually generates is
h_n = sqrt(n) {pi(eta_0 + C g) - pi_0} ('exact'); the theorem covers h_n -> h, and using h_n itself is
the natural finite-n refinement. Both are reported.
For the quadratic-form tests two discretisations are given: the M = 3000 reference sample (the limit, as
in 01_are_irrelevant.py) and M = n ('n-surrogate': the Gaussian-error version of the statistic at the
study's own sample size, averaged over reference samples), which keeps the diagonal share of the moment
reference that vanishes only as n -> infinity.
Simulated: raw rejection at p <= .05 with binomial SE (an NA counts as p = 1, the study's rule).

Run:  python -u 02_predict_vs_sim.py
"""
import os, sys, json, glob, time
import numpy as np
import pandas as pd
from ir_core import *

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = os.path.join(os.path.dirname(os.path.dirname(HERE)), "niche_paper", "NICHE_ROWS.csv")
OUT = os.path.join(HERE, "PRED_VS_SIM.json")
ARMS = ["HL", "Stukel", "proj", "DG", "DG_all"]
QARMS = ("proj",)
KK, CIRR = 8, C_IRR

# ---- simulation rows ----
fs = [ROWS]
D = pd.read_csv(ROWS, keep_default_na=False, na_values=["NA", ""])  # 'null' is a scenario name, not NA
D = D.rename(columns={"Stukel_joint": "Stukel", "proj_superseded": "proj"})
assert not D.duplicated(["scen", "rep"]).any()
na = {a: pd.to_numeric(D[a], errors="coerce").isna() for a in ARMS}
for a in ARMS:
    D[a] = pd.to_numeric(D[a], errors="coerce").fillna(1.0)
print("rows: %d replicates in %d files (%s)" % (len(D), len(fs), time.strftime("%Y-%m-%d %H:%M")))
sim = {}
for (fam, C, n), g in D.groupby(["fam", "C", "n"]):
    key = "%s_C%.1f_n%d" % (fam, C, n)
    rej = {a: float((g[a] <= ALPHA).mean()) for a in ARMS}
    sim[key] = dict(fam=fam, C=float(C), n=int(n), R=int(len(g)),
                    raw=rej, se={a: float(np.sqrt(rej[a] * (1 - rej[a]) / len(g))) for a in ARMS},
                    sel12=float(g["sel12"].mean()), na={a: float(na[a][g.index].mean()) for a in ARMS})
print("cells:", ", ".join("%s (R=%d)" % (k, s["R"]) for k, s in sim.items()))
cells = [(key, s) for key, s in sim.items() if s["fam"] != "null"]
ns = sorted(set(int(n) for n in D["n"]))

# ---- predictions ----
t0 = time.time()
des = GridDesign(KK, cirr=CIRR); dgp = DGPower(des); hlp = HLPower(des)
print("design built %.0f s" % (time.time() - t0)); sys.stdout.flush()
lin = {dep: des.drift(h_local(dep)) for dep in DEPS}

def build(arm, M, seed):
    x, X, _, v = ref_sample(KK, M, seed, cirr=CIRR)
    R = proj_kernel(x)
    q = QuadForm(R, X, v, seed=seed + 7)
    q.prep(np.zeros(M), "null")
    for dep in DEPS:
        q.prep(v * gfun(dep, x[:, 0], x[:, 1]), dep)
    for key, s in cells:
        q.prep(hexact(s["fam"], x, CIRR, s["C"], s["n"]), key)
    q.drop()
    return q

# asymptotic reference samples (seed 1 of 01_are_irrelevant.py, same draw as there)
Q3 = {"proj": build("proj", 1000, 2000 + KK)}
print("asymptotic kernels built %.0f s" % (time.time() - t0)); sys.stdout.flush()
NSUR = {"proj": 2}
Qn = {(arm, n): [build(arm, n, 7000 + 10 * n + s) for s in range(NSUR[arm])] for n in ns for arm in QARMS}
print("n-surrogate kernels built %.0f s" % (time.time() - t0)); sys.stdout.flush()

def qf_power(q, arm, tau, key):
    return q.power_bootstrap(tau, key)

pred = {"size": {"DG": dgp.power(np.zeros(K * K)), "HL": hlp.power(0, np.zeros(10)), "Stukel": ALPHA}}
for arm in QARMS:
    pred["size"][arm] = qf_power(Q3[arm], arm, 0, "null")
    for n in ns:
        pred["size"]["%s_n%d" % (arm, n)] = float(np.mean([qf_power(q, arm, 0, "null") for q in Qn[(arm, n)]]))
for key, s in cells:
    fam, C, n = s["fam"], s["C"], s["n"]; tau = C * np.sqrt(n)
    a, ahl, zh, b = lin[fam]
    ae, ahle, zhe, be = des.drift(h_finite(fam, C, n, KK, CIRR))
    P = dict(tau=tau)
    P["DG_lin"] = dgp.power(tau * des.delta_dg(a, b)); P["DG_exact"] = dgp.power(des.delta_dg(ae, be))
    P["HL_lin"] = hlp.power(tau, des.drift_hl(ahl, b)); P["HL_exact"] = hlp.power(1.0, des.drift_hl(ahle, be))
    P["Stukel_lin"] = power_stukel(tau, des.ncp_stukel(zh, b)); P["Stukel_exact"] = power_stukel(1.0, des.ncp_stukel(zhe, be))
    for arm in QARMS:
        P[arm + "_lin"] = qf_power(Q3[arm], arm, tau, fam)
        P[arm + "_exact"] = qf_power(Q3[arm], arm, 1.0, key)
        P[arm + "_nsur_lin"] = float(np.mean([qf_power(q, arm, tau, fam) for q in Qn[(arm, n)]]))
        P[arm + "_nsur_exact"] = float(np.mean([qf_power(q, arm, 1.0, key) for q in Qn[(arm, n)]]))
    pred[key] = P
    print("%-16s tau %5.2f | " % (key, tau) + " ".join("%s %.3f/%.3f" % (arm, P[arm + "_lin"], P[arm + "_exact"]) for arm in ("DG", "HL", "Stukel", "proj")))
    sys.stdout.flush()

json.dump({"sim": sim, "pred": pred, "note": __doc__}, open(OUT, "w"), indent=1)

# ---- table ----
print("\npredicted (exact drift; quadratic forms: limit / n-surrogate) vs simulated raw power, binomial SE, z = (sim - ref)/SE, ref = n-sur if available else limit")
print("%-16s %-8s %7s %7s %7s %7s %7s" % ("cell", "arm", "lin", "limit", "n-sur", "sim", "SE  z"))
for key, s in cells:
    P = pred[key]
    for arm in ("DG", "HL", "Stukel", "proj"):
        lim = P[arm + "_exact"]; nsur = P.get(arm + "_nsur_exact", float("nan"))
        ref = nsur if np.isfinite(nsur) else lim
        z = (s["raw"][arm] - ref) / max(s["se"][arm], 1e-9)
        print("%-16s %-8s %7.3f %7.3f %7.3f %7.3f %7.3f %6.1f" % (key, arm, P[arm + "_lin"], lim, nsur, s["raw"][arm], s["se"][arm], z))
print("\nnull size (simulated raw, SE) vs limit size and n-surrogate:")
for key, s in sim.items():
    if s["fam"] == "null":
        for arm in ("DG", "HL", "Stukel", "proj"):
            ks = "%s_n%d" % (arm, s["n"])
            print("  %-10s %-8s sim %.4f (SE %.4f)  limit %.4f  n-sur %s" % (key, arm, s["raw"][arm], s["se"][arm], pred["size"][arm],
                  ("%.4f" % pred["size"][ks]) if ks in pred["size"] else "-"))
print("wrote", OUT)
