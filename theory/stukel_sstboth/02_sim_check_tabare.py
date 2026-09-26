"""02_sim_check_tabare.py -- SstBoth on the stored datasets of theory/local_power/03_convergence.py, CONV_RAW.jsonl (null and
the five departures at DeepGOF-1's tau_50, n = 2,000, 8,000 and 32,000, 1,000 replicates each; SIM_RAW.json of
02_simulate.py stores no seeds, so its datasets cannot be regenerated), against the limit of 01_tab_are.py.

Each dataset is regenerated from its stored seed exactly as lp_core.one_rep draws it (x, then y; the bootstrap
draws come after and are not needed). Pairing check: the joint-form Stukel p-value recomputed here must equal the
stored one. SstBoth undefined (no fitted risk on one side of .5) counts as p = 1 and is counted.

Run:  python -u 02_sim_check_tabare.py   -> SIM_CHECK_TABARE.json
"""
import os, sys, json
import numpy as np
from scipy.special import expit
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
LP = os.path.join(os.path.dirname(HERE), "local_power")
sys.path.insert(0, HERE)
from sst_core import sstboth_stat, ALPHA

BETA0 = np.array([-0.5, 1.0, 0.6])
FORM = json.load(open(os.path.join(LP, "FORMULA.json")))
LIM = json.load(open(os.path.join(HERE, "TAB_ARE_SSTBOTH.json")))
RAW = [json.loads(l) for l in open(os.path.join(LP, "CONV_RAW.jsonl"))]

def gfun(name, x1, x2):
    g = {"quadratic": x1 ** 2 - 1, "interaction": x1 * x2, "cubic": x1 ** 3 - 3 * x1,
         "bump": 3 * np.exp(-((x1 - 1) ** 2 + (x2 + 0.5) ** 2) / 0.5), "threshold": 2.0 * (x2 > 1)}[name]
    return np.clip(g, -6, 6)

def newton(X, y, iters=25):                               # lp_core.newton for a single dataset
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        P = expit(X @ b)
        b = b + np.linalg.solve((X * (P * (1 - P))[:, None]).T @ X, X.T @ (y - P))
    return b

res = {}; maxdiff = 0.0; nund = 0
for (name, n, p_dg, p_hl, p_st, nclip, seed) in RAW:
    tau = 0.0 if name == "null" else FORM["dep"][name]["tau50"]["DG"]
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2)); X = np.column_stack([np.ones(n), x])
    pi0 = expit(X @ BETA0)
    h = 0.0 if name == "null" else tau * pi0 * (1 - pi0) * gfun(name, x[:, 0], x[:, 1])
    pin = np.clip(pi0 + h / np.sqrt(n), 1e-9, 1 - 1e-9)
    y = (rng.uniform(size=n) < pin).astype(float)
    b = newton(X, y)
    # joint form, verbatim from lp_core.one_rep (pairing check)
    eh = X @ b; ph = expit(eh); vh = ph * (1 - ph); r = y - ph
    Z = np.column_stack([0.5 * eh ** 2 * (eh >= 0), -0.5 * eh ** 2 * (eh < 0)])
    s = Z.T @ r; A = Z.T @ (vh[:, None] * X)
    Cv = Z.T @ (vh[:, None] * Z) - A @ np.linalg.solve(X.T @ (vh[:, None] * X), A.T)
    pj = chi2.sf(float(s @ np.linalg.solve(Cv, s)), 2)
    maxdiff = max(maxdiff, abs(pj - p_st))
    T, pm = sstboth_stat(X, y, b)
    if not np.isfinite(pm):
        nund += 1; pm = 1.0
    res.setdefault((name, n), []).append((pm, p_st, p_dg))
print("pairing: max |joint p recomputed - stored| = %.2e over %d datasets; SstBoth undefined in %d" % (maxdiff, len(RAW), nund))

out = {"pairing_maxdiff": maxdiff, "n_datasets": len(RAW), "n_undefined": nund, "cells": {}}
print("%-12s %5s | %8s %8s %6s %6s | joint: formula sim" % ("departure", "n", "limit", "sim", "SE", "z"))
for name in ["null", "quadratic", "interaction", "cubic", "bump", "threshold"]:
    lim = LIM["size_sstboth"] if name == "null" else LIM["dep"][name]["power_at_tauDG_sstboth"]
    limj = 0.05 if name == "null" else LIM["dep"][name]["power_at_tauDG_joint"]
    for n in (2000, 8000, 32000):
        P = np.array(res[(name, n)])
        sim = float((P[:, 0] <= ALPHA).mean()); simj = float((P[:, 1] <= ALPHA).mean())
        se = float(np.sqrt(max(sim * (1 - sim), 1e-4) / len(P)))
        z = (sim - lim) / se
        out["cells"]["%s|%d" % (name, n)] = dict(limit=lim, sim=sim, se=se, z=z, reps=len(P), joint_limit=limj, joint_sim=simj)
        print("%-12s %5d | %8.3f %8.3f %6.3f %6.2f | %.3f %.3f" % (name, n, lim, sim, se, z, limj, simj))
json.dump(out, open(os.path.join(HERE, "SIM_CHECK_TABARE.json"), "w"), indent=1)
print("wrote SIM_CHECK_TABARE.json")
