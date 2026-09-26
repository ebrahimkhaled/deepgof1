"""03_irrelevant.py -- the Stukel curve of Figure fig:irrelevant in the SstBoth form, and the predicted-vs-simulated
check against the niche rows (theory/addon_arms/niche_rows, column Stukel_marginal).

Population quantities: ir_core.GridDesign (unchanged; G = 1200, 24 Gauss-Hermite nodes), Omega = inverse of its
Oinv, s(h) = E[z h] - E[v z x'] I^-1 b(h). DeepGOF-1's tau_50 is read from ARE_IRRELEVANT.json (unchanged). The
joint form is recomputed as a check against ARE_IRRELEVANT.json.

Resumable: one record per (variant, k, dep) appended to IRR_RAW.jsonl.
Run:  python -u 03_irrelevant.py       -> IRR_SSTBOTH.json, PRED_VS_SIM_SSTBOTH.json
"""
import os, sys, json, time, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
IRD = os.path.join(os.path.dirname(HERE), "local_power", "irrelevant")
sys.path.insert(0, IRD); sys.path.insert(0, HERE)
from ir_core import GridDesign, h_local, h_finite, power_stukel, C_IRR, DEPS
from ir_core import tau_for as tau_for_ir
from sst_core import SstBoth, tau_for, ALPHA

ARE = json.load(open(os.path.join(IRD, "ARE_IRRELEVANT.json")))
KS = ARE["k"]
RAW = os.path.join(HERE, "IRR_RAW.jsonl")
done = {}
if os.path.exists(RAW):
    for line in open(RAW):
        try:
            r = json.loads(line); done[(r["variant"], r["k"], r["dep"])] = r
        except ValueError:
            pass

VARS = {"study": C_IRR, "pure": 0.0}
des8 = None
for vname, cirr in VARS.items():
    for k in KS:
        need = [dep for dep in DEPS if (vname, k, dep) not in done and str(k) in ARE["variants"][vname][dep]]
        if not need and not (vname == "study" and k == 8):
            continue
        t0 = time.time()
        des = GridDesign(k, cirr=cirr)
        Omega = np.linalg.inv(des.Oinv); S = SstBoth(Omega)
        for dep in need:
            a, ahl, zh, b = des.drift(h_local(dep))
            s = zh - des.Ezx @ des.Iinv @ b
            ncp1 = des.ncp_stukel(zh, b)
            st = ARE["variants"][vname][dep][str(k)]
            tdg = st["tau50"]["DG"]
            tj = tau_for_ir(lambda t: power_stukel(t, ncp1))
            tm = tau_for(lambda t: S.power(t, s))
            rec = dict(variant=vname, k=k, dep=dep, rho=S.rho, weights=S.lam.tolist(), size_sstboth=S.size(),
                       s=s.tolist(), tau50_DG=tdg, tau50_joint=tj, tau50_joint_stored=st["tau50"]["Stukel"],
                       ratio_joint=(tj / tdg) ** 2, ratio_joint_stored=st["ratio_vs_DG"]["Stukel"],
                       tau50_sstboth=tm, ratio_sstboth=(tm / tdg) ** 2, secs=time.time() - t0)
            with open(RAW, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            done[(vname, k, dep)] = rec
            print("[%s k=%d %-5s] rho %.4f size %.4f | ratio joint %.3f (stored %.3f) | SstBoth %.3f" %
                  (vname, k, dep, S.rho, rec["size_sstboth"], rec["ratio_joint"], rec["ratio_joint_stored"], rec["ratio_sstboth"]))
            sys.stdout.flush()
        if vname == "study" and k == 8:
            des8, S8 = des, S
        else:
            del des

out = {"k": KS, "variants": {vn: {dep: {str(k): done[(vn, k, dep)] for k in KS if (vn, k, dep) in done} for dep in DEPS}
                             for vn in VARS}}
json.dump(out, open(os.path.join(HERE, "IRR_SSTBOTH.json"), "w"), indent=1)
print("wrote IRR_SSTBOTH.json")

# ---- predicted vs simulated at k = 8 (the niche design) ----
rates = pd.read_csv(os.path.join(os.path.dirname(HERE), "addon_arms", "niche_addon_rates.csv"))
rates = rates[rates.arm.isin(["Stukel_marginal", "Stukel_joint"])]
pv = {"size_limit_sstboth": S8.size(), "rho_k8": S8.rho, "cells": {}}
print("\nk = 8 niche design: predicted (linear drift / exact finite-n drift) vs simulated SstBoth")
for _, g in rates.groupby("scen"):
    scen = g.scen.iloc[0]
    sm = g[g.arm == "Stukel_marginal"].iloc[0]; sj = g[g.arm == "Stukel_joint"].iloc[0]
    if scen.startswith("null"):
        pv["cells"][scen] = dict(pred=S8.size(), sim=float(sm.rate), se=float(sm.se), na=float(sm.na), N=int(sm.N),
                                 z=float((sm.rate - S8.size()) / sm.se), joint_sim=float(sj.rate))
        print("  %-16s limit %.4f  sim %.4f (SE %.4f)  z %.1f   [joint sim %.3f]" % (scen, S8.size(), sm.rate, sm.se, pv["cells"][scen]["z"], sj.rate))
        continue
    fam, C, n = scen.split("_")[0], float(scen.split("_")[1][1:]), int(scen.split("_")[2][1:])
    tau = C * np.sqrt(n)
    a, ahl, zh, b = des8.drift(h_local(fam)); s_lin = zh - des8.Ezx @ des8.Iinv @ b
    ae, ahle, zhe, be = des8.drift(h_finite(fam, C, n, 8, C_IRR)); s_ex = zhe - des8.Ezx @ des8.Iinv @ be
    p_lin = S8.power(tau, s_lin); p_ex = S8.power(1.0, s_ex)
    z = (sm.rate - p_ex) / sm.se
    pv["cells"][scen] = dict(fam=fam, C=C, n=n, tau=tau, pred_lin=p_lin, pred_exact=p_ex, sim=float(sm.rate), se=float(sm.se),
                             na=float(sm.na), N=int(sm.N), z_exact=float(z), z_lin=float((sm.rate - p_lin) / sm.se),
                             joint_sim=float(sj.rate))
    print("  %-16s tau %5.2f  pred lin %.3f exact %.3f  sim %.3f (SE %.3f)  z(exact) %+.1f   [joint sim %.3f]"
          % (scen, tau, p_lin, p_ex, sm.rate, sm.se, z, sj.rate))
json.dump(pv, open(os.path.join(HERE, "PRED_VS_SIM_SSTBOTH.json"), "w"), indent=1)
print("wrote PRED_VS_SIM_SSTBOTH.json")
