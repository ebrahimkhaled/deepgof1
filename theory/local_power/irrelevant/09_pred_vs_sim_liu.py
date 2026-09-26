"""09_pred_vs_sim_liu.py -- predicted vs simulated power of the projection test with the printed weight on the raw
covariates (arm proj_liu_raw), niche design, k = 8.

Predicted: PROJ_LIU_RAW.json 'pred' (08_proj_liu_raw.py: six reference samples of M = 1000, the exact finite-n drift
sqrt(n){pi(eta_0 + C g) - pi_0}, the B = 199 rule at its limit), with its Monte Carlo SE over the samples.
Simulated: raw rejection at p <= .05 of proj_liu_raw in theory/addon_arms/proj_rows_niche/ (READ ONLY), binomial SE;
checked against theory/addon_arms/proj_rates_niche.csv. z = (sim - pred) / SE_binomial; z2 also adds the
prediction's own SE. Paired with the stored DeepGOF-1 p-values (theory/niche/rows/, READ ONLY) on (scen, rep):
difference in rejection rate DG - proj_liu_raw with its paired SE and the discordant counts.
Writes PRED_VS_SIM_LIU.json and prints the table.   Run:  python 09_pred_vs_sim_liu.py
"""
import os, glob, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.path.dirname(os.path.dirname(HERE))
PR = json.load(open(os.path.join(HERE, "PROJ_LIU_RAW.json")))["pred"]
OLD = json.load(open(os.path.join(HERE, "PRED_VS_SIM.json")))
rd = lambda pat: pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA", ""]) for f in glob.glob(pat)],
                           ignore_index=True)
L = rd(os.path.join(TH, "addon_arms", "proj_rows_niche", "*.csv"))
S = rd(os.path.join(TH, "niche_paper", "NICHE_ROWS.csv")).rename(columns={"proj_superseded": "proj"})
assert not L.duplicated(["scen", "rep"]).any() and not S.duplicated(["scen", "rep"]).any()
J = L.merge(S[["scen", "rep", "DG", "proj"]], on=["scen", "rep"], how="inner", validate="1:1")
assert len(J) == len(L) == len(S), (len(J), len(L), len(S))
assert (np.abs(J["proj_old"] - J["proj"]).dropna() <= 1e-12).all()        # pairing check where proj_old exists
CSV = pd.read_csv(os.path.join(TH, "addon_arms", "proj_rates_niche.csv")).set_index("scen")

out = {"cells": {}, "null": {}}
print("%-16s %6s | %-22s %8s %6s | %s" % ("cell", "tau", "pred (MC SE)", "sim", "SE", "z  z2 | DG sim  DG-proj (paired SE)  disc DG-only/proj-only"))
for scen, g in J.groupby("scen"):
    m = len(g); pl = g["proj_liu_raw"].fillna(1.0) <= 0.05; dg = g["DG"].fillna(1.0) <= 0.05
    r = float(pl.mean()); se = float(np.sqrt(r * (1 - r) / m))
    assert abs(r - CSV.loc[scen, "proj_liu_raw"]) < 1e-12
    d = dg.astype(float) - pl.astype(float)
    rec = dict(m=m, sim=r, se=se, na=int(g["proj_liu_raw"].isna().sum()), DG_sim=float(dg.mean()),
               diff_DG_minus_proj=float(d.mean()), diff_se=float(d.std(ddof=1) / np.sqrt(m)),
               DG_only=int((dg & ~pl).sum()), proj_only=int((pl & ~dg).sum()))
    if scen.startswith("null"):
        out["null"][scen] = rec
        print("%-16s %6s | %-22s %8.3f %6.3f | (limit size %.3f) | DG %.3f" % (scen, "0", "", r, se, PR["size"]["mean"], rec["DG_sim"]))
        continue
    p = PR[scen]; C, n = float(scen.split("_C")[1].split("_")[0]), int(scen.split("_n")[1])
    rec.update(tau=C * np.sqrt(n), pred=p["mean"], pred_se=p["se"], z=(r - p["mean"]) / se,
               z2=(r - p["mean"]) / np.sqrt(se ** 2 + p["se"] ** 2),
               DG_pred=OLD["pred"][scen]["DG_exact"])
    out["cells"][scen] = rec
    print("%-16s %6.1f | %.3f (%.3f)            %8.3f %6.3f | %+5.1f %+5.1f | %.3f  %+.3f (%.3f)  %d/%d"
          % (scen, rec["tau"], p["mean"], p["se"], r, se, rec["z"], rec["z2"], rec["DG_sim"],
             rec["diff_DG_minus_proj"], rec["diff_se"], rec["DG_only"], rec["proj_only"]))
zs = [c["z"] for c in out["cells"].values()]
out["summary"] = dict(max_abs_z=float(np.max(np.abs(zs))), n_within_2=int(np.sum(np.abs(zs) <= 2)),
                      n_cells=len(zs))
print("max |z| %.2f; %d of %d cells within 2 SE" % (out["summary"]["max_abs_z"], out["summary"]["n_within_2"], len(zs)))
json.dump(out, open(os.path.join(HERE, "PRED_VS_SIM_LIU.json"), "w"), indent=1)
