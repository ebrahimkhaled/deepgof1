"""liu_summary_aug.py -- rebuild the per-block LIU summary WITH the AMENDMENT arms.

The blocks' own summary.csv predates the Stute-Zhu / proj augmentation, so Figure 8 was plotting
an eleven-arm world that no longer exists. This regenerates the same schema
(setting,n,seed,arm,C,raw,matched,external,own) over all 24 blocks with the two new arms folded in,
into results/liu/SUMMARY_AUG.csv. Nothing is recomputed from scratch: the base p-values, the
augment p-values and the deployed scores are merged on (kind, rep, C).
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probit_common import randomized_matched_power, narobust_cct, ALPHA

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
EXT8 = ["HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","Stukel","Osius-Rojek","le-Cessie"]
NEW  = ["Stute-Zhu","proj"]
OWN  = ["EF","DEF.poly2","DEF.poly3","DEF.stukel","A.poly","A.spline"]

rows = []
for d in sorted(glob.glob(os.path.join(PAP, "results/liu/S*_n*_s*"))):
    tag = os.path.basename(d)
    ap = os.path.join(d, "augment_pvalues.csv")
    if not os.path.exists(ap):
        print("  SKIP (no augment):", tag); continue
    b = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    a = pd.read_csv(ap).rename(columns={"SZ": "Stute-Zhu"})
    p = pd.read_csv(os.path.join(d, "deep_pvalues.csv"))
    for f in (b, a, p): f["_C"] = f["C"].fillna(-1.0)
    m = (b.merge(a[["kind","rep","_C"]+NEW], on=["kind","rep","_C"])
           .merge(p[["kind","rep","_C","p_deep"]], on=["kind","rep","_C"])
           .rename(columns={"p_deep": "DEEP"}))
    m["CCT-8"]  = narobust_cct(m[EXT8].to_numpy())
    m["CCT-10"] = narobust_cct(m[EXT8 + NEW].to_numpy())
    setting = int(tag.split("_")[0][1:]); n = int(tag.split("_")[1][1:]); seed = int(tag.split("_s")[1])
    nul = m[m.kind == 0]
    for arm in EXT8 + NEW + OWN + ["CCT-8","CCT-10","DEEP"]:
        if arm not in m.columns: continue
        pn = nul[arm].to_numpy()
        for c in sorted(m.loc[m.kind == 2, "_C"].unique()):
            pa = m[(m.kind == 2) & (np.isclose(m["_C"], c))][arm].to_numpy()
            if not len(pa): continue
            rows.append(dict(setting=setting, n=n, seed=seed, arm=arm, C=float(c),
                             raw=float(np.mean(np.nan_to_num(pa, nan=1.0) <= ALPHA)),
                             matched=randomized_matched_power(pn, pa)[0],
                             external=(arm in EXT8 + NEW), own=(arm in OWN)))
        # the C = 0 row is the size, reported on the same curve as Liu et al. plot it
        rows.append(dict(setting=setting, n=n, seed=seed, arm=arm, C=0.0,
                         raw=float(np.mean(np.nan_to_num(pn, nan=1.0) <= ALPHA)),
                         matched=np.nan, external=(arm in EXT8 + NEW), own=(arm in OWN)))

R = pd.DataFrame(rows)
out = os.path.join(PAP, "results/liu/SUMMARY_AUG.csv")
R.to_csv(out, index=False)
print(f"\nblocks: {R.groupby(['setting','n','seed']).ngroups}  rows: {len(R):,}  arms: {R.arm.nunique()}")
print("arms:", sorted(R.arm.unique()))
print("->", out)
