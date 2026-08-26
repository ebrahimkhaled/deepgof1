"""liu_augscore.py -- interim scoring of GATE LIU with Stute-Zhu and proj added.

INTERIM: runs over whatever blocks have an augment_pvalues.csv, so it can be read before the
last n=500 block lands. Reports how many blocks fed each number; the final table must be
regenerated once all 24 are present.

Arms: the eight external published tests + Stute-Zhu + proj (the two the AMENDMENT added),
their NA-robust Cauchy combinations, and DEEP. The author's own EF/DEF/EDGE family is present
in the deposit but is NEVER a comparator here (speckit/26 SS2).
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probit_common import randomized_matched_power, narobust_cct, ALPHA

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
EXT8 = ["HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","Stukel","Osius-Rojek","le-Cessie"]
NEW  = ["Stute-Zhu","proj"]
CS   = [0.25, 0.5, 0.75, 1.0]

rows, used, skipped = [], [], []
for d in sorted(glob.glob(os.path.join(PAP, "results/liu/S*_n*_s*"))):
    ap = os.path.join(d, "augment_pvalues.csv")
    if not os.path.exists(ap):
        skipped.append(os.path.basename(d)); continue
    b = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    a = pd.read_csv(ap).rename(columns={"SZ": "Stute-Zhu"})
    p = pd.read_csv(os.path.join(d, "deep_pvalues.csv"))
    for f in (b, a, p): f["_C"] = f["C"].fillna(-1.0)
    m = b.merge(a[["kind","rep","_C","Stute-Zhu","proj"]], on=["kind","rep","_C"], how="inner")
    m = m.merge(p[["kind","rep","_C","p_deep"]], on=["kind","rep","_C"], how="inner")
    m = m.rename(columns={"p_deep": "DEEP"})
    # the two Cauchy combinations: the original 8, and the 8 + the two new arms
    m["CCT-8"]  = narobust_cct(m[EXT8].to_numpy())
    m["CCT-10"] = narobust_cct(m[EXT8 + NEW].to_numpy())
    arms = EXT8 + NEW + ["CCT-8", "CCT-10", "DEEP"]
    nul = m[m.kind == 0]
    tag = os.path.basename(d); S, N = tag.split("_")[0], int(tag.split("_")[1][1:])
    for arm in arms:
        pn = nul[arm].to_numpy()
        size = float(np.mean(np.nan_to_num(pn, nan=1.0) <= ALPHA))
        for c in CS:
            pa = m[(m.kind == 2) & (np.isclose(m["_C"], c))][arm].to_numpy()
            if len(pa) == 0: continue
            rows.append(dict(block=tag, setting=S, n=N, arm=arm, C=c, size=size,
                             matched=randomized_matched_power(pn, pa)[0]))
    used.append(tag)

R = pd.DataFrame(rows)
R.to_csv(os.path.join(PAP, "results/liu/AUGSCORE_INTERIM.csv"), index=False)
print(f"blocks scored: {len(used)}   | not yet available: {skipped}\n")
print("=== SIZE at nominal .05 (mean over blocks, nulls) ===")
for arm, v in R.groupby("arm")["size"].mean().sort_values().items():
    flag = "  <-- inflated" if v > .075 else ("  <-- conservative" if v < .035 else "")
    print(f"  {arm:15s} {v:.4f}{flag}")
print("\n=== MATCHED POWER (mean over all alternative cells) ===")
pw = R.groupby("arm")["matched"].mean().sort_values(ascending=False)
for i, (arm, v) in enumerate(pw.items(), 1):
    print(f"  {i:2d}. {arm:15s} {v:.4f}")
print("\n=== DEEP vs each arm, by sample size (matched power) ===")
piv = R.pivot_table(index="arm", columns="n", values="matched", aggfunc="mean")
print(piv.round(4).to_string())
