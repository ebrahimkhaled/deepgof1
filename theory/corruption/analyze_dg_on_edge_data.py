# analyze_dg_on_edge_data.py -- false-alarm rates at .05 on the stored corrupted data sets: DeepGOF-1 (this run)
# beside the stored p-values of the other tests on the SAME replicates. Writes CORRUPTION_FALSE_ALARMS_MT.csv.
import glob, gzip, os
import numpy as np, pandas as pd

BAT = r"C:\Users\ebrah\.gemini\Projects\PDFs\Paper_Ebrahim_Frangiton\simulations\battery"
HERE = os.path.dirname(os.path.abspath(__file__))
dg = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(HERE, "rows_edge_mt", "*.csv"))])
dg["block"] = dg["block"].astype(str)
print("fingerprint mismatches:", int((~dg["match"]).sum()), "of", len(dg))
TESTS = ["HL.G10", "Stk.joint", "GiViTI", "Cubic.LR", "EDGE.poly3.u.Grule", "proj", "BAGofT"]
rows = []
for (blk, fl), g in dg.groupby(["block", "file"], sort=False):
    S = pd.read_csv(os.path.join(BAT, blk, fl + "_pvalues.csv.gz"))
    m = g.merge(S, on=["rep", "seed"], how="left")
    r = dict(block=blk, cell=fl, reps=len(m), DeepGOF1=(m["DG"] <= .05).mean())
    for t in TESTS:
        if t in m:
            v = m[t].dropna()
            r[t] = (v <= .05).mean() if len(v) else np.nan
    r["se_at_.05"] = np.sqrt(.05 * .95 / len(m))
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(os.path.join(HERE, "CORRUPTION_FALSE_ALARMS_MT.csv"), index=False)
pd.set_option("display.width", 200)
print(R.round(3).to_string(index=False))
