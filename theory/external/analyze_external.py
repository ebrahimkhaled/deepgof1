# analyze_external.py -- level and power of run_external.R. Writes EXTERNAL_LEVEL.csv and EXTERNAL_POWER.csv.
# Power is raw and at a matched level of .05 (randomized threshold on the same test's nulls at n = 500).
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"]) for f in glob.glob(os.path.join(HERE, "rows", "*.csv"))])
d = d.drop_duplicates(["fam", "C", "n", "rep"])
T = ["DG_axes", "DG_allpairs", "DG_combined", "proj", "lecessie", "pearson", "HL", "stukel", "cox", "spiegelhalter", "giviti"]

def matched(p0, p1, a=.05):
    p0, p1 = p0[np.isfinite(p0)], p1[np.isfinite(p1)]
    for c in np.unique(p0):
        if (p0 <= c).mean() >= a:
            eq = (p0 == c).mean()
            g = (a - (p0 < c).mean()) / eq if eq > 0 else 0
            return (p1 < c).mean() + g * (p1 == c).mean()
    return 1.0

nul = d[d.fam == "null"]
L = nul.groupby("n")[T].agg(lambda s: (s.dropna() <= .05).mean()).T
L.to_csv(os.path.join(HERE, "EXTERNAL_LEVEL.csv"))
print("level (reps", nul.groupby("n").size().to_dict(), ")"); print(L.round(3).to_string())
n0 = nul[nul.n == 500]
rows = []
for (fam, C), g in d[d.fam != "null"].groupby(["fam", "C"], sort=False):
    for t in T:
        rows.append(dict(fam=fam, C=C, reps=len(g), test=t, raw=(g[t].dropna() <= .05).mean(),
                         matched=matched(n0[t].to_numpy(float), g[t].to_numpy(float))))
P = pd.DataFrame(rows); P.to_csv(os.path.join(HERE, "EXTERNAL_POWER.csv"), index=False)
pd.set_option("display.width", 200)
print("matched-level power, n = 500 (reps", P.groupby(["fam", "C"]).reps.first().min(), ")")
print(P.pivot_table(index=["fam", "C"], columns="test", values="matched", sort=False)[T].round(2).to_string())
