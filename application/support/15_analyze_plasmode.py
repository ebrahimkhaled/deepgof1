# analyze_plasmode2.py -- level and power of plasmode2_SUPPORT.R. Writes PLASMODE2_LEVEL.csv and PLASMODE2_POWER.csv.
# Level: rejection rate at .05 with its Monte Carlo SE. Power: raw, and at a matched level of .05 (a randomized
# threshold set on the same test's nulls at the same n, so a liberal test earns no power from its excess size).
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=["NA"])
               for f in glob.glob(os.path.join(HERE, "plasmode2_rows", "*.csv"))])
T = ["DG_axes", "DG_all", "DG_comb", "proj", "StukelJ", "StukelS", "HL", "GiViTI"]

def matched(p0, p1, a=.05):
    p0, p1 = p0[np.isfinite(p0)], p1[np.isfinite(p1)]
    if len(p0) == 0 or len(p1) == 0:
        return np.nan
    for c in np.unique(p0):
        if (p0 <= c).mean() >= a:
            eq = (p0 == c).mean()
            g = (a - (p0 < c).mean()) / eq if eq > 0 else 0
            return (p1 < c).mean() + g * (p1 == c).mean()
    return (p1 <= 1).mean()

lev, pw = [], []
for n, g in d.groupby("n"):
    nul, alt = g[g.truth == "null"], g[g.truth == "alt"]
    for t in T:
        x = nul[t].to_numpy(float); x = x[np.isfinite(x)]
        if len(x):
            r = (x <= .05).mean()
            lev.append(dict(n=n, test=t, reps=len(x), size=r, se=np.sqrt(.05 * .95 / len(x))))
        if len(alt):
            y = alt[t].to_numpy(float); y = y[np.isfinite(y)]
            if len(y):
                pw.append(dict(n=n, test=t, reps=len(y), raw=(y <= .05).mean(), matched=matched(x, y)))
L, P = pd.DataFrame(lev), pd.DataFrame(pw)
L.to_csv(os.path.join(HERE, "PLASMODE2_LEVEL.csv"), index=False)
P.to_csv(os.path.join(HERE, "PLASMODE2_POWER.csv"), index=False)
print(L.pivot(index="test", columns="n", values="size").reindex(T).round(3).to_string())
print("reps per n:", L.groupby("n").reps.max().to_dict(), " SE at .05:", L.groupby("n").se.max().round(3).to_dict())
print(P.pivot(index="test", columns="n", values="raw").reindex(T).round(3).to_string())
print(P.pivot(index="test", columns="n", values="matched").reindex(T).round(3).to_string())
print("alt reps per n:", P.groupby("n").reps.max().to_dict())
print("median seconds per replicate:", d.groupby("n").secs.median().round(1).to_dict())
