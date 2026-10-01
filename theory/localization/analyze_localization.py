# analyze_localization.py -- summarise rows/*.csv of run_localization.R into LOCALIZATION.csv.
# Per cell (family, C), over all datasets:
#   rej        the tool rejects at .05 (DeepGOF-1: its p-value; LWY and spline screen: any covariate flagged at
#              the Bonferroni level .05/5, the tool's own global rejection)
#   right      DeepGOF-1: the displayed pair contains the true covariate(s) (both x1 and x3 for inter);
#              LWY / spline: the flagged set is exactly the true covariate(s)
#   located    DeepGOF-1: rejects, flags at least one cell, and every flagged cell is inside the true region;
#              LWY / spline: rejects and flags exactly the true covariate(s)
#   false      at least one flag outside the truth (any flag at all under the null)
import glob, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
d = pd.concat([pd.read_csv(f, keep_default_na=False) for f in glob.glob(os.path.join(HERE, "rows", "*.csv"))])
TRUTH = {"null": set(), "ushape": {"x1"}, "thresh": {"x2"}, "inter": {"x1", "x3"}}

def split(s):
    return set(s.split("|")) if s else set()

out = []
for (fam, C), g in d.groupby(["fam", "C"], sort=False):
    t = TRUTH[fam]
    for tool in ["DG axis rule", "DG all-pairs", "Lin-Wei-Ying", "spline screen"]:
        if tool.startswith("DG"):
            p, pair, fl, fin = (g["p_axes"], g["axes"], g["flagA"], g["flagA_in"]) if tool == "DG axis rule" else \
                               (g["p_all"], g["top"], g["flagT"], g["flagT_in"])
            rej = p <= .05
            right = pair.map(lambda s: bool(t) and t <= split(s))
            located = rej & (fl > 0) & (fin == fl)
            false = rej & (fl > fin)
        else:
            fs = (g["lwy_flags"] if tool == "Lin-Wei-Ying" else g["spl_flags"]).map(split)
            rej = fs.map(bool)
            right = fs.map(lambda s: bool(t) and s == t)
            located = rej & right
            false = fs.map(lambda s: bool(s - t))
        n = len(g)
        row = dict(family=fam, C=C, tool=tool, n=n, rej=rej.mean(), right=right.mean(),
                   located=located.mean(), false=false.mean())
        row["se_max"] = np.sqrt(.25 / n)
        out.append(row)
res = pd.DataFrame(out)
res.to_csv(os.path.join(HERE, "LOCALIZATION.csv"), index=False)
pd.set_option("display.width", 160)
print(res.round(3).to_string(index=False))
