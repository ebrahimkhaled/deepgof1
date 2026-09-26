"""analyze_size_by_n.py -- Addendum 1: does the size of BAGofT (and of DeepGOF-1) drift with n?

Pre-declared in ADDENDUM_1_size_by_n.md before any data existed. Reads the 2,400 null rows
(replicates 13-212 of each setting x n block), pairs them with DeepGOF-1 on the same datasets, and
writes SIZE_BY_N.csv and SIZE_BY_N.md. Refuses to run until the run is complete.

    python analyze_size_by_n.py [--rows DIR]
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.proportion import proportion_confint

PAP = r"C:\Users\ebrah\.gemini\Projects\PDFs\Paper_Ebrahim_Frangiton\paper_deepgof"
HERE = os.path.join(PAP, "confirm_bagoft")
ALPHA = 0.05
EXPECTED_NULLS = 4 * 3 * 200

ap = argparse.ArgumentParser()
ap.add_argument("--rows", default=os.path.join(HERE, "rows"))
args = ap.parse_args()

R = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(args.rows, "*.csv"))],
              ignore_index=True)
if (R.label != "CONFIRMATORY").any():
    sys.exit("rows other than CONFIRMATORY found -- refusing")
N0 = R[R.C == 0]
if len(N0) < EXPECTED_NULLS:
    sys.exit("only %d of %d null datasets on disk -- runs once, on the complete run"
             % (len(N0), EXPECTED_NULLS))
N0 = N0[N0.status != "dropped_degenerate_y"].copy()

deep = []
for SET in (1, 2, 3, 4):
    for n in (50, 100, 500):
        d = pd.read_csv(os.path.join(PAP, "results", "liu", "S%d_n%d_s101000" % (SET, n),
                                     "deep_pvalues.csv"))
        d = d[d.kind == 0][["rep", "p_deep"]]
        d["setting"], d["n"] = SET, n
        deep.append(d)
N0 = N0.merge(pd.concat(deep), on=["setting", "n", "rep"], how="left")
if N0.p_deep.isna().any():
    sys.exit("null datasets without a DeepGOF-1 p-value -- pairing broken")

rows, lines = [], ["# Addendum 1 -- size against n (ADDENDUM_1_size_by_n.md)", ""]
for arm, col in (("BAGofT", "p"), ("DeepGOF-1", "p_deep")):
    N0["rej"] = (N0[col].fillna(1.0) <= ALPHA).astype(int)     # NA cannot reject
    lines += ["## %s" % arm, "", "| n | nulls | size | Wilson 95% CI |", "|---|---|---|---|"]
    for n, g in N0.groupby("n"):
        k, m = int(g.rej.sum()), len(g)
        lo, hi = proportion_confint(k, m, alpha=0.05, method="wilson")
        rows.append(dict(arm=arm, setting="all", n=n, nulls=m, size=k / m, lo=lo, hi=hi))
        lines.append("| %d | %d | %.4f | [%.4f, %.4f] |" % (n, m, k / m, lo, hi))
    for (s, n), g in N0.groupby(["setting", "n"]):
        k, m = int(g.rej.sum()), len(g)
        lo, hi = proportion_confint(k, m, alpha=0.05, method="wilson")
        rows.append(dict(arm=arm, setting=s, n=n, nulls=m, size=k / m, lo=lo, hi=hi))
    # setting as a factor under its own name: the rows carry a column called C (the departure
    # strength), which would shadow the formula language's C() and break the fit
    N0["setting_f"] = N0.setting.astype(str)
    fit = smf.logit("rej ~ np.log(n) + setting_f", data=N0).fit(disp=0)
    b = fit.params["np.log(n)"]
    ci = fit.conf_int().loc["np.log(n)"]
    p = fit.pvalues["np.log(n)"]
    rows.append(dict(arm=arm, setting="trend", n=np.nan, nulls=len(N0), size=np.nan,
                     lo=ci[0], hi=ci[1], slope_logn=b, p_trend=p))
    lines += ["", "Trend (logit of rejection on log n, setting fixed effects): slope %.3f "
              "[%.3f, %.3f], p = %.4g -> %s" % (b, ci[0], ci[1], p,
                                                 "TREND PRESENT" if p < 0.05 else "no trend detected"),
              ""]

pd.DataFrame(rows).to_csv(os.path.join(HERE, "SIZE_BY_N.csv"), index=False)
with open(os.path.join(HERE, "SIZE_BY_N.md"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
