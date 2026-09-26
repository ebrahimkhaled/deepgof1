"""analyze_confirm.py -- the pre-declared analysis of the BAGofT confirmatory arm.

Implements O1-O4 of PREDECLARATION.md and nothing else. Written, and its hash recorded, before
any row of the run existed. It refuses to run on an incomplete run unless called with --partial,
in which case every output is labelled PARTIAL.

    python analyze_confirm.py [--partial] [--rows DIR]

Writes CONFIRM_ROW.csv (the Table 2 row), CONFIRM_PAIRED.csv (O2 + O3) and CONFIRM_REPORT.md.
"""
import argparse
import glob
import os
import sys
from math import comb, sqrt

import numpy as np
import pandas as pd

PAP = r"C:\Users\ebrah\.gemini\Projects\PDFs\Paper_Ebrahim_Frangiton\paper_deepgof"
HERE = os.path.join(PAP, "confirm_bagoft")
sys.path.insert(0, os.path.join(PAP, "code"))
from probit_common import randomized_matched_power, clean, ALPHA  # noqa: E402  same as Table 2

SIZE_CAP = 0.075                                   # make_tables.py's rule for "(inflated)"
EXPECTED = 4 * 3 * 4 * 100 + 4 * 3 * 200           # 4,800 alternatives + 2,400 nulls

ap = argparse.ArgumentParser()
ap.add_argument("--partial", action="store_true")
ap.add_argument("--rows", default=os.path.join(HERE, "rows"))
args = ap.parse_args()

files = glob.glob(os.path.join(args.rows, "*.csv"))
R = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
if (R.label != "CONFIRMATORY").any():
    sys.exit("rows other than CONFIRMATORY found -- refusing (smoke-test rows mixed in?)")
TAG = ""
if len(R) < EXPECTED:
    if not args.partial:
        sys.exit("only %d of %d datasets on disk -- the analysis runs once, on the complete run "
                 "(use --partial for a labelled look at coverage)" % (len(R), EXPECTED))
    TAG = "PARTIAL (%d of %d datasets) -- " % (len(R), EXPECTED)

dropped = int((R.status == "dropped_degenerate_y").sum())
errors = int((R.status == "error").sum())
R = R[R.status != "dropped_degenerate_y"].copy()      # dropped as liu_grid.R drops them

# ---- DeepGOF-1 on the same datasets ---------------------------------------------------------
deep = []
for SET in (1, 2, 3, 4):
    for n in (50, 100, 500):
        d = pd.read_csv(os.path.join(PAP, "results", "liu", "S%d_n%d_s101000" % (SET, n),
                                     "deep_pvalues.csv"))
        d = d[d.kind.isin([0, 2])].copy()
        d["C"] = d["C"].fillna(0.0)
        d["setting"], d["n"] = SET, n
        deep.append(d[["setting", "n", "C", "rep", "p_deep"]])
deep = pd.concat(deep, ignore_index=True)
R = R.merge(deep, on=["setting", "n", "C", "rep"], how="left")
if R.p_deep.isna().any():
    sys.exit("%d datasets have no DeepGOF-1 p-value -- the pairing is broken" % R.p_deep.isna().sum())

nul = R[R.C == 0]
alt = R[R.C > 0]

# ---- O1: the Table 2 row --------------------------------------------------------------------
blocks, cells = [], []
for (SET, n), g0 in nul.groupby(["setting", "n"]):
    pn, _ = clean(g0.p)
    blocks.append(dict(setting=SET, n=n, size=float((pn <= ALPHA).mean()), nulls=len(g0)))
    for C, ga in alt[(alt.setting == SET) & (alt.n == n)].groupby("C"):
        pw = randomized_matched_power(g0.p.to_numpy(), ga.p.to_numpy())[0]
        cells.append(dict(setting=SET, n=n, C=C, matched=pw, alts=len(ga)))
B = pd.DataFrame(blocks)
Cc = pd.DataFrame(cells)
Cc["var"] = Cc.matched * (1 - Cc.matched) / Cc.alts      # binomial MC variance, per cell


def mean_se(frame):
    return frame.matched.mean(), sqrt(frame["var"].sum()) / len(frame)


row = dict(arm="BAGofT", size=B["size"].mean())
for n in (50, 100, 500):
    row["pw%d" % n], row["se%d" % n] = mean_se(Cc[Cc.n == n])
row["overall"], row["se_overall"] = mean_se(Cc)
row["inflated"] = row["size"] > SIZE_CAP
row.update(null_datasets=len(nul), alt_datasets=len(alt), errors=errors, dropped=dropped,
           error_rate=errors / max(len(R), 1))
pd.DataFrame([row]).to_csv(os.path.join(HERE, "CONFIRM_ROW.csv"), index=False)

# ---- O2 + O3: the paired test ---------------------------------------------------------------
def mcnemar(b, c):
    n = b + c
    return 1.0 if n == 0 else min(1.0, 2 * sum(comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n)


def paired(frame, label):
    f = frame[frame.p.notna() & frame.p_deep.notna()]
    rb, rd = f.p <= ALPHA, f.p_deep <= ALPHA
    b, c = int((rb & ~rd).sum()), int((~rb & rd).sum())
    return dict(group=label, pairs=len(f), excluded_na=len(frame) - len(f),
                bagoft_rate=rb.mean(), deepgof_rate=rd.mean(),
                only_bagoft=b, only_deepgof=c, p_exact=mcnemar(b, c))


P = [paired(alt, "all 48 power cells (O2, primary)")]
sub = [paired(alt[alt.n == n], "n = %d" % n) for n in (50, 100, 500)]
sub += [paired(alt[alt.setting == s], "Setting %d" % s) for s in (1, 2, 3, 4)]
sub += [paired(alt[alt.C == C], "C = %g" % C) for C in (0.25, 0.5, 0.75, 1.0)]
# Holm across the eleven secondary tests
ps = np.array([s["p_exact"] for s in sub])
order = np.argsort(ps)
adj = np.empty_like(ps)
running = 0.0
for k, i in enumerate(order):
    running = max(running, min(1.0, (len(ps) - k) * ps[i]))
    adj[i] = running
for s, a in zip(sub, adj):
    s["p_holm"] = a
P[0]["p_holm"] = np.nan
PD = pd.DataFrame(P + sub)
PD.to_csv(os.path.join(HERE, "CONFIRM_PAIRED.csv"), index=False)

# ---- O4: both tests' sizes on this run's own nulls --------------------------------------------
pb, _ = clean(nul.p)
pdp, _ = clean(nul.p_deep)
size_b, size_d = float((pb <= ALPHA).mean()), float((pdp <= ALPHA).mean())

# ---- report ---------------------------------------------------------------------------------
o2 = P[0]
verdict = ("no difference detected" if o2["p_exact"] >= 0.05 else
           ("DeepGOF-1 AHEAD" if o2["only_deepgof"] > o2["only_bagoft"] else "BAGofT AHEAD"))
L = ["# BAGofT confirmatory arm -- %sreport" % TAG, "",
     "Pre-declared analysis (PREDECLARATION.md), `analyze_confirm.py`. Implementation `bagoft_fast()` "
     "(identical to BAGofT 1.0.0, VERIFY_FAST.log).", "",
     "Datasets: %d nulls, %d alternatives, %d dropped as degenerate, %d errors (%.2f%%)."
     % (len(nul), len(alt), dropped, errors, 100 * row["error_rate"]), "",
     "## O1 -- the Table 2 row", "",
     "| size | n=50 | n=100 | n=500 | overall | MC SE (overall) | inflated? |", "|---|---|---|---|---|---|---|",
     "| %.3f | %.3f | %.3f | %.3f | **%.3f** | %.4f | %s |" % (row["size"], row["pw50"], row["pw100"],
                                                          row["pw500"], row["overall"],
                                                          row["se_overall"], row["inflated"]), "",
     "## O2 -- primary paired test: **%s**" % verdict, "",
     "| group | pairs | BAGofT | DeepGOF-1 | only BAGofT | only DeepGOF-1 | exact p | Holm p |",
     "|---|---|---|---|---|---|---|---|"]
for _, r in PD.iterrows():
    L.append("| %s | %d | %.3f | %.3f | %d | %d | %.4g | %s |" % (
        r.group, r.pairs, r.bagoft_rate, r.deepgof_rate, r.only_bagoft, r.only_deepgof,
        r.p_exact, "--" if np.isnan(r.p_holm) else "%.4g" % r.p_holm))
L += ["", "## O4 -- sizes on this run's own nulls (replicates 13-212)", "",
      "BAGofT %.3f, DeepGOF-1 %.3f. O2 compares raw decisions at p <= .05; the more liberal test is "
      "favoured by at most roughly the difference." % (size_b, size_d), ""]
io_path = os.path.join(HERE, "CONFIRM_REPORT.md")
with open(io_path + ".part", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
os.replace(io_path + ".part", io_path)
print("\n".join(L))
