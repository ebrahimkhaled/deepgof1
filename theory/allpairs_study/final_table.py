"""final_table.py -- the two DeepGOF rows and the shipped-vs-design-trained gap table.

Consumes phd/rows/ (Settings 2-4, four readings) and phd/rows_s1/ (Setting 1, shipped weights) and
emits, for make_tables.py to pick up:

  * DEEP_SHIPPED.csv     the deployed arm: shipped weights, ALL-PAIRS reading.
                         All-pairs IS v1 wherever p <= 2, by construction -- Setting 1 has no pairs
                         and Settings 3 and 4 have exactly one, which is the pair v1 selects. Only
                         Setting 2 differs. That is why the arm is assembled rather than read off a
                         single column, and why the assembly is written down here rather than done
                         by hand.
  * DEEP_SHIPPED_V1.csv  the same weights with the deployed axis rule, for the honest comparison.
  * SHIPGAP.csv          shipped minus design-trained, cell by cell, all four settings.
  * MACROS.txt           the numbers EDIT_PLAN_2026-09-21_protocol.md names.

Matched-level power uses the same estimator and the same held-out-null split as every other arm in
Table 2 (code/liu_score.py: the second half of each block's nulls calibrates).

Run:  python final_table.py
"""
import os, sys, glob
import numpy as np
import pandas as pd

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
sys.path.insert(0, os.path.join(PAP, "code"))
from probit_common import randomized_matched_power, ALPHA

HERE = os.path.join(PAP, "phd")
NNULL = 2000


def load(pattern):
    fs = sorted(glob.glob(pattern))
    if not fs:
        return None
    out = []
    for f in fs:
        d = pd.read_csv(f)
        if "allpairs" not in d.columns:
            d["allpairs"] = np.nan
        d["block"] = os.path.basename(f).replace(".csv", "")
        out.append(d)
    return pd.concat(out, ignore_index=True)


R24 = load(os.path.join(HERE, "rows", "*.csv"))
R1 = load(os.path.join(HERE, "rows_s1", "*.csv"))
if R1 is None:
    sys.exit("Setting 1 has not finished: phd/rows_s1/ is empty. Nothing is written.")
if len(R1.block.unique()) < 6:
    sys.exit("Setting 1 incomplete: %d of 6 blocks. Refusing to write a partial table."
             % len(R1.block.unique()))
R = pd.concat([R24, R1], ignore_index=True)
print("blocks %d   rows %d   settings %s"
      % (R.block.nunique(), len(R), sorted(R.setting.unique())))

R["heldout"] = (R.kind == 0) & (R.rep > NNULL // 2)

# ---- the deployed arm: all-pairs where it exists, v1 where all-pairs IS v1 --------------------
R["deployed"] = np.where(R.allpairs.notna(), R.allpairs, R.v1)
assert R.loc[R.setting == 2, "allpairs"].notna().all(), "Setting 2 must carry all-pairs"
assert R.loc[R.setting != 2, "allpairs"].isna().all(), "only Setting 2 should carry all-pairs"

cells, sizes = [], []
for (blk, arm) in [(b, a) for b in sorted(R.block.unique()) for a in ("deployed", "v1")]:
    g = R[R.block == blk]
    pn = g.loc[g.heldout, arm].to_numpy(float)
    sizes.append(dict(block=blk, setting=int(g.setting.iloc[0]), n=int(g.n.iloc[0]), arm=arm,
                      size=float(np.mean(pn <= ALPHA)), nulls=len(pn)))
    for C, ga in g[g.kind == 2].groupby("C"):
        mp = randomized_matched_power(pn, ga[arm].to_numpy(float))[0]
        cells.append(dict(block=blk, setting=int(g.setting.iloc[0]), n=int(g.n.iloc[0]), arm=arm,
                          C=float(C), matched=mp, raw=float(np.mean(ga[arm] <= ALPHA))))
CE = pd.DataFrame(cells); SZ = pd.DataFrame(sizes)

by_set = CE.pivot_table(index="arm", columns="setting", values="matched", aggfunc="mean")
by_n = CE.pivot_table(index="arm", columns="n", values="matched", aggfunc="mean")
overall = CE.groupby("arm")["matched"].mean()
size_m = SZ.groupby("arm")["size"].mean()

# ---- the design-trained column, from the deposit ----------------------------------------------
dep_rows = []
for blk, g in R.groupby("block"):
    s, n, sd = int(g.setting.iloc[0]), int(g.n.iloc[0]), int(g.seed.iloc[0])
    f = os.path.join(PAP, "results", "liu", "S%d_n%d_s%d" % (s, n, sd), "deep_pvalues.csv")
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f); d["C"] = d["C"].fillna(0.0)
    half = d[(d.kind == 0) & (d.rep > NNULL // 2)].p_deep.to_numpy(float)
    for C, ga in d[d.kind == 2].groupby("C"):
        dep_rows.append(dict(block=blk, setting=s, n=n, C=float(C), arm="design-trained",
                             matched=randomized_matched_power(half, ga.p_deep.to_numpy(float))[0],
                             raw=float(np.mean(ga.p_deep <= ALPHA))))
DP = pd.DataFrame(dep_rows)

print("\n=== matched-level power by setting ===")
tab = pd.concat([by_set,
                 DP.pivot_table(index="arm", columns="setting", values="matched", aggfunc="mean")])
tab["overall"] = tab.mean(axis=1)
tab = tab.rename(index={"deployed": "DeepGOF-1 (shipped, all-pairs)",
                        "v1": "DeepGOF-1 (shipped, deployed axis rule)"})
print(tab.round(3).to_string())

print("\n=== by sample size, deployed arm ===")
both = pd.concat([CE, DP])
print(both.pivot_table(index="arm", columns="n", values="matched", aggfunc="mean").round(3).to_string())

print("\n=== size on held-out nulls ===")
for a in ("deployed", "v1"):
    s = SZ[SZ.arm == a]
    print("  %-9s mean %.4f   sd %.4f   range [%.3f, %.3f]   blocks %d"
          % (a, s["size"].mean(), s["size"].std(), s["size"].min(), s["size"].max(), len(s)))

# ---- the gap table ----------------------------------------------------------------------------
gap = (CE[CE.arm == "v1"].merge(DP, on=["block", "setting", "n", "C"], suffixes=("_ship", "_dep")))
gap["diff_raw"] = gap.raw_ship - gap.raw_dep
G = gap.groupby(["setting", "n", "C"])[["raw_ship", "raw_dep", "diff_raw"]].mean().reset_index()
print("\n=== shipped vs design-trained, raw rejection, %d cells ===" % len(G))
print("  mean difference %+.4f   worst %+.4f   weaker in %d of %d"
      % (G.diff_raw.mean(), G.diff_raw.min(), (G.diff_raw < -0.005).sum(), len(G)))

CE.to_csv(os.path.join(HERE, "DEEP_SHIPPED.csv"), index=False)
DP.to_csv(os.path.join(HERE, "DEEP_DESIGNTRAINED.csv"), index=False)
G.to_csv(os.path.join(HERE, "SHIPGAP.csv"), index=False)

with open(os.path.join(HERE, "MACROS.txt"), "w", encoding="ascii") as fh:
    def emit(k, v):
        fh.write("%s\t%s\n" % (k, v)); print("  %-24s %s" % (k, v))
    print("\n=== macros for make_tables.py ===")
    emit("powdeepshipped", "%.3f" % overall["deployed"])
    emit("powdeepshippedv1", "%.3f" % overall["v1"])
    emit("powdeepdesigntrained", "%.3f" % DP.groupby("setting")["matched"].mean().mean())
    emit("sizedeepshipped", "%.4f" % size_m["deployed"])
    emit("sdsizedeepshipped", "%.4f" % SZ[SZ.arm == "deployed"]["size"].std())
    emit("shipgappower", "%+.3f" % G[G.C > 0].diff_raw.mean())
    emit("shipgapworst", "%+.3f" % G[G.C > 0].diff_raw.min())
    emit("shipgapcells", "%d" % len(G[G.C > 0]))
print("\nwrote DEEP_SHIPPED.csv, DEEP_DESIGNTRAINED.csv, SHIPGAP.csv, MACROS.txt")
