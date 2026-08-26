"""s3_verdict.py -- GATE S3 stage 3: the main table and the verdict against frozen criteria.

Judged against speckit/21_S3_PREDECLARATION.md section 6. The criteria are printed alongside the
result so the verdict and the frozen text can be read together, and every number is rebuilt from the
deposited per-cell summaries rather than from anything held in memory.

Runs on a PARTIAL grid: it reports what exists and says plainly how much is missing, so it can be
used to watch the run without pretending an incomplete grid is a verdict.
"""
import os, sys, glob
import numpy as np
import pandas as pd

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
NS = [50, 75, 100, 150, 200]
DESIGNS = ["D1", "D2", "D3", "D4"]
SEEDS = [101000, 102000, 103000]
DESC = {"D1": "textbook p=2", "D2": "p=3, 3-level factor", "D3": "p=3, AR(1) rho=.7",
        "D4": "p=2, 15% events"}

f = sorted(glob.glob(os.path.join(PAP, "results", "s3", "*", "summary.csv")))
if not f:
    sys.exit("no S3 summaries yet")
V = pd.concat([pd.read_csv(x) for x in f], ignore_index=True)
V = V[V.rule == "pooled_D1"]                     # THE pre-declared rule (21 section 5)
V.to_csv(os.path.join(PAP, "results", "s3", "VERDICT.csv"), index=False)

done = V.groupby(["n", "design"]).seed.nunique()
tot = len(NS) * len(DESIGNS) * len(SEEDS)
have = int(V.groupby(["n", "design", "seed"]).ngroups)
print(f"\ncells present: {have}/{tot} (n x design x seed)")
if have < tot:
    print("*** PARTIAL GRID -- this is not the final verdict ***")

# ---- the main table: gain of the learned rule over the pre-declared comparator ------------------
def gain_table(comparator):
    rows = []
    for n in NS:
        for d in DESIGNS:
            s = V[(V.n == n) & (V.design == d)]
            if not len(s): continue
            deep = s[s.arm == "DEEP"].set_index("seed")
            if comparator == "CCT-11":
                riv = s[s.arm == "CCT-11"].set_index("seed").matched_power
            else:
                riv = s[~s.arm.isin(["DEEP", "CCT-11"])].groupby("seed").matched_power.max()
            g = deep.matched_power - riv.reindex(deep.index)
            rows.append(dict(n=n, design=d, seeds=len(g), deep=deep.matched_power.mean(),
                             rival=riv.mean(), gain=g.mean(), gain_lo=g.min(), gain_hi=g.max(),
                             dsize=deep.raw_size.mean(), size_min=deep.raw_size.min(),
                             size_max=deep.raw_size.max(), auc=deep.score_auc.mean()))
    return pd.DataFrame(rows)

for comp in ["CCT-11", "best-of-11"]:
    T = gain_table(comp)
    if not len(T): continue
    print("\n" + "=" * 92)
    print(f"MAIN TABLE -- matched-level gain over {comp}   (mean of available seeds; range in brackets)")
    print("=" * 92)
    print(f"{'n':>4} {'design':>22} {'DEEP':>7} {comp:>10} {'gain':>8} {'[range]':>18} {'size':>7} {'seeds':>5}")
    for _, r in T.iterrows():
        print(f"{int(r.n):>4} {DESC[r.design]:>22} {r.deep:>7.4f} {r.rival:>10.4f} {r.gain:>+8.4f} "
              f"{f'[{r.gain_lo:+.3f},{r.gain_hi:+.3f}]':>18} {r.dsize:>7.4f} {int(r.seeds):>5}")
    print(f"{'':>4} {'MEAN':>22} {T.deep.mean():>7.4f} {T.rival.mean():>10.4f} {T.gain.mean():>+8.4f}")

T = gain_table("CCT-11")                          # criteria are judged against the pre-declared comparator

# ---- the frozen criteria ------------------------------------------------------------------------
print("\n" + "=" * 92)
print("CRITERIA (speckit/21 section 6, frozen before any data)")
print("=" * 92)
bad = V[(V.arm == "DEEP") & ((V.raw_size < .035) | (V.raw_size > .065))]
a = len(bad) == 0
print(f"  [{'PASS' if a else 'FAIL'}] S3a  deployed size in [.035,.065] in EVERY cell")
if not a:
    print(f"         {len(bad)} violating cell(s):")
    for _, r in bad.iterrows():
        print(f"           n={int(r.n)} {r.design} seed={int(r.seed)}: size {r.raw_size:.4f}")
pos = (T.gain > 0).mean() if len(T) else 0
b = pos >= .80
print(f"  [{'PASS' if b else 'FAIL'}] S3b  gain > 0 in >= 80% of cells        observed {100*pos:.1f}% "
      f"({int((T.gain>0).sum())}/{len(T)})")
c = T.gain.mean() >= .03 if len(T) else False
print(f"  [{'PASS' if c else 'FAIL'}] S3c  mean gain >= +.03                  observed {T.gain.mean():+.4f}")
if 50 in set(T.n) and 200 in set(T.n):
    g50, g200 = T[T.n == 50].gain.mean(), T[T.n == 200].gain.mean()
    d_ = g50 > g200
    print(f"  [{'PASS' if d_ else 'FAIL'}] S3d  gain(n=50) > gain(n=200)          {g50:+.4f} vs {g200:+.4f}")
else:
    print("  [ -- ] S3d  needs both n=50 and n=200")

# ---- the pooling cost, a declared diagnostic (21 section 5) -------------------------------------
W = pd.concat([pd.read_csv(x) for x in f], ignore_index=True)
if "per_design" in set(W.rule):
    print("\n" + "-" * 92)
    print("POOLING COST (diagnostic only, never the headline -- 21 section 8)")
    p = W[(W.arm == "DEEP")].pivot_table(index=["n", "design", "seed"], columns="rule",
                                         values="matched_power")
    if {"pooled_D1", "per_design"} <= set(p.columns):
        diff = (p["per_design"] - p["pooled_D1"]).dropna()
        print(f"  per-design rule minus pooled-D1 rule: mean {diff.mean():+.4f} "
              f"(range {diff.min():+.4f} to {diff.max():+.4f}, {len(diff)} cells)")
        print("  positive => retraining per design would help; the shipped test does NOT offer that.")

print(f"\nwrote results/s3/VERDICT.csv")
