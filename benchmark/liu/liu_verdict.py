"""liu_verdict.py -- GATE LIU stage 3: the benchmark verdict against speckit/26 §4.

Rebuilt from the deposited per-block summaries. Reports the RAW curves (the quantity Liu et al.
plot and the one their "power functions are not even increasing with C" claim refers to) and the
matched-level ranking, across both corpora, with ranges.

The author's own EF/DEF/EDGE family is printed for reference and never used as a comparator.
`proj` and Stute-Zhu are absent; no "best available" claim is permitted (26 §6).
"""
import os, sys, glob
import numpy as np
import pandas as pd

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
EXT = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "Stukel", "Osius-Rojek", "le-Cessie"]
OWN = ["EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "A.poly", "A.spline"]
LIU_SET = ["HL", "Pigeon-Heyse", "Stukel", "Osius-Rojek"]   # the tests Liu et al. actually included
DESC = {1: "S1 quadratic, p=1", 2: "S2 quadratic OFF-INDEX, p=3",
        3: "S3 interaction, binary x2", 4: "S4 interaction, both continuous"}

f = sorted(glob.glob(os.path.join(PAP, "results", "liu", "S*_n*_s*", "summary.csv")))
if not f:
    sys.exit("no LIU summaries yet")
T = pd.concat([pd.read_csv(x) for x in f], ignore_index=True)
T.to_csv(os.path.join(PAP, "results", "liu", "VERDICT.csv"), index=False)
ncell = T.groupby(["setting", "n", "seed"]).ngroups
print(f"\nblocks scored: {ncell}/24" + ("" if ncell == 24 else "   *** PARTIAL ***"))

rise = (T.sort_values("C").groupby(["setting", "n", "seed", "arm"])
          .raw.agg(lambda s: s.iloc[-1] - s.iloc[0]).rename("rise").reset_index())
size = T[T.C == 0][["setting", "n", "seed", "arm", "raw"]].rename(columns={"raw": "size"})

print("\n" + "=" * 96)
print("RAW rejection: RISE from C=0 to C=1, mean over corpora (range in brackets)")
print("=" * 96)
for S in sorted(T.setting.unique()):
    print(f"\n-- {DESC[S]} --")
    print(f"{'n':>5} {'OUR TEST':>18} {'best external':>22} {'best of LIU-SET':>22} {'our rank':>9}")
    for n in sorted(T[T.setting == S].n.unique()):
        r = rise[(rise.setting == S) & (rise.n == n)]
        g = r.groupby("arm").rise
        m, lo, hi = g.mean(), g.min(), g.max()
        if "DEEP" not in m.index: continue
        ext = m[[t for t in EXT if t in m.index]]
        liu = m[[t for t in LIU_SET if t in m.index]]
        rk = int((ext > m["DEEP"]).sum()) + 1
        print(f"{n:>5} {m['DEEP']:>+8.3f} [{lo['DEEP']:+.3f},{hi['DEEP']:+.3f}]"
              f" {ext.idxmax():>12s} {ext.max():>+8.3f}"
              f" {liu.idxmax():>12s} {liu.max():>+8.3f}   {rk}/{len(ext)+1}")

print("\n" + "=" * 96)
print("CRITERIA (speckit/26 §4, frozen before any data)")
print("=" * 96)
d = size[size.arm == "DEEP"]
bad = d[(d["size"] < .035) | (d["size"] > .065)]
print(f"  [{'PASS' if len(bad)==0 else 'FAIL'}] LIU1  size in [.035,.065] every cell   "
      f"({len(bad)} violation(s) of {len(d)}; range {d['size'].min():.4f}-{d['size'].max():.4f})")
for _, r in bad.iterrows():
    print(f"          S{int(r.setting)} n={int(r.n)} seed={int(r.seed)}: {r['size']:.4f}")
s2 = rise[(rise.setting == 2) & (rise.arm == "DEEP")].rise
ok2 = s2.mean() > .15
print(f"  [{'PASS' if ok2 else 'FAIL'}] LIU2  Setting 2 rise > +.15            "
      f"mean {s2.mean():+.4f} (range {s2.min():+.4f}..{s2.max():+.4f})")
mp = T[T.C > 0]
dp = mp[mp.arm == "DEEP"].groupby(["setting", "n", "seed", "C"]).matched.mean()
c8 = mp[mp.arm == "CCT-8"].groupby(["setting", "n", "seed", "C"]).matched.mean().reindex(dp.index)
g3 = (dp - c8).mean()
print(f"  [{'PASS' if g3 > .03 else 'FAIL'}] LIU3  matched gain over CCT-8 > +.03   observed {g3:+.4f}")
be = mp[mp.arm.isin(EXT)].groupby(["setting", "n", "seed", "C"]).matched.max().reindex(dp.index)
print(f"  [ rep] LIU4  vs best external, per setting (matched level):")
for S in sorted(T.setting.unique()):
    i = [k for k in dp.index if k[0] == S]
    print(f"          {DESC[S]:32s} {(dp[i]-be[i]).mean():+.4f}")

print("\n" + "-" * 96)
print("THE EXTERNAL FINDING: is Setting 2 a blind spot of the LITERATURE, or of LIU'S TEST SET?")
r2 = rise[rise.setting == 2].groupby("arm").rise.mean()
print("  tests Liu et al. INCLUDED:")
for t in LIU_SET:
    if t in r2.index: print(f"      {t:16s} {r2[t]:+.3f}")
print("  published tests they OMITTED:")
for t in ["Tsiatis", "Xie", "le-Cessie", "HL-equalwidth"]:
    if t in r2.index: print(f"      {t:16s} {r2[t]:+.3f}")
print(f"      {'OUR TEST':16s} {r2['DEEP']:+.3f}")
print("\n  the author's own family (reference only, never a comparator):")
for t in OWN:
    if t in r2.index: print(f"      {t:16s} {r2[t]:+.3f}")
print("\nwrote results/liu/VERDICT.csv")
