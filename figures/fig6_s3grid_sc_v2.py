"""AUDIT VERSION (v2, 2026-09-25): left panel against the best published partition test.
Figure 6 at Statistics and Computing print size (174 mm, Arial 8/9 pt, fonts embedded).

Same data (results/s3/VERDICT.csv, results/s3/EXTERNAL_COMPARATOR.csv), same two panels, same
reference objects as fig6_s3grid.py. Changes for the journal's artwork rules only:
- canvas = printed width, lettering 8 pt (ticks, direct labels, annotations) / 9 pt (axis labels,
  panel titles), Arial embedded;
- each design gets its own marker shape as well as its colour (circle, square, triangle, diamond),
  in both panels, so the four series stay distinct in greyscale;
- the direct labels of the two light series (tan, light green) are set in a darker shade of the
  same hue so the text reaches 4.5:1 contrast; the lines and dots keep the original colours.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()

PAP = str(S.PAP)
V = pd.read_csv(os.path.join(PAP, "results", "s3", "VERDICT.csv"))
d = V[V.arm == "DEEP"].set_index(["n", "design", "seed"])
# audit version (2026-09-25): the comparator is the best of the three published partition tests
# (HL, Pigeon-Heyse, HL equal-width) in each cell, chosen in hindsight, from VERDICT.csv itself
P = V.pivot_table(index=["n", "design", "seed"], columns="arm", values="matched_power")
G = (P["DEEP"] - P[["HL", "HL-equalwidth", "Pigeon-Heyse"]].max(axis=1)).rename("gain").reset_index()
G = G.merge(d.raw_size.rename("size").reset_index(), on=["n", "design", "seed"])
assert len(G) == 60

DESC = {"D1": "p=2, textbook", "D2": "p=3, 3-level factor",
        "D3": "p=3, correlated $\\rho$=.7", "D4": "p=2, 15% events"}
COL = {"D1": S.DEEP, "D2": S.CLASSIC, "D3": S.TAN, "D4": S.DEEP_2}
MRK = {"D1": "o", "D2": "s", "D3": "^", "D4": "D"}
MS = {"o": 4.2, "s": 3.9, "^": 4.6, "D": 3.6}
NS = sorted(G.n.unique())

fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.FULL_IN, 2.85),
                             gridspec_kw={"width_ratios": [1.35, 1]})

# ---- LEFT: the gain that never arrived -----------------------------------------------------
ax.axhspan(0, .06, color=S.BAND, alpha=.06, lw=0)
ax.axhline(0, color=S.INK, lw=1.0)
ax.annotate("parity", (NS[0], -.006), fontsize=S.TICK, color=S.INK, alpha=.85, va="top")
for des in ["D4", "D1", "D3", "D2"]:
    s = G[G.design == des].groupby("n").gain
    m, lo, hi = s.mean(), s.min(), s.max()
    ax.fill_between(NS, lo.reindex(NS), hi.reindex(NS), color=COL[des], alpha=.13, lw=0)
    ax.plot(NS, m.reindex(NS), color=COL[des], lw=1.7, marker=MRK[des], ms=MS[MRK[des]],
            markeredgecolor="white", markeredgewidth=.8)
    S.endlab(ax, NS[-1], m.reindex(NS).iloc[-1], DESC[des], COL[des], dx=5)
ax.set_xscale("log"); S.plain_log_ticks(ax, NS)
ax.set_xlim(45, 380)
ax.set_xlabel("sample size $n$")
ax.set_ylabel("matched power,\nlearned $-$ best partition test", linespacing=1.15)
ax.set_title("Ahead of the Hosmer–Lemeshow family in every design", fontsize=S.LABEL,
             pad=6, loc="left")

# ---- RIGHT: the claim that survived ---------------------------------------------------------
bx.axhspan(.035, .065, color=S.BAND, alpha=.10, lw=0)
bx.axhline(.05, color=S.INK, lw=.8, ls=(0, (4, 3)), alpha=.6)
rng = np.random.default_rng(5)                     # same jitter as the original figure
G = G.copy()
G["xj"] = np.nan
for des in ["D1", "D2", "D3", "D4"]:
    m = G.design == des
    G.loc[m, "xj"] = G.loc[m, "n"].to_numpy() * np.exp(rng.uniform(-.05, .05, int(m.sum())))
    bx.scatter(G.loc[m, "xj"], G.loc[m, "size"], marker=MRK[des],
               s=(MS[MRK[des]] * 1.25) ** 2, color=COL[des], alpha=.85,
               edgecolor="white", linewidth=.6, zorder=3)

out = G[(G["size"] < .035) | (G["size"] > .065)]
bx.scatter(out["xj"], out["size"], s=110, facecolor="none", edgecolor=S.CLASSIC,
           linewidth=1.4, zorder=4)
for _, r in out.iterrows():
    bx.annotate(f"{r['size']:.3f}", (r["xj"], r["size"]), xytext=(8, -3),
                textcoords="offset points", fontsize=S.TICK, color=S.CLASSIC, fontweight="semibold")
bx.set_xscale("log"); S.plain_log_ticks(bx, NS); bx.set_xlim(42, 260)
bx.set_ylim(.02, .08)
bx.set_xlabel("sample size $n$"); bx.set_ylabel("deployed type-I error")
bx.set_title("58 of 60 cells hold the nominal level", fontsize=S.LABEL, pad=6, loc="left")
# starts at n = 57, clear of the n = 50 dots at .062 (at 8 pt the label from 43.5 ran over them)
bx.annotate("pre-declared band", (57, .0644), fontsize=S.TICK, color=S.DEEP, alpha=.9, va="top")

fig.tight_layout(pad=0.25, w_pad=1.2)
S.save(fig, "fig6_s3grid")

print("\nby design (mean gain):")
for des in ["D1", "D2", "D3", "D4"]:
    print(f"  {DESC[des]:26s} {G[G.design==des].gain.mean():+.4f}")
print(f"  {'OVERALL':26s} {G.gain.mean():+.4f}")
print(f"size: mean {G['size'].mean():.4f}, in band {((G['size']>=.035)&(G['size']<=.065)).sum()}/60")
for des in ["D3", "D4"]:
    c = S.ink_for(COL[des])
    print(f"label ink {des}: {COL[des]} -> {c}  contrast {S.contrast(COL[des]):.2f} -> {S.contrast(c):.2f}")
