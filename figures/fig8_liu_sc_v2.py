"""Figure 8, audit version (2026-09-25): the Liu et al. (2024) benchmark with the arm set of Table 2.

Same canvas, lattice and lettering as fig8_liu_sc.py (kept unchanged). Data: the per-block scores
written by manuscript/sc_audit/make_tables_audit.py (liu_table2_scores.csv): raw rejection at .05
for each departure strength, and the null rejection rate at C = 0, averaged over the two corpora.
Classes: Liu et al.'s tests (rust: HL, mHL, Osius-Rojek, Stukel, Stute-Zhu), their projection test
(ink, dash-dot; Liu et al.'s statistic on the raw covariates), the Hosmer-Lemeshow variants and
GiViTI (tan dashes), and DeepGOF-1 (green, ribbon = range over the two training corpora).
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()

src = os.path.join(str(S.PAP), "results", "liu", "TABLE2_SCORES.csv")   # written by tables/make_tables.py
R = pd.read_csv(src)
R["seed"] = R.block.str.split("_s").str[1]
R["setting"] = R.setting.str[1:].astype(int)
# one row per (block, arm, C) with the raw rate; C = 0 carries the null rate of the block
null = R.drop_duplicates(["block", "arm"]).assign(C=0.0, raw=lambda z: z["size"])
T = pd.concat([null, R])[["setting", "n", "seed", "arm", "C", "raw"]]

LIU = ["HL", "mHL", "Osius-Rojek", "Stukel", "Stute-Zhu"]
OTHER = ["Pigeon-Heyse", "HL-equalwidth", "GiViTI"]
TITLE = {1: "Setting 1\nquadratic,\n$p=1$",
         2: "Setting 2\nquadratic OFF-INDEX,\n$p=3$",
         3: "Setting 3\ninteraction,\nbinary $x_2$",
         4: "Setting 4\ninteraction,\nboth continuous"}
SETS = sorted(T.setting.unique()); NS = sorted(T.n.unique())

fig, AX = plt.subplots(len(NS), len(SETS), figsize=(S.FULL_IN, 5.30),
                       sharex=True, sharey=True,
                       gridspec_kw=dict(wspace=.10, hspace=.14))
AX = np.atleast_2d(AX)

for r, n in enumerate(NS):
    for c, Sx in enumerate(SETS):
        ax = AX[r, c]
        s = T[(T.setting == Sx) & (T.n == n)]
        CS = sorted(s.C.unique())
        cur = lambda arm: [s[(s.arm == arm) & (s.C == v)].raw.mean() for v in CS]
        ax.axhspan(.035, .065, color=S.BAND, alpha=.07, lw=0)
        ax.axhline(.05, color=S.INK, lw=.6, ls=(0, (3, 3)), alpha=.45)
        for t in LIU:
            ax.plot(CS, cur(t), color=S.CLASSIC, lw=1.0, alpha=.85, zorder=3)
        for t in OTHER:
            ax.plot(CS, cur(t), color=S.TAN, lw=1.2, ls=(0, (4, 2)), alpha=.95, zorder=4)
        ax.plot(CS, cur("proj"), color=S.INK, lw=1.4, ls=(0, (5, 1.6, 1, 1.6)),
                alpha=.9, zorder=5)
        d = s[s.arm == "DEEP"]
        lo = [d[d.C == v].raw.min() for v in CS]; hi = [d[d.C == v].raw.max() for v in CS]
        ax.fill_between(CS, lo, hi, color=S.DEEP, alpha=.20, lw=0, zorder=6)
        ax.plot(CS, cur("DEEP"), color=S.DEEP, lw=1.9, zorder=7, marker="o", ms=3.2,
                markeredgecolor="white", markeredgewidth=.7)
        ax.set_ylim(-.03, 1.03); ax.set_xlim(-.05, 1.05)
        ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .25, .5, .75, 1])
        ax.grid(True, color=S.GRID, lw=.5, alpha=.6)
        if r == 0:
            ax.set_title(TITLE[Sx], fontsize=S.LABEL, pad=5, linespacing=1.15,
                         color=S.DEEP if Sx == 2 else S.INK,
                         fontweight="semibold" if Sx == 2 else "normal")
        if c == 0: ax.set_ylabel(f"$n={n}$\nrejection at $\\alpha$=.05", fontsize=S.LABEL,
                                 linespacing=1.15)
        if r == len(NS)-1: ax.set_xlabel("departure parameter  $C$", fontsize=S.LABEL)
        if Sx == 2:
            for sp in ax.spines.values(): sp.set_color(S.DEEP); sp.set_linewidth(1.2)

leg = [Line2D([], [], color=S.DEEP, lw=1.9, marker="o", ms=3.4, label="DeepGOF-1"),
       Line2D([], [], color=S.CLASSIC, lw=1.0, label="other tests of Liu et al."),
       Line2D([], [], color=S.INK, lw=1.4, ls=(0, (5, 1.6, 1, 1.6)),
              label="their projection test"),
       Line2D([], [], color=S.TAN, lw=1.2, ls=(0, (4, 2)), label="HL variants and GiViTI")]
fig.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=S.TICK,
           bbox_to_anchor=(.5, 0.0), columnspacing=3.0, handlelength=2.8, labelspacing=.45)
fig.subplots_adjust(top=.885, bottom=.165, left=.115, right=.99)
S.save(fig, "fig8_liu")
