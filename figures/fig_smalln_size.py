"""fig_smalln_size.py -- size at nominal .05 against n (20 to 100) for five tests, two covariates, in the house
style (ekml_style_sc). Data: theory/smalln_link/LEVEL.csv (2,000 nulls at n = 20 and 30, 1,000 at 50 and 100).
Colour conventions follow fig8_liu_sc_v2.py: DeepGOF-1 deep green, the projection test ink dash-dot, the classical
tests rust (Stukel and Hosmer-Lemeshow told apart by marker and line), GiViTI tan dashes; the band is the
paper's nominal band [.035, .065].
Output: figures/fig_smalln_size.pdf (+ _colour.png, _grey.png), copied to the SiM figures folder.
"""
import os, sys, csv, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()

rows = list(csv.DictReader(open(S.PAP / "theory" / "smalln_link" / "LEVEL.csv")))
NS = [int(r["n"]) for r in rows]
val = lambda r, k: float(r[k].split()[0])
SERIES = {"DG": "DG", "proj": "proj", "StukelJ": "StukelJ", "GiViTI": "GiViTI", "HL": "HL"}
Y = {k: [val(r, c) for r in rows] for k, c in SERIES.items()}

PROJ_LS = (0, (5, 1.6, 1, 1.6))
STYLE = {"HL":      dict(color=S.CLASSIC, lw=1.1, ls=(0, (1.2, 1.6)), marker="o", ms=3.6, mfc="white", mec=S.CLASSIC, mew=0.9),
         "StukelJ": dict(color=S.CLASSIC, lw=1.1, ls="-", marker="^", ms=4.0, mfc=S.CLASSIC, mec="white", mew=0.6),
         "GiViTI":  dict(color=S.TAN, lw=1.2, ls=(0, (4, 2)), marker="D", ms=3.2, mfc=S.TAN, mec="white", mew=0.6),
         "proj":    dict(color=S.INK, lw=1.3, ls=PROJ_LS, marker="s", ms=3.4, mfc=S.INK, mec="white", mew=0.6),
         "DG":      dict(color=S.DEEP, lw=2.0, ls="-", marker="o", ms=4.2, mfc=S.DEEP, mec="white", mew=0.7)}
LABEL = {"DG": "DeepGOF-1", "proj": "projection (Liu et al.)", "StukelJ": "Stukel (joint)", "GiViTI": "GiViTI",
         "HL": "Hosmer–Lemeshow"}

fig, ax = plt.subplots(figsize=(S.HALF_IN, 2.9))
fig.subplots_adjust(left=0.19, right=0.97, top=0.96, bottom=0.40)
ax.axhspan(.035, .065, color=S.BAND, alpha=.10, lw=0, zorder=0)
ax.axhline(.05, color=S.GRID, lw=0.8, zorder=1)
for k in STYLE:
    ax.plot(NS, Y[k], zorder=5 if k == "DG" else 4, **STYLE[k])
ax.set_xscale("log"); ax.xaxis.set_major_locator(FixedLocator(NS)); ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d" % v)); ax.set_xlim(18, 110)
ax.set_ylim(0.0, 0.10); ax.yaxis.set_major_locator(FixedLocator([0, .025, .05, .075, .10]))
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: ("%.3f" % v).lstrip("0").rstrip("0").rstrip(".") if v else "0"))
ax.set_xlabel("sample size $n$"); ax.set_ylabel("size at nominal .05")
ax.text(105, .0362, "nominal band", ha="right", va="bottom", fontsize=S.TICK - 0.5, color=S.ink_for(S.BAND))
leg = [Line2D([], [], label=LABEL[k], **STYLE[k]) for k in ["DG", "proj", "StukelJ", "GiViTI", "HL"]]
fig.legend(handles=leg, loc="lower center", ncol=2, frameon=False, fontsize=S.TICK - 0.5,
           bbox_to_anchor=(0.55, 0.0), handlelength=2.6, columnspacing=1.2)
S.check_layout(fig)
out = S.PAP / "figures" / "fig_smalln_size"
fig.savefig(str(out) + ".pdf")
fig.savefig(str(out) + "_colour.png", dpi=220)
from PIL import Image
Image.open(str(out) + "_colour.png").convert("L").save(str(out) + "_grey.png")
shutil.copy(str(out) + ".pdf", S.PAP / "manuscript" / "sim" / "figures" / "fig_smalln_size.pdf")
print("WROTE", out.with_suffix(".pdf"), "%.1f x %.1f mm" % tuple(np.array(fig.get_size_inches()) * 25.4))
