"""fig_are_irrelevant_sc.py -- the local sample-size ratio figure (SiM Fig. 3) in the house style of the other
figures (ekml_style_sc: Arial, 8/9 pt, grid, semibold titles, 174 mm canvas).

Same data and the same curves as theory/stukel_sstboth/fig_are_irrelevant_sstboth.py (which drew them in grey):
local sample-size ratio n*_test / n*_DeepGOF-1 at limit power one half against the number k of covariates that
carry no signal, in the niche design. Colour conventions follow fig8_liu_sc_v2.py: DeepGOF-1 deep green (here the
reference line at one), the projection test ink dash-dot, the classical tests rust; Stukel and Hosmer-Lemeshow are
told apart by marker and line style, so the figure also reads in greyscale.
Output: figures/fig_are_irrelevant_sc.pdf (+ _colour.png, _grey.png previews)
"""
import os, sys, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullLocator, NullFormatter, FuncFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()

TH = S.PAP / "theory"
J = json.load(open(TH / "local_power" / "irrelevant" / "ARE_IRRELEVANT.json"))
P = json.load(open(TH / "local_power" / "irrelevant" / "PROJ_LIU_RAW.json"))["curve"]
SB = json.load(open(TH / "stukel_sstboth" / "IRR_SSTBOTH.json"))["variants"]["study"]
KS = J["k"]; V = J["variants"]["study"]


def series(dep, test):
    ks = [k for k in KS if str(k) in P[dep]]
    if test == "proj":
        return ks, [P[dep][str(k)]["ratio_of_mean_tau"] for k in ks], [P[dep][str(k)]["ratio_of_mean_tau_se"] for k in ks]
    if test == "Stukel":                                   # SstBoth form of Liu et al. (2024)
        return ks, [SB[dep][str(k)]["ratio_sstboth"] for k in ks], None
    return ks, [V[dep][str(k)]["ratio_vs_DG"][test] for k in ks], None


PROJ_LS = (0, (5, 1.6, 1, 1.6))
STYLE = {"proj":   dict(color=S.INK, lw=1.4, ls=PROJ_LS, marker="s", ms=3.6, mfc=S.INK, mec="white", mew=0.6),
         "Stukel": dict(color=S.CLASSIC, lw=1.1, ls="-", marker="^", ms=4.0, mfc=S.CLASSIC, mec="white", mew=0.6),
         "HL":     dict(color=S.CLASSIC, lw=1.1, ls=(0, (1.2, 1.6)), marker="o", ms=3.6, mfc="white",
                        mec=S.CLASSIC, mew=0.9)}
LABEL = {"proj": "projection test (Liu et al.)", "Stukel": "Stukel (SstBoth)", "HL": "Hosmer–Lemeshow"}
TITLES = {"inter": "interaction  $x_1x_2$", "quad": "quadratic  $x_1^2-1$"}

allv = [v for d in ("inter", "quad") for t in STYLE for v in series(d, t)[1]]
lo, hi = min(min(allv) / 1.2, 0.45), max(max(allv) * 1.2, 4.6)

fig, axes = plt.subplots(1, 2, figsize=(S.FULL_IN, 2.55), sharey=True)
fig.subplots_adjust(left=0.085, right=0.985, top=0.88, bottom=0.30, wspace=0.08)
for ax, dep in zip(axes, ["inter", "quad"]):
    ax.axhline(1.0, color=S.DEEP, lw=1.9, zorder=2)
    for t, st in STYLE.items():
        ks, y, e = series(dep, t)
        ax.plot(ks, y, zorder=4, **st)
        if e is not None:
            for k, yy, ee in zip(ks, y, e):
                ax.plot([k, k], [yy - ee, yy + ee], color=st["color"], lw=0.7, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(lo, hi)
    ax.set_title(TITLES[dep], loc="left", pad=5)
    ax.set_xticks(KS); ax.set_xlim(-0.4, KS[-1] + 0.4)
    ticks = [t for t in (0.25, 0.5, 1, 2, 4, 8) if lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.set_xlabel("added covariates without misfit, $q$")
    # above the reference line at the left edge: every curve is below it there, so nothing collides
    ax.text(-0.25, 1.06, "DeepGOF-1", ha="left", va="bottom", fontsize=S.TICK,
            color=S.ink_for(S.DEEP), fontweight="semibold")
axes[0].set_ylabel("sample-size ratio,\ntest / DeepGOF-1")

leg = [Line2D([], [], color=S.DEEP, lw=1.9, label="DeepGOF-1 (reference)")] + \
      [Line2D([], [], label=LABEL[t], **STYLE[t]) for t in STYLE]
fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=S.TICK,
           bbox_to_anchor=(0.53, 0.0), handlelength=2.8, columnspacing=1.6)

S.check_layout(fig)
out = S.PAP / "figures" / "fig_are_irrelevant_sc"
fig.savefig(str(out) + ".pdf")
fig.savefig(str(out) + "_colour.png", dpi=200)
from PIL import Image
Image.open(str(out) + "_colour.png").convert("L").save(str(out) + "_grey.png")
print("WROTE", out.with_suffix(".pdf"), "%.1f x %.1f mm" % tuple(np.array(fig.get_size_inches()) * 25.4))
