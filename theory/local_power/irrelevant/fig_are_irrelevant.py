"""fig_are_irrelevant.py -- local sample-size ratio n*_arm / n*_DeepGOF-1 at limit power one half against the
number k of irrelevant covariates, in the niche design (c_irr = 0.15). Paper arms only: DeepGOF-1 (the reference
line at one), the projection test of Liu et al. (2024) with its printed angle weight on the raw covariates
(PROJ_LIU_RAW.json, 08_proj_liu_raw.py: twelve reference samples of 1,000 draws, off-diagonal discretisation), Stukel and
Hosmer-Lemeshow (ARE_IRRELEVANT.json, exact quadratures). Two panels (interaction, quadratic), exactly 84 mm wide
(no tight bounding box), Arial 8 pt, greyscale-legible (line style and marker carry the arm, not colour). For the
projection test the thin vertical bar is plus/minus one Monte Carlo SE over the reference samples.
Run:  python fig_are_irrelevant.py   -> fig_are_irrelevant.pdf (+ .png preview)
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(os.path.join(HERE, "ARE_IRRELEVANT.json")))
P = json.load(open(os.path.join(HERE, "PROJ_LIU_RAW.json")))["curve"]
KS = J["k"]; V = J["variants"]["study"]

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold"})
ARMS = [("proj", "projection", dict(color="0.0", ls="-", marker="s", mfc="0.0")),
        ("Stukel", "Stukel", dict(color="0.35", ls="--", marker="^", mfc="white")),
        ("HL", "Hosmer–Lemeshow", dict(color="0.55", ls=":", marker="o", mfc="0.55"))]
TITLES = {"inter": "interaction  $x_1x_2$", "quad": "quadratic  $x_1^2-1$"}


def series(dep, arm):
    ks = [k for k in KS if str(k) in P[dep]]
    if arm == "proj":
        y = [P[dep][str(k)]["ratio_of_mean_tau"] for k in ks]
        e = [P[dep][str(k)]["ratio_of_mean_tau_se"] for k in ks]
        return ks, y, e
    return ks, [V[dep][str(k)]["ratio_vs_DG"][arm] for k in ks], None


allv = [v for d in ("inter", "quad") for a, _, _ in ARMS for v in series(d, a)[1]]
lo, hi = min(min(allv) / 1.2, 0.45), max(max(allv) * 1.2, 4.6)
MM = 1 / 25.4
fig, axes = plt.subplots(1, 2, figsize=(84 * MM, 58 * MM), sharey=True)
for ax, dep in zip(axes, ["inter", "quad"]):
    ax.axhline(1.0, color="0.7", lw=0.8, zorder=0)
    for arm, lab, st in ARMS:
        ks, y, e = series(dep, arm)
        ax.plot(ks, y, lw=0.9, ms=3.2, mew=0.7, label=lab, zorder=3, **st)
        if e is not None:
            for k, yy, ee in zip(ks, y, e):
                ax.plot([k, k], [yy - ee, yy + ee], color=st["color"], lw=0.6, zorder=2)
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.set_title(TITLES[dep], loc="left", pad=3)
    ax.set_xticks(KS)
    ax.set_xlim(-0.5, KS[-1] + 0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    ticks = [t for t in (0.25, 0.5, 1, 2, 4, 8) if lo <= t <= hi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
axes[0].set_ylabel("sample-size ratio, test / DeepGOF-1")
axes[1].text(KS[-1] + 0.4, 0.95, "DeepGOF-1", ha="right", va="top", fontsize=8, color="0.4")
fig.supxlabel("irrelevant covariates $k$", fontsize=8, y=0.115)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.55, 0.0), handlelength=2.4,
           columnspacing=1.0, handletextpad=0.4)
fig.subplots_adjust(left=0.13, right=0.985, top=0.92, bottom=0.26, wspace=0.1)
fig.savefig(os.path.join(HERE, "fig_are_irrelevant.pdf"))
fig.savefig(os.path.join(HERE, "fig_are_irrelevant.png"), dpi=300)
print("wrote fig_are_irrelevant.pdf / .png, %.1f x %.1f mm" % tuple(np.array(fig.get_size_inches()) * 25.4))
