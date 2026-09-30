"""Application figure: the DeepGOF-1 residual map of the SUPPORT in-hospital mortality model over two
continuous covariates, mean arterial pressure x respiratory rate (the pair on which the all-pairs reading of
M0 reaches its largest score), along the modelling ladder M0 -> M1 -> M2. Same colour scale, sign marks
(threshold 1.0) and canvas as Fig. 1 (fig3_maps_sc_signs.py); axes labelled in clinical units.

Input : medical_example/screen_SUPPORT/maps_allpairs.csv, maps_allpairs_bins.csv, maps_allpairs_stats.csv
        (export_maps_allpairs_SUPPORT.R)
Output: figures/fig_support_maps.pdf (+ _colour.png, _grey.png for the visual check), copied to the SiM
        manuscript's figures folder
"""
import sys, os, csv, shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()
K, THRESH, MARK_PT = 6, 1.0, 6.5
SRC = S.PAP / "medical_example" / "screen_SUPPORT"
OUT = S.PAP / "figures"
SIM = S.PAP / "manuscript" / "sim" / "figures"

maps = {m: np.zeros((K, K)) for m in ("M0", "M1", "M2")}
with open(SRC / "maps_allpairs.csv") as f:
    for r in csv.DictReader(f):
        # rows of the image = respiratory-rate sixth, columns = mean arterial pressure sixth
        maps[r["model"]][int(r["resp_bin"]) - 1, int(r["meanbp_bin"]) - 1] = float(r["value"])
bins = list(csv.DictReader(open(SRC / "maps_allpairs_bins.csv")))
st = next(csv.DictReader(open(SRC / "maps_allpairs_stats.csv")))


def lab(s, i, two_lines):
    a, b = (int(float(v)) for v in s.split("|"))
    if i == K - 1:
        return f"≥{a}"
    if i == 0:
        return f"≤{b}"
    return f"{a}–\n{b}" if two_lines else f"{a}–{b}"   # two lines, so neighbouring ranges never touch


xlab = [lab(b["meanbp"], i, True) for i, b in enumerate(bins)]
ylab = [lab(b["resp"], i, False) for i, b in enumerate(bins)]


def fmt_p(p):
    return f"{float(p):.3f}".lstrip("0")


panels = [("M0", "M0: linear in every covariate", st["S0"], st["p0"]),
          ("M1", "M1: spline in blood pressure", st["S1"], st["p1"]),
          ("M2", "M2: splines in seven variables", st["S2"], st["p2"])]

cmap = LinearSegmentedColormap.from_list("ek_div", [S.CLASSIC, "#F4EFE7", S.DEEP])
norm = Normalize(vmin=-4.0, vmax=4.0)
fig, axes = plt.subplots(1, 3, figsize=(S.FULL_IN, 2.75))
fig.subplots_adjust(left=0.085, right=0.885, bottom=0.25, top=0.80, wspace=0.14)
for ax, (key, ttl, sc, pv) in zip(axes, panels):
    M = maps[key]
    im = ax.imshow(M, cmap=cmap, norm=norm, origin="lower", interpolation="nearest", aspect="auto")
    for i in range(K):
        for j in range(K):
            v = M[i, j]
            if abs(v) > THRESH:
                ink = "black" if S._lum(cmap(norm(v))) > 0.18 else "white"
                ax.text(j, i - 0.04, "+" if v > 0 else "−", ha="center", va="center",
                        fontsize=MARK_PT, color=ink, fontweight="bold")
    s = float(sc)
    ax.set_title(f"{ttl}\nscore {s:.2f}".replace("-", "−") + f", $p$ = {fmt_p(pv)}",
                 fontsize=S.LABEL, pad=5, linespacing=1.2)
    ax.set_xticks(range(K)); ax.set_xticklabels(xlab, fontsize=S.TICK)
    ax.set_yticks(range(K)); ax.set_yticklabels(ylab if ax is axes[0] else [], fontsize=S.TICK)
    ax.tick_params(length=0, pad=2); ax.grid(False)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(S.GRID)
    ax.set_xlabel("mean arterial pressure (mmHg)", fontsize=S.LABEL, labelpad=3)
axes[0].set_ylabel("respiratory rate (per min)", fontsize=S.LABEL, labelpad=3)
pos = [a.get_position() for a in axes]
cax = fig.add_axes([pos[-1].x1 + 0.012, pos[-1].y0, 0.011, pos[-1].height])
cb = fig.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=S.TICK, length=2.5, width=0.6)
cb.set_label("standardized residual\nsum per cell", fontsize=S.LABEL, linespacing=1.15)
cb.outline.set_visible(False)

fig.savefig(OUT / "fig_support_maps.pdf")
fig.savefig(OUT / "fig_support_maps_colour.png", dpi=200)
from PIL import Image
Image.open(OUT / "fig_support_maps_colour.png").convert("L").save(OUT / "fig_support_maps_grey.png")
shutil.copy(OUT / "fig_support_maps.pdf", SIM / "fig_support_maps.pdf")
print("WROTE", OUT / "fig_support_maps.pdf", "and copied to", SIM)
