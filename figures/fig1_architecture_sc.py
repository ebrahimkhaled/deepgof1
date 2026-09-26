r"""Figure 1 at Statistics and Computing print size (174 mm, Arial 8/9 pt, fonts embedded).

Same three-lane diagram as fig1_architecture.py: same boxes, arrows, colours and wording. What
changes, for the journal's artwork rules:
- canvas = 6.85 in (174 mm), placed at that width, so sizes below are printed sizes;
- lane titles and lane-box labels 9 pt, everything else (sub-labels, lane subtitles, the network
  blocks, the arrow labels) 8 pt;
- Arial embedded, including the math. Arial has no circled-plus or tilde operator, so
  "GlobalMax (+) Mean" is set with a plain "+" (the sub-label already reads 32+32 = 64) and the
  model formula uses an ASCII "~";
- grey and light sub-label text darkened/lightened just enough for 4.5:1 contrast on its box.
The data-unit grid is unchanged (136 x 88), and the canvas keeps the units square.
`check_boxes()` confirms every label sits inside its own box.
"""
import os, sys
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()

L9, L8 = S.LABEL, S.TICK
W_IN = S.FULL_IN
fig, ax = plt.subplots(figsize=(W_IN, W_IN * 88 / 136))
fig.subplots_adjust(0, 0, 1, 1)
ax.set_xlim(0, 136); ax.set_ylim(0, 88); ax.axis("off")
DARK = {S.DEEP, S.CLASSIC}
GREY_T = S.ink_for(S.GREY, "#EAF1EF")          # sub-label grey, 4.5:1 on the lightest-tinted box
BOXES = []                                      # (patch-extent in data units, [text artists])


def box(x, y, w, h, label, sub="", fc="white", ec=S.INK, lw=0.9, fs=L9,
        bold=False, dy_lab=3.5, dy_sub=-2.9, fs_sub=L8):
    dark = fc in DARK
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.6",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3))
    t1 = ax.text(x+w/2, y+h/2 + (dy_lab if sub else 0), label, ha="center", va="center",
                 fontsize=fs, color="white" if dark else S.INK, zorder=4,
                 fontweight="semibold" if bold else "normal")
    ts = [t1]
    if sub:
        sc = S.ink_for("#DCE9E4", fc) if dark else S.ink_for(GREY_T, fc)
        ts.append(ax.text(x+w/2, y+h/2+dy_sub, sub, ha="center", va="center", fontsize=fs_sub,
                          color=sc, zorder=4, linespacing=1.25))
    BOXES.append(((x - 0.6, y - 0.6, x + w + 0.6, y + h + 0.6), ts))


def arrow(x1, y1, x2, y2, color=S.INK, lw=1.0, ls="-", z=5):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9,
                                 color=color, linewidth=lw, linestyle=ls, zorder=z,
                                 shrinkA=1, shrinkB=1))


def lane(y, h, color, num, title, subtitle):
    ax.add_patch(Rectangle((1.5, y), 133, h, facecolor=color, alpha=.045, lw=0, zorder=0))
    ax.text(3.4, y+h-3.4, f"{num} · {title}", fontsize=L9, fontweight="semibold",
            color=color, zorder=1, va="center")
    ax.text(3.4, y+h-7.6, subtitle, fontsize=L8, color=S.ink_for(S.GREY), zorder=1, va="center")


# ------------------------------------------------------------------ 1. PRETRAINING lane
LY, LH = 61.0, 27.0
BY, BH = 62.5, 15.0
lane(LY, LH, S.CLASSIC, "1", "PRETRAINING — once, offline, on simulated data",
     "the analyst never runs this step")
box(4.0,   BY, 21, BH, "simulator", "logistic DGP\n~100 random\ndeparture families",
    fc="#FBF7F2", ec=S.CLASSIC)
box(29.0,  BY, 21, BH, "fit  $y$ ~ $x_1+x_2$", "the WRONG model,\non purpose",
    fc="#FBF7F2", ec=S.CLASSIC)
box(54.0,  BY, 23, BH, "residual maps", "null class =\nbootstrap replicates\n(deployment-matched)",
    fc="#FBF7F2", ec=S.CLASSIC)
box(81.0,  BY, 20, BH, "train CNN", "25 epochs\nAdam · BCE", fc="#FBF7F2", ec=S.CLASSIC)
box(105.0, BY, 26, BH, "FROZEN WEIGHTS", "18,273 parameters\nshipped frozen",
    fc=S.CLASSIC, ec=S.CLASSIC, bold=True)
for a, b in ((25, 29), (50, 54), (77, 81), (101, 105)):
    arrow(a, BY+BH/2, b, BY+BH/2, color=S.CLASSIC)

# ------------------------------------------------------------------ 2. THE NETWORK module
ax.text(3.4, 57.4, "2 · THE NETWORK — trained above, applied below",
        fontsize=L9, fontweight="semibold", color=S.DEEP, va="center")
ax.text(3.4, 53.6, "one module, shared by both lanes", fontsize=L8, color=S.ink_for(S.GREY),
        va="center")

MX, MY, MW, MH = 2.5, 31.0, 131.0, 20.0
ax.add_patch(FancyBboxPatch((MX, MY), MW, MH, boxstyle="round,pad=0.4,rounding_size=1.4",
                            facecolor="none", edgecolor=S.DEEP, linewidth=0.9,
                            linestyle=(0, (5, 3)), alpha=.5, zorder=1))
NY, NH = 35.2, 13.0
blocks = [
    (5.0,   15.0, "residual map",       "$1\\times6\\times6$",   "#EAF1EF"),
    (23.0,  14.0, "Conv $3\\times3$",   "$16\\times6\\times6$",  "#DCE9E5"),
    (40.0,  14.0, "Conv $3\\times3$",   "$32\\times6\\times6$",  "#CFE1DC"),
    (57.0,  11.5, "MaxPool 2",          "$32\\times3\\times3$",  "#BFD8D1"),
    (71.5,  14.0, "Conv $3\\times3$",   "$32\\times3\\times3$",  "#CFE1DC"),
    (88.5,  19.0, "GlobalMax + Mean",   "$32+32 = 64$",          "#AECFC6"),
    (110.5, 9.5,  "FC 64",              "ReLU",                  "#9CC5BA"),
    (123.0, 10.5, "score  $S$",         "one scalar",            S.DEEP),
]
for i, (x, w, lab, sub, fc) in enumerate(blocks):
    box(x, NY, w, NH, lab, sub, fc=fc, ec=S.DEEP, lw=1.0, fs=L8, bold=(fc == S.DEEP),
        dy_lab=2.5, dy_sub=-2.5, fs_sub=L8)
    if i:
        px, pw = blocks[i-1][0], blocks[i-1][1]
        arrow(px+pw, NY+NH/2, x, NY+NH/2, color=S.DEEP, lw=1.0)

# ------------------------------------------------------------------ 3. DEPLOYMENT lane
DY, DHt = 1.0, 27.0
DBY, DBH = 2.5, 15.0
lane(DY, DHt, S.DEEP, "3", "DEPLOYMENT — per dataset, weights frozen",
     "everything the analyst actually runs — no retraining")
box(4.0,   DBY, 21, DBH, "your data", "$(y, X)$\n$n = 50$–$200$", fc="#EAF1EF", ec=S.DEEP)
box(29.0,  DBY, 22, DBH, "fit + residual map", "grid over\ncovariate RANKS",
    fc="#EAF1EF", ec=S.DEEP)
box(55.0,  DBY, 17, DBH, "frozen CNN", "$S_{obs}$", fc="#EAF1EF", ec=S.DEEP)
box(76.0,  DBY, 27, DBH, "$B=199$ bootstrap refits",
    "simulate from $\\hat\\pi$\n$\\to S_1 \\ldots S_B$", fc="#EAF1EF", ec=S.DEEP)
box(106.0, DBY, 27.5, DBH, "exact rank $p$-value",
    "$p=(1+\\#\\{S_b \\geq S_{obs}\\})/(B{+}1)$", fc=S.DEEP, ec=S.DEEP, bold=True)
for a, b in ((25, 29), (51, 55), (72, 76), (103, 106)):
    arrow(a, DBY+DBH/2, b, DBY+DBH/2, color=S.DEEP)

# ------------------------------------------------------------------ the two module arrows
arrow(118.0, BY, 118.0, MY+MH+0.4, color=S.CLASSIC, lw=1.2, ls=(0, (4, 2.5)))
ax.text(119.7, (BY + MY+MH)/2, "weights", fontsize=L8, color=S.CLASSIC,
        fontweight="semibold", va="center")
arrow(63.5, MY-0.4, 63.5, DBY+DBH+0.4, color=S.DEEP, lw=1.2, ls=(0, (4, 2.5)))
ax.text(65.2, 29.5, "applied unchanged", fontsize=L8, color=S.DEEP,
        fontweight="semibold", va="center")


def check_boxes():
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = []
    for (x0, y0, x1, y1), ts in BOXES:
        for t in ts:
            e = t.get_window_extent(r)
            (a0, b0), (a1, b1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
            if a0 < x0 or a1 > x1 or b0 < y0 or b1 > y1:
                bad.append(f"OUTSIDE BOX  {t.get_text()!r}  text x[{a0:.1f},{a1:.1f}] "
                           f"y[{b0:.1f},{b1:.1f}]  box x[{x0:.1f},{x1:.1f}] y[{y0:.1f},{y1:.1f}]")
    return bad


for p in check_boxes():
    print("   ", p)
S.save(fig, "fig1_architecture")
for fc in ("#FBF7F2", "#EAF1EF", "#AECFC6", "#9CC5BA", S.CLASSIC, S.DEEP):
    sc = S.ink_for("#DCE9E4", fc) if fc in DARK else S.ink_for(GREY_T, fc)
    print(f"sub-label on {fc}: {sc}  contrast {S.contrast(sc, fc):.2f}")
