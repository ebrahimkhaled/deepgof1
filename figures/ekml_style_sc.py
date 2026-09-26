"""ekml_style_sc.py -- print-size variant of ekml_style for Statistics and Computing (2026-09-24).

Same palette, same spines, same grid as ekml_style. What changes is only the lettering and the
canvas, to meet the journal's artwork rules:

- the canvas IS the printed size: 174 mm (6.85 in) for a full-width figure, 84 mm (3.31 in) for a
  single column. LaTeX places it at that width, so nothing is rescaled and a font size written here
  is the size the reader gets. savefig.bbox is NOT "tight" -- tight-bbox changes the page size.
- Arial throughout, including mathtext (custom fontset, no fallback font), TrueType-embedded
  (pdf.fonttype 42) so the PDF carries real Arial rather than Type 3 outlines.
- 8 pt for ticks, legends and in-panel annotations; 9 pt for axis labels and panel titles.
- text in colour is checked against its background for a 4.5:1 contrast ratio; `ink_for()`
  darkens (or lightens) a series colour just enough for its text label, the line keeps the colour.

`check_layout(fig)` fails loudly if any text leaves the canvas or two text items overlap.
"""
import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

from ekml_style import DEEP, DEEP_2, CLASSIC, TAN, GREY, GRID, INK, BAND  # noqa: F401  (same palette)

FULL_IN = 174 / 25.4        # 6.850 in
HALF_IN = 84 / 25.4         # 3.307 in
TICK = 8.0                  # ticks, legend, annotations
LABEL = 9.0                 # axis labels, panel titles

PAP = Path(__file__).resolve().parents[1]                     # the archive root
OUTDIR = PAP / "figures" / "out"                               # the PDFs the manuscript includes
PNGDIR = PAP / "figures" / "out"                               # PNG previews for checking
PNGDIR.mkdir(parents=True, exist_ok=True)


def use():
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300,
        "font.family": "sans-serif", "font.sans-serif": ["Arial"],
        "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.sf": "Arial",
        "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
        "mathtext.fallback": None,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.size": TICK, "axes.titlesize": LABEL, "axes.labelsize": LABEL,
        "xtick.labelsize": TICK, "ytick.labelsize": TICK, "legend.fontsize": TICK,
        "axes.titleweight": "semibold",
        "axes.edgecolor": INK, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.alpha": 0.8,
        "axes.axisbelow": True,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK, "axes.labelcolor": INK,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.bbox": None, "savefig.pad_inches": 0.0,
    })


# ------------------------------------------------------------------ contrast
def _lum(c):
    r, g, b = mpl.colors.to_rgb(c)
    f = lambda v: v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(fg, bg="white"):
    a, b = sorted((_lum(fg), _lum(bg)), reverse=True)
    return (a + 0.05) / (b + 0.05)


def ink_for(color, bg="white", target=4.5):
    """The same hue, moved in lightness just far enough to reach `target` contrast on `bg`."""
    if contrast(color, bg) >= target:
        return color
    h, l, s = colorsys.rgb_to_hls(*mpl.colors.to_rgb(color))
    step = -0.01 if _lum(bg) > 0.18 else 0.01
    while 0 < l < 1 and contrast(colorsys.hls_to_rgb(h, l, s), bg) < target:
        l += step
    return mpl.colors.to_hex(colorsys.hls_to_rgb(h, min(max(l, 0), 1), s))


# ------------------------------------------------------------------ helpers kept from ekml_style
def plain_log_ticks(ax, ticks, axis="x"):
    import matplotlib.ticker as mt
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(mt.FixedLocator(list(ticks)))
    a.set_major_formatter(mt.FixedFormatter([f"{t:g}" for t in ticks]))
    a.set_minor_locator(mt.NullLocator())
    a.set_minor_formatter(mt.NullFormatter())


def endlab(ax, x, y, text, color, dx=5, fontsize=TICK, weight="semibold", va="center"):
    return ax.annotate(text, (x, y), xytext=(dx, 0), textcoords="offset points",
                       color=ink_for(color), fontsize=fontsize, fontweight=weight, va=va)


# ------------------------------------------------------------------ checks + save
def _texts(fig):
    """Every text that is actually drawn. Tick labels count only for ticks inside the view limits
    (matplotlib keeps Text objects for off-axis ticks but never draws them)."""
    tick_labels, drawn = set(), []
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                tick_labels.update({id(tk.label1), id(tk.label2)})
            if not (ax.axison and axis.get_visible()):
                continue                         # axis("off"): tick labels are never drawn
            for tk in axis._update_ticks():
                drawn += [lb for lb in (tk.label1, tk.label2) if lb.get_visible()]
    out = [t for t in drawn if t.get_text().strip()]
    for t in fig.findobj(mpl.text.Text):
        if id(t) in tick_labels:
            continue
        if t.get_visible() and t.get_text().strip() and t.get_figure() is fig:
            out.append(t)
    return out


def check_layout(fig, allow=()):
    """Every text inside the canvas; no two texts overlapping (pairs listed in `allow` excepted)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    boxes = [(t, t.get_window_extent(r)) for t in _texts(fig)]
    bad = []
    for t, b in boxes:
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5:
            bad.append(f"CLIPPED  {t.get_text()!r}")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            (t1, b1), (t2, b2) = boxes[i], boxes[j]
            if b1.overlaps(b2):
                ov = (min(b1.x1, b2.x1) - max(b1.x0, b2.x0)) * (min(b1.y1, b2.y1) - max(b1.y0, b2.y0))
                if ov > 1.0 and (t1.get_text(), t2.get_text()) not in allow:
                    bad.append(f"OVERLAP  {t1.get_text()!r} x {t2.get_text()!r}")
    return bad


def save(fig, name, width_in=FULL_IN, allow=()):
    w, h = fig.get_size_inches()
    assert abs(w - width_in) < 1e-6, f"canvas {w:.3f} in, printed width {width_in:.3f} in"
    problems = check_layout(fig, allow)
    fig.savefig(OUTDIR / f"{name}.pdf")
    fig.savefig(PNGDIR / f"{name}.png")
    sizes = sorted({round(t.get_fontsize(), 2) for t in _texts(fig)})
    print(f"WROTE {name}.pdf  {w:.3f} x {h:.3f} in  font sizes {sizes}")
    for p in problems:
        print("   ", p)
    return problems
