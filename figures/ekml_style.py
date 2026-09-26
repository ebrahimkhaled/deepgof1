"""ekml_style.py -- the paper's shared figure language (2026-08-21)

ML-conference clarity (clean spines, direct labels, generous whitespace) carried in the EK palette,
so the figures are distinctive rather than default-matplotlib. Standing rules baked in:
- every figure contains the object that would falsify its claim (nominal band, SE bar, boundary)
- one coordinate system or a lattice of identical ones; no page furniture
- colorblind-robust: the two protagonist colors differ in lightness, not only hue
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# --- the palette -------------------------------------------------------------------------
DEEP    = "#0B3B33"   # the learned test (protagonist)
DEEP_2  = "#4E9A8A"   # learned, secondary arm
CLASSIC = "#A9522F"   # the classical battery (antagonist, rust)
TAN     = "#B8A88A"
GREY    = "#8A8175"
GRID    = "#E8E2D8"
INK     = "#1F2430"
BAND    = "#0B3B33"   # nominal/acceptance bands at low alpha

FIGDIR = Path(__file__).resolve().parent / "out"          # archive: figures/out/
FIGDIR.mkdir(parents=True, exist_ok=True)

def use():
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
        "font.size": 10.5, "axes.titlesize": 12, "axes.labelsize": 10.5,
        "axes.titleweight": "semibold",
        "axes.edgecolor": INK, "axes.linewidth": 0.9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7, "grid.alpha": 0.8,
        "axes.axisbelow": True,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "axes.labelcolor": INK,
        "legend.frameon": False, "legend.fontsize": 9.5,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.pad_inches": 0.08,
    })

TEXTWIDTH_IN = 6.5          # article, 12pt, margin=1in -> \textwidth = 6.5in

_K = 1.0                    # current canvas-to-page scale factor, set by fit_to_textwidth()


def fit_to_textwidth(figwidth, base=8.4, textwidth=TEXTWIDTH_IN):
    r"""Make font sizes mean POINTS ON THE PRINTED PAGE, not points on the canvas.

    THE BUG THIS EXISTS TO KILL (found 2026-08-26, after it had shipped).
    Every figure here is placed with `\includegraphics[width=\textwidth]`, so a canvas wider than
    \textwidth is scaled DOWN. Font sizes are baked into the PDF in absolute points and do not
    scale with the drawing, so the printed size is

        effective pt = fontsize x (textwidth / figwidth)

    At the old figsize of 17.6in that multiplier was 0.37 and a 7pt tick label reached the reader
    at **2.6pt**. Four of the nine figures were below 6pt; the architecture diagram was the worst
    thing in the paper. Nothing in the build warns about it -- the PNG looks perfect on screen,
    because on screen nobody is scaling it by 0.37.

    Call this immediately after use() and before creating the figure, passing the SAME width you
    will pass to figsize. It rescales the rcParams defaults, and `pt()` rescales the explicit
    `fontsize=` arguments, so every number you write in a figure script is what the reader gets.

    `base` is the body size for tick labels and axis titles. 8.4pt against 12pt body text is the
    usual journal relation -- caption-sized, clearly subordinate, still readable. Do not go below
    7pt effective for anything a referee has to read.
    """
    global _K
    _K = figwidth / textwidth
    mpl.rcParams.update({
        "font.size":      base * _K,
        "axes.labelsize": base * _K,
        "axes.titlesize": (base + 1.4) * _K,
        "xtick.labelsize": (base - 0.5) * _K,
        "ytick.labelsize": (base - 0.5) * _K,
        "legend.fontsize": (base - 0.5) * _K,
    })
    return _K


def pt(size_on_page):
    """Convert a printed-page point size into the canvas point size to pass to `fontsize=`.

    Write `fontsize=S.pt(8.6)` and the reader gets 8.6pt, whatever the canvas width is."""
    return size_on_page * _K


def caption(fig, text, y=None, fontsize=8.6, x=0.008, chars_per_inch=15.5):
    """Deck line under the title, WRAPPED to the figure's own width.

    WHY THIS EXISTS (2026-08-25). Every figure carried its argument in a long one-line
    `fig.text(...)`. With `savefig.bbox = "tight"` matplotlib grows the saved canvas until it
    contains every artist -- including text running far past the right edge -- so a 12.4x8.2in
    figure was being written out at aspect 3.9 with the panels crushed into the left third.
    Wrapping the deck to the figure width fixes the aspect for all nine figures at once.

    `chars_per_inch` is calibrated for the ~8.6pt deck in this style; drop it if you raise the
    font size. Returns the wrapped string as well as drawing it.
    """
    import textwrap
    width = max(40, int(fig.get_figwidth() * chars_per_inch * (8.6 / fontsize)))
    wrapped = "\n".join(textwrap.wrap(text, width=width))
    n = wrapped.count("\n") + 1
    H = fig.get_figheight()
    if y is None:                       # TOP of the deck, a fixed step below the suptitle
        y = 1 - 0.050 * (8.2 / H)
    fig.text(x, y, wrapped, fontsize=fontsize, color=GREY, va="top", linespacing=1.45)
    # The deck is drawn va="top", so it grows DOWNWARD from y. Return where it ends so the caller
    # can start the axes below it -- computing y as if it were the bottom (the first version of
    # this helper) pushed four-line decks straight through the top row of panels.
    line_h = (fontsize * 1.45 / 72.0) / H
    return y - n * line_h


def legend(fig, handles, ncol=None, y=-0.004, fontsize=9.6):
    """One figure-level legend, below the lattice. Prefer this over end-of-line labels whenever
    three or more series converge -- direct labels are cleaner until they collide, and then they
    are worse than a legend (fig5 shipped with three labels stacked on one point)."""
    return fig.legend(handles=handles, loc="lower center", ncol=ncol or len(handles),
                      frameon=False, fontsize=fontsize, bbox_to_anchor=(.5, y))


def plain_log_ticks(ax, ticks, axis="x"):
    """Log axis showing ONLY the ticks asked for, as plain integers.

    matplotlib keeps its own minor locator on a log axis, so a hand-set major list still gets
    "6 x 10^1" / "3 x 10^2" drawn underneath it -- they collided with the real labels in fig6.
    """
    import matplotlib.ticker as mt
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(mt.FixedLocator(list(ticks)))
    a.set_major_formatter(mt.FixedFormatter([f"{t:g}" for t in ticks]))
    a.set_minor_locator(mt.NullLocator())
    a.set_minor_formatter(mt.NullFormatter())


def save(fig, name, expect_aspect=None):
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"{name}.{ext}")
    try:                                # report the SAVED aspect: tight-bbox can differ from figsize
        from PIL import Image
        w, h = Image.open(FIGDIR / f"{name}.png").size
        want = fig.get_figwidth() / fig.get_figheight()
        note = f"  {w}x{h}  aspect {w/h:.2f} (declared {want:.2f})"
        if w / h > want * 1.18:         # tight-bbox grew the canvas well past the declared shape
            note += "   <-- OVERFLOW: an artist runs past the figure edge; wrap it or move it in"
    except Exception:
        note = ""
    print(f"WROTE figures/smalln/{name}.png|.pdf{note}")

def endlab(ax, x, y, text, color, dx=6, fontsize=9.5, weight="semibold", va="center"):
    """direct end-of-line labeling instead of a legend box (the ML-paper convention)"""
    ax.annotate(text, (x, y), xytext=(dx, 0), textcoords="offset points",
                color=color, fontsize=fontsize, fontweight=weight, va=va)
