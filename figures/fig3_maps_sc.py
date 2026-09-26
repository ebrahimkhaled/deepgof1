"""Figure 3 at Statistics and Computing print size (174 mm, Arial 8/9 pt, fonts embedded).

Same maps as fig3_maps.py -- same simulator, seeds, grid, colour scale -- only the canvas and the
lettering change. Two titles and the colour-bar label are wrapped onto one more line, because at
9 pt on a 1.4 in panel the single-line versions ran into the neighbouring panel.
"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ekml_style_sc as S
S.use()
K, N = 12, 8000


def fit_logit(X, y, iters=60):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1/(1+np.exp(-X@b)); W = p*(1-p)
        b += np.linalg.solve((X*W[:, None]).T@X + 1e-9*np.eye(X.shape[1]), X.T@(y-p))
    return b


def make_map(delta_fn, seed):
    r = np.random.default_rng(seed)
    x1, x2 = r.uniform(-3, 3, N), r.normal(0, 1, N)
    eta = 0.3 + 0.8*x1 - 0.5*x2 + delta_fn(x1, x2)
    y = (r.uniform(size=N) < 1/(1+np.exp(-eta))).astype(float)
    X = np.column_stack([np.ones(N), x1, x2])
    b = fit_logit(X, y)
    ph = 1/(1+np.exp(-X@b)); res = y - ph
    i1 = np.minimum(K-1, (K * (np.argsort(np.argsort(x1))) / N).astype(int))
    i2 = np.minimum(K-1, (K * (np.argsort(np.argsort(x2))) / N).astype(int))
    M = np.zeros((K, K)); V = np.zeros((K, K))
    for a, bb, rr, v in zip(i1, i2, res, ph*(1-ph)):
        M[bb, a] += rr; V[bb, a] += v
    return M / np.sqrt(np.maximum(V, 1e-8))


panels = [
    ("correctly specified\n(what 'healthy'\nlooks like)",       lambda x1, x2: 0*x1),
    ("omitted $x_1^2$\n(curvature stripe)",                     lambda x1, x2: 1.15*(x1**2 - 3)/3),
    ("local bump in $(x_1,x_2)$\n(invisible to\nrisk-ordering)",
     lambda x1, x2: 2.6*np.exp(-((x1-1.4)**2+(x2-0.9)**2)/0.9)),
    ("omitted interaction\n$x_1 x_2$ (saddle)",                 lambda x1, x2: 0.65*x1*x2),
]
cmap = LinearSegmentedColormap.from_list("ek_div", [S.CLASSIC, "#F4EFE7", S.DEEP])

fig, axes = plt.subplots(1, 4, figsize=(S.FULL_IN, 2.30))
fig.subplots_adjust(left=0.045, right=0.875, bottom=0.10, top=0.74, wspace=0.10)
for ax, (ttl, fn), sd in zip(axes, panels, (11, 12, 13, 14)):
    M = make_map(fn, sd)
    im = ax.imshow(M, cmap=cmap, vmin=-3.2, vmax=3.2, origin="lower", interpolation="nearest")
    ax.set_title(ttl, fontsize=S.LABEL, pad=5, linespacing=1.15)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color(S.GRID)
    if ax is axes[0]:
        ax.set_xlabel("rank of $x_1$  $\\rightarrow$", fontsize=S.LABEL, labelpad=3)
        ax.set_ylabel("rank of $x_2$  $\\rightarrow$", fontsize=S.LABEL, labelpad=3)
pos = [a.get_position() for a in axes]
cax = fig.add_axes([pos[-1].x1 + 0.012, pos[-1].y0, 0.011, pos[-1].height])
cb = fig.colorbar(im, cax=cax)
cb.ax.tick_params(labelsize=S.TICK, length=2.5, width=0.6)
cb.set_label("standardized residual\nsum per cell", fontsize=S.LABEL, linespacing=1.15)
cb.outline.set_visible(False)
S.save(fig, "fig3_maps")
