"""09_alpha_full.py -- probabilistic CROWN with alpha-OPTIMISED intermediate bounds.

In 07/08 the intermediate bounds of conv5 and the head are 40-124x wider than the sampled ranges. Every
lower slope in [0, 1] of an unstable ReLU gives a valid linear lower bound, so for each intermediate
bound we may choose the slopes that make THAT bound tightest (alpha-CROWN, Xu et al. 2021), with
separate slopes per bounded unit. Soundness is unchanged: each bound is still 'linear function of e
minus z_int times its Gaussian sd', valid on the event that all earlier bands hold, and the union bound
over all ReLU units still covers every band.

Optimised: every layer from `OPT_FROM` on (default 2: all but the two first conv layers, whose bounds are
cheap and nearly exact already). Validation as in 07: sampled band violations <= delta_int, and the
certified power at the certified point below the Monte Carlo power.

Run:  python -u 09_alpha_full.py [direction]
"""
import os, sys, json, time
import numpy as np
import torch
from scipy.stats import norm
from crown_core import NET, fwd, K

torch.set_default_dtype(torch.float64)
torch.set_num_threads(8)
HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "CERTIFIED.json")))
t_bar = C["t_bar"]
N_RELU = int(sum(l["relu"].sum() for l in NET))
Z_INT = norm.ppf(1 - 0.025 / (2 * N_RELU)); Z_F = norm.ppf(0.975)
OPT_FROM = int(os.environ.get("OPT_FROM", "2")); STEPS = int(os.environ.get("STEPS", "40"))
Wt = [torch.tensor(l["W"]) for l in NET]; bt = [torch.tensor(l["b"]) for l in NET]
RL = [torch.tensor(l["relu"]) for l in NET]

def relax(l, u, rl):
    dl = torch.ones_like(l); du = torch.ones_like(l); cu = torch.zeros_like(l)
    off = rl & (u <= 0); dl = torch.where(off, 0.0, dl); du = torch.where(off, 0.0, du)
    un = rl & (l < 0) & (u > 0)
    du = torch.where(un, u / (u - l).clamp(min=1e-12), du)
    cu = torch.where(un, -u * l / (u - l).clamp(min=1e-12), cu)
    dl = torch.where(un, (u >= -l).double(), dl)
    return dl, du, cu, un

def backsub(t, sign, bnds, c0, logits=None):
    """Rows = units of layer t. logits[s]: (rows x n_s) free lower slopes for layer s (or None)."""
    Lam = sign * Wt[t]; const = sign * bt[t].clone()
    for s in range(t - 1, -1, -1):
        dl, du, cu, un = relax(*bnds[s], RL[s])
        Lp, Ln = Lam.clamp(min=0), Lam.clamp(max=0)
        const = const + Ln @ cu
        if logits is not None and logits[s] is not None:
            dlm = torch.where(un[None, :], torch.sigmoid(logits[s]), dl[None, :])
            Lam = Lp * dlm + Ln * du[None, :]
        else:
            Lam = Lp * dl[None, :] + Ln * du[None, :]
        const = const + Lam @ bt[s]
        Lam = Lam @ Wt[s]
    return Lam @ c0 + const, Lam.norm(dim=1)

def bound_layer(t, sign, bnds, c0, z, optimise):
    a, s = backsub(t, sign, bnds, c0)
    best = (a - z * s).detach()
    if not optimise:
        return best
    rows = Wt[t].shape[0]
    logits = [None] * t
    for s_ in range(t):
        _, _, _, un = relax(*bnds[s_], RL[s_])
        if un.any():
            init = torch.where((bnds[s_][1] >= -bnds[s_][0]), 3.0, -3.0)
            logits[s_] = init[None, :].repeat(rows, 1).clone().requires_grad_(True)
    params = [p for p in logits if p is not None]
    if not params:
        return best
    opt = torch.optim.Adam(params, lr=0.25)
    for _ in range(STEPS):
        opt.zero_grad()
        a, s = backsub(t, sign, bnds, c0, logits)
        val = a - z * s
        (-val.sum()).backward(); opt.step()
        best = torch.maximum(best, val.detach())
    return best

def pcert(c0_np, verbose=False):
    c0 = torch.tensor(c0_np); bnds = []
    for t in range(len(NET) - 1):
        opt = t >= OPT_FROM and bool(RL[t].any())
        lo = bound_layer(t, 1.0, bnds, c0, Z_INT, opt)
        up = -bound_layer(t, -1.0, bnds, c0, Z_INT, opt)
        bnds.append((lo, up))
        if verbose:
            un = int((RL[t] & (lo < 0) & (up > 0)).sum())
            print("   layer %2d  median width %9.3f  unstable %d" % (t, float((up - lo).median()), un), flush=True)
    fin = bound_layer(len(NET) - 1, 1.0, bnds, c0, Z_F, True)
    return float(fin[0]), bnds

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}

if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "stripe_x1"
    v = DIRS[name] / np.linalg.norm(DIRS[name])
    for t in (15.0, 20.0, 25.0):
        t0 = time.time()
        q, _ = pcert(t * v, verbose=(t == 15.0))
        print("%s t=%5.1f  certified 2.5%%-quantile lower bound %9.2f   (need > %.3f)  [%.0fs]"
              % (name, t, q, t_bar, time.time() - t0), flush=True)
