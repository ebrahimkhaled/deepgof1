"""08_slopes.py -- does the lower-slope rule of the unstable ReLUs matter for the probabilistic
certificate? Compares, at the same point, CROWN-adaptive (0/1), parallel (lower slope = upper slope
u/(u-l)), and an alpha-optimised final bound. Any slope in [0, 1] is valid, so every variant is sound."""
import os, json, time
import numpy as np
import torch
from scipy.stats import norm
from crown_core import NET, K

torch.set_default_dtype(torch.float64)
HERE = os.path.dirname(os.path.abspath(__file__))
t_bar = json.load(open(os.path.join(HERE, "CERTIFIED.json")))["t_bar"]
N_RELU = int(sum(l["relu"].sum() for l in NET))
Z_INT = norm.ppf(1 - 0.025 / (2 * N_RELU))

def relax(l, u, rl, mode):
    dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
    off = rl & (u <= 0); dl[off] = 0; du[off] = 0
    un = rl & (l < 0) & (u > 0)
    du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
    dl[un] = du[un] if mode == "parallel" else (u[un] >= -l[un]).astype(float)
    return dl, du, cu, un

def pcrown(c0, mode):
    bnds = []
    def backsub(t, sign):
        Lam = sign * NET[t]["W"]; const = sign * NET[t]["b"].copy()
        for s in range(t - 1, -1, -1):
            dl, du, cu, _ = relax(*bnds[s], NET[s]["relu"], mode)
            Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
            const = const + Ln @ cu
            Lam = Lp * dl[None, :] + Ln * du[None, :]
            const = const + Lam @ NET[s]["b"]
            Lam = Lam @ NET[s]["W"]
        return Lam @ c0 + const, np.linalg.norm(Lam, axis=1)
    for t in range(len(NET) - 1):
        a_lo, s_lo = backsub(t, 1.0); a_up, s_up = backsub(t, -1.0)
        bnds.append((a_lo - Z_INT * s_lo, -a_up + Z_INT * s_up))
    return bnds

def final_bound(c0, bnds, opt_steps=0, mode="adaptive"):
    """Final lower function a + lam'e; score = a - z_f * ||lam|| (the certified 97.5% lower quantile)."""
    Wt = [torch.tensor(l["W"]) for l in NET]; bt = [torch.tensor(l["b"]) for l in NET]
    parts = [relax(*bnds[s], NET[s]["relu"], mode) for s in range(len(NET) - 1)]
    logit = [torch.where(torch.tensor(p[0]) > 0.5, 3.0, -3.0).requires_grad_(True) for p in parts]
    zf = norm.ppf(0.975); c0t = torch.tensor(c0)
    def val():
        T = len(NET) - 1
        Lam = Wt[T].clone(); const = bt[T].clone()
        for s in range(T - 1, -1, -1):
            dl, du, cu, un = [torch.tensor(q) for q in parts[s]]
            if opt_steps: dl = torch.where(un, torch.sigmoid(logit[s]), dl)
            Lp, Ln = Lam.clamp(min=0), Lam.clamp(max=0)
            const = const + Ln @ cu
            Lam = Lp * dl + Ln * du
            const = const + Lam @ bt[s]
            Lam = Lam @ Wt[s]
        return (Lam @ c0t + const)[0] - zf * Lam.norm()
    best = float(val())
    if opt_steps:
        opt = torch.optim.Adam(logit, lr=0.2)
        for _ in range(opt_steps):
            opt.zero_grad(); v = val(); (-v).backward(); opt.step(); best = max(best, float(v))
    return best

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]
v = np.repeat(stripe, K); v /= np.linalg.norm(v)
print("certified 2.5%%-quantile lower bound of S(t v + e); need > t_bar = %.3f" % t_bar)
print("%6s %12s %12s %14s %14s" % ("t", "adaptive", "parallel", "adapt+alphaF", "par+alphaF"))
for t in (15, 20, 25, 30, 40):
    b_a = pcrown(t * v, "adaptive"); b_p = pcrown(t * v, "parallel")
    print("%6.0f %12.2f %12.2f %14.2f %14.2f" % (t, final_bound(t * v, b_a), final_bound(t * v, b_p, mode="parallel"),
          final_bound(t * v, b_a, 150), final_bound(t * v, b_p, 150, "parallel")), flush=True)
