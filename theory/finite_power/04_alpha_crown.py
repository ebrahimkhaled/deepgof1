"""04_alpha_crown.py -- tighten the certified detection threshold of 03_certify.py.

CROWN fixes the lower slope of every unstable ReLU at 0 or 1. Any slope in [0, 1] gives a valid lower
bound, so the slopes can be OPTIMISED for the one quantity we certify -- the minimum of S over the
noise ball around t v -- by gradient ascent (alpha-CROWN; Xu et al., 2021). The intermediate bounds
stay the CROWN ones, so the result is still rigorous: it is CROWN's bound with better slopes.

Reuses the network, the pure-ReLU rewrite and CROWN from 03_certify.py (everything up to its
validation block), so the two files certify the same object.

Run:  python -u 04_alpha_crown.py
"""
import os, json, time
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "03_certify.py"), encoding="utf-8").read()
cut = src.index("T0 = time.time()\ncrown(NET")            # definitions only, no validation or search
exec(compile(src[:cut], "03_certify.py", "exec"))

def alpha_lb(Ls, c0, r, bnds, steps=80, lr=0.1):
    """Optimised lower bound of the output over the l2 ball; intermediate bounds bnds are fixed."""
    Wt_ = [torch.tensor(l["W"]) for l in Ls]; bt_ = [torch.tensor(l["b"]) for l in Ls]
    fixed, unst = [], []
    for s in range(len(Ls) - 1):
        l, u = bnds[s]; rl = Ls[s]["relu"]
        dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
        off = rl & (u <= 0); dl[off] = 0; du[off] = 0
        un = rl & (l < 0) & (u > 0)
        du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
        dl[un] = (u[un] >= -l[un]).astype(float)
        fixed.append((torch.tensor(dl), torch.tensor(du), torch.tensor(cu)))
        unst.append(torch.tensor(un))
    # logits for the lower slopes, started at CROWN's choice
    logit = [torch.where(fx[0] > 0.5, 3.0, -3.0).clone().requires_grad_(True) for fx in fixed]
    c0t = torch.tensor(c0)
    def bound():
        T = len(Ls) - 1
        Lam = Wt_[T].clone(); const = bt_[T].clone()
        for s in range(T - 1, -1, -1):
            dl, du, cu = fixed[s]
            dl = torch.where(unst[s], torch.sigmoid(logit[s]), dl)
            Lp, Ln = Lam.clamp(min=0), Lam.clamp(max=0)
            const = const + Ln @ cu
            Lam = Lp * dl + Ln * du
            const = const + Lam @ bt_[s]
            Lam = Lam @ Wt_[s]
        return (Lam @ c0t + const - r * Lam.norm(dim=1))[0]
    opt = torch.optim.Adam(logit, lr=lr)
    best = float(bound())
    for it in range(steps):
        opt.zero_grad(); v = bound(); (-v).backward(); opt.step()
        best = max(best, float(v))
    return max(best, float(bound()))

# ---- soundness: the optimised bound must still sit below every sample and attack ------------
rng2 = np.random.default_rng(777)
for c0, r in [(np.zeros(36), 7.14), (rng2.normal(size=36) * 5, 3.0), (rng2.normal(size=36) * 12, 7.14)]:
    lb0, ub0, bnds = crown(NET, c0, r, want_pre=True)
    lba = alpha_lb(NET, c0, r, bnds)
    pts = sample_ball(c0, r, 20000)
    smin = min(fwd(NET, pts).min(), attack_min(c0, r, iters=60, starts=8))
    print("ball r=%.2f |c0|=%.1f : CROWN lb %.3f  alpha-CROWN lb %.3f  sampled/attack min %.3f"
          % (r, np.linalg.norm(c0), lb0, lba, smin))
    assert lba <= smin + 1e-7, "alpha-CROWN unsound -- bug"

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}          # as in 03_certify.py

C = json.load(open(os.path.join(HERE, "CERTIFIED.json")))
t_bar, R2, rho = C["t_bar"], C["R2"], C["rho"]
print("\ncritical value %.4f, noise radius %.2f, rho %.2f" % (t_bar, R2, rho))

def cert(t, v):
    _, _, bnds = crown(NET, t * v, R2, want_pre=True)
    return alpha_lb(NET, t * v, R2, bnds) > t_bar

rows = []
print("%-10s %9s %9s %9s %10s %10s %8s" % ("direction", "t_crown", "t_alpha", "t_emp",
                                           "n*_crown", "n*_alpha", "n*_emp"))
for row in C["rows"]:
    name = row["direction"]
    v = DIRS[name] / np.linalg.norm(DIRS[name])
    hi_t = row["t_cert"]                          # CROWN certifies here, so alpha-CROWN does too
    grid = hi_t * np.exp(np.linspace(np.log(0.05), 0, 30))
    ok = np.array([cert(t, v) for t in grid])
    bad = np.where(~ok)[0]
    i0 = bad.max() + 1 if len(bad) else 0
    lo_t, hi_t = (grid[i0 - 1] if i0 > 0 else 0.0), grid[i0]
    for _ in range(20):
        mid = 0.5 * (lo_t + hi_t)
        if cert(mid, v): hi_t = mid
        else: lo_t = mid
    te = row["t_emp"]
    rows.append(dict(row, t_alpha=hi_t, n_alpha=(hi_t / rho) ** 2))
    print("%-10s %9.2f %9.2f %9.2f %10.0f %10.0f %8.0f" % (name, row["t_cert"], hi_t, te,
          row["n_cert"], (hi_t / rho) ** 2, row["n_emp"]))

C["rows"] = rows
json.dump(C, open(os.path.join(HERE, "CERTIFIED.json"), "w"), indent=2)
print("\nupdated theory/finite_power/CERTIFIED.json with t_alpha / n_alpha")
