"""05_ball_limit.py -- how good can ANY certificate over the 95% noise ball be?

03/04 certify min_{||e|| <= r} S(t v + e) > t_bar and need t 10-12x the Monte Carlo truth. That gap has
two possible sources: the RELAXATION (CROWN is loose) or the CERTIFICATE TYPE (the worst case over a ball
that holds 95% of the noise is itself far below the typical case). This separates them.

A projected-gradient attack finds a point in the ball with a low score, so it gives an UPPER bound on the
exact minimum. t_attack := the smallest t at which even the attack cannot push S below t_bar is therefore
a LOWER bound on the t that ANY exact ball certificate needs. If t_attack / t_emp is already far above
sqrt(10) ~ 3.16, then no ball-type certificate, however tight its relaxation, meets the 10x-in-n target,
and that is a proof, not an opinion.

Run:  python -u 05_ball_limit.py
"""
import os, json
import numpy as np
import torch
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "local_power"))
from lp_core import score, W, MU, SD, K
import torch.nn.functional as F

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "CERTIFIED.json")))
t_bar, R2 = C["t_bar"], C["R2"]
torch.set_default_dtype(torch.float64)

def S_t(m):
    z = ((m - MU) / SD).view(-1, 1, K, K)
    h = F.relu(F.conv2d(z, W["conv.0.weight"], W["conv.0.bias"], padding=1))
    h = F.relu(F.conv2d(h, W["conv.2.weight"], W["conv.2.bias"], padding=1))
    h = F.max_pool2d(h, 2)
    h = F.relu(F.conv2d(h, W["conv.5.weight"], W["conv.5.bias"], padding=1))
    e = torch.cat([h.amax(dim=(2, 3)), h.mean(dim=(2, 3))], 1)
    h = F.relu(e @ W["head.0.weight"].T + W["head.0.bias"])
    return (h @ W["head.2.weight"].T + W["head.2.bias"]).squeeze(-1)

def attack_min(c0, r, starts=64, iters=300, seed=0):
    """Batched PGD on the l2 sphere/ball: returns the lowest score found (upper bound on the min)."""
    g = torch.Generator().manual_seed(seed)
    e = torch.randn(starts, 36, generator=g)
    e = e / e.norm(dim=1, keepdim=True) * r * torch.rand(starts, 1, generator=g) ** (1 / 36)
    c0 = torch.as_tensor(c0)
    best = float("inf")
    for it in range(iters):
        e.requires_grad_(True)
        v = S_t(c0[None, :] + e)
        best = min(best, float(v.min()))
        gr, = torch.autograd.grad(v.sum(), e)
        step = 0.05 * r * (1 - it / iters) + 0.005 * r
        with torch.no_grad():
            e = e - step * gr / (gr.norm(dim=1, keepdim=True) + 1e-12)
            nr = e.norm(dim=1, keepdim=True)
            e = torch.where(nr > r, e * r / nr, e)
    with torch.no_grad():
        best = min(best, float(S_t(c0[None, :] + e).min()))
    return best

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}

rows = []
print("%-10s %8s %8s %8s | %10s %10s" % ("direction", "t_emp", "t_attack", "t_cert", "t_att/t_emp", "n-ratio"))
for row in C["rows"]:
    v = DIRS[row["direction"]]; v = v / np.linalg.norm(v)
    lo, hi = 0.0, row.get("t_alpha", row["t_cert"])
    for _ in range(18):
        mid = (lo + hi) / 2
        if attack_min(mid * v, R2) > t_bar: hi = mid
        else: lo = mid
    # hi: the attack can no longer push below t_bar -> a lower bound on any exact ball certificate
    ta = hi
    rows.append(dict(direction=row["direction"], t_emp=row["t_emp"], t_attack=ta,
                     t_cert=row.get("t_alpha", row["t_cert"]), n_ratio_floor=(ta / row["t_emp"]) ** 2))
    print("%-10s %8.2f %8.2f %8.2f | %10.2f %10.1f" % (row["direction"], row["t_emp"], ta,
          rows[-1]["t_cert"], ta / row["t_emp"], (ta / row["t_emp"]) ** 2))
json.dump(rows, open(os.path.join(HERE, "BALL_LIMIT.json"), "w"), indent=2)
print("\nn-ratio = the smallest n*_cert / n*_emp that ANY exact certificate over the 95% l2 ball could reach")
print("(an upper bound on the exact minimum is used, so the true floor is at least this -- up to attack strength)")
