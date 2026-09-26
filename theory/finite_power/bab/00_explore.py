"""00_explore.py -- timing and structure before building the branch-and-bound (diagnostic only)."""
import os, sys, time, json
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scipy.stats import norm
from crown_core import NET, fwd, K, S_ref
import importlib.util
spec = importlib.util.spec_from_file_location("pc", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "07_pcrown.py"))

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]
v = np.repeat(stripe, K); v /= np.linalg.norm(v)
t_bar = 2.4865397005924965
rng = np.random.default_rng(5)

# active subspace: E[grad grad'] at t v + e
for t in (17.4, 40.0):
    E = rng.normal(size=(4000, 36))
    x = torch.tensor(t * v + E, requires_grad=True)
    S_ref(x).sum().backward()
    G = x.grad.numpy()
    C = G.T @ G / len(G)
    ev = np.linalg.eigvalsh(C)[::-1]
    print("t=%.1f grad-cov eigen share, cumulative:" % t, np.round(np.cumsum(ev)[:10] / ev.sum(), 3))
    val = fwd(NET, t * v + E)
    print("   S quantiles 1/2.5/50%%: %s   mean grad norm %.3f" % (np.round(np.quantile(val, [.01, .025, .5]), 3), np.linalg.norm(G, axis=1).mean()))
    # first-layer stability at z = 4.47
    W0, b0 = NET[0]["W"], NET[0]["b"]
    mu = W0 @ (t * v) + b0; sd = np.linalg.norm(W0, axis=1)
    print("   first-layer |mu|/sd > 4.47 : %d of %d; >2: %d" % ((np.abs(mu / sd) > 4.47).sum(), len(mu), (np.abs(mu / sd) > 2).sum()))
    # empirical unstable counts per layer
    _, pres = fwd(NET, t * v + E, keep=True)
    print("   empirically unstable ReLUs per layer:", [int(((p.min(1) < 0) & (p.max(1) > 0) & l["relu"]).sum()) for p, l in zip(pres[:-1], NET[:-1])])
