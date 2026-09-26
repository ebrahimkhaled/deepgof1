"""02_pinning_limit.py -- DIAGNOSTIC: the infinite-refinement limit of k-dimensional stratification.

A point cell (lo = hi = s0) is the limit of splitting the k pinned coordinates forever.  If even point
cells at typical s0 fail at a given t, NO amount of k-dim input splitting can certify that t.
Bases: 'grad' = eigenvectors of E[grad S grad S'] at t v + e;  'w0' = eigenvectors of W0'W0
(first-layer weight energy).  Reports, for each k, the conditional fail bound (z optimised) at s0 = 0
and at s0 = -1.645 e_1 (the unfavourable side of the dominant direction)."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
from bab_core import NET, S_ref, DIRS, cell_pcrown, cell_fail_bound, T_BAR

v = DIRS["stripe_x1"]
rng = np.random.default_rng(1)
ZG = np.arange(2.5, 5.51, 0.5)

def grad_basis(t):
    E = rng.normal(size=(4000, 36))
    x = torch.tensor(t * v + E, requires_grad=True)
    S_ref(x).sum().backward()
    G = x.grad.numpy()
    w, V = np.linalg.eigh(G.T @ G)
    return V[:, ::-1]

W0 = NET[0]["W"]
w, V0 = np.linalg.eigh(W0.T @ W0); V0 = V0[:, ::-1]
out = {}
for t in (17.4, 40.0):
    Vg = grad_basis(t)
    for bname, B in (("grad", Vg), ("w0", V0)):
        for k in (0, 1, 2, 4, 8, 16, 24, 36):
            U = B[:, :k]
            res = []
            for s0v in (0.0, -1.645):
                s0 = np.zeros(k)
                if k: s0[0] = s0v
                best = (9, None)
                for z in ZG:
                    a, sd = cell_pcrown(t * v, U, s0, s0, z)
                    fb = cell_fail_bound(a, sd, z)
                    if fb < best[0]: best = (fb, z, a, sd)
                res.append(best)
            print("t=%5.1f basis=%-4s k=%2d | s0=0: bound %.3g (z=%.2f a=%.1f sd=%.2f) | s0=-1.645: bound %.3g (a=%.1f)"
                  % (t, bname, k, res[0][0], res[0][1], res[0][2], res[0][3], res[1][0], res[1][2]), flush=True)
            out["%g_%s_%d" % (t, bname, k)] = [r[0] for r in res]
json.dump(out, open("PINNING_LIMIT.json", "w"), indent=1)
