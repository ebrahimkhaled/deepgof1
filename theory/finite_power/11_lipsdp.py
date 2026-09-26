"""11_lipsdp.py -- W4, the Lipschitz route: a CERTIFIED global Lipschitz constant of the score.

Why: Gaussian concentration (Tsirelson-Ibragimov-Sudakov) needs no relaxation of any layer:
P{S(tv + e) <= E S(tv + e) - L u} <= exp(-u^2/2) for e ~ N(0, I) and L = Lip(S). The norm product gives
L <= 50.65, far too loose (largest sampled gradient 3.82).

LipSDP-Neuron (Fazlyab, Robey, Hassani, Morari & Pappas, NeurIPS 2019) for a one-hidden-layer block
g(x) = W1 relu(W0 x + b): Lip(g) <= sqrt(rho) whenever, for some diagonal T >= 0,
        [[ -rho I ,  W0' T          ],
         [  T W0  , -2T + W1' W1    ]]   <=  0 .
Written this way every entry is linear in ONE multiplier, so the SDP compiles fast (the first version
built A' diag(t) A and never finished compiling).

Blocks of the score in map units (scaler folded into conv0):
  A = conv2 o relu o conv0            one-hidden-layer LipSDP, 36 -> (576) -> 1152
  relu, 2x2 max-pool                  1-Lipschitz exactly (disjoint windows)
  conv5                               exact operator norm
  relu, global max + mean per channel sqrt(1 + 1/9) exactly (max 1-Lip, mean of 9 is 1/3-Lip, per channel)
  H = head2 o relu o head0            one-hidden-layer LipSDP, 64 -> (64) -> 1
Lip(S) <= Lip(A) * ||conv5|| * sqrt(10/9) * Lip(H). Each LipSDP value is checked against the largest
sampled gradient of its own block (a lower bound).

Run:  python -u 11_lipsdp.py
"""
import os, json, time
import numpy as np
import cvxpy as cp
import torch
from crown_core import A1, A2, A3, b1, H0, H2

HERE = os.path.dirname(os.path.abspath(__file__))
torch.set_default_dtype(torch.float64)

def lipsdp_1hidden(W0, W1, solver="CLARABEL"):
    n0, h = W0.shape[1], W0.shape[0]
    rho = cp.Variable(nonneg=True); t = cp.Variable(h, nonneg=True)
    TW0 = cp.multiply(t[:, None] @ np.ones((1, n0)), W0)          # diag(t) W0, entrywise
    Q = W1.T @ W1
    M = cp.bmat([[-rho * np.eye(n0), TW0.T], [TW0, -2 * cp.diag(t) + Q]])
    prob = cp.Problem(cp.Minimize(rho), [(M + M.T) / 2 << 0])
    t0 = time.time()
    if solver == "SCS":
        prob.solve(solver="SCS", max_iters=50000, eps=1e-7)
    else:
        prob.solve(solver=solver)
    return float(np.sqrt(rho.value)), prob.status, time.time() - t0

def sampled_grad(W0, b0, W1, n=40000, seed=1):
    g = torch.Generator().manual_seed(seed)
    W0t, b0t, W1t = torch.tensor(W0), torch.tensor(b0), torch.tensor(W1)
    best = 0.0
    for s in (0.3, 1.0, 3.0, 10.0):
        x = torch.randn(n, W0.shape[1], generator=g) * s
        # Jacobian norm of W1 relu(W0 x + b): W1 diag(mask) W0, spectral norm per sample (batched)
        mask = ((x @ W0t.T + b0t) > 0).double()
        J = torch.einsum("oh,nh,hi->noi", W1t, mask, W0t)
        best = max(best, float(torch.linalg.matrix_norm(J, 2).max()))
    return best

if __name__ == "__main__":
    from crown_core import hb0
    out = {}
    # the small head block first: a fast check that the SDP is set up right
    LH, st, sec = lipsdp_1hidden(H0, H2.reshape(1, -1))
    pH = np.linalg.norm(H0, 2) * np.linalg.norm(H2)
    eH = sampled_grad(H0, hb0, H2.reshape(1, -1))
    print("block H  head0-relu-head2 : LipSDP %.4f [%s, %.0fs]  norm product %.4f  sampled >= %.4f"
          % (LH, st, sec, pH, eH), flush=True)
    assert LH >= eH - 1e-6, "LipSDP below a sampled gradient -- unsound"
    LA, st, sec = lipsdp_1hidden(A1, A2)
    pA = np.linalg.norm(A1, 2) * np.linalg.norm(A2, 2)
    eA = sampled_grad(A1, b1, A2, n=800)
    print("block A  conv0-relu-conv2 : LipSDP %.4f [%s, %.0fs]  norm product %.4f  sampled >= %.4f"
          % (LA, st, sec, pA, eA), flush=True)
    assert LA >= eA - 1e-6, "LipSDP below a sampled gradient -- unsound"
    n5 = np.linalg.norm(A3, 2); rd = np.sqrt(1 + 1 / 9)
    L = LA * n5 * rd * LH
    print()
    print("conv5 %.3f  readout %.3f" % (n5, rd))
    print("CERTIFIED Lip(S) <= %.3f   (norm product 50.65; largest sampled gradient 3.82)" % L)
    out.update(A=LA, A_product=pA, A_sampled=eA, H=LH, H_product=pH, H_sampled=eH, conv5=n5, readout=rd, L=L)
    json.dump(out, open(os.path.join(HERE, "LIPSDP.json"), "w"), indent=2)
