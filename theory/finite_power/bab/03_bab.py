"""03_bab.py -- probabilistic branch-and-bound by STRATIFYING the Gaussian input along k fixed directions.

Certificate (see README.md).  U (36 x k) orthonormal, s = U'e ~ N(0, I_k), cells are boxes in s whose
edges are Gaussian quantiles, so P(C) = prod_i (ph_i - pl_i) EXACTLY.  For each cell a band level z_C is
chosen from a fixed grid (a deterministic function of the cell, never of the data), and
    P(S <= t_bar) = sum_C P(C) P(S <= t_bar | C) <= sum_C P(C) * min(1, 2 N_RELU Phi(-z_C) + tail_C)
where tail_C = Phi((t_bar - a_C)/sd_C) from the cell's final CROWN lower function (bab_core.py).
The band charge used is the USAGE-AWARE one (bab_core.cell_pcrown_ua), which is <= 2 N_RELU Phi(-z_C).
Power >= .95 is certified when the right-hand side is <= .05.

Branching: pop the cell with the largest contribution P(C) b_C, bisect it at the conditional median
(equal Gaussian mass halves) along the pinned coordinate with the largest |g_fin,i| x (box width in
probability units), where g_fin = U' lambda_fin of the cell's final lower function.

Usage: python -u 03_bab.py <direction> <basis grad|w0> <k> <t> <max_evals>
"""
import os, sys, json, time, heapq
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import torch
torch.set_num_threads(4)
from scipy.stats import norm
from bab_core import NET, S_ref, DIRS, cell_pcrown, cell_pcrown_ua, cell_fail_bound, T_BAR, N_RELU

ZG = (4.2, 4.25, 4.3, 4.35, 4.4, 4.5, 4.6)
DELTA = 0.05


def basis(name, t, v, seed=1):
    if name == "w0":
        W0 = NET[0]["W"]; w, V = np.linalg.eigh(W0.T @ W0); return V[:, ::-1]
    rng = np.random.default_rng(seed)
    E = rng.normal(size=(4000, 36))
    x = torch.tensor(t * v + E, requires_grad=True)
    S_ref(x).sum().backward()
    G = x.grad.numpy(); w, V = np.linalg.eigh(G.T @ G)
    return V[:, ::-1]


def eval_cell(c0, U, pl, ph):
    """A coordinate still at its full range [0, 1] is left Gaussian (moved into the perp part): only the
    restricted coordinates are pinned.  z is scanned upward and the scan stops once the bound worsens
    (a deterministic rule, so the choice of z depends on the cell only)."""
    act = (pl > 0) | (ph < 1)
    Ua = U[:, act]
    lo, hi = norm.ppf(pl[act]), norm.ppf(ph[act])
    best = (2.0, None, None, None)
    for z in ZG:
        a, sd, ch = cell_pcrown_ua(c0, Ua, lo, hi, z)
        if not np.isfinite(a) or not np.isfinite(ch): b = 1.0
        else: b = min(1.0, ch + (norm.cdf((T_BAR - a) / sd) if sd > 0 else float(a <= T_BAR)))
        if b < best[0]: best = (b, z, a, sd)
        elif best[1] is not None and b > best[0]: break
    return best


def final_g(c0, U, pl, ph, z):
    """g_fin = U' lambda_fin for the branching heuristic (one extra pass)."""
    act = (pl > 0) | (ph < 1)
    lo, hi = norm.ppf(pl[act]), norm.ppf(ph[act])
    a, sd, bnds = cell_pcrown(c0, U[:, act], lo, hi, z, want_bnds=True)
    T = len(NET) - 1
    Lam = NET[T]["W"].copy()
    for s in range(T - 1, -1, -1):
        l, u = bnds[s]; rl = NET[s]["relu"]
        dl = np.ones_like(l); du = np.ones_like(l)
        off = rl & (u <= 0); dl[off] = 0; du[off] = 0
        un = rl & (l < 0) & (u > 0)
        du[un] = u[un] / (u[un] - l[un]); dl[un] = (u[un] >= -l[un]).astype(float)
        Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
        Lam = (Lp * dl[None, :] + Ln * du[None, :]) @ NET[s]["W"]
    return (Lam @ U)[0]


def run(direction, bname, k, t, max_evals, log=print):
    v = DIRS[direction]; c0 = t * v
    U = basis(bname, t, v)[:, :k]
    t0 = time.time()
    pl, ph = np.zeros(k), np.ones(k)
    b, z, a, sd = eval_cell(c0, U, pl, ph)
    cells = [(-(1.0 * b), 0, pl, ph, b, z, a, sd, 1.0)]
    total = b; nev = 1; uid = 1
    hist = [(nev, total)]
    while total > DELTA and nev < max_evals and cells:
        negc, _, pl, ph, b, z, a, sd, pc = heapq.heappop(cells)
        if b <= 0: heapq.heappush(cells, (negc, _, pl, ph, b, z, a, sd, pc)); break
        if k == 0: break
        # branching coordinate
        if z is not None and np.isfinite(a):
            g = np.abs(final_g(c0, U, pl, ph, z))
        else:
            g = np.ones(k)
        wid = ph - pl
        i = int(np.argmax(g * wid + 1e-12 * wid))
        mid = 0.5 * (pl[i] + ph[i])
        total -= pc * b
        for side in (0, 1):
            npl, nph = pl.copy(), ph.copy()
            if side == 0: nph[i] = mid
            else: npl[i] = mid
            npc = float(np.prod(nph - npl))
            nb, nz, na, nsd = eval_cell(c0, U, npl, nph); nev += 1
            total += npc * nb
            heapq.heappush(cells, (-npc * nb, uid, npl, nph, nb, nz, na, nsd, npc)); uid += 1
        hist.append((nev, total))
        if nev % 20 < 2:
            log("   evals %4d  cells %4d  bound %.4f  (%.0fs)" % (nev, len(cells), total, time.time() - t0))
    return dict(direction=direction, basis=bname, k=k, t=t, bound=float(total), certified=bool(total <= DELTA),
                n_cells=len(cells), n_evals=nev, wall_s=time.time() - t0, history=hist,
                cells=[dict(pl=c[2].tolist(), ph=c[3].tolist(), b=c[4], z=c[5], a=c[6], sd=c[7], p=c[8]) for c in cells])


if __name__ == "__main__":
    d, bn, k, t, me = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
    r = run(d, bn, k, t, me)
    print("RESULT %s %s k=%d t=%.2f bound=%.4f certified=%s cells=%d evals=%d wall=%.0fs"
          % (d, bn, k, t, r["bound"], r["certified"], r["n_cells"], r["n_evals"], r["wall_s"]))
    json.dump(r, open("run_%s_%s_k%d_t%g.json" % (d, bn, k, t), "w"), indent=1)
