"""Referee check: the projection rate with Escanciano's printed (angle-only, degree-zero) weight
F0(al, be, ga) = pi - arccos(ga / sqrt(al be)).  Expected: p^2 sigma_k -> Q0 = F0ab (trM)^2 + F0gg/2 ||M||^2,
p nu_k in [sqrt(2 (Ev)^2 C0), sqrt(C0 / 8)], C0 = 4 F0ab^2 + F0gg^2 / 2."""
import os, sys, json, time
os.environ["IR_THREADS"] = "4"
H = r"/path/to/project\paper_deepgof\theory\local_power\irrelevant"
sys.path.insert(0, H)
import numpy as np
from ir_core import GridDesign, gfun, BASE, DEPS, project_kernel, ref_sample, expit

R = json.load(open(os.path.join(H, "PROJ_RATE.json")))
c = R["constants"]; Ev = c["Ev"]
MM = {d: np.array(c["M"][d]) for d in DEPS}

def F0(al, be, ga):
    S = np.sqrt(al * be)
    rho = np.clip(np.where(S > 1e-12, ga / np.where(S > 1e-12, S, 1.0), 0.0), -1.0, 1.0)
    return np.pi - np.arccos(rho)

eps = 1e-3
fab = (F0(2 + eps, 2 + eps, 1) - F0(2 + eps, 2 - eps, 1) - F0(2 - eps, 2 + eps, 1) + F0(2 - eps, 2 - eps, 1)) / (4 * eps * eps)
fgg = (F0(2, 2, 1 + eps) - 2 * F0(2, 2, 1) + F0(2, 2, 1 - eps)) / (eps * eps)
s3 = np.sqrt(3)
print("F0ab fd %.6f analytic 1/(12 sqrt3) %.6f | F0gg fd %.6f analytic 1/(3 sqrt3) %.6f" % (fab, 1 / (12 * s3), fgg, 1 / (3 * s3)))
fab, fgg = 1 / (12 * s3), 1 / (3 * s3)
C0 = 4 * fab ** 2 + 0.5 * fgg ** 2
Q0 = {d: fab * np.trace(MM[d]) ** 2 + 0.5 * fgg * (MM[d] ** 2).sum() for d in DEPS}
lo, hi = np.sqrt(2 * Ev ** 2 * C0), np.sqrt(C0 / 8)
print("C0 %.5f  Q0 %s  p nu band [%.4f, %.4f]  p e band: %s" % (C0, Q0, lo, hi,
      {d: (round(Q0[d] / hi, 4), round(Q0[d] / lo, 4)) for d in DEPS}))

# drift, as in 04_proj_rate.signal with F0 in place of F
des0 = GridDesign(0, cirr=0.0)
BB = {}
for dep in DEPS:
    a, ahl, zh, b = des0.drift(lambda x1, x2, pi0, v, u, dep=dep: v * gfun(dep, x1, x2))
    BB[dep] = des0.Iinv @ b

def htilde(x):
    X = np.column_stack([np.ones(len(x)), x]); v = expit(X @ BASE); v = v * (1 - v)
    P = X @ np.linalg.solve(X.T @ X, X.T)
    out = {}
    for dep in DEPS:
        h = v * gfun(dep, x[:, 0], x[:, 1]) - v * (X @ BB[dep]); out[dep] = h - P @ h
    return out

def signal(p, seed, MS=1500, NW=16):
    k = p - 2; rng = np.random.default_rng(seed)
    xa, xb = rng.normal(size=(MS, 2)), rng.normal(size=(MS, 2)); Ha, Hb = htilde(xa), htilde(xb)
    per = {d: [] for d in DEPS}
    for w in range(NW):
        zR = rng.normal(size=2); xI, xIp, zI = rng.normal(size=(3, k)); aI, bI = xI - zI, xIp - zI
        da, db = xa - zR, xb - zR
        F = F0(aI @ aI + (da * da).sum(1)[:, None], bI @ bI + (db * db).sum(1)[None, :], aI @ bI + da @ db.T)
        for d in DEPS:
            per[d].append(float(Ha[d] @ F @ Hb[d]) / MS ** 2)
    return {d: np.mean(per[d]) for d in DEPS}

t0 = time.time()
for p in (20, 50, 200):
    vals = {d: [] for d in DEPS}
    for seed in range(100, 108):
        s = signal(p, seed)
        for d in DEPS: vals[d].append(p * p * s[d])
    print("p %d: p^2 sigma  " % p + "  ".join("%s %.5f (SE %.5f) vs Q0 %.5f" % (d, np.mean(vals[d]), np.std(vals[d], ddof=1) / np.sqrt(8), Q0[d]) for d in DEPS)
          + "  (%.0f s)" % (time.time() - t0)); sys.stdout.flush()

def kern0(x):
    X = (x - x.mean(0)) / x.std(0, ddof=1); n = len(X); A = np.zeros((n, n))
    for r in range(n):
        D = X - X[r]; nm = np.sqrt((D * D).sum(1)); G = D @ D.T; NN = np.outer(nm, nm)
        cs = np.where(NN > 1e-12, G / np.where(NN > 1e-12, NN, 1), 0.0)
        A += np.pi - np.arccos(np.clip(cs, -1, 1))
    return A

MN = 1000
for p in (10, 20, 50):
    x, X, pi0, v = ref_sample(p - 2, MN, 5000 + p, cirr=0.0)
    Rc, EQ, VarQ = project_kernel(kern0(x), X, v); sv = np.sqrt(v)
    Bm = sv[:, None] * Rc * sv[None, :]; Bm[np.diag_indices(MN)] = 0.0
    nu = np.sqrt(2 * (Bm ** 2).sum()) / MN ** 2
    print("p %d: p nu (U) %.4f  band [%.4f, %.4f]  (%.0f s)" % (p, p * nu, lo, hi, time.time() - t0)); sys.stdout.flush()
