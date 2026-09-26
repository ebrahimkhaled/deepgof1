"""10_oracle_bounds.py -- DIAGNOSTIC ONLY (not a certificate): where is the probabilistic CROWN's slack?

Replace every intermediate bound by an 'oracle' band -- the empirical min/max of that pre-activation over
400,000 Gaussian draws -- and keep everything else (the ReLU relaxations and the Gaussian tail of the final
linear lower function). If even these near-perfect bands leave the certificate far above the truth, the
slack is in the triangle relaxation itself (only branch-and-bound on ReLUs helps). If the oracle version
is close to the truth, the fight is to certify tighter intermediate bands."""
import json, os
import numpy as np
from scipy.stats import norm
from crown_core import NET, fwd, K

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "CERTIFIED.json"))); t_bar = C["t_bar"]
rng = np.random.default_rng(99)

def final_lower(c0, bnds):
    T = len(NET) - 1
    Lam = NET[T]["W"].copy(); const = NET[T]["b"].copy()
    for s in range(T - 1, -1, -1):
        l, u = bnds[s]; rl = NET[s]["relu"]
        dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
        off = rl & (u <= 0); dl[off] = 0; du[off] = 0
        un = rl & (l < 0) & (u > 0)
        du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
        dl[un] = (u[un] >= -l[un]).astype(float)
        Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
        const = const + Ln @ cu
        Lam = Lp * dl[None, :] + Ln * du[None, :]
        const = const + Lam @ NET[s]["b"]
        Lam = Lam @ NET[s]["W"]
    return float((Lam @ c0 + const)[0]), float(np.linalg.norm(Lam))

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]
v = np.repeat(stripe, K); v /= np.linalg.norm(v)
print("t    oracle-band cert. 2.5%% quantile   true 2.5%% quantile   unstable units (oracle)")
for t in (6.0, 9.0, 12.0, 17.4, 25.0):
    lo = [np.full(len(l["b"]), np.inf) for l in NET[:-1]]; hi = [np.full(len(l["b"]), -np.inf) for l in NET[:-1]]
    vals = []
    for _ in range(20):
        val, pres = fwd(NET, t * v + rng.normal(size=(20000, 36)), keep=True)
        vals.append(val)
        for s in range(len(NET) - 1):
            lo[s] = np.minimum(lo[s], pres[s].min(1)); hi[s] = np.maximum(hi[s], pres[s].max(1))
    bnds = list(zip(lo, hi))
    a, sd = final_lower(t * v, bnds)
    unst = sum(int((NET[s]["relu"] & (lo[s] < 0) & (hi[s] > 0)).sum()) for s in range(len(NET) - 1))
    print("%5.1f %20.2f %20.2f %12d" % (t, a - norm.ppf(0.975) * sd, np.quantile(np.concatenate(vals), 0.025), unst))
print("(need > t_bar = %.3f; truth crosses it at t_emp = 5.52 for this direction)" % t_bar)
