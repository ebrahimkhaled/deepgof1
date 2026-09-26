"""01_zsens.py -- diagnostic: how the probabilistic CROWN final lower function depends on the band z."""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scipy.stats import norm
from crown_core import NET, fwd, K

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]
v = np.repeat(stripe, K); v /= np.linalg.norm(v)
t_bar = 2.4865397005924965

def pcrown(c0, z):
    bnds = []
    def backsub(t, sign):
        Lam = sign * NET[t]["W"]; const = sign * NET[t]["b"].copy()
        for s in range(t - 1, -1, -1):
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
        return Lam @ c0 + const, np.linalg.norm(Lam, axis=1)
    nun = []
    for t in range(len(NET) - 1):
        a_lo, s_lo = backsub(t, 1.0); a_up, s_up = backsub(t, -1.0)
        l, u = a_lo - z * s_lo, -a_up + z * s_up
        bnds.append((l, u)); nun.append(int((NET[t]["relu"] & (l < 0) & (u > 0)).sum()))
    a, sd = backsub(len(NET) - 1, 1.0)
    return float(a[0]), float(sd[0]), nun

for t in (17.4, 40.0):
    print("t = %.1f   S(t v) = %.3f" % (t, fwd(NET, (t * v)[None, :])[0]))
    for z in (0.0, 1.0, 2.0, 3.0, 4.47):
        t0 = time.time(); a, sd, nun = pcrown(t * v, z)
        print("  z=%.2f  a=%9.3f sd=%7.3f  a-1.96sd=%9.3f  unstable=%s  (%.2fs)" % (z, a, sd, a - 1.96 * sd, nun, time.time() - t0))
