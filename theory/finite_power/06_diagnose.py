"""06_diagnose.py -- where does CROWN lose its tightness? Per-layer bound width against the sampled
range, at the point that matters: c0 = t v with t at the exact-ball floor, r = 7.14."""
import numpy as np
from crown_core import NET, crown, fwd, sample_ball, K

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]
v = np.repeat(stripe, K); v = v / np.linalg.norm(v)
c0, r = 12.73 * v, 7.14
lb, ub, bnds = crown(NET, c0, r, want_pre=True)
val, pres = fwd(NET, sample_ball(c0, r, 40000), keep=True)
names = ["conv0", "conv2", "pool1", "pool2", "conv5", "gmaxA", "gmaxB", "gmaxC", "gmaxD", "head0", "out"]
print("output: CROWN [%.2f, %.2f]  sampled [%.2f, %.2f]" % (lb, ub, val.min(), val.max()))
print("%-7s %6s %9s %9s %10s %10s %7s" % ("layer", "units", "relu", "unstable", "CROWN w", "sample w", "ratio"))
for t, nm in enumerate(names):
    l, u = bnds[t]; rl = NET[t]["relu"]
    w_c = np.median(u - l); w_s = np.median(pres[t].max(1) - pres[t].min(1))
    unst = int((rl & (l < 0) & (u > 0)).sum())
    print("%-7s %6d %9d %9d %10.3f %10.3f %7.1f" % (nm, len(l), int(rl.sum()), unst, w_c, w_s, w_c / max(w_s, 1e-9)))
