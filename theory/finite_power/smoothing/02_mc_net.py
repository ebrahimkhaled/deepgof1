"""02_mc_net.py -- Monte Carlo at a fixed net of signal strengths along one direction.  Resumable:
every net point is appended to net_<direction>.jsonl as soon as it is done; points already present are
skipped.  Run: python -u 02_mc_net.py <direction>

The net is a deterministic function of the direction (fixed before any draw), so the Bonferroni split of
ETA_NET over its J points is legitimate.  At each point t_j, M_j independent maps t_j v + e are scored and
the exceedance counts over every threshold (T_BAR for the fixed test; the layer-cake critical values
c_w for the bootstrap test) are recorded.  Nothing else is stored.
"""
import os, sys, json, time
import numpy as np
from smooth_core import (forward, DIRS, T_BAR, T_HAT, T_EMP_PAPER, LEVELS, HERE, append_jsonl, read_jsonl)

def design(direction):
    te = T_EMP_PAPER[direction]
    pts = []
    t = round(0.3 * te, 3)
    while t < 1.5 * te:                              # dense: spacing .1, M = 100,000
        pts.append((round(t, 3), 100_000)); t += 0.1
    while t < 4.0 * te:                              # middle: spacing .5, M = 20,000
        pts.append((round(t, 3), 20_000)); t += 0.5
    while t <= 250.0 + 1e-9:                         # far: spacing 1.25, M = 4,000
        pts.append((round(t, 3), 4_000)); t += 1.25
    return pts

if __name__ == "__main__":
    d = sys.argv[1]; v = DIRS[d]
    ucb = json.load(open(os.path.join(HERE, "NULL_UCB.json")))["ucb"]
    cvals = np.array([ucb[str(w)] for w in LEVELS])
    path = os.path.join(HERE, "net_%s.jsonl" % d)
    done = {r["t"] for r in read_jsonl(path)}
    pts = design(d)
    print("%s: %d net points, %d done" % (d, len(pts), len(done)), flush=True)
    t0 = time.time()
    for i, (t, M) in enumerate(pts):
        if t in done:
            continue
        seed = int(round(t * 1000))                  # deterministic (Python's hash of str is salted)
        rng = np.random.default_rng([seed, sum(map(ord, d)), 20260924])
        s = forward(t * v + rng.normal(size=(M, 36)))
        rec = dict(t=t, M=M, n_fixed=int((s > T_BAR).sum()), n_hat=int((s > T_HAT).sum()),
                   n_levels=[int((s > c).sum()) for c in cvals], mean_s=float(s.mean()), seed=seed)
        append_jsonl(path, rec)
        if i % 10 == 0:
            print("  t=%7.3f M=%6d  power(fixed) %.4f   [%d/%d, %.0f s]"
                  % (t, M, rec["n_fixed"] / M, i + 1, len(pts), time.time() - t0), flush=True)
    print("%s net done in %.0f s" % (d, time.time() - t0))
