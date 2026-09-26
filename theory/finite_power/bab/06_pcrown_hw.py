"""06_pcrown_hw.py -- probabilistic CROWN (no splitting) with the band budget charged by Hunter-Worsley
instead of the union bound, and z chosen on a grid.  Sound: for each (t, z) the band forms are a
deterministic function of (t, z); P(E^c) <= HW(z); P(S <= t_bar) <= HW(z) + Phi((t_bar - a)/sd).

Usage: python -u 06_pcrown_hw.py <direction> <t1> <t2> ...
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
from bab_core import NET, DIRS, cell_pcrown, T_BAR, N_RELU
from hw import hw_bound

E0 = np.zeros((36, 0)); Z0 = np.zeros(0)


def evaluate(direction, t, zgrid):
    v = DIRS[direction]
    rows = []
    for z in zgrid:
        t0 = time.time()
        a, sd, bnds, F = cell_pcrown(t * v, E0, Z0, Z0, z, want_forms=True)
        tail = norm.cdf((T_BAR - a) / sd) if np.isfinite(a) else 1.0
        if tail > 0.05:
            rows.append(dict(z=z, a=a, sd=sd, tail=tail, hw=None, union=2 * N_RELU * norm.sf(z), bound=1.0)); continue
        hwb, un = hw_bound(F, z)
        rows.append(dict(z=z, a=a, sd=sd, tail=tail, hw=hwb, union=un, bound=hwb + tail, secs=time.time() - t0))
    return rows


if __name__ == "__main__":
    d = sys.argv[1]
    out = {}
    for t in map(float, sys.argv[2:]):
        rows = evaluate(d, t, (3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3))
        for r in rows:
            print("%s t=%.2f z=%.2f a=%8.2f tail=%.4f union=%.4f HW=%s bound=%.4f"
                  % (d, t, r["z"], r["a"], r["tail"], r["union"], "%.4f" % r["hw"] if r["hw"] is not None else "  -   ", r["bound"]), flush=True)
        out[t] = rows
    json.dump(out, open("hw_%s.json" % d, "w"), indent=1)
