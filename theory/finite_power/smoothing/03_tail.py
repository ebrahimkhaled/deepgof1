"""03_tail.py -- the deterministic tail of the certificate along one direction.  Resumable (tail_<d>.jsonl).

For t beyond the Monte Carlo net, a chain of intervals [t_a, t_b] is certified by CROWN over the ball
with centre ((t_a + t_b)/2) v and radius (t_b - t_a)/2 + R_TAIL, which contains every t v + e with
t in [t_a, t_b] and ||e|| <= R_TAIL.  When the CROWN lower bound exceeds c_max (the largest critical
value in play), P{S(tv + e) <= c_max} <= P{||e|| > R_TAIL} = EPS_TAIL for every t in the interval.
The chain runs from the last net point to t_global(c_max, kappa), beyond which Lemma 2 of the paper
with the certified constants (c_cert, L_cert) gives the same bound for all larger t.
Run: python -u 03_tail.py <direction>
"""
import os, sys, json, time
import numpy as np
from smooth_core import (NET, crown, DIRS, KAPPA, T_BAR, R_TAIL, EPS_TAIL, LEVELS, HERE,
                         tail_t_global, append_jsonl, read_jsonl)
sys.path.insert(0, HERE)
from importlib import import_module
design = import_module("02_mc_net").design

if __name__ == "__main__":
    d = sys.argv[1]; v = DIRS[d]
    ucb = json.load(open(os.path.join(HERE, "NULL_UCB.json")))["ucb"]
    c_max = max([T_BAR] + [ucb[str(w)] for w in LEVELS])
    t_start = max(t for t, M in design(d))
    t_glob = tail_t_global(c_max, KAPPA[d])
    path = os.path.join(HERE, "tail_%s.jsonl" % d)
    ivs = read_jsonl(path)
    t_a = ivs[-1]["t_b"] if ivs else t_start
    print("%s: tail from t=%.2f to t_global=%.1f (c_max %.4f, R_TAIL %.3f, kappa %.4f); %d intervals done"
          % (d, t_start, t_glob, c_max, R_TAIL, KAPPA[d], len(ivs)), flush=True)
    t0 = time.time(); ncalls = 0
    while t_a < t_glob:
        w = 0.05 * t_a; best = None
        # grow while it certifies (at most 3 doublings), shrink while it fails (down to 0.25)
        lb = crown(NET, (t_a + w) * v, w + R_TAIL)[0]; ncalls += 1
        if lb > c_max:
            best = (w, lb)
            for _ in range(3):
                w2 = 2 * w
                lb2 = crown(NET, (t_a + w2) * v, w2 + R_TAIL)[0]; ncalls += 1
                if lb2 > c_max: w, best = w2, (w2, lb2)
                else: break
        else:
            while w > 0.25:
                w = w / 2
                lb = crown(NET, (t_a + w) * v, w + R_TAIL)[0]; ncalls += 1
                if lb > c_max:
                    best = (w, lb); break
        if best is None:
            append_jsonl(path, dict(t_a=t_a, t_b=None, failed=True))
            print("  FAILED to certify any interval at t_a=%.3f -- the tail is open from here" % t_a)
            break
        w, lb = best
        t_b = min(t_a + 2 * w, t_glob)
        append_jsonl(path, dict(t_a=t_a, t_b=t_b, half_width=w, radius=w + R_TAIL, crown_lb=lb, c_max=c_max))
        if ncalls % 20 < 4:
            print("  [%.2f, %.2f]  lb %.1f  (%d CROWN calls, %.0f s)" % (t_a, t_b, lb, ncalls, time.time() - t0), flush=True)
        t_a = t_b
    ivs = read_jsonl(path)
    closed = bool(ivs) and ivs[-1].get("t_b") is not None and ivs[-1]["t_b"] >= t_glob - 1e-9
    json.dump(dict(direction=d, t_start=t_start, t_global=t_glob, c_max=c_max, r_tail=R_TAIL, eps_tail=EPS_TAIL,
                   n_intervals=len(ivs), closed=closed, crown_calls=ncalls, wall_s=time.time() - t0),
              open(os.path.join(HERE, "tail_%s_summary.json" % d), "w"), indent=1)
    print("%s tail %s: %d intervals, %d CROWN calls, %.0f s" % (d, "CLOSED" if closed else "OPEN", len(ivs), ncalls, time.time() - t0))
