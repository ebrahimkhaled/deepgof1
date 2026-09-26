"""04_truth.py -- the Monte Carlo truth for one direction: the signal t at which the power of the Gaussian
map model reaches the target, for the fixed-threshold test (S > T_BAR) and the deployed bootstrap test
(B = 199, alpha = .05; power = E[Bin(9; 199, W)] with W = 1 - F_0(S) taken from the 2,000,000 null draws),
at targets .95 and .50.  Bisection with M = 200,000 fresh draws per evaluation (seeds independent of
the net), then a standard error by the delta method from the local slope.  Resumable: truth_<d>.jsonl.
Run: python -u 04_truth.py <direction>
"""
import os, sys, json, time
import numpy as np
from scipy.stats import binom
from smooth_core import forward, DIRS, T_BAR, T_HAT, T_EMP_PAPER, HERE, B_BOOT, K_REJ, append_jsonl, read_jsonl

M = 200_000
null_sorted = np.load(os.path.join(HERE, "null_sorted.npy")); N0 = len(null_sorted)

def powers(t, v, rng):
    s = forward(t * v + rng.normal(size=(M, 36)))
    W = 1.0 - np.searchsorted(null_sorted, s, side="right") / N0        # P(S* > s)
    return dict(fixed=float((s > T_BAR).mean()), hat=float((s > T_HAT).mean()),
                b199=float(binom.cdf(K_REJ - 1, B_BOOT, W).mean()))

def solve(version, target, v, seed, log):
    rng = np.random.default_rng([seed, 4])
    lo, hi = 0.0, 2.0
    while powers(hi, v, rng)[version] < target:
        hi *= 2
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if powers(mid, v, rng)[version] >= target: hi = mid
        else: lo = mid
    t_hat = 0.5 * (lo + hi)
    # standard error: SE(t) = SE(q) / slope, slope from +-0.25 around t_hat with fresh draws
    h = 0.25
    q0 = powers(t_hat, v, rng)[version]; qp = powers(t_hat + h, v, rng)[version]; qm = powers(t_hat - h, v, rng)[version]
    slope = (qp - qm) / (2 * h)
    se_q = np.sqrt(max(q0 * (1 - q0), 1e-12) / M)
    se_t = se_q / slope if slope > 0 else float("nan")
    log("  %s target %.2f: t = %.4f  (power there %.4f, slope %.4f/unit t, SE(t) = %.4f)" % (version, target, t_hat, q0, slope, se_t))
    return dict(version=version, target=target, t=t_hat, power_at_t=q0, slope=slope, se_t=se_t, M=M)

if __name__ == "__main__":
    d = sys.argv[1]; v = DIRS[d]
    path = os.path.join(HERE, "truth_%s.jsonl" % d)
    done = {(r["version"], r["target"]) for r in read_jsonl(path)}
    t0 = time.time()
    print("%s: paper's t_emp (threshold T_HAT, power .95) = %.4f" % (d, T_EMP_PAPER[d]), flush=True)
    for i, (version, target) in enumerate([("fixed", .95), ("fixed", .50), ("b199", .95), ("b199", .50), ("hat", .95)]):
        if (version, target) in done:
            continue
        rec = solve(version, target, v, seed=sum(map(ord, d)) * 100 + i, log=lambda s: print(s, flush=True))
        rec["direction"] = d
        append_jsonl(path, rec)
    print("%s truth done in %.0f s" % (d, time.time() - t0))
