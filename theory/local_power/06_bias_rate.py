"""06_bias_rate.py -- W1: how fast does the finite-n bias of the threshold departure vanish?

Declared before running (2026-09-23): threshold departure at its DeepGOF-1 tau_50, n in
{1000, 4000, 8000, 16000}, 4,000 fresh-seed replicates each (clipping recorded), plus the existing n = 2000
(RECHECK_RAW.jsonl, 4,000 reps) and n = 32000 (CONV_RAW.jsonl, 1,000 reps). Fit
bias(n) = c * n^(-gamma) by weighted nonlinear least squares (weights 1/SE^2) over the cells with no
clipping, and report gamma with its SE and a bootstrap interval (resampling replicates within cells).
Reading: gamma ~ 0.5 is the usual O(n^-1/2) finite-sample term; gamma ~ 0 would mean no convergence.
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import curve_fit
import lp_core

HERE = os.path.dirname(os.path.abspath(__file__))
NW = int(os.environ.get("LP_WORKERS", "12"))
NEW = (1000, 4000, 8000, 16000); REPS = 4000

def rows(path, name, n):
    return [json.loads(l) for l in open(path) if json.loads(l)[0] == name and json.loads(l)[1] == n]

if __name__ == "__main__":
    F = json.load(open(os.path.join(HERE, "FORMULA.json")))
    tau = F["dep"]["threshold"]["tau50"]["DG"]; form = F["dep"]["threshold"]["power_at_tauDG"]["DG"]
    RAW = os.path.join(HERE, "BIAS_RAW.jsonl")
    done = set()
    if os.path.exists(RAW):
        for line in open(RAW): done.add(json.loads(line)[6])
    jobs = [("threshold", tau, n, 60_000_000 + n * 10 + r) for n in NEW for r in range(REPS)]
    todo = [j for j in jobs if j[3] not in done]
    print("bias-rate: %d to run (%d on disk)" % (len(todo), len(done)), flush=True)
    t0 = time.time()
    with Pool(NW, maxtasksperchild=200) as pool, open(RAW, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(lp_core.one_rep, todo, chunksize=10)):
            fh.write(json.dumps(list(r)) + "\n"); fh.flush()
            if (i + 1) % 2000 == 0:
                print("  %d / %d (%.1f min)" % (i + 1, len(todo), (time.time() - t0) / 60), flush=True)
    cells = {}
    for n in NEW: cells[n] = rows(RAW, "threshold", n)
    cells[2000] = rows(os.path.join(HERE, "RECHECK_RAW.jsonl"), "threshold", 2000)
    cells[32000] = rows(os.path.join(HERE, "CONV_RAW.jsonl"), "threshold", 32000)
    ns, bias, se, rej_by_n = [], [], [], {}
    print("\n%6s %6s %8s %8s %7s %7s %8s" % ("n", "reps", "sim", "bias", "SE", "z", "clipped"))
    for n in sorted(cells):
        rr = cells[n]; rej = np.array([r[2] <= 0.05 + 1e-12 for r in rr], float); clip = np.mean([r[5] > 0 for r in rr])
        s = rej.mean(); e = np.sqrt(s * (1 - s) / len(rej)); b = s - form
        print("%6d %6d %8.4f %+8.4f %7.4f %+7.2f %7.0f%%" % (n, len(rej), s, b, e, b / e, 100 * clip))
        if clip == 0:
            ns.append(n); bias.append(b); se.append(e); rej_by_n[n] = rej
    ns, bias, se = np.array(ns, float), np.array(bias), np.array(se)
    f = lambda n, c, g: c * n ** (-g)
    (c, g), cov = curve_fit(f, ns, bias, p0=(1.0, 0.5), sigma=se, absolute_sigma=True, maxfev=20000)
    rng = np.random.default_rng(7); gb = []
    for _ in range(2000):
        bb = np.array([rng.choice(rej_by_n[n], len(rej_by_n[n])).mean() - form for n in ns])
        try:
            gb.append(curve_fit(f, ns, bb, p0=(c, g), sigma=se, absolute_sigma=True, maxfev=20000)[0][1])
        except Exception:
            pass
    lo, hi = np.quantile(gb, [0.025, 0.975])
    print("\nfit bias(n) = c n^-gamma:  c = %.3f,  gamma = %.3f (SE %.3f; bootstrap 95%% [%.3f, %.3f]; %d fits)"
          % (c, g, np.sqrt(cov[1, 1]), lo, hi, len(gb)))
    json.dump(dict(formula=form, n=ns.tolist(), bias=bias.tolist(), se=se.tolist(), c=c, gamma=g,
                   gamma_se=float(np.sqrt(cov[1, 1])), gamma_boot=[float(lo), float(hi)]),
              open(os.path.join(HERE, "BIAS_RATE.json"), "w"), indent=2)
