"""05_recheck.py -- W1: the three cells that missed the 2-SE rule, rerun on FRESH seeds, 4,000 replicates.

Rule fixed before running (2026-09-23): the new replicates are judged ALONE (not pooled with the first
run) against the Theorem L formula. |z| <= 2 means the first run's miss was noise; |z| > 2 means a real
finite-n bias, to be reported as such. Every replicate is appended to disk (RECHECK_RAW.jsonl).

Cells: (null, 8000), (interaction, 8000), (threshold, 2000).
Run:  python -u 05_recheck.py        Env: LP_WORKERS (12)
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
import lp_core

HERE = os.path.dirname(os.path.abspath(__file__))
NW = int(os.environ.get("LP_WORKERS", "12"))
REPS = 4000
CELLS = [("null", 8000), ("interaction", 8000), ("threshold", 2000)]

if __name__ == "__main__":
    F = json.load(open(os.path.join(HERE, "FORMULA.json")))
    jobs = []
    for k, (name, n) in enumerate(CELLS):
        tau = 0.0 if name == "null" else F["dep"][name]["tau50"]["DG"]
        jobs += [(name, tau, n, 40_000_000 + k * 1_000_000 + r) for r in range(REPS)]
    RAW = os.path.join(HERE, "RECHECK_RAW.jsonl")
    done = set()
    if os.path.exists(RAW):
        for line in open(RAW):
            r = json.loads(line); done.add(r[6])
    todo = [j for j in jobs if j[3] not in done]
    print("recheck: %d to run (%d on disk), %d workers" % (len(todo), len(done), NW), flush=True)
    t0 = time.time()
    with Pool(NW, maxtasksperchild=200) as pool, open(RAW, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(lp_core.one_rep, todo, chunksize=10)):
            fh.write(json.dumps(list(r)) + "\n"); fh.flush()
            if (i + 1) % 1000 == 0:
                print("  %d / %d (%.1f min)" % (i + 1, len(todo), (time.time() - t0) / 60), flush=True)
    R = [json.loads(l) for l in open(RAW)]
    out = {}
    for name, n in CELLS:
        rr = [r for r in R if r[0] == name and r[1] == n]
        p = np.array([r[2] for r in rr]); sim = float((p <= 0.05 + 1e-12).mean())
        se = float(np.sqrt(sim * (1 - sim) / len(p)))
        form = F["size"]["DG"] if name == "null" else F["dep"][name]["power_at_tauDG"]["DG"]
        z = (sim - form) / se; clip = float(np.mean([r[5] > 0 for r in rr]))
        out["%s|%d" % (name, n)] = dict(sim=sim, formula=form, se=se, z=z, reps=len(p), clip=clip)
        print("%-12s n=%5d  sim %.4f  formula %.4f  SE %.4f  z %+.2f  (%d reps, clipped %.0f%%)  -> %s"
              % (name, n, sim, form, se, z, len(p), 100 * clip, "noise (within 2 SE)" if abs(z) <= 2 else "REAL finite-n bias"))
    json.dump(out, open(os.path.join(HERE, "RECHECK.json"), "w"), indent=2)
