"""02_simulate.py -- W1, part 2: does Theorem L predict the finite-sample power?

Simulates the three tests in full (fit, map, score, 199 bootstrap refits for DeepGOF-1) at the DeepGOF-1
tau_50 of each departure and under the null, at n = 500 and n = 2000, and compares each rejection rate
with the limit formula in FORMULA.json. The exit test of W1: every DeepGOF-1 cell within 2 Monte Carlo
standard errors of its formula. The rivals' agreement is reported too, because the efficiency table is
only as good as their formulas.

Run:  python -u 02_simulate.py        Env: LP_REPS (1000), LP_WORKERS (16)
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
import lp_core

HERE = os.path.dirname(os.path.abspath(__file__))
REPS = int(os.environ.get("LP_REPS", "1000"))
NW = int(os.environ.get("LP_WORKERS", "16"))
TESTS = ["DG", "HL", "Stukel"]

if __name__ == "__main__":
    FORM = json.load(open(os.path.join(HERE, "FORMULA.json")))
    jobs = []
    for n in (500, 2000):
        jobs += [("null", 0.0, n, 7_000_000 + n * 10 + r) for r in range(REPS)]
        for k, name in enumerate(lp_core.DEPS):
            t = FORM["dep"][name]["tau50"]["DG"]
            jobs += [(name, t, n, 8_000_000 + k * 100_000 + n * 10 + r) for r in range(REPS)]
    print("simulating %d datasets on %d workers" % (len(jobs), NW), flush=True)
    t0 = time.time()
    res = []
    with Pool(NW) as pool:
        for i, r in enumerate(pool.imap_unordered(lp_core.one_rep, jobs, chunksize=10)):
            res.append(r)
            if (i + 1) % 1000 == 0:
                print("  %d / %d  (%.1f min)" % (i + 1, len(jobs), (time.time() - t0) / 60), flush=True)
    json.dump(res, open(os.path.join(HERE, "SIM_RAW.json"), "w"))       # every p-value, per rep
    out = {}
    print("\n%-12s %5s  %-9s %8s %8s %7s %6s" % ("departure", "n", "test", "formula", "sim", "SE", "z"))
    for name in ["null"] + lp_core.DEPS:
        form = FORM["size"] if name == "null" else FORM["dep"][name]["power_at_tauDG"]
        for n in (500, 2000):
            P = np.array([r[2:5] for r in res if r[0] == name and r[1] == n])
            for j, t in enumerate(TESTS):
                sim = float((P[:, j] <= lp_core.ALPHA + 1e-12).mean())
                se = float(np.sqrt(max(sim * (1 - sim), 1e-4) / len(P)))
                z = (sim - form[t]) / se
                out["%s|%d|%s" % (name, n, t)] = dict(formula=form[t], sim=sim, se=se, z=z, reps=len(P))
                flag = "" if abs(z) <= 2 else "   <-- outside 2 SE"
                print("%-12s %5d  %-9s %8.3f %8.3f %7.3f %6.2f%s" % (name, n, t, form[t], sim, se, z, flag))
    dg = [v for k, v in out.items() if k.endswith("|DG")]
    print("\nW1 exit test: DeepGOF-1 cells within 2 SE: %d of %d" % (sum(abs(v["z"]) <= 2 for v in dg), len(dg)))
    json.dump(out, open(os.path.join(HERE, "SIMULATION.json"), "w"), indent=2)
    print("wrote theory/local_power/SIMULATION.json")
