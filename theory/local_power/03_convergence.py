"""03_convergence.py -- W1 exit test where the theorem applies: does finite-n power converge to the
Theorem L formula as n grows, once the local alternative is realisable without clipping?

02_simulate.py showed that at n = 500 the DeepGOF-1 tau_50 of several departures pushes
pi_0 + h/sqrt(n) outside (0, 1) (the cubic needs n >~ (tau max|g|)^2 ~ 6000), so the simulated law was
not the local alternative at all. Here n in {2000, 8000, 32000}; every replicate records how many
records had to be clipped. DeepGOF-1, HL and Stukel are simulated.

Exit (amended in GOAL_THEORY.md, with the reason): in every departure, the DeepGOF-1 cells whose
clipping share is zero are within 2 SE of the formula, and |sim - formula| does not grow with n.

Run:  python -u 03_convergence.py      Env: LP_REPS (1000), LP_WORKERS (16)
"""
import os, json, time
import numpy as np
from multiprocessing import Pool
import lp_core

HERE = os.path.dirname(os.path.abspath(__file__))
REPS = int(os.environ.get("LP_REPS", "1000"))
NW = int(os.environ.get("LP_WORKERS", "16"))
NS = (2000, 8000, 32000)

if __name__ == "__main__":
    FORM = json.load(open(os.path.join(HERE, "FORMULA.json")))
    jobs = []
    for n in NS:
        jobs += [("null", 0.0, n, 5_000_000 + n * 10 + r) for r in range(REPS)]
        for k, name in enumerate(lp_core.DEPS):
            t = FORM["dep"][name]["tau50"]["DG"]
            jobs += [(name, t, n, 6_000_000 + k * 1_000_000 + n * 10 + r) for r in range(REPS)]
    # every finished replicate is appended to disk at once, so a crash loses nothing (the first run died
    # of out-of-memory at 12,000 of 18,000 with its results held in RAM); a rerun skips what is done
    RAW = os.path.join(HERE, "CONV_RAW.jsonl")
    done = set()
    if os.path.exists(RAW):
        for line in open(RAW):
            r = json.loads(line); done.add((r[0], r[1], r[6]))
    todo = [j for j in jobs if (j[0], j[2], j[3]) not in done]
    print("simulating %d datasets (%d already on disk) on %d workers" % (len(todo), len(done), NW), flush=True)
    t0 = time.time()
    with Pool(NW) as pool, open(RAW, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(lp_core.one_rep, todo, chunksize=5)):
            fh.write(json.dumps(list(r)) + "\n"); fh.flush()
            if (i + 1) % 1000 == 0:
                print("  %d / %d  (%.1f min)" % (i + 1, len(todo), (time.time() - t0) / 60), flush=True)
    res = [json.loads(line) for line in open(RAW)]
    out = {}
    print("\n%-12s %6s  %-7s %8s %8s %7s %6s %9s" % ("departure", "n", "test", "formula", "sim", "SE", "z", "clipped"))
    for name in ["null"] + lp_core.DEPS:
        form = FORM["size"] if name == "null" else FORM["dep"][name]["power_at_tauDG"]
        for n in NS:
            rr = [r for r in res if r[0] == name and r[1] == n]
            P = np.array([r[2:5] for r in rr]); clip = np.mean([r[5] > 0 for r in rr])
            for j, t in enumerate(["DG", "HL", "Stukel"]):
                sim = float((P[:, j] <= lp_core.ALPHA + 1e-12).mean())
                se = float(np.sqrt(max(sim * (1 - sim), 1e-4) / len(P)))
                z = (sim - form[t]) / se
                out["%s|%d|%s" % (name, n, t)] = dict(formula=form[t], sim=sim, se=se, z=z, clip_share=float(clip))
                print("%-12s %6d  %-7s %8.3f %8.3f %7.3f %6.2f %8.1f%%%s" % (name, n, t, form[t], sim, se, z,
                      100 * clip, "" if abs(z) <= 2 else "   <-- outside 2 SE"))
    ok = [v for k, v in out.items() if k.endswith("|DG") and v["clip_share"] == 0]
    print("\nW1 exit (DeepGOF-1 cells with no clipping within 2 SE): %d of %d" % (sum(abs(v["z"]) <= 2 for v in ok), len(ok)))
    json.dump(out, open(os.path.join(HERE, "CONVERGENCE.json"), "w"), indent=2)
