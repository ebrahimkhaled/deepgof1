"""06_proj_offdiag.py -- the projection-test limit recomputed with the U-statistic (off-diagonal) discretisation.

Why: 01_are_irrelevant.py discretises the projection test's limit at a reference sample of M = 1000 draws and
keeps the diagonal of the M x M matrix. The deployed weight function at x = x' is about 2 pi p, against a
typical off-diagonal value of (4 pi / 3) p, so the diagonal adds O(p/M) to the drift and O(p^2/M) to
sum lambda^2 -- at p = 10 and M = 1000 comparable to the true values (checked in 04_proj_rate.py). In the
n -> infinity limit the diagonal only shifts the trace, which a bootstrap-calibrated test ignores. Here the
diagonal of the projected matrix is set to zero (ir_core.QuadForm(offdiag=True)); everything else, the
reference samples (same seeds as 01) and the B = 199 rule, is unchanged. The with-diagonal value is
recomputed alongside as a control (it must reproduce ARE_RAW.jsonl).

Tasks (argument list; each 'variant:k' or 'M2000:variant:k' or 'pred'):
  variant:k      tau50 of proj per seed and departure, off-diagonal and with diagonal   -> PROJ_OFFDIAG.jsonl
  M2000:v:k      the same at M = 2000 (seed 1 only), a discretisation check
  extra:v:k      four more reference samples (seeds 3-6), off-diagonal only, to reduce the Monte Carlo spread
  pred           k = 8 study design: limit power of proj (off-diagonal) at the niche cells' exact drift
Resumable; merge into ARE_IRRELEVANT.json with:  python 06_proj_offdiag.py assemble
"""
import os, sys, json, time
import numpy as np
from ir_core import *

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "PROJ_OFFDIAG.jsonl")
VARS = {"study": dict(cirr=C_IRR, seeds=[1, 2]), "pure": dict(cirr=0.0, seeds=[1])}
CELLS = [(f, C, n) for f in DEPS for C in (0.3, 0.5, 0.8) for n in (200, 500)]

def load():
    out = {}
    if os.path.exists(RAW):
        for line in open(RAW):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            out[r["key"]] = r
    return out

def log(rec):
    with open(RAW, "a") as fh:
        fh.write(json.dumps(rec) + "\n")

def assemble():
    D = load()
    A = json.load(open(os.path.join(HERE, "ARE_IRRELEVANT.json")))
    A["proj_offdiag_note"] = ("arm 'proj' = projection test, limit by U-statistic (off-diagonal) discretisation at M = 1000 "
                              "(06_proj_offdiag.py); 'proj_diag' = the first computation with the diagonal kept (finite-M "
                              "bias O(p/M) in the drift and O(p^2/M) in sum lambda^2), kept for the archive")
    for vname, V in A["variants"].items():
        for dep in DEPS:
            for ks, rec in V[dep].items():
                key = "%s:%s:%s" % (vname, ks, dep)
                if key not in D or "proj_diag" in rec["tau50"]:
                    continue
                r = D[key]; t = rec["tau50"]["DG"]
                rec["arms"]["proj_diag"] = rec["arms"]["proj"]
                rec["tau50"]["proj_diag"] = rec["tau50"]["proj"]
                rec["ratio_vs_DG"]["proj_diag"] = rec["ratio_vs_DG"]["proj"]
                ex = D.get("extra:" + key)
                rec["arms"]["proj"] = r["per_seed"] + (ex["per_seed"] if ex else [])
                rec["tau50"]["proj"] = float(np.mean([s["tau50"] for s in rec["arms"]["proj"]]))
                rec["ratio_vs_DG"]["proj"] = (rec["tau50"]["proj"] / t) ** 2
                rec["size"]["proj"] = float(np.mean([s["size"] for s in rec["arms"]["proj"]]))
    A["proj_M2000"] = {k: v for k, v in D.items() if k.startswith("M2000")}
    if "pred" in D:
        per = {}
        for kk, r in D.items():
            if kk == "pred" or kk.startswith("predx:"):
                for cell in r["cells"]:
                    per.setdefault(cell, []).extend(r.get("cells_per_seed", {}).get(cell, [r["cells"][cell]] * len(r.get("seeds", [1, 2]))))
        A["proj_pred_k8"] = {"cells": {c: float(np.mean(v)) for c, v in per.items()},
                             "cells_se": {c: float(np.std(v, ddof=1) / np.sqrt(len(v))) for c, v in per.items()},
                             "n_ref_samples": {c: len(v) for c, v in per.items()}, "size": D["pred"]["size"],
                             "note": "reference samples 1-2 stored as a mean only (both given that mean); 3-6 per seed"}
    json.dump(A, open(os.path.join(HERE, "ARE_IRRELEVANT.json"), "w"), indent=1)
    print("merged %d off-diagonal records into ARE_IRRELEVANT.json" % len(D))

def per_seed_record(q, dep, sd):
    pw = lambda t: q.power_bootstrap(t, dep)
    return dict(seed=sd, tau50=tau_for(pw), size=pw(0.0), signal=q.signal(dep), noise=q.noise(),
                snr=q.signal(dep) / q.noise(), max_share=q.max_share, n_neg=q.n_neg, n_keep=q.n_keep)

args = sys.argv[1:]
if args and args[0] == "assemble":
    assemble(); sys.exit(0)
done = load()
print("bench %.1e | %d records done | tasks %s" % (check_bench(), len(done), args)); sys.stdout.flush()
t0 = time.time()
for task in args:
    parts = task.split(":")
    if parts[0] in ("pred", "predx"):
        # 'pred' = reference samples 1, 2; 'predx:a,b' = further samples, stored per seed and averaged at assembly
        pkey = "pred" if parts[0] == "pred" else "predx:" + parts[1]
        if pkey in done:
            continue
        k = 8; out = {"cells": {}, "cells_per_seed": {}}
        qs = []
        pseeds = VARS["study"]["seeds"] if parts[0] == "pred" else [int(a) for a in parts[1].split(",")]
        for sd in pseeds:
            xp, Xp, pi0p, vp = ref_sample(k, 1000, 2000 * sd + k, cirr=C_IRR)
            q = QuadForm(proj_kernel(xp), Xp, vp, seed=sd + 100, offdiag=True)
            for (f, C, n) in CELLS:
                q.prep(hexact(f, xp, C_IRR, C, n), "%s_C%.1f_n%d" % (f, C, n))
            q.drop(); qs.append(q)
            print("[pred] seed %d kernel %.0f s" % (sd, time.time() - t0)); sys.stdout.flush()
        for (f, C, n) in CELLS:
            key = "%s_C%.1f_n%d" % (f, C, n)
            per = [q.power_bootstrap(1.0, key) for q in qs]
            out["cells_per_seed"][key] = per
            out["cells"][key] = float(np.mean(per))
        out["size"] = float(np.mean([q.power_bootstrap(0.0, CELLS[0][0] + "_C0.3_n200") for q in qs]))
        out["seeds"] = pseeds
        rec = dict(key=pkey, **out); log(rec); done[pkey] = rec
        print("[pred]", json.dumps(out)); sys.stdout.flush()
        continue
    big = parts[0] == "M2000"
    extra = parts[0] == "extra"
    vname, k = (parts[1], int(parts[2])) if (big or extra) else (parts[0], int(parts[1]))
    V = VARS[vname]
    pre = "M2000:" if big else ("extra:" if extra else "")
    keys = ["%s%s:%d:%s" % (pre, vname, k, dep) for dep in DEPS]
    if all(kk in done for kk in keys):
        continue
    M = 2000 if big else 1000
    seeds = [1] if big else ([3, 4, 5, 6] if extra else V["seeds"])
    res = {dep: {"offdiag": [], "diag": []} for dep in DEPS}
    for sd in seeds:
        xp, Xp, pi0p, vp = ref_sample(k, M, 2000 * sd + k, cirr=V["cirr"])
        R = proj_kernel(xp)
        for mode in (("offdiag",) if extra else ("offdiag", "diag")):
            q = QuadForm(R, Xp, vp, seed=sd + 100, offdiag=(mode == "offdiag"))
            for dep in DEPS:
                q.prep(vp * gfun(dep, xp[:, 0], xp[:, 1]), dep)
            q.drop()
            for dep in DEPS:
                res[dep][mode].append(per_seed_record(q, dep, sd))
            del q
        del R
        print("[%s] seed %d done %.0f s" % (task, sd, time.time() - t0)); sys.stdout.flush()
    for dep, kk in zip(DEPS, keys):
        rec = dict(key=kk, variant=vname, k=k, dep=dep, M=M, per_seed=res[dep]["offdiag"], per_seed_diag=res[dep]["diag"])
        log(rec); done[kk] = rec
        print("  %s tau50 offdiag %s | with diagonal %s" % (kk, [round(s["tau50"], 2) for s in res[dep]["offdiag"]],
                                                          [round(s["tau50"], 2) for s in res[dep]["diag"]])); sys.stdout.flush()
