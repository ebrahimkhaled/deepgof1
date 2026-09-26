"""05_summary.py -- summary tables from ARE_IRRELEVANT.json, PRED_VS_SIM.json and PROJ_RATE.json
(reads only). Prints markdown.   Run:  python 05_summary.py > SUMMARY.log
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
A = json.load(open(os.path.join(HERE, "ARE_IRRELEVANT.json")))
KS = A["k"]
ARMS = ["proj", "proj_diag", "Stukel", "HL"]

print("## curve: ratio n*_arm / n*_DG (range over reference samples for the quadratic-form arms)\n")
for vname in ("study", "pure"):
    V = A["variants"][vname]
    print("### variant %s\n" % vname)
    print("| dep | arm | " + " | ".join("k=%d" % k for k in KS) + " |")
    print("|---|---|" + "---|" * len(KS))
    for dep in ("inter", "quad"):
        for arm in ARMS + ["tau50_DG"]:
            cells = []
            for k in KS:
                r = V[dep].get(str(k))
                if r is None:
                    cells.append("–"); continue
                if arm == "tau50_DG":
                    cells.append("%.2f" % r["tau50"]["DG"]); continue
                val = r["ratio_vs_DG"][arm]
                if arm in r.get("arms", {}) and len(r["arms"][arm]) > 1:
                    ts = np.array([s["tau50"] for s in r["arms"][arm]])
                    rr = (ts / r["tau50"]["DG"]) ** 2
                    se = 2 * val * ts.std(ddof=1) / np.sqrt(len(ts)) / ts.mean()
                    cells.append("%.2f ±%.2f [%.2f-%.2f; %d]" % (val, se, rr.min(), rr.max(), len(ts)))
                else:
                    cells.append("%.2f" % val)
            print("| %s | %s | %s |" % (dep, arm, " | ".join(cells)))
    print()

print("## projection test: theta* = tau50^2 sigma/nu, p e_k (per reference sample; study variant)\n")
V = A["variants"]["study"]
for dep in ("inter", "quad"):
    for k in KS:
        r = V[dep].get(str(k))
        if r is None:
            continue
        s = r["arms"]["proj"]
        print("%s k=%d: theta* %s | p e_k %s | df %s | max share %s" % (
            dep, k, ", ".join("%.2f" % (x["tau50"] ** 2 * x["snr"]) for x in s),
            ", ".join("%.4f" % ((2 + k) * x["snr"]) for x in s), ", ".join("%.0f" % x.get("df", float("nan")) for x in s),
            ", ".join("%.2f" % x["max_share"] for x in s)))
print()

P = json.load(open(os.path.join(HERE, "PRED_VS_SIM.json")))
print("## predicted vs simulated (paper arms; k = 8, study design)\n")
print("| cell | tau | DG pred/sim (z) | proj limit (U) / limit with diagonal / n-sur / sim (z U-limit, z nsur) | Stukel pred/sim (z) | HL pred/sim (z) | sel12 |")
print("|---|---|---|---|---|---|---|")
for key in sorted(k for k in P["sim"] if not k.startswith("null")):
    s, pr = P["sim"][key], P["pred"][key]
    n = s["n"]
    def z(sim, pred, se):
        return (sim - pred) / se if se > 0 else float("nan")
    row = ["%s" % key, "%.2f" % pr["tau"]]
    row.append("%.3f / %.3f±%.3f (%+.1f)" % (pr["DG_exact"], s["raw"]["DG"], s["se"]["DG"], z(s["raw"]["DG"], pr["DG_exact"], s["se"]["DG"])))
    nsur = pr.get("proj_nsur_exact", float("nan"))
    lim = A.get("proj_pred_k8", {}).get("cells", {}).get(key, float("nan"))
    row.append("%.3f / %.3f / %.3f / %.3f±%.3f (%+.1f, %+.1f)" % (lim, pr["proj_exact"], nsur, s["raw"]["proj"], s["se"]["proj"],
               z(s["raw"]["proj"], lim, s["se"]["proj"]), z(s["raw"]["proj"], nsur, s["se"]["proj"])))
    for a in ("Stukel", "HL"):
        row.append("%.3f / %.3f±%.3f (%+.1f)" % (pr[a + "_exact"], s["raw"][a], s["se"][a], z(s["raw"][a], pr[a + "_exact"], s["se"][a])))
    row.append("%.3f" % s["sel12"])
    print("| " + " | ".join(row) + " |")
print()
print("null sizes (sim ± SE):")
for key in sorted(k for k in P["sim"] if k.startswith("null")):
    s = P["sim"][key]
    print("  %s: " % key + ", ".join("%s %.4f±%.4f" % (a, s["raw"][a], s["se"][a]) for a in ("DG", "proj", "Stukel", "HL")))
print("limits:", {a: round(P["pred"]["size"][a], 4) for a in ("DG", "proj", "Stukel", "HL")})

f = os.path.join(HERE, "PROJ_RATE.json")
if os.path.exists(f):
    R = json.load(open(f))
    c = R["constants"]
    print("\n## projection rate check (pure design)\nQ_P", c["Q_P"], "C_P %.4f Ev %.4f" % (c["C_P"], c["Ev"]))
    print("| p | inter p*sigma (SE), 2 samples | quad p*sigma (SE) |")
    print("|---|---|---|")
    ps = sorted({r["p"] for r in R["signal"]})
    for p in ps:
        rs = [r for r in R["signal"] if r["p"] == p]
        print("| %d | %s | %s |" % (p, "; ".join("%.4f (%.4f)" % (p * r["sigma_inter"], p * r["se_inter"]) for r in rs),
                                    "; ".join("%.4f (%.4f)" % (p * r["sigma_quad"], p * r["se_quad"]) for r in rs)))
    for p in (50, 200):
        rs = [r for r in R.get("signal_rep", []) if r["p"] == p]
        for d in ("inter", "quad"):
            if rs:
                v = np.array([p * r["sigma_" + d] for r in rs])
                print("signal_rep p=%d %s: %.5f (SE %.5f, %d samples) vs Q_P %.5f" % (p, d, v.mean(), v.std(ddof=1) / np.sqrt(len(v)), len(v), c["Q_P"][d]))
    print("\n| p | seed | nu (U) | nu with diagonal | lower (finite k) | lower asym | upper asym | max share | q95 std |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in R["noise"]:
        print("| %d | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f | %.2f |" % (r["p"], r["seed"], r["nu"], r["nu_with_diag"],
              np.sqrt(r["lower_finite"]), np.sqrt(r["lower_asym"]), np.sqrt(r["upper_asym"]), r["max_share"], r["q95_std"]))
