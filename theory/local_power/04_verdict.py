"""04_verdict.py -- W1 exit verdict from CONV_RAW.jsonl (n = 2000, 8000, 32000) and FORMULA.json.

Reports, for DeepGOF-1:
  * every cell: simulated power vs the Theorem L formula, z, clipping share;
  * the per-cell rule of the amended exit (cells with zero clipping within 2 SE);
  * a JOINT test: sum of z^2 over the zero-clipping cells against chi-square(k). With k ~ 15 cells, about
    one cell outside 2 SE is expected by chance alone, so the joint test is the fair summary; both are
    reported;
  * the trend of |sim - formula| with n per departure.
And re-checks the rivals' formulas against the n = 500 / 2000
simulations saved by 02_simulate.py (SIM_RAW.json), zero-clipping cells only.
"""
import os, json
import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
F = json.load(open(os.path.join(HERE, "FORMULA.json")))
R = [json.loads(l) for l in open(os.path.join(HERE, "CONV_RAW.jsonl"))]
NAMES = ["null", "quadratic", "interaction", "cubic", "bump", "threshold"]
A = 0.05

def form(name, t):
    return F["size"][t] if name == "null" else F["dep"][name]["power_at_tauDG"][t]

print("DeepGOF-1, Theorem L vs simulation")
print("%-12s %6s %8s %8s %7s %7s %8s" % ("departure", "n", "formula", "sim", "SE", "z", "clipped"))
zs, cells = [], []
for name in NAMES:
    for n in (2000, 8000, 32000):
        rr = [r for r in R if r[0] == name and r[1] == n]
        if not rr: continue
        p = np.array([r[2] for r in rr]); sim = float((p <= A + 1e-12).mean())
        se = float(np.sqrt(max(sim * (1 - sim), 1e-4) / len(p)))
        clip = float(np.mean([r[5] > 0 for r in rr])); z = (sim - form(name, "DG")) / se
        cells.append(dict(name=name, n=n, sim=sim, formula=form(name, "DG"), se=se, z=z, clip=clip, reps=len(p)))
        if clip == 0: zs.append(z)
        print("%-12s %6d %8.3f %8.3f %7.3f %7.2f %7.0f%%%s" % (name, n, form(name, "DG"), sim, se, z, 100 * clip,
              "" if abs(z) <= 2 or clip > 0 else "   <- outside 2 SE"))
zs = np.array(zs); X2 = float((zs ** 2).sum()); k = len(zs)
print("\nper-cell rule: %d of %d zero-clipping cells within 2 SE (expected by chance alone: %.1f outside)"
      % (int((np.abs(zs) <= 2).sum()), k, 0.0455 * k))
print("JOINT test: sum z^2 = %.2f on %d df, p = %.3f  (a small p would reject the formula)"
      % (X2, k, chi2.sf(X2, k)))
print("\n|sim - formula| by n (zero-clipping cells):")
for name in NAMES:
    row = "  %-12s" % name
    for n in (2000, 8000, 32000):
        c = [c for c in cells if c["name"] == name and c["n"] == n]
        row += "  n=%d: %s" % (n, ("%.3f" % abs(c[0]["sim"] - c[0]["formula"])) if c and c[0]["clip"] == 0 else "  -  ")
    print(row)

out = dict(cells=cells, joint_X2=X2, joint_df=k, joint_p=float(chi2.sf(X2, k)))
if os.path.exists(os.path.join(HERE, "SIM_RAW.json")):
    S = json.load(open(os.path.join(HERE, "SIM_RAW.json")))
    print("\nrivals at n = 500 / 2000 vs the corrected formulas (02_simulate.py data; clipping not recorded there,")
    print("so only departures that cannot clip at these n are shown: null, bump, threshold at n = 2000)")
    for name, n in [("null", 500), ("null", 2000), ("bump", 2000), ("threshold", 2000)]:
        rr = [r for r in S if r[0] == name and r[1] == n]
        P = np.array([r[2:5] for r in rr])
        row = "  %-10s n=%4d" % (name, n)
        for j, t in enumerate(["DG", "HL", "Stukel"]):
            sim = (P[:, j] <= A + 1e-12).mean(); se = np.sqrt(max(sim * (1 - sim), 1e-4) / len(P))
            row += "   %s %.3f/%.3f z=%+.1f" % (t, sim, form(name, t), (sim - form(name, t)) / se)
        print(row)
json.dump(out, open(os.path.join(HERE, "CONVERGENCE.json"), "w"), indent=2)
