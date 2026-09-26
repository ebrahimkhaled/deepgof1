"""02_analyse.py -- the small-n study, analysed exactly as PREDECLARATION.md says.

Q1  size of each base test under one corrupted record (C = 0); misbehaves if > .075.
Q2  size of the LOO-unanimity rule under each contamination; must not exceed the base test's CLEAN size
    by more than 2 SE (Theorem 1 of THEORY_LOO.md bounds it by the clean size at n - 1).
Power cost of the rule: base vs LOO power at C in {0.5, 1.0, 1.5}, clean and with a decimal slip.
POST HOC (declared as such, not part of the pre-declaration): the calibrated variant of Proposition 8,
alpha' chosen so the LOO rule's clean size is .05, with its contaminated size and its power.

Run:  python 02_analyse.py
"""
import os, glob, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(HERE, "rows", "*.csv")))],
              ignore_index=True)
TESTS = ["HL", "Stukel", "DG"]
A = 0.05
print("rows: %d  (expected 2*4*1000 + 2*3*2*300 = %d)" % (len(D), 2 * 4 * 1000 + 2 * 3 * 2 * 300))

def rate(x):
    x = np.asarray(x, bool); r = x.mean(); return r, np.sqrt(max(r * (1 - r), 1e-4) / len(x))

out = {"Q1": {}, "Q2": {}, "power": {}, "posthoc_alpha_prime": {}}
print("\n=== Q1: size of the base test with one corrupted record (C = 0); misbehave if > .075 ===")
print("%-4s %-8s" % ("n", "type") + "".join("%16s" % t for t in TESTS))
null = D[D.C == 0]
for n in (30, 50):
    for ty in ("none", "decimal", "sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        row = "%-4d %-8s" % (n, ty)
        for t in TESTS:
            r, se = rate(s["p_" + t] <= A)
            out["Q1"]["%d|%s|%s" % (n, ty, t)] = dict(size=r, se=se, reps=len(s), misbehaves=bool(r > .075))
            row += "   %.3f (%.3f)%s" % (r, se, "*" if r > .075 else " ")
        print(row)
print("  (* = size above .075; the SE is in brackets)")

print("\n=== Q2: size of the LOO rule; must stay within 2 SE of the base test's CLEAN size ===")
print("%-4s %-8s" % ("n", "type") + "".join("%22s" % t for t in TESTS))
for n in (30, 50):
    clean = null[(null.n == n) & (null.type == "none")]
    for ty in ("none", "decimal", "sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        row = "%-4d %-8s" % (n, ty)
        for t in TESTS:
            c0, c0se = rate(clean["p_" + t] <= A)
            r, se = rate(s["loo_" + t] <= A)
            lim = c0 + 2 * np.sqrt(c0se ** 2 + se ** 2)
            ok = r <= lim
            out["Q2"]["%d|%s|%s" % (n, ty, t)] = dict(loo_size=r, se=se, clean_base=c0, limit=lim, holds=bool(ok))
            row += "   %.3f vs %.3f %s" % (r, c0, "ok " if ok else "FAIL")
        print(row)

print("\n=== power: base test vs LOO rule (300 reps; SE <= .029) ===")
alt = D[D.C > 0]
print("%-4s %-4s %-8s" % ("n", "C", "type") + "".join("%16s" % t for t in TESTS) + "     (base / LOO)")
for n in (30, 50):
    for C in (0.5, 1.0, 1.5):
        for ty in ("none", "decimal"):
            s = alt[(alt.n == n) & (alt.C == C) & (alt.type == ty)]
            row = "%-4d %-4.1f %-8s" % (n, C, ty)
            for t in TESTS:
                b, _ = rate(s["p_" + t] <= A); l, _ = rate(s["loo_" + t] <= A)
                out["power"]["%d|%.1f|%s|%s" % (n, C, ty, t)] = dict(base=b, loo=l, reps=len(s))
                row += "   %.3f / %.3f" % (b, l)
            print(row)

print("\n=== POST HOC: calibrated alpha' (clean LOO size = .05), Proposition 8 ===")
for n in (30, 50):
    clean = null[(null.n == n) & (null.type == "none")]
    for t in TESTS:
        ap = float(np.quantile(clean["loo_" + t], A))           # P(loo <= alpha') = .05 on clean nulls
        sizes = {ty: rate(null[(null.n == n) & (null.type == ty)]["loo_" + t] <= ap)[0]
                 for ty in ("decimal", "sign", "outcome")}
        base_at = {ty: rate(null[(null.n == n) & (null.type == ty)]["p_" + t] <= A)[0]
                   for ty in ("decimal", "sign", "outcome")}
        pw = {C: (rate(alt[(alt.n == n) & (alt.C == C) & (alt.type == "none")]["p_" + t] <= A)[0],
                  rate(alt[(alt.n == n) & (alt.C == C) & (alt.type == "none")]["loo_" + t] <= ap)[0])
              for C in (0.5, 1.0, 1.5)}
        out["posthoc_alpha_prime"]["%d|%s" % (n, t)] = dict(alpha_prime=ap, contaminated_size=sizes,
                                                           base_contaminated_size=base_at, power_base_vs_cal=pw)
        print("n=%d %-6s alpha'=%.3f  size decimal/sign/outcome %.3f/%.3f/%.3f (base %.3f/%.3f/%.3f)  "
              "power C=.5/1/1.5 base %.2f/%.2f/%.2f  cal %.2f/%.2f/%.2f"
              % (n, t, ap, sizes["decimal"], sizes["sign"], sizes["outcome"],
                 base_at["decimal"], base_at["sign"], base_at["outcome"],
                 pw[0.5][0], pw[1.0][0], pw[1.5][0], pw[0.5][1], pw[1.0][1], pw[1.5][1]))

print("\n=== diagnostics: how often does deleting the CORRUPTED record rescue a base rejection? ===")
for n in (30, 50):
    for ty in ("decimal", "sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        row = "n=%d %-8s" % (n, ty)
        for t in TESTS:
            rej = s[s["p_" + t] <= A]
            hit = (rej["amax_" + t] == rej["corrupted"]).mean() if len(rej) else float("nan")
            row += "  %s: %3d base rejections, most-rescuing deletion = corrupted record %.0f%%" % (t, len(rej), 100 * hit)
        print(row)
json.dump(out, open(os.path.join(HERE, "RESULTS.json"), "w"), indent=2)
print("\nwrote theory/small_n/RESULTS.json")
