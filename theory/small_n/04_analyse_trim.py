"""04_analyse_trim.py -- the confirmation study of covariate-only trimming, exactly as
PREDECLARATION_TRIM.md declares: (1) protection under decimal corruption, (2) no harm when clean,
(3) mechanism, (4) sign / outcome reported without claim."""
import os, glob, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(HERE, "rows_trim", "*.csv")))], ignore_index=True)
T = ["HL", "Stukel", "DG"]; A = 0.05
print("rows %d (expected 11600)" % len(D))
def rate(x):
    x = np.asarray(x, bool); r = x.mean(); return r, np.sqrt(max(r * (1 - r), 1e-4) / len(x))
out = {"check1": {}, "check2": {}, "check3": {}, "check4": {}}
null = D[D.C == 0]

print("\n(1) PROTECTION: trimmed size under decimal corruption vs the base test's CLEAN size (2 SE)")
for n in (30, 50):
    cl = null[(null.n == n) & (null.type == "none")]; de = null[(null.n == n) & (null.type == "decimal")]
    row = "n=%d" % n
    for t in T:
        c0, s0 = rate(cl["p_" + t] <= A); b, _ = rate(de["p_" + t] <= A); r, s = rate(de["t_" + t] <= A)
        ok = r <= c0 + 2 * np.sqrt(s0 ** 2 + s ** 2)
        out["check1"]["%d|%s" % (n, t)] = dict(clean_base=c0, decimal_base=b, decimal_trim=r, se=s, holds=bool(ok))
        row += "   %s: base %.3f -> trimmed %.3f (clean %.3f) %s" % (t, b, r, c0, "ok" if ok else "FAIL")
    print(row)

print("\n(2) NO HARM WHEN CLEAN: trimmed vs base rejection rate on clean data (null and C > 0)")
for n in (30, 50):
    for C in (0.0, 0.5, 1.0, 1.5):
        s = D[(D.n == n) & (D.C == C) & (D.type == "none")]
        row = "n=%d C=%.1f" % (n, C)
        for t in T:
            b, sb = rate(s["p_" + t] <= A); r, sr = rate(s["t_" + t] <= A)
            diff = r - b; se = np.sqrt(max(np.mean((s["t_" + t] <= A) != (s["p_" + t] <= A)), 1e-4) / len(s))
            ok = abs(diff) <= 2 * max(se, 1e-3)
            out["check2"]["%d|%.1f|%s" % (n, C, t)] = dict(base=b, trim=r, diff=diff, paired_se=se, holds=bool(ok))
            row += "   %s %.3f/%.3f%s" % (t, b, r, "" if ok else "!")
        print(row)
print("   (base/trimmed; '!' = differs by more than 2 paired SE)")

print("\n(3) MECHANISM under decimal corruption")
for n in (30, 50):
    de = D[(D.n == n) & (D.type == "decimal")]
    hit = (de.trimmed == de.corrupted).mean(); none = (de.trimmed == 0).mean()
    missed = de[de.trimmed != de.corrupted]
    m0 = null[(null.n == n) & (null.type == "decimal")]; mm = m0[m0.trimmed != m0.corrupted]
    out["check3"][str(n)] = dict(corrupted_deleted=float(hit), nothing_deleted=float(none),
                                 size_when_missed={t: float((mm["t_" + t] <= A).mean()) for t in T},
                                 n_missed_null=int(len(mm)))
    print("n=%d  corrupted record deleted in %.3f; nothing deleted in %.3f" % (n, hit, none))
    print("      null datasets where the slip was NOT deleted (%d): size " % len(mm) +
          "  ".join("%s %.3f" % (t, (mm["t_" + t] <= A).mean()) for t in T))

print("\n(4) sign / outcome corruption (reported, no claim)")
for n in (30, 50):
    for ty in ("sign", "outcome"):
        s = null[(null.n == n) & (null.type == ty)]
        print("n=%d %-8s " % (n, ty) + "  ".join("%s %.3f->%.3f" % (t, (s["p_" + t] <= A).mean(), (s["t_" + t] <= A).mean()) for t in T))
        out["check4"]["%d|%s" % (n, ty)] = {t: [float((s["p_" + t] <= A).mean()), float((s["t_" + t] <= A).mean())] for t in T}
json.dump(out, open(os.path.join(HERE, "RESULTS_TRIM.json"), "w"), indent=2)
print("\nwrote RESULTS_TRIM.json")
