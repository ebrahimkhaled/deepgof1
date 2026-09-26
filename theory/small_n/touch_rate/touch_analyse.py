"""touch_analyse.py -- W5(ii) analysis: reads cells/*.npz (touch_sim.py) and rows_power/*.csv
(power_check_d10.R), computes
  (1) the clean touch rate of the declared rule (c = 4.5) with exact binomial (Clopper-Pearson)
      95% intervals and the binomial SE, per law x d x n;
  (2) three analytic approximations for normal covariates and their error against the simulation:
      A1  naive union / extreme-value form   1 - (1 - 2 Phi(-c))^(n d)
      A2  A1 with the Croux-Rousseeuw finite-sample MAD factor: E[1.4826 MAD_n] ~ (n - 0.8)/n, so the
          robust z is inflated by n/(n - 0.8):  1 - (1 - 2 Phi(-c (n - 0.8)/n))^(n d)
      A3  the MAD-integrated form: within one covariate the n robust z's share the same sample MAD s,
          so touch = 1 - ( E_s[(1 - 2 Phi(-c s / tau))^n] )^d, with s ~ N((n-.8)/n, 1.361/n) (the
          scaled MAD's asymptotic law under normality, efficiency .367) and tau^2 = 1 + pi/(2n) (the
          sample median's variance added to each deviation); Gauss-Hermite quadrature over s
  (3) the miss probability of the x10 decimal slip at c = 4.5 (overall, gross |x| > .45, mild);
  (4) a dimension-adjusted threshold t(n, d): the Bonferroni form
          t_B(n, d) = (n/(n - 0.8)) Phi^{-1}(1 - .025/(n d))
      and the A3-calibrated t_5(n, d) solving touch = .05; their clean touch rates and slip-miss
      probabilities from the SAME stored replicates (any threshold is evaluable after the fact);
  (5) the paired power cost of one covariate-selected deletion at d = 10 (McNemar-type paired SE).
Writes TOUCH_RATE.json and prints the tables that RESULT.md quotes.
"""
import os, glob, json
import numpy as np
from scipy import stats, optimize

HERE = os.path.dirname(os.path.abspath(__file__))
C0 = 4.5
NS = (30, 50, 100, 200, 500); DS = (2, 5, 10)
LAWS = ("normal", "t4", "unif", "lognormal", "bern", "liu2")


def cp_int(k, n, conf=0.95):
    lo = 0.0 if k == 0 else stats.beta.ppf((1 - conf) / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - (1 - conf) / 2, k + 1, n - k)
    return float(lo), float(hi)


def rate(mask):
    k = int(mask.sum()); n = int(mask.size); r = k / n
    return dict(rate=r, se=float(np.sqrt(r * (1 - r) / n)), ci=cp_int(k, n), k=k, R=n)


# ---------- analytic approximations (normal covariates) ----------
def A1(c, n, d):
    return 1 - (1 - 2 * stats.norm.cdf(-c)) ** (n * d)


def A2(c, n, d):
    return 1 - (1 - 2 * stats.norm.cdf(-c * (n - 0.8) / n)) ** (n * d)


_GH = np.polynomial.hermite_e.hermegauss(60)          # probabilists' Hermite: weight exp(-x^2/2)


def A3(c, n, d):
    mu, sd = (n - 0.8) / n, np.sqrt(1.361 / n)
    tau = np.sqrt(1 + np.pi / (2 * n))
    x, w = _GH
    s = mu + sd * x
    s = np.maximum(s, 1e-3)
    inner = (1 - 2 * stats.norm.cdf(-c * s / tau)) ** n
    Es = float(np.sum(w * inner) / np.sqrt(2 * np.pi))
    return 1 - Es ** d


def A4(c, n, d):
    """A3 with a LOGNORMAL law for the scaled MAD (same mean (n-.8)/n and variance 1.361/n): the
    finite-sample law of the MAD is right-skewed, and the normal form of A3 puts too much mass on
    small MADs, each of which inflates all n robust z's of that covariate."""
    mu, var = (n - 0.8) / n, 1.361 / n
    v = np.log(1 + var / mu ** 2); m = np.log(mu) - v / 2
    tau = np.sqrt(1 + np.pi / (2 * n))
    x, w = _GH
    s = np.exp(m + np.sqrt(v) * x)
    inner = (1 - 2 * stats.norm.cdf(-c * s / tau)) ** n
    Es = float(np.sum(w * inner) / np.sqrt(2 * np.pi))
    return 1 - Es ** d


def t_bonf(n, d, alpha=0.05):
    return float(n / (n - 0.8) * stats.norm.ppf(1 - alpha / 2 / (n * d)))


def t_A3(n, d, target=0.05):
    return float(optimize.brentq(lambda c: A3(c, n, d) - target, 2.0, 12.0))


# ---------- load cells ----------
cells = {}
for f in sorted(glob.glob(os.path.join(HERE, "cells", "*.npz"))):
    z = np.load(f)
    cells[(str(z["law"]), int(z["d"]), int(z["n"]))] = {k: z[k] for k in ("omax_clean", "o0_slip", "oth_slip", "x01")}
print("cells loaded: %d" % len(cells))


def touch(cell, t):
    return rate(cell["omax_clean"] > t)


def miss(cell, t):
    caught = (cell["o0_slip"] > t) & (cell["o0_slip"] > cell["oth_slip"])
    gross = np.abs(cell["x01"]) > 0.45
    out = dict(overall=rate(~caught))
    if gross.any() and (~gross).any():
        out["gross"] = rate(~caught[gross]); out["mild"] = rate(~caught[~gross])
    return out


out = dict(rule="delete argmax_i max_j |x_ij - med_j| / (1.4826 MAD_j) if it exceeds t; declared t = 4.5",
           R=20000, touch_c45={}, analytic_normal={}, miss_c45={}, adjusted={}, quantiles={}, power_cost={})

# (1) touch rate at the declared threshold, all laws
print("\n(1) CLEAN TOUCH RATE at c = 4.5 (rate, binomial SE)")
for law in LAWS:
    dlist = (3,) if law == "liu2" else DS
    for d in dlist:
        row = "%-9s d=%2d " % (law, d)
        for n in NS:
            if (law, d, n) not in cells: continue
            r = touch(cells[(law, d, n)], C0)
            out["touch_c45"]["%s|%d|%d" % (law, d, n)] = r
            row += "  n=%3d %.4f (%.4f)" % (n, r["rate"], r["se"])
        print(row)

# (2) analytic vs simulated, normal
print("\n(2) NORMAL: simulated vs A1 (naive union) / A2 (+MAD factor) / A3 (MAD-integrated)")
for d in DS:
    for n in NS:
        r = touch(cells[("normal", d, n)], C0)
        a = dict(sim=r["rate"], se=r["se"], A1=A1(C0, n, d), A2=A2(C0, n, d), A3=A3(C0, n, d), A4=A4(C0, n, d))
        a["z_A3"] = (a["A3"] - a["sim"]) / max(r["se"], 1e-6); a["z_A4"] = (a["A4"] - a["sim"]) / max(r["se"], 1e-6)
        out["analytic_normal"]["%d|%d" % (d, n)] = a
        print("d=%2d n=%3d  sim %.4f (%.4f)   A1 %.4f  A2 %.4f  A3 %.4f  A4 %.4f   (A3-sim)/SE %+.1f  (A4-sim)/SE %+.1f"
              % (d, n, a["sim"], a["se"], a["A1"], a["A2"], a["A3"], a["A4"], a["z_A3"], a["z_A4"]))

# (3) miss probability at 4.5, normal (and the other laws, overall only)
print("\n(3) x10 DECIMAL SLIP MISSED at c = 4.5, normal: overall / gross |x|>.45 / mild")
for d in DS:
    row = "d=%2d " % d
    for n in NS:
        m = miss(cells[("normal", d, n)], C0)
        out["miss_c45"]["normal|%d|%d" % (d, n)] = m
        row += "  n=%3d %.3f/%.3f/%.3f" % (n, m["overall"]["rate"], m["gross"]["rate"], m["mild"]["rate"])
    print(row)
for law in ("t4", "lognormal", "liu2"):
    dlist = (3,) if law == "liu2" else DS
    for d in dlist:
        row = "%-9s d=%2d overall miss" % (law, d)
        for n in NS:
            m = miss(cells[(law, d, n)], C0); out["miss_c45"]["%s|%d|%d" % (law, d, n)] = m
            row += "  n=%3d %.3f" % (n, m["overall"]["rate"])
        print(row)

# (4) dimension-adjusted thresholds
print("\n(4) ADJUSTED THRESHOLDS (normal): t_B Bonferroni, t_5 A3-calibrated to 5%%, empirical q95 of max o")
for d in DS:
    for n in NS:
        cell = cells[("normal", d, n)]
        tb, t5 = t_bonf(n, d), t_A3(n, d)
        q95 = float(np.quantile(cell["omax_clean"], 0.95))
        tq = float(np.ceil(q95 * 10) / 10)          # the empirical q95 rounded UP to one decimal
        rec = dict(t_bonf=tb, t_A3=t5, q95_empirical=q95, t_q95_rounded=tq,
                   touch_tB=touch(cell, tb), touch_t5=touch(cell, t5), touch_tq=touch(cell, tq),
                   miss_tB=miss(cell, tb), miss_t5=miss(cell, t5), miss_tq=miss(cell, tq))
        out["adjusted"]["normal|%d|%d" % (d, n)] = rec
        print("d=%2d n=%3d  t_B %.2f (touch %.3f, miss %.3f gross %.3f)   t_5 %.2f (touch %.3f, miss %.3f gross %.3f)   q95 %.2f -> %.1f (touch %.3f, miss %.3f gross %.3f)"
              % (d, n, tb, rec["touch_tB"]["rate"], rec["miss_tB"]["overall"]["rate"], rec["miss_tB"]["gross"]["rate"],
                 t5, rec["touch_t5"]["rate"], rec["miss_t5"]["overall"]["rate"], rec["miss_t5"]["gross"]["rate"],
                 q95, tq, rec["touch_tq"]["rate"], rec["miss_tq"]["overall"]["rate"], rec["miss_tq"]["gross"]["rate"]))
# the adjusted threshold applied to the heavy-tailed and skewed laws (sensitivity)
for law in ("t4", "lognormal"):
    for d in DS:
        row = "%-9s d=%2d touch at t_B" % (law, d)
        for n in NS:
            r = touch(cells[(law, d, n)], t_bonf(n, d)); out["adjusted"]["%s|%d|%d" % (law, d, n)] = dict(t_bonf=t_bonf(n, d), touch_tB=r)
            row += "  n=%3d %.3f" % (n, r["rate"])
        print(row)
for key, cell in cells.items():
    fin = cell["omax_clean"][np.isfinite(cell["omax_clean"])]
    q = np.quantile(fin, [0.5, 0.9, 0.95, 0.99]) if fin.size else [float("nan")] * 4
    out["quantiles"]["%s|%d|%d" % key] = dict(q50=float(q[0]), q90=float(q[1]), q95=float(q[2]), q99=float(q[3]),
                                             frac_inf=float(np.isinf(cell["omax_clean"]).mean()))
print("\nbern (MAD = 0 -> o = inf) share of clean designs with an infinite max o:",
      {"%d|%d" % (d, n): round(out["quantiles"]["bern|%d|%d" % (d, n)]["frac_inf"], 3) for d in DS for n in NS})

# (5) paired power cost at d = 10
rows = sorted(glob.glob(os.path.join(HERE, "rows_power", "r*.csv")))
if rows:
    import pandas as pd
    P = pd.concat([pd.read_csv(f) for f in rows], ignore_index=True)
    print("\n(5) PAIRED POWER COST of deleting the most outlying record, d=10 n=100 C=1 (%d replicates)" % len(P))
    for t in ("DG",):
        ok = P[["p_" + t, "t_" + t]].notna().all(1); Q = P[ok]
        a = (Q["p_" + t] <= 0.05).values; b = (Q["t_" + t] <= 0.05).values
        n10, n01 = int((a & ~b).sum()), int((~a & b).sum()); m = len(Q)
        diff = b.mean() - a.mean(); se = np.sqrt((n10 + n01) / m - diff ** 2) / np.sqrt(m)
        pmc = float(stats.binomtest(n01, n10 + n01, 0.5).pvalue) if n10 + n01 > 0 else 1.0
        rec = dict(test=t, reps=m, power_n=float(a.mean()), power_n_minus_1=float(b.mean()), diff=float(diff),
                   paired_se=float(se), ci95=[float(diff - 1.96 * se), float(diff + 1.96 * se)],
                   discordant=[n10, n01], mcnemar_exact_p=pmc, omax_mean=float(P["omax"].mean()),
                   share_omax_gt_45=float((P["omax"] > 4.5).mean()))
        out["power_cost"][t] = rec
        print("%-4s power n=100 %.3f -> n=99 %.3f  diff %+.3f (paired SE %.3f, 95%% CI %+.3f..%+.3f)  discordant %d:%d  McNemar p=%.3f"
              % (t, rec["power_n"], rec["power_n_minus_1"], diff, se, rec["ci95"][0], rec["ci95"][1], n10, n01, pmc))
    print("     share of these clean d=10 designs the declared rule would touch: %.3f" % ((P["omax"] > 4.5).mean()))

json.dump(out, open(os.path.join(HERE, "TOUCH_RATE.json"), "w"), indent=1)
print("\nwrote TOUCH_RATE.json")
