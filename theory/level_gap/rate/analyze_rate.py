"""analyze_rate.py -- level-rate study: sizes, |size - alpha| with SEs, and the log-log slope.

Reads RATE_RAW/*.csv (one row per null replicate, p = released deepgof1 p-value, B = 199),
writes RATE.json. The slope of log|size - alpha| on log n is fitted by weighted least squares
with delta-method weights, and its sampling distribution is checked by a replicate-level
bootstrap (resampling p-values within each cell), because |error| is not a smooth statistic
when the error is inside its own standard error.
"""
import glob, json, math, os, csv, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHAS = [0.05, 0.10]
NS = [50, 100, 200, 400, 800]
DES = ["D1", "D4", "D3"]
rng = np.random.default_rng(20260924)

cells = collections.defaultdict(list)
meta = collections.defaultdict(lambda: {"sep": 0, "ninf": 0, "n": 0})
for f in glob.glob(os.path.join(HERE, "RATE_RAW", "*.csv")):
    with open(f) as fh:
        for r in csv.DictReader(fh):
            k = (r["design"], int(r["n"]))
            cells[k].append(float(r["p"]))
            meta[k]["sep"] += int(r["sep"]); meta[k]["ninf"] += int(float(r["ninf"])); meta[k]["n"] += 1
P = {k: np.array(v) for k, v in cells.items()}


def size_row(p, a):
    R = len(p); s = float(np.mean(p <= a + 1e-12))
    se = math.sqrt(s * (1 - s) / R)
    return {"reps": R, "size": s, "se": se, "err": s - a, "abs_err": abs(s - a), "z": (s - a) / se}


def wls_slope(ns, errs, ses):
    """log|err| = a + b log n, weights 1/Var(log|err|) with Var = (se/|err|)^2 (delta method).
    Returns b, se_b, and how many points have |err| < 2 se (those carry almost no information)."""
    x = np.log(np.array(ns, float)); e = np.array(errs, float); s = np.array(ses, float)
    ae = np.abs(e); w = (ae / s) ** 2          # inverse variance of log|err|
    y = np.log(np.maximum(ae, 1e-12))
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    XtWX = X.T @ W @ X
    beta = np.linalg.solve(XtWX, X.T @ W @ y)
    cov = np.linalg.inv(XtWX)
    return float(beta[1]), float(math.sqrt(cov[1, 1])), int(np.sum(ae < 2 * s))


def ols_slope(ns, errs):
    x = np.log(np.array(ns, float)); y = np.log(np.maximum(np.abs(np.array(errs, float)), 1e-12))
    X = np.column_stack([np.ones_like(x), x])
    beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(x) - 2
    if dof <= 0:
        return float(beta[1]), float("nan")
    s2 = float(np.sum((y - X @ beta) ** 2) / dof)
    cov = s2 * np.linalg.inv(X.T @ X)
    return float(beta[1]), float(math.sqrt(cov[1, 1]))


out = {"design": {}, "pooled": {}, "slopes": {}, "notes": []}
for d in DES:
    out["design"][d] = {}
    for n in NS:
        p = P[(d, n)]
        m = meta[(d, n)]
        out["design"][d][str(n)] = {str(a): size_row(p, a) for a in ALPHAS}
        out["design"][d][str(n)]["sep_rate"] = m["sep"] / m["n"]
        out["design"][d][str(n)]["boot_refit_failures_per_rep"] = m["ninf"] / m["n"]

for n in NS:
    p = np.concatenate([P[(d, n)] for d in DES])
    out["pooled"][str(n)] = {str(a): size_row(p, a) for a in ALPHAS}

# slopes: per design and pooled, at each alpha; WLS (delta) + OLS + replicate bootstrap
def boot_slopes(getp, a, Bb=2000):
    sl = []
    for _ in range(Bb):
        errs = []
        for n in NS:
            p = getp(n); pb = rng.choice(p, size=len(p), replace=True)
            errs.append(float(np.mean(pb <= a + 1e-12)) - a)
        sl.append(ols_slope(NS, errs)[0])
    sl = np.array(sl)
    return {"mean": float(sl.mean()), "sd": float(sl.std(ddof=1)),
            "q025": float(np.quantile(sl, .025)), "q975": float(np.quantile(sl, .975))}

for label, getp in [("pooled", lambda n: np.concatenate([P[(d, n)] for d in DES]))] + \
                   [(d, (lambda d: (lambda n: P[(d, n)]))(d)) for d in DES]:
    out["slopes"][label] = {}
    for a in ALPHAS:
        rows = [size_row(getp(n), a) for n in NS]
        errs = [r["err"] for r in rows]; ses = [r["se"] for r in rows]
        b, seb, nweak = wls_slope(NS, errs, ses)
        bo, sebo = ols_slope(NS, errs)
        out["slopes"][label][str(a)] = {
            "wls_delta": {"slope": b, "se": seb, "points_within_2se_of_zero": nweak},
            "ols": {"slope": bo, "se": sebo},
            "replicate_bootstrap_of_ols_slope": boot_slopes(getp, a),
            "abs_errors": [r["abs_err"] for r in rows], "ses": ses, "signed_errors": errs}

# a direct test of "no detectable n-trend": does the signed error depend on n? (chi-square, pooled)
for a in ALPHAS:
    rows = [out["pooled"][str(n)][str(a)] for n in NS]
    zs = np.array([r["z"] for r in rows])
    out["pooled"]["heterogeneity_" + str(a)] = {"chi2_4df": float(np.sum(zs ** 2)),
                                                "max_abs_z": float(np.max(np.abs(zs)))}
out["notes"].append("size = P(p <= alpha); p is on the grid k/200 so alpha = .05 and .10 are attainable (k = 10, 20).")
out["notes"].append("SE = sqrt(s(1-s)/R); at s = .05, R = 7,600 gives .0025; pooled R = 22,800 gives .0014.")
out["notes"].append("A log|error| slope is only meaningful when every |error| is well above its SE; the field points_within_2se_of_zero counts the cells where it is not.")

with open(os.path.join(HERE, "RATE.json"), "w") as fh:
    json.dump(out, fh, indent=1)

# console summary
print("cell           size.05    SE     |err|    z   |  size.10    SE     |err|    z   | sep  ninf")
for d in DES:
    for n in NS:
        c = out["design"][d][str(n)]; a, b = c["0.05"], c["0.1"]
        print(f"{d} n={n:4d}   {a['size']:.4f}  {a['se']:.4f}  {a['abs_err']:.4f} {a['z']:5.2f} | "
              f"{b['size']:.4f}  {b['se']:.4f}  {b['abs_err']:.4f} {b['z']:5.2f} | {c['sep_rate']:.4f} {c['boot_refit_failures_per_rep']:.4f}")
for n in NS:
    a, b = out["pooled"][str(n)]["0.05"], out["pooled"][str(n)]["0.1"]
    print(f"pooled n={n:4d} {a['size']:.4f}  {a['se']:.4f}  {a['abs_err']:.4f} {a['z']:5.2f} | "
          f"{b['size']:.4f}  {b['se']:.4f}  {b['abs_err']:.4f} {b['z']:5.2f}")
for lab in out["slopes"]:
    for a in ALPHAS:
        s = out["slopes"][lab][str(a)]
        print(f"slope {lab} a={a}: WLS {s['wls_delta']['slope']:+.2f} (SE {s['wls_delta']['se']:.2f}, "
              f"{s['wls_delta']['points_within_2se_of_zero']}/5 pts within 2SE of 0) | OLS {s['ols']['slope']:+.2f} (SE {s['ols']['se']:.2f}) | "
              f"boot sd {s['replicate_bootstrap_of_ols_slope']['sd']:.2f} 95% [{s['replicate_bootstrap_of_ols_slope']['q025']:+.2f}, {s['replicate_bootstrap_of_ols_slope']['q975']:+.2f}]")
for a in ALPHAS:
    h = out["pooled"]["heterogeneity_" + str(a)]
    print(f"pooled heterogeneity a={a}: chi2(4)={h['chi2_4df']:.2f}, max|z|={h['max_abs_z']:.2f}")
