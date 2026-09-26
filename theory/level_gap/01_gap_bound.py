"""01_gap_bound.py -- W3: a data-computable confidence bound on the calibration gap of Proposition 2.

Proposition 2: size <= alpha + eps + eta, with Delta+ = sup_t (F_thetahat(t) - F_theta0(t))_+ the gap
between the bootstrap law and the true null law of the score, conditional on the covariates, and
eta = P(Delta+ > eps). The AE's objection: eps is unknown. THEOREM G (THEOREM_G.md): for ANY confidence
set C with P(theta0 in C) >= 1 - gamma,
     eps_hat = sup_{theta in C} Dhat+(thetahat, theta) + (Monte Carlo slack)
satisfies P(Delta+ > eps_hat) <= gamma + delta_MC. And the gap-adjusted p-value
p_tilde = p + eps_hat is valid at level alpha + gamma + delta_MC at EVERY sample size.

This script, per outer dataset (covariates held fixed, as Proposition 2 does):
  * F_thetahat and F_theta0 by N simulated datasets each (batched Newton refits; a failed refit scores
    +inf, the deployed convention);
  * the TRUE gap Delta+ (full sup, and the tail-localised sup over t >= the 80% point);
  * eps_hat over the 95% Wald ellipsoid (30 points: 6 principal-axis extremes, 16 boundary, 8
    interior) plus the DKW slack for the Monte Carlo cdfs, union-bounded;
  * coverage: is Delta+ <= eps_hat?  (W3 exit: >= 95% of datasets)
  * the deployed p (B = infinity limit from F_thetahat) and p_tilde, for size and power.

Run:  python -u 01_gap_bound.py      Env: G_N (20000), G_OUTER (100), G_WORKERS (8), G_NS (50,100,200)
"""
import os, sys, json, time
import numpy as np
from scipy.special import expit
from scipy.stats import chi2
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "local_power"))

HERE = os.path.dirname(os.path.abspath(__file__))
NSIM = int(os.environ.get("G_N", "20000"))
OUTER = int(os.environ.get("G_OUTER", "100"))
NW = int(os.environ.get("G_WORKERS", "8"))
NS = tuple(int(s) for s in os.environ.get("G_NS", "50,100,200").split(","))
BETA0 = np.array([-0.5, 1.0, 0.6]); K = 6; ALPHA = 0.05
D_MC = 0.01                                     # total Monte Carlo failure probability
TAIL_Q = 0.85                                   # tail localisation: t above the 85% point of F_thetahat

def cells_of(x):
    n = len(x)
    r1 = np.empty(n, int); r1[np.argsort(x[:, 0], kind="stable")] = np.arange(n)
    r2 = np.empty(n, int); r2[np.argsort(x[:, 1], kind="stable")] = np.arange(n)
    c = np.minimum(K * r1 // n, K - 1) * K + np.minimum(K * r2 // n, K - 1)
    O = np.zeros((n, K * K)); O[np.arange(n), c] = 1
    return O

def fit_batch(X, Y, beta, iters=25):
    """Batched Newton; returns beta and a convergence flag per row (finite, small score)."""
    for _ in range(iters):
        P = expit(beta @ X.T); Wt = P * (1 - P)
        grad = (Y - P) @ X
        H = np.einsum("bn,ni,nj->bij", Wt, X, X) + 1e-10 * np.eye(X.shape[1])
        step = np.linalg.solve(H, grad[:, :, None])[:, :, 0]
        beta = beta + np.clip(step, -5, 5)
    P = expit(beta @ X.T)
    g = np.abs((Y - P) @ X).max(1)
    ok = np.isfinite(beta).all(1) & (g < 1e-6 * X.shape[0]) & (np.abs(beta).max(1) < 15)
    return beta, ok

def scores_at(theta, X, O, N, rng, score, chunk=5000, want_u=False):
    """N draws of the score when y ~ Bern(pi_theta(X)), each refitted and mapped; failed refit = +inf.
    want_u: also return the data score U = X'(y - pi_theta) of every draw (unsorted order kept)."""
    out = np.empty(N); U = np.empty((N, X.shape[1]))
    pi = expit(X @ theta)
    for a in range(0, N, chunk):
        b = min(N, a + chunk)
        Y = (rng.uniform(size=(b - a, len(pi))) < pi[None, :]).astype(float)
        U[a:b] = (Y - pi[None, :]) @ X
        be, ok = fit_batch(X, Y, np.tile(theta, (b - a, 1)))
        P = expit(be @ X.T)
        M = ((Y - P) @ O) / np.sqrt(np.maximum((P * (1 - P)) @ O, 1e-8))
        s = np.full(b - a, np.inf)
        if ok.any():
            s[ok] = score(M[ok])
        out[a:b] = s
    return (out, U) if want_u else np.sort(out)

def eps_gradient(s, U, Iinv, rng, qs=(0.95, 0.99), M=20000):
    """The gradient-quantile bound. dF_theta(t)/dtheta at thetahat = E[1{T<=t} U] (likelihood-ratio
    derivative, estimated from the SAME draws at thetahat). Delta+ ~ sup_t (g(t)'Z)_+ with
    Z ~ N(0, I^-1) the sampling law of thetahat - theta0; return its quantiles."""
    fin = np.isfinite(s)
    # tail-localised (Proposition S.8): only the upper tail of the score law enters the size at alpha
    grid = np.quantile(s[fin], np.linspace(TAIL_Q, 0.995, 40))
    Uc = U - U.mean(0)                                                    # E[U] = 0 exactly: centre it
    G = np.stack([(Uc * (s <= t)[:, None]).mean(0) for t in grid])       # 40 x p
    Z = rng.multivariate_normal(np.zeros(Iinv.shape[0]), Iinv, size=M)      # M x p
    sup = np.clip((Z @ G.T).max(1), 0, None)
    return {q: float(np.quantile(sup, q)) for q in qs}

def dplus(Fa, Fb, tmin=-np.inf):
    """sup_t (F_a(t) - F_b(t))_+ over t >= tmin, for sorted samples Fa, Fb."""
    grid = np.concatenate([Fa, Fb]); grid = grid[np.isfinite(grid) & (grid >= tmin)]
    if len(grid) == 0: return 0.0
    ca = np.searchsorted(Fa, grid, side="right") / len(Fa)
    cb = np.searchsorted(Fb, grid, side="right") / len(Fb)
    return float(max(0.0, (ca - cb).max()))

def one_dataset(args):
    n, rep, C, seed = args
    from lp_core import score
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2)); X = np.column_stack([np.ones(n), x]); O = cells_of(x)
    eta = X @ BETA0 + C * (x[:, 0] ** 2 - 1)
    y = (rng.uniform(size=n) < expit(eta)).astype(float)
    bh, ok = fit_batch(X, y[None, :], np.zeros((1, 3)))
    if not ok[0]:
        return dict(n=n, rep=rep, C=C, fitted=False)
    bh = bh[0]; ph = expit(X @ bh)
    m = ((y - ph) @ O) / np.sqrt(np.maximum((ph * (1 - ph)) @ O, 1e-8))
    T = float(score(m[None, :])[0])
    s_raw, U = scores_at(bh, X, O, NSIM, rng, score, want_u=True)
    Fh = np.sort(s_raw)
    p = float((1 + np.sum(Fh >= T)) / (NSIM + 1))
    Iinf = (X * (ph * (1 - ph))[:, None]).T @ X
    Iinv = np.linalg.inv(Iinf)
    # (1) the gradient-quantile bound: no extra simulation at all
    eq = eps_gradient(s_raw, U, Iinv, rng)
    # (2) the worst case over the 95% Wald ellipsoid (the MMC-type bound), principal-axis extremes only
    r = np.sqrt(chi2.ppf(0.95, 3))
    ev, evec = np.linalg.eigh(Iinv)
    pts = [bh + s * r * np.sqrt(ev[j]) * evec[:, j] for j in range(3) for s in (-1, 1)]
    tmin = np.quantile(Fh[np.isfinite(Fh)], TAIL_Q)
    d_hat = 0.0
    if os.environ.get("G_WC", "0") == "1":               # the worst-case (MMC-type) arm: 6 extra cdfs
        for th in pts:
            d_hat = max(d_hat, dplus(Fh, scores_at(th, X, O, NSIM, rng, score), tmin))
    else:
        d_hat = float("nan")
    eps_dkw = np.sqrt(np.log(2 * 7 / D_MC) / (2 * NSIM))
    out = dict(n=n, rep=rep, C=C, fitted=True, T=T, p=p,
               eps_q95=eq[0.95], eps_q99=eq[0.99], eps_wc=d_hat + 2 * eps_dkw, d_hat=d_hat, eps_dkw=eps_dkw,
               theta0_in_wald=bool((bh - BETA0) @ Iinf @ (bh - BETA0) <= r ** 2))
    if C == 0:
        F0 = scores_at(BETA0, X, O, NSIM, rng, score)
        out["gap"] = dplus(Fh, F0, tmin)                  # the tail-localised gap of Proposition S.8
    return out

if __name__ == "__main__":
    jobs = []
    for n in NS:
        jobs += [(n, r, 0.0, 9_100_000 + n * 1000 + r) for r in range(OUTER)]
        jobs += [(n, r, 1.0, 9_300_000 + n * 1000 + r) for r in range(OUTER // 2)]
    # every finished dataset goes to disk at once and a rerun skips it: the first full run hung at
    # 1325/1350 when a worker was killed under memory pressure, with all results still in RAM
    RAW = os.path.join(HERE, "GAP_RAW.jsonl")
    done = set()
    if os.path.exists(RAW):
        for line in open(RAW):
            r = json.loads(line); done.add((r["n"], r["rep"], r["C"]))
    todo = [j for j in jobs if (j[0], j[1], j[2]) not in done]
    print("W3: %d outer datasets to run (%d on disk), N = %d per cdf, %d workers"
          % (len(todo), len(done), NSIM, NW), flush=True)
    t0 = time.time()
    with Pool(NW, maxtasksperchild=20) as pool, open(RAW, "a") as fh:
        for i, r in enumerate(pool.imap_unordered(one_dataset, todo)):
            fh.write(json.dumps(r) + "\n"); fh.flush()
            if (i + 1) % 25 == 0:
                print("  %d / %d (%.1f min)" % (i + 1, len(todo), (time.time() - t0) / 60), flush=True)
    res = [json.loads(line) for line in open(RAW)]
    for n in NS:
        R0 = [r for r in res if r["n"] == n and r["C"] == 0 and r["fitted"]]
        R1 = [r for r in res if r["n"] == n and r["C"] == 1 and r["fitted"]]
        g = np.array([r["gap"] for r in R0])
        print("\nn=%d  null datasets %d (unfitted %d)" % (n, len(R0), sum(1 for r in res if r["n"] == n and not r["fitted"])))
        print("   true gap Delta+: median %.4f  q90 %.4f  max %.4f" % (np.median(g), np.quantile(g, .9), g.max()))
        for key in ("eps_q95", "eps_q99", "eps_wc"):
            e = np.array([r[key] for r in R0])
            if not np.isfinite(e).all(): continue
            cv = np.mean(g <= e)
            se = np.sqrt(max(cv * (1 - cv), 1e-4) / len(g))
            print("   %-8s median %.4f q90 %.4f | COVERAGE P(Delta+ <= eps) = %.3f (SE %.3f)"
                  % (key, np.median(e), np.quantile(e, .9), cv, se))
        print("   size: deployed %.3f | p + eps_q99 %.3f | p + eps_wc %.3f"
              % (np.mean([r["p"] <= ALPHA for r in R0]), np.mean([r["p"] + r["eps_q99"] <= ALPHA for r in R0]),
                 np.mean([r["p"] + r["eps_wc"] <= ALPHA for r in R0])))
        if R1:
            print("   power (C=1): deployed %.3f | p + eps_q99 %.3f | p + eps_wc %.3f" % (
                np.mean([r["p"] <= ALPHA for r in R1]), np.mean([r["p"] + r["eps_q99"] <= ALPHA for r in R1]),
                np.mean([r["p"] + r["eps_wc"] <= ALPHA for r in R1])))
