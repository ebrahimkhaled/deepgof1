"""07_certify.py -- the best SOUND certificate found by this prototype, with validation.

Method: probabilistic CROWN with the usage-aware band charge (bab_core.cell_pcrown_ua, k = 0) and the band
level z chosen from a fixed grid; smallest certified t found by bisection.  Then at the certified t, 200,000
Gaussian draws check
  (i)   the share of draws on which ANY relaxation actually used is invalid (a stable-active unit < 0, a
        stable-inactive unit > 0, or a used unstable unit outside [l, u]) is <= the charged budget,
  (ii)  on draws where every used relaxation is valid, S >= the final linear lower function (numerically),
  (iii) Monte Carlo power >= .95.
The script stops (assert) if any check fails.

Usage: python -u 07_certify.py <direction> [<direction> ...]
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
from bab_core import NET, fwd, DIRS, cell_pcrown_ua, T_BAR, N_RELU

ZG = np.round(np.arange(4.0, 4.71, 0.02), 2)
E0 = np.zeros((36, 0)); Z0 = np.zeros(0)
PC = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "PCERT.json")))
PCROWS = {r["direction"]: r for r in PC["rows"]}


def best_bound(c0):
    """Scan z upward; stop after the bound has risen 3 grid steps past its minimum (deterministic)."""
    best = (np.inf, None, None, None, None); worse = 0
    for z in ZG:
        a, sd, ch = cell_pcrown_ua(c0, E0, Z0, Z0, z)
        tail = norm.cdf((T_BAR - a) / sd) if np.isfinite(a) else 1.0
        b = ch + tail
        if b < best[0]: best = (b, z, a, sd, ch); worse = 0
        elif np.isfinite(best[0]):
            worse += 1
            if worse >= 3: break
    return best


def certify(d, log=print):
    v = DIRS[d]; t0 = time.time()
    hi = PCROWS[d]["t_pcert"] * 1.02
    assert best_bound(hi * v)[0] <= 0.05
    lo = hi * 0.8
    while best_bound(lo * v)[0] <= 0.05: lo *= 0.9
    for _ in range(9):
        mid = 0.5 * (lo + hi)
        if best_bound(mid * v)[0] <= 0.05: hi = mid
        else: lo = mid
    b, z, a, sd, ch = best_bound(hi * v)
    wall = time.time() - t0
    log("%s: certified t = %.3f (07_pcrown %.3f)  z=%.2f charge=%.4f (union %.4f) tail=%.4f  [%.0fs]"
        % (d, hi, PCROWS[d]["t_pcert"], z, ch, 2 * N_RELU * norm.sf(z), b - ch, wall))
    return dict(direction=d, t_cert=hi, t_pcrown07=PCROWS[d]["t_pcert"], t_emp=PCROWS[d]["t_emp"], z=float(z),
                charge=ch, union_charge=float(2 * N_RELU * norm.sf(z)), tail=float(b - ch), bound=float(b),
                a=a, sd=sd, wall_s=wall, n_regions=1)


def validate(d, t, z, n=200000, seed=777):
    v = DIRS[d]
    a, sd, ch, det = cell_pcrown_ua(t * v, E0, Z0, Z0, z, details=True)
    rng = np.random.default_rng(seed)
    viol = 0; below = 0; pw = 0; worst_gap = np.inf
    for s0 in range(0, n, 20000):
        E = rng.normal(size=(20000, 36))
        val, pres = fwd(NET, t * v + E, keep=True)
        bad = np.zeros(len(E), bool)
        for s, (l, u) in enumerate(det["bnds"]):
            rl = NET[s]["relu"]
            act = rl & (l >= 0); ina = rl & (u <= 0) & ~act; un = rl & (l < 0) & (u > 0) & det["used"][s]
            P = pres[s]
            bad |= (P[act] < 0).any(0) | (P[ina] > 0).any(0)
            bad |= ((P[un] < l[un, None] - 1e-9) | (P[un] > u[un, None] + 1e-9)).any(0)
        F = (t * v + E) @ det["lam"] + det["const"]
        gap = (val - F)[~bad]
        worst_gap = min(worst_gap, gap.min() if len(gap) else np.inf)
        below += (gap < -1e-7).sum()
        viol += bad.sum(); pw += (val > T_BAR).sum()
    return dict(violation_rate=viol / n, charge=ch, lower_fn_breaches=int(below), min_gap=float(worst_gap),
                mc_power=pw / n, ok=bool(viol / n <= ch and below == 0 and pw / n >= 0.95))


if __name__ == "__main__":
    out = []
    for d in sys.argv[1:]:
        r = certify(d)
        val = validate(d, r["t_cert"], r["z"])
        r["validation"] = val
        r["n_ratio_vs_emp"] = (r["t_cert"] / r["t_emp"]) ** 2
        print("   validation: used-relaxation violations %.5f (budget %.5f), lower-fn breaches %d (min gap %.3f), "
              "MC power %.4f -> %s" % (val["violation_rate"], val["charge"], val["lower_fn_breaches"],
                                       val["min_gap"], val["mc_power"], "ok" if val["ok"] else "FAILED"), flush=True)
        assert val["ok"], "certificate failed validation -- stop"
        out.append(r)
        json.dump(out, open("CERT_UA_%s.json" % "_".join(sys.argv[1:]), "w"), indent=1)
