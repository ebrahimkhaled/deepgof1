"""07_pcrown.py -- W4: a PROBABILISTIC CROWN certificate for random (not adversarial) noise.

Why the ball certificate cannot work (06_diagnose.py): over the 95% l2 ball (r = 7.14) the noise may
point along any single neuron's weight vector with 7.14 standard deviations, while Gaussian noise
projects onto a fixed direction with sd 1. At the certifying point EVERY ReLU is undecided and the head
bounds are 124x too wide.

THE CERTIFICATE. Map m = c0 + e, e ~ N(0, Sigma). CROWN produces, for every unit j of every layer, linear
functions with z_j >= L_j(e) = lam_j' e + a_j and z_j <= U_j(e), valid WHENEVER all earlier-layer bounds
hold at e. Because lam_j' e is exactly N(0, lam_j' Sigma lam_j), set
      l_j = a_j - z_int * sd_j,     u_j = a'_j + z_int * sd'_j,
and let E = the event that every such Gaussian term stays inside its z_int band, over all ReLU units.
By induction over layers, on E every intermediate bound holds, so every relaxation is valid, so
S(c0 + e) >= lam' e + a (the final CROWN lower function). Hence
      P{S(c0 + e) <= t_bar} <= P(E^c) + P{lam' e + a <= t_bar}
                             <= delta_int + P{N(0,1) <= (t_bar - a)/sd}.
With delta_int = .025 split over all ReLU units (union bound) and the final term forced below .025,
power >= .95 is CERTIFIED for the Gaussian map model -- the same model the Monte Carlo truth t_emp uses.
No worst case over a ball is taken anywhere.

VALIDATION (the script stops on failure): at each certified t, draw 200,000 Gaussian maps; the share
violating ANY intermediate bound must be <= delta_int, and the Monte Carlo power must be >= the certified
power.

Run:  python -u 07_pcrown.py
"""
import os, json
import numpy as np
from scipy.stats import norm
from crown_core import NET, fwd, K

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "CERTIFIED.json")))
t_bar = C["t_bar"]
SIG_ROOT = np.eye(36)          # the idealised null map N(0, I), the model t_emp and t_bar are computed in
D_INT, D_FIN = 0.025, 0.025
N_RELU = int(sum(l["relu"].sum() for l in NET))
Z_INT = norm.ppf(1 - D_INT / (2 * N_RELU))
print("ReLU units %d -> intermediate band z = %.3f (ball certificate used r = %.2f)" % (N_RELU, Z_INT, C["R2"]))

def pcrown(Ls, c0, want_pre=False):
    """Probabilistic CROWN: returns (a, sd) of the final lower linear function and the bounds used."""
    bnds = []
    def backsub(t, sign):
        Lam = sign * Ls[t]["W"]; const = sign * Ls[t]["b"].copy()
        for s in range(t - 1, -1, -1):
            l, u = bnds[s]; rl = Ls[s]["relu"]
            dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
            off = rl & (u <= 0); dl[off] = 0; du[off] = 0
            un = rl & (l < 0) & (u > 0)
            du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
            dl[un] = (u[un] >= -l[un]).astype(float)
            Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
            const = const + Ln @ cu
            Lam = Lp * dl[None, :] + Ln * du[None, :]
            const = const + Lam @ Ls[s]["b"]
            Lam = Lam @ Ls[s]["W"]
        return Lam @ c0 + const, np.linalg.norm(Lam @ SIG_ROOT, axis=1)
    for t in range(len(Ls) - 1):
        a_lo, s_lo = backsub(t, 1.0); a_up, s_up = backsub(t, -1.0)
        bnds.append((a_lo - Z_INT * s_lo, -(a_up) + Z_INT * s_up))
    a, sd = backsub(len(Ls) - 1, 1.0)
    return (float(a[0]), float(sd[0]), bnds) if want_pre else (float(a[0]), float(sd[0]))

def cert_power(t, v):
    a, sd = pcrown(NET, t * v)
    return 1 - D_INT - norm.cdf((t_bar - a) / sd)          # certified lower bound on P(S > t_bar)

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}
rng = np.random.default_rng(777)
rows = []
print("\n%-10s %8s %8s %8s | %8s %8s %8s | validation" % ("direction", "t_emp", "t_pcert", "t_ball",
                                                         "n*_pc", "n*_emp", "ratio"))
for row in C["rows"]:
    name = row["direction"]; v = DIRS[name] / np.linalg.norm(DIRS[name])
    hi = 1.0
    while cert_power(hi, v) < 0.95:
        hi *= 1.5
        if hi > 1e5: break
    lo = hi / 1.5
    for _ in range(25):
        mid = (lo + hi) / 2
        if cert_power(mid, v) >= 0.95: hi = mid
        else: lo = mid
    tp = hi
    # validation at the certified point: bounds violated rarely enough, and true power above the claim
    a, sd, bnds = pcrown(NET, tp * v, want_pre=True)
    E = rng.normal(size=(200000, 36)) @ SIG_ROOT.T
    viol = np.zeros(len(E), bool); pw = []
    for s0 in range(0, len(E), 20000):
        val, pres = fwd(NET, tp * v + E[s0:s0 + 20000], keep=True)
        bad = np.zeros(pres[0].shape[1], bool)
        for t, (l, u) in enumerate(bnds):
            rl = NET[t]["relu"]
            bad |= ((pres[t][rl] < l[rl, None] - 1e-9) | (pres[t][rl] > u[rl, None] + 1e-9)).any(0)
        viol[s0:s0 + 20000] = bad; pw.append(val > t_bar)
    vr, mc = viol.mean(), np.concatenate(pw).mean()
    ok = vr <= D_INT and mc >= 0.95
    rows.append(dict(direction=name, t_emp=row["t_emp"], t_pcert=tp, t_ball=row.get("t_alpha", row["t_cert"]),
                     n_pcert=(tp / 0.5) ** 2, n_emp=row["n_emp"], ratio=(tp / row["t_emp"]) ** 2,
                     bound_violation_rate=float(vr), mc_power_at_cert=float(mc)))
    print("%-10s %8.2f %8.2f %8.2f | %8.0f %8.0f %8.1f | violations %.5f, MC power %.4f %s"
          % (name, row["t_emp"], tp, rows[-1]["t_ball"], (tp / 0.5) ** 2, row["n_emp"],
             (tp / row["t_emp"]) ** 2, vr, mc, "ok" if ok else "FAILED"))
    assert ok, "probabilistic certificate unsound -- bug"
json.dump(dict(t_bar=t_bar, z_int=Z_INT, delta_int=D_INT, delta_fin=D_FIN, rows=rows),
          open(os.path.join(HERE, "PCERT.json"), "w"), indent=2)
print("\nW4 exit test (n*_cert / n*_emp <= 10 for every direction): %s"
      % ("PASS" if all(r["ratio"] <= 10 for r in rows) else "not yet"))
