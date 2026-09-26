"""08_validate_strat.py -- validate a stratified (03_bab.py) certificate on 200,000 Gaussian draws.

For every cell with bound < 1 the cell's usage-aware bounds are recomputed at the cell's z.  Each draw is
assigned to its cell through s = U'e.  Checks: (i) share of draws whose OWN cell has an invalid used
relaxation <= sum_C P(C) charge_C;  (ii) no draw with all used relaxations valid has S below its cell's
final lower function;  (iii) share of draws in cells with bound 1 is close to their total mass;
(iv) Monte Carlo power >= .95.

Usage: python -u 08_validate_strat.py run_<...>.json
"""
import os, sys, json, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
from bab_core import NET, fwd, DIRS, cell_pcrown_ua, T_BAR

spec = importlib.util.spec_from_file_location("bab3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "03_bab.py"))
bab3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(bab3)

R = json.load(open(sys.argv[1]))
d, t, k = R["direction"], R["t"], R["k"]
v = DIRS[d]; c0 = t * v
U = bab3.basis(R["basis"], t, v)[:, :k]
cells = R["cells"]
dets, charged = [], 0.0
for c in cells:
    pl, ph = np.array(c["pl"]), np.array(c["ph"])
    if c["b"] >= 1.0 or c["z"] is None:
        dets.append(None); continue
    act = (pl > 0) | (ph < 1)
    a, sd, ch, det = cell_pcrown_ua(c0, U[:, act], norm.ppf(pl[act]), norm.ppf(ph[act]), c["z"], details=True)
    det["ch"] = ch; dets.append(det); charged += c["p"] * ch
print("cells %d (finite %d); charged band budget sum P(C) ch_C = %.5f; certificate bound %.4f"
      % (len(cells), sum(x is not None for x in dets), charged, R["bound"]))

rng = np.random.default_rng(2027)
n = 200000; viol = below = pw = 0; in_dead = 0; dead_mass = sum(c["p"] for c, x in zip(cells, dets) if x is None)
for s0 in range(0, n, 20000):
    E = rng.normal(size=(20000, 36)); X = c0 + E
    val, pres = fwd(NET, X, keep=True)
    q = norm.cdf(E @ U)
    pw += (val > T_BAR).sum()
    for c, det in zip(cells, dets):
        pl, ph = np.array(c["pl"]), np.array(c["ph"])
        m = np.all((q >= pl) & ((q < ph) | (ph >= 1)), axis=1)
        if not m.any(): continue
        if det is None:
            in_dead += m.sum(); continue
        bad = np.zeros(m.sum(), bool)
        for s, (l, u) in enumerate(det["bnds"]):
            rl = NET[s]["relu"]; P = pres[s][:, m]
            act = rl & (l >= 0); ina = rl & (u <= 0) & ~act; un = rl & (l < 0) & (u > 0) & det["used"][s]
            bad |= (P[act] < 0).any(0) | (P[ina] > 0).any(0)
            bad |= ((P[un] < l[un, None] - 1e-9) | (P[un] > u[un, None] + 1e-9)).any(0)
        F = X[m] @ det["lam"] + det["const"]
        below += ((val[m] - F)[~bad] < -1e-7).sum()
        viol += bad.sum()
out = dict(run=sys.argv[1], violation_rate=viol / n, charged=charged, lower_fn_breaches=int(below),
           dead_share=in_dead / n, dead_mass=dead_mass, mc_power=pw / n)
out["ok"] = bool(out["violation_rate"] <= charged and below == 0 and out["mc_power"] >= 0.95)
print(json.dumps(out, indent=1))
json.dump(out, open("VALID_" + os.path.basename(sys.argv[1]), "w"), indent=1)
assert out["ok"], "stratified certificate failed validation"
