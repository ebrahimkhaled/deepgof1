"""05_accounting_ceiling.py -- DIAGNOSTIC (not a certificate): the best any sharper bound on P(E^c) could do.

probabilistic CROWN's band event E_z (all ReLU bands at z sd) is charged 2 N_RELU Phi(-z) by the union
bound.  Replace that charge by the TRUE P(E_z^c), estimated from 40,000 draws (an oracle for any
Bonferroni-improving bound such as Hunter-Worsley).  For each t, report the certificate
    true P(E_z^c) + Phi((t_bar - a_z)/sd_z)
at the best z.  If this oracle-accounted certificate fails at t, no refinement of the band accounting
alone can certify t."""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
from bab_core import NET, fwd, DIRS, cell_pcrown, T_BAR, N_RELU

v = DIRS["stripe_x1"]
rng = np.random.default_rng(3)
E = rng.normal(size=(40000, 36))
out = []
for t in (20.0, 25.0, 30.0, 35.0, 38.0):
    pres_all = [fwd(NET, t * v + E[i:i + 10000], keep=True)[1] for i in range(0, len(E), 10000)]
    best = None
    for z in (1.5, 2.0, 2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5):
        a, sd, bnds = cell_pcrown(t * v, np.zeros((36, 0)), np.zeros(0), np.zeros(0), z, want_bnds=True)
        viol = 0
        for pres in pres_all:
            bad = np.zeros(pres[0].shape[1], bool)
            for s, (l, u) in enumerate(bnds):
                rl = NET[s]["relu"]
                bad |= ((pres[s][rl] < l[rl, None]) | (pres[s][rl] > u[rl, None])).any(0)
            viol += bad.sum()
        pv = viol / len(E)
        tot = pv + (norm.cdf((T_BAR - a) / sd) if np.isfinite(a) else 1.0)
        print("t=%5.1f z=%.2f  true P(E^c)=%.4f  union=%.4f  a=%9.2f sd=%.2f  oracle-accounted bound=%.4f"
              % (t, z, pv, min(1, 2 * N_RELU * norm.sf(z)), a, sd, tot), flush=True)
        out.append(dict(t=t, z=z, p_viol=pv, union=2 * N_RELU * norm.sf(z), a=a, sd=sd, bound=tot))
json.dump(out, open("ACCOUNTING_CEILING.json", "w"), indent=1)
