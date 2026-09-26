"""01_null_quantiles.py -- null draws of the Gaussian map model and upper confidence bounds on the
critical values used by the bootstrap layer cake.  Writes null_sorted.npy and NULL_UCB.json.

For each level w in LEVELS the order statistic S_(k), k = binom.ppf(1 - a, N0, 1 - w) + 1, exceeds the
true (1 - w)-quantile of F_0 with probability >= 1 - a (proof in THEOREM_SMOOTH.md, Lemma 3), with
a = ETA_CRIT / #levels so that all levels hold together with probability >= 1 - ETA_CRIT.
Run: python -u 01_null_quantiles.py
"""
import os, json, time
import numpy as np
from smooth_core import forward, ucb_quantile, LEVELS, ETA_CRIT, T_BAR, T_HAT, HERE

N0 = 2_000_000
out_npy = os.path.join(HERE, "null_sorted.npy")
t0 = time.time()
if os.path.exists(out_npy):
    s0 = np.load(out_npy)
    print("loaded %d null scores" % len(s0))
else:
    rng = np.random.default_rng(2026_09_24_01)
    s0 = np.sort(forward(rng.normal(size=(N0, 36))))
    np.save(out_npy, s0)
    print("drew and sorted %d null scores in %.0f s" % (N0, time.time() - t0))
a = ETA_CRIT / len(LEVELS)
ucb = {str(w): ucb_quantile(s0, 1 - w, a) for w in LEVELS}
chk = ucb_quantile(s0, 0.95, 0.001)
print("MC .95 quantile %.4f (paper's MC %.4f); own 99.9%% UCB %.4f (paper's T_BAR %.4f)"
      % (np.quantile(s0, .95), T_HAT, chk, T_BAR))
for w in LEVELS:
    print("  level w=%.3f  UCB c_w = %.4f  (MC quantile %.4f)" % (w, ucb[str(w)], np.quantile(s0, 1 - w)))
json.dump(dict(N0=int(len(s0)), a_per_level=a, eta_crit=ETA_CRIT, levels=LEVELS, ucb=ucb,
               own_ucb_95=chk, t_bar_paper=T_BAR, mc_q95=float(np.quantile(s0, .95))),
          open(os.path.join(HERE, "NULL_UCB.json"), "w"), indent=1)
print("wrote NULL_UCB.json")
