"""smooth_core.py -- shared pieces of the Gaussian-smoothing power certificate (THEOREM_SMOOTH.md).

Model (the paper's, Online Resource 1 section 'A power certificate under a Gaussian model of the map'):
map m = t v + e, e ~ N(0, I_36), v a unit direction, S the frozen network (crown_core.S_ref), fixed
threshold T_BAR = 2.4865 (a 99.9% upper confidence bound on the .95 null quantile, from 400,000 null
draws in 03_certify.py).  Signal to sample size: n = 4 t^2.

The certificate rests on the Gaussian-shift lemma: q_c(mu) = P{S(mu + e) > c} has Phi^{-1}(q_c) 1-Lipschitz
in mu.  Along the ray mu = t v this gives, for net points t_j with Clopper-Pearson lower bounds q_lo,j,
    Phi^{-1} q_c(t) >= max_j [ Phi^{-1} q_lo,j - |t - t_j| ]      for every t,
on the event (probability >= 1 - eta) that every lower bound holds.

Two test versions are certified:
  FIXED     reject iff S > T_BAR                                   (the paper's table)
  B199      the deployed rank test, B = 199, alpha = .05 (k = 10): power = E[Bin(9; 199, W)],
            W = 1 - F_0(S(tv + e)), bounded from below by a layer cake over levels w_j with
            critical values c_j = F_0^{-1}(1 - w_j), each replaced by an upper confidence bound.
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
import torch
torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
from scipy.stats import norm, beta, binom, chi2
import contextlib, io

HERE = os.path.dirname(os.path.abspath(__file__))
FP = os.path.dirname(HERE)
sys.path.insert(0, FP)
with contextlib.redirect_stdout(io.StringIO()):
    from crown_core import NET, S_ref, crown, K

CERT = json.load(open(os.path.join(FP, "CERTIFIED.json")))
PCERT = json.load(open(os.path.join(FP, "PCERT.json")))
T_BAR = float(CERT["t_bar"])                 # 2.4865397 (99.9% UCB of the .95 null quantile)
T_HAT = float(CERT["t_hat"])                 # 2.4615683 (the MC estimate the paper's t_emp used)
C_CERT, L_CERT = float(CERT["c_cert"]), float(CERT["L_cert"])
KAPPA = {r["direction"]: float(r["kappa"]) for r in CERT["rows"]}
T_EMP_PAPER = {r["direction"]: float(r["t_emp"]) for r in CERT["rows"]}
RHO = 0.5                                     # ||mu|| per sqrt(n): n = (t / RHO)^2 = 4 t^2

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}
DIRS = {k: v / np.linalg.norm(v) for k, v in DIRS.items()}
ORDER = ["step_x1", "diffuse", "stripe_x1", "saddle", "random"]

# ---- the bootstrap layer (Theorem L item 4, B = 199, alpha = .05 -> k = 10, reject iff N <= 9) ------
B_BOOT, K_REJ = 199, 10
LEVELS = [.05, .045, .04, .035, .03, .025, .02, .015, .01, .007, .005, .003, .002, .001]
def g_bin(w):
    """P(Bin(199, w) <= 9): the rejection probability given W = w (decreasing in w)."""
    return binom.cdf(K_REJ - 1, B_BOOT, w)
LAYER_W = np.array([g_bin(LEVELS[0])] + [g_bin(LEVELS[j]) - g_bin(LEVELS[j - 1]) for j in range(1, len(LEVELS))])
assert (LAYER_W > 0).all() and LAYER_W.sum() <= 1.0 + 1e-12

# ---- probability budgets --------------------------------------------------------------------------
ETA_NET = 0.001          # Monte Carlo at the net points (Bonferroni over the net, per direction)
ETA_CRIT = 0.001         # the critical values: T_BAR carries .001 from the paper's null draws; the
                         # layer-cake levels share .001 by Bonferroni (01_null_quantiles.py)
EPS_TAIL = 1e-4          # the deterministic tail: P(||e|| > R_TAIL) = EPS_TAIL
R_TAIL = float(np.sqrt(chi2.ppf(1 - EPS_TAIL, 36)))

def forward(m, batch=25000):
    """S at raw maps m (rows), float64, the shipped network."""
    out = []
    with torch.no_grad():
        for a in range(0, len(m), batch):
            out.append(S_ref(torch.tensor(m[a:a + batch])).numpy())
    return np.concatenate(out)

def cp_lower(x, m, a):
    """Clopper-Pearson lower confidence bound at level 1 - a for x successes in m."""
    return 0.0 if x == 0 else float(beta.ppf(a, x, m - x + 1))

def cp_upper(x, m, a):
    return 1.0 if x == m else float(beta.ppf(1 - a, x + 1, m - x))

def ucb_quantile(sorted_null, p, a):
    """Order statistic that exceeds the p-quantile of the null law with probability >= 1 - a."""
    n0 = len(sorted_null)
    k = int(binom.ppf(1 - a, n0, p)) + 1
    return float(sorted_null[min(k, n0) - 1])

def tail_t_global(c_max, kappa):
    """Beyond this t, S(tv+e) >= t kappa - L ||e|| - c > c_max whenever ||e|| <= R_TAIL (Lemma 2 + the
    certified constants), so P(S <= c_max) <= EPS_TAIL."""
    return (c_max + C_CERT + L_CERT * R_TAIL) / kappa

def append_jsonl(path, rec):
    with open(path, "a") as fh:
        fh.write(json.dumps(rec) + "\n")

def read_jsonl(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass                      # a partial last line after a crash
    return out

def rscript_count():
    try:
        import subprocess
        s = subprocess.run(["tasklist"], capture_output=True, text=True).stdout.lower()
        return s.count("rscript")
    except Exception:
        return 0
