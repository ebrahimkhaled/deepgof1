"""08_proj_liu_raw.py -- the projection-test limits recomputed with the CORRECTED weight of Liu et al. (2024),
after Escanciano (2006), on the RAW covariates (theory/addon_arms/proj_liu.R, proj_liu_c.c):

    T = e'Ae / n^2,  A_ij = sum_l A_ijl,  A_ijl = (pi - angle(x_i - x_l, x_j - x_l)) / (2 pi),
    A_ijl = 1/2 if exactly one of the two differences is zero, 1 if both are; i = j: 1/2 (angle 0).

No length factor, no standardisation. Everything else is as in 06_proj_offdiag.py: Gaussian quadratic form with
drift at reference samples of M covariate draws (same seeds as 01/06), the diagonal of the projected M x M matrix
removed (U-statistic discretisation; the weight is 1/2 on the diagonal against about 1/3 off it, an O(1/M) share of
the drift but an O(1) share of sum lambda^2 at M = 1000), the deployed B = 199 rule at its limit (Dufour).
The with-diagonal value is recomputed alongside as a control.

Tasks (argument list):
  curve:k:seed    study design (c_irr = 0.15), k irrelevant covariates, M = 1000, seeds 1-6
  M2000:k:seed    the same at M = 2000 (discretisation check)
  pred:seed       k = 8 study design: limit power at the niche cells' exact finite-n drift, M = 1000
  are:seed        the design of Table tab:are (FORMULA.json: p = 2, eta_0 = -0.5 + x1 + 0.6 x2, five departures)
  areM2000:seed   the same at M = 2000
  rate:p:seed     pure design (c_irr = 0), p = 2 + k: nu, p nu, q_k, largest share (proposition check)
One record per task, appended to PROJLIU_<tag>.jsonl (tag = env PL_TAG); all PROJLIU_*.jsonl are read on start,
so the run resumes. Assemble with:  python 08_proj_liu_raw.py assemble   (-> PROJ_LIU_RAW.json)
"""
import os, sys, json, time, glob, types
import numpy as np

# The projection arm needs no network. lp_core loads torch and the shipped weights at import, which commits
# 1-2 GB per process on Windows; a light stand-in with the same constants and departures is registered instead
# (values copied from lp_core.py: K = 6, B = 199, k = 10, alpha = .05, BETA0, the five bounded departures).
_lp = types.ModuleType("lp_core")
_lp.K, _lp.B, _lp.KREJ, _lp.ALPHA = 6, 199, 10, 0.05
_lp.BETA0 = np.array([-0.5, 1.0, 0.6])
_lp.DEPS = ["quadratic", "interaction", "cubic", "bump", "threshold"]
def _g5(name, x1, x2):
    g = {"quadratic": x1 ** 2 - 1, "interaction": x1 * x2, "cubic": x1 ** 3 - 3 * x1,
         "bump": 3 * np.exp(-((x1 - 1) ** 2 + (x2 + 0.5) ** 2) / 0.5), "threshold": 2.0 * (x2 > 1)}[name]
    return np.clip(g, -6, 6)
_lp.gfun = _g5
_lp.score = None
_lp.check_bench = lambda: 0.0
sys.modules.setdefault("lp_core", _lp)
_torch = types.ModuleType("torch"); _torch.set_num_threads = lambda n: None
sys.modules.setdefault("torch", _torch)
from ir_core import (QuadForm, ref_sample, tau_for, hexact, gfun, C_IRR, DEPS, expit, check_bench)
import lp_core

HERE = os.path.dirname(os.path.abspath(__file__))
TAG = os.environ.get("PL_TAG", "0")
RAW = os.path.join(HERE, "PROJLIU_%s.jsonl" % TAG)
CELLS = [(f, C, n) for f in DEPS for C in (0.3, 0.5, 0.8) for n in (200, 500)]
TWO_PI = 2 * np.pi


def liu_kernel(x, tol=1e-10):
    """A of proj_A_liu(X, std = FALSE) in proj_liu.R: raw covariates, angle-only weight, zero rule."""
    x = np.asarray(x, float); n = len(x)
    eps = tol * max(1.0, float(np.sqrt((x * x).sum(1).max())))
    A = np.zeros((n, n))
    for l in range(n):
        D = x - x[l]; nm = np.sqrt((D * D).sum(1)); z = nm <= eps
        U = D / np.where(z, 1.0, nm)[:, None]
        Kl = (np.pi - np.arccos(np.clip(U @ U.T, -1.0, 1.0))) / TWO_PI
        np.fill_diagonal(Kl, 0.5)
        if z.any():
            Kl[z, :] = 0.5; Kl[:, z] = 0.5; Kl[np.ix_(z, z)] = 1.0
        A += Kl
    return A


def liu_kernel_base(x, zb):
    """Sum over external base points zb (independent draws of the covariate law) of the same angle weight: a Monte
    Carlo of n K(x_i, x_j) whose base-point noise is independent of the reference sample (check task 'base')."""
    x = np.asarray(x, float); A = np.zeros((len(x), len(x)))
    for l in range(len(zb)):
        D = x - zb[l]; nm = np.sqrt((D * D).sum(1))
        U = D / nm[:, None]
        Kl = (np.pi - np.arccos(np.clip(U @ U.T, -1.0, 1.0))) / TWO_PI
        np.fill_diagonal(Kl, 0.5)
        A += Kl
    return A * (len(x) / len(zb))


def load():
    out = {}
    for f in glob.glob(os.path.join(HERE, "PROJLIU_*.jsonl")):
        for line in open(f):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            out[r["key"]] = r
    return out


def log(rec):
    with open(RAW, "a") as fh:
        fh.write(json.dumps(rec) + "\n"); fh.flush(); os.fsync(fh.fileno())


def qstats(q, M):
    """Scale-free and operator-scale summaries of the null law; operator scale = matrix / M^2, reported in units of
    F0 = pi - arccos (i.e. times 2 pi) so that they compare with PRINTED_WEIGHT_CHECK.log."""
    q95 = float((np.quantile(q.A0, 0.95) - q.sum_lam) / q.noise())
    return dict(noise_raw=q.noise(), nu_F0=q.noise() / M ** 2 * TWO_PI, q95=q95, max_share=q.max_share,
                n_neg=q.n_neg, n_keep=q.n_keep)


def arm_record(q, key, M):
    pw = lambda t: q.power_bootstrap(t, key)
    t50 = tau_for(pw)
    sig, nz = q.signal(key), q.noise()
    return dict(tau50=t50, size=pw(0.0), signal_raw=sig, sigma_F0=sig / M ** 2 * TWO_PI, snr=sig / nz,
                theta_star=t50 ** 2 * sig / nz)


class LightQuadForm(QuadForm):
    """ir_core.QuadForm with the null draws A0 = (Z^2) lam accumulated chunk by chunk, so that the 20000 x M float64
    copy of Z is never formed (at M = 2000 that copy is 320 MB, which the loaded machine could not commit). Same
    random stream (same chunks, same float32 rounding), same formulas; used for the M = 2000 checks only."""
    def __init__(self, R, X, v, NZ=20000, seed=1234, chunk=4000, offdiag=False):
        from ir_core import project_kernel, chi2, ALPHA
        self.M = len(v); self.v = v
        Rc, self.EQ, self.VarQ = project_kernel(R, X, v)
        Bm = np.sqrt(v)[:, None] * Rc * np.sqrt(v)[None, :]
        del Rc
        if offdiag:
            Bm[np.diag_indices(len(v))] = 0.0
        lam, U = np.linalg.eigh(Bm)
        del Bm
        keep = np.abs(lam) > 1e-12 * np.abs(lam).max()
        self.lam, self.U = lam[keep], U[:, keep]
        self.n_neg = int((self.lam < 0).sum()); self.n_keep = int(keep.sum())
        self.sum_lam = float(self.lam.sum()); self.sum_lam2 = float((self.lam ** 2).sum())
        self.max_share = float((self.lam ** 2).max() / self.sum_lam2)
        rng = np.random.default_rng(seed)
        self.Z = np.empty((NZ, len(self.lam)), dtype=np.float32); self.A0 = np.empty(NZ)
        for s in range(0, NZ, chunk):
            z = rng.normal(size=(min(chunk, NZ - s), len(self.lam))).astype(np.float32)
            self.Z[s:s + len(z)] = z
            self.A0[s:s + len(z)] = z.astype(np.float64) ** 2 @ self.lam
            del z
        self.F0 = np.sort(self.A0)
        self.df = 2 * self.EQ ** 2 / self.VarQ
        self.crit_mm = chi2.ppf(1 - ALPHA, self.df) * self.VarQ / (2 * self.EQ)
        self._cache = {}


def build(x, X, v, M, seed, offdiag, R=None):
    QF = LightQuadForm if M > 1000 else QuadForm
    return QF(liu_kernel(x) if R is None else R, X, v, seed=seed + 100, offdiag=offdiag)


def run_task(task):
    parts = task.split(":"); kind = parts[0]; t0 = time.time()
    if kind in ("curve", "M2000", "base"):
        # base:L:k:seed -- the curve's reference sample (M = 1000) with L external base points instead of the
        # sample's own points (base-point Monte Carlo check); 'diag' mode kept for a uniform record
        if kind == "base":
            L, k, sd = int(parts[1]), int(parts[2]), int(parts[3]); M = 1000
        else:
            k, sd = int(parts[1]), int(parts[2]); M = 2000 if kind == "M2000" else 1000
        x, X, pi0, v = ref_sample(k, M, 2000 * sd + k, cirr=C_IRR)
        if kind == "base":
            zb = np.random.default_rng(777000 + 100 * sd + k).normal(size=(L, 2 + k))
            R = liu_kernel_base(x, zb)
        else:
            R = liu_kernel(x)
        rec = dict(key=task, kind=kind, k=k, p=2 + k, seed=sd, M=M, L=(L if kind == "base" else M))
        for mode in ("offdiag", "diag"):
            q = build(x, X, v, M, sd, mode == "offdiag", R)
            for dep in DEPS:
                q.prep(v * gfun(dep, x[:, 0], x[:, 1]), dep)
            q.drop()
            rec[mode] = dict(null=qstats(q, M), **{dep: arm_record(q, dep, M) for dep in DEPS})
            del q
    elif kind == "pred":
        sd = int(parts[1]); k, M = 8, 1000
        x, X, pi0, v = ref_sample(k, M, 2000 * sd + k, cirr=C_IRR)
        q = build(x, X, v, M, sd, True)
        keys = []
        for (f, C, n) in CELLS:
            kk = "%s_C%.1f_n%d" % (f, C, n); q.prep(hexact(f, x, C_IRR, C, n), kk); keys.append(kk)
        q.prep(np.zeros(M), "null"); q.drop()
        rec = dict(key=task, kind=kind, seed=sd, M=M, cells={kk: q.power_bootstrap(1.0, kk) for kk in keys},
                   size=q.power_bootstrap(0.0, "null"))
    elif kind == "nsur":
        # nsur:n:seed -- finite-n surrogate: a reference sample of the study's own size n, its own points as base
        # points (as the statistic does), Gaussian errors; with the diagonal removed and kept
        n, sd = int(parts[1]), int(parts[2]); k = 8
        x, X, pi0, v = ref_sample(k, n, 9000 + 10 * n + sd, cirr=C_IRR)
        R = liu_kernel(x); rec = dict(key=task, kind=kind, n=n, seed=sd, M=n)
        for mode in ("offdiag", "diag"):
            q = build(x, X, v, n, sd, mode == "offdiag", R)
            keys = []
            for (f, C, nn) in CELLS:
                if nn == n:
                    kk = "%s_C%.1f_n%d" % (f, C, n); q.prep(hexact(f, x, C_IRR, C, n), kk); keys.append(kk)
            q.prep(np.zeros(n), "null"); q.drop()
            rec[mode] = dict(cells={kk: q.power_bootstrap(1.0, kk) for kk in keys}, size=q.power_bootstrap(0.0, "null"))
            del q
    elif kind in ("are", "areM2000"):
        sd = int(parts[1]); M = 2000 if kind == "areM2000" else 1000
        rng = np.random.default_rng(3000 * sd + 17)
        x = rng.normal(size=(M, 2)); X = np.column_stack([np.ones(M), x])
        pi0 = expit(X @ lp_core.BETA0); v = pi0 * (1 - pi0)
        R = liu_kernel(x); rec = dict(key=task, kind=kind, seed=sd, M=M)
        for mode in ("offdiag", "diag"):
            q = build(x, X, v, M, sd, mode == "offdiag", R)
            for dep in lp_core.DEPS:
                q.prep(v * lp_core.gfun(dep, x[:, 0], x[:, 1]), dep)
            q.drop()
            rec[mode] = dict(null=qstats(q, M), **{dep: arm_record(q, dep, M) for dep in lp_core.DEPS})
            del q
    elif kind == "rate":
        p, sd = int(parts[1]), int(parts[2]); M = 1000
        x, X, pi0, v = ref_sample(p - 2, M, 5000 + 10 * sd + p, cirr=0.0)
        q = build(x, X, v, M, sd, True)
        for dep in DEPS:
            q.prep(v * gfun(dep, x[:, 0], x[:, 1]), dep)
        q.drop()
        st = qstats(q, M)
        rec = dict(key=task, kind=kind, p=p, seed=sd, M=M, null=st, p_nu_F0=p * st["nu_F0"],
                   **{dep: dict(p_e=p * q.signal(dep) / q.noise()) for dep in DEPS})
    else:
        raise ValueError(task)
    rec["secs"] = time.time() - t0
    return rec


def assemble():
    D = load(); A = json.load(open(os.path.join(HERE, "ARE_IRRELEVANT.json")))
    F = json.load(open(os.path.join(os.path.dirname(HERE), "FORMULA.json")))
    se = lambda a: float(np.std(a, ddof=1) / np.sqrt(len(a))) if len(a) > 1 else float("nan")

    def rat_tau(ts, tdg):
        # headline convention of 06_proj_offdiag.py: (mean tau50 / tau50_DG)^2, SE by the delta method
        ts = np.asarray(ts, float); val = float((ts.mean() / tdg) ** 2)
        return val, (float(2 * val * ts.std(ddof=1) / np.sqrt(len(ts)) / ts.mean()) if len(ts) > 1 else float("nan"))
    out = {"weight": "A_ijl = (pi - angle(x_i - x_l, x_j - x_l)) / (2 pi), 1/2 if one difference is zero, 1 if both; "
                     "raw covariates (theory/addon_arms/proj_liu.R proj_A_liu(std = FALSE))",
           "method": "Gaussian quadratic form at reference samples, diagonal of the projected matrix removed; "
                     "B = 199 rule at its limit; ratio = (tau50 proj / tau50 DG)^2 averaged over reference samples, "
                     "SE = sd / sqrt(#samples)", "curve": {}, "M2000": {}, "pred": {}, "are": {}, "rate": {}}
    for dep in DEPS:
        out["curve"][dep] = {}
        for k in (0, 1, 2, 4, 6, 8):
            recs = [r for r in D.values() if r["kind"] == "curve" and r["k"] == k]
            if not recs:
                continue
            tdg = A["variants"]["study"][dep][str(k)]["tau50"]["DG"]
            rat = [(r["offdiag"][dep]["tau50"] / tdg) ** 2 for r in recs]
            ratd = [(r["diag"][dep]["tau50"] / tdg) ** 2 for r in recs]
            pe = [(2 + k) * r["offdiag"][dep]["snr"] for r in recs]
            out["curve"][dep][str(k)] = dict(
                n=len(recs), tau50_DG=tdg, tau50=[r["offdiag"][dep]["tau50"] for r in recs],
                ratio=float(np.mean(rat)), ratio_se=se(rat), ratio_range=[float(min(rat)), float(max(rat))],
                ratio_of_mean_tau=rat_tau([r["offdiag"][dep]["tau50"] for r in recs], tdg)[0],
                ratio_of_mean_tau_se=rat_tau([r["offdiag"][dep]["tau50"] for r in recs], tdg)[1],
                ratio_diag=float(np.mean(ratd)), ratio_diag_se=se(ratd),
                p_e=float(np.mean(pe)), p_e_se=se(pe),
                q95=[r["offdiag"]["null"]["q95"] for r in recs],
                theta_star=[r["offdiag"][dep]["theta_star"] for r in recs],
                size=float(np.mean([r["offdiag"][dep]["size"] for r in recs])),
                ratio_HL=A["variants"]["study"][dep][str(k)]["ratio_vs_DG"]["HL"],
                ratio_Stukel=A["variants"]["study"][dep][str(k)]["ratio_vs_DG"]["Stukel"])
    for r in D.values():
        if r["kind"] == "M2000":
            dd = {}
            for dep in DEPS:
                tdg = A["variants"]["study"][dep][str(r["k"])]["tau50"]["DG"]
                dd[dep] = dict(ratio=(r["offdiag"][dep]["tau50"] / tdg) ** 2, ratio_diag=(r["diag"][dep]["tau50"] / tdg) ** 2)
            out["M2000"]["k%d_seed%d" % (r["k"], r["seed"])] = dd
    out["base"] = {}
    for r in D.values():
        if r["kind"] == "base":
            own = D.get("curve:%d:%d" % (r["k"], r["seed"]))
            dd = {}
            for dep in DEPS:
                tdg = A["variants"]["study"][dep][str(r["k"])]["tau50"]["DG"]
                dd[dep] = dict(ratio=(r["offdiag"][dep]["tau50"] / tdg) ** 2,
                               ratio_own=(own["offdiag"][dep]["tau50"] / tdg) ** 2 if own else None,
                               sigma=r["offdiag"][dep]["sigma_F0"], sigma_own=own["offdiag"][dep]["sigma_F0"] if own else None)
            dd["nu"] = r["offdiag"]["null"]["nu_F0"]; dd["nu_own"] = own["offdiag"]["null"]["nu_F0"] if own else None
            out["base"]["L%d_k%d_seed%d" % (r["L"], r["k"], r["seed"])] = dd
    out["nsur"] = {}
    for n in (200, 500):
        recs = [r for r in D.values() if r["kind"] == "nsur" and r["n"] == n]
        if not recs:
            continue
        for mode in ("offdiag", "diag"):
            for kk in recs[0][mode]["cells"]:
                vals = [r[mode]["cells"][kk] for r in recs]
                out["nsur"].setdefault(kk, {})[mode] = dict(mean=float(np.mean(vals)), se=se(vals), n=len(vals))
            sz = [r[mode]["size"] for r in recs]
            out["nsur"].setdefault("null_n%d" % n, {})[mode] = dict(mean=float(np.mean(sz)), se=se(sz), n=len(sz))
    preds = [r for r in D.values() if r["kind"] == "pred"]
    if preds:
        for kk in preds[0]["cells"]:
            vals = [r["cells"][kk] for r in preds]
            out["pred"][kk] = dict(mean=float(np.mean(vals)), se=se(vals), n=len(vals))
        out["pred"]["size"] = dict(mean=float(np.mean([r["size"] for r in preds])), n=len(preds))
    for kind in ("are", "areM2000"):
        recs = [r for r in D.values() if r["kind"] == kind]
        if not recs:
            continue
        o = {}
        for dep in lp_core.DEPS:
            tdg = F["dep"][dep]["tau50"]["DG"]
            rat = [(r["offdiag"][dep]["tau50"] / tdg) ** 2 for r in recs]
            ratd = [(r["diag"][dep]["tau50"] / tdg) ** 2 for r in recs]
            o[dep] = dict(n=len(recs), tau50_DG=tdg, tau50=[r["offdiag"][dep]["tau50"] for r in recs],
                          ratio_of_mean_tau=rat_tau([r["offdiag"][dep]["tau50"] for r in recs], tdg)[0],
                          ratio_of_mean_tau_se=rat_tau([r["offdiag"][dep]["tau50"] for r in recs], tdg)[1],
                          ratio=float(np.mean(rat)), ratio_se=se(rat), ratio_range=[float(min(rat)), float(max(rat))],
                          ratio_diag=float(np.mean(ratd)), size=float(np.mean([r["offdiag"][dep]["size"] for r in recs])),
                          ratio_HL=F["dep"][dep]["ARE_vs_DG"]["HL"], ratio_Stukel=F["dep"][dep]["ARE_vs_DG"]["Stukel"])
        out["are"][kind] = o
    for r in D.values():
        if r["kind"] == "rate":
            out["rate"]["p%d_seed%d" % (r["p"], r["seed"])] = dict(p_nu_F0=r["p_nu_F0"], q95=r["null"]["q95"],
                                                                    max_share=r["null"]["max_share"],
                                                                    p_e={d: r[d]["p_e"] for d in DEPS})
    json.dump(out, open(os.path.join(HERE, "PROJ_LIU_RAW.json"), "w"), indent=1)
    print(json.dumps(out, indent=1)[:200000])


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "assemble":
        assemble(); sys.exit(0)
    done = load()
    print("%d records on disk | tasks %s" % (len(done), args)); sys.stdout.flush()
    for task in args:
        if task in done:
            continue
        rec = run_task(task); log(rec); done[task] = rec
        print("%s done in %.0f s" % (task, rec["secs"])); sys.stdout.flush()
