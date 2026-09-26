"""ir_core.py -- local-power machinery of Theorem L (lp_core.py / 01_local_power.py) extended to the
niche design: two relevant covariates plus k irrelevant ones.

Design (theory/niche/run_niche.R):  x ~ N(0, I_p), p = 2 + k,
    eta_0 = -0.5 + 1.0 x1 + 0.8 x2 + c_irr * sum_{j>=3} x_j          (c_irr = 0.15 in the study, 0 = 'pure')
    local departure on the logit scale: eta_0 + tau g(x1, x2)/sqrt(n), g = x1 x2 ('inter') or x1^2 - 1 ('quad');
    h = tau v g with v = pi_0 (1 - pi_0).   The study uses C f with C = tau / sqrt(n).

Population integrals for DeepGOF-1, Hosmer-Lemeshow and Stukel are exact quadratures on a product grid in
probability space for (x1, x2) and Gauss-Hermite nodes for u = k^{-1/2} sum_{j>=3} x_j, the only function of
the irrelevant block that enters eta_0. The irrelevant block of the information matrix and of G, b follows
from x_j = u/sqrt(k) + w_j with w independent of (x1, x2, u), mean zero, Cov(w_j, w_l) = delta_jl - 1/k.

Liu's projection test is a Gaussian quadratic form with drift,
discretised at a reference sample of M covariate draws (Nystrom), ALL eigenvalues kept with their signs.
"""
import os, sys
NT = os.environ.get("IR_THREADS", "1")                         # BLAS / torch threads per process
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, NT)
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.special import ndtri, expit
from scipy.stats import chi2, binom, ncx2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from lp_core import score, check_bench, K, B, KREJ, ALPHA      # the shipped network, K = 6, B = 199, k = 10
import torch
torch.set_num_threads(int(NT))

BASE = np.array([-0.5, 1.0, 0.8])
C_IRR = 0.15
DEPS = ["inter", "quad"]

def gfun(name, x1, x2):
    return {"inter": x1 * x2, "quad": x1 ** 2 - 1}[name]

def eta0_of(x, cirr):
    """eta_0 at covariate rows x (n x p)."""
    e = BASE[0] + BASE[1] * x[:, 0] + BASE[2] * x[:, 1]
    if x.shape[1] > 2:
        e = e + cirr * x[:, 2:].sum(1)
    return e

def hexact(name, x, cirr, C, n):
    """The finite-n drift the study actually generates: sqrt(n) {pi(eta_0 + C g) - pi_0}."""
    e0 = eta0_of(x, cirr)
    return np.sqrt(n) * (expit(e0 + C * gfun(name, x[:, 0], x[:, 1])) - expit(e0))

# =====================================================================================================
# 1. DeepGOF-1 / HL / Stukel: exact population quantities on the (x1, x2, u) grid
# =====================================================================================================
class GridDesign:
    def __init__(self, k, cirr=C_IRR, G=1200, nGH=24):
        self.k, self.cirr, self.p = k, cirr, 2 + k
        u = (np.arange(G) + 0.5) / G
        xg = ndtri(u)
        bins = np.minimum((K * u).astype(int), K - 1)
        X1, X2 = np.meshgrid(xg, xg, indexing="ij")
        self.X1, self.X2 = X1.ravel(), X2.ravel()
        self.CELL = (bins[:, None] * K + bins[None, :]).ravel()
        self.N = len(self.X1)
        if k > 0:
            nodes, wts = hermegauss(nGH)
            self.U, self.WU = nodes, wts / wts.sum()
        else:
            self.U, self.WU = np.array([0.0]), np.array([1.0])
        self.var_eta = 1.0 + 0.64 + cirr ** 2 * k
        self._accumulate()

    def _eta(self, u):
        return BASE[0] + BASE[1] * self.X1 + BASE[2] * self.X2 + self.cirr * np.sqrt(self.k) * u

    def _accumulate(self):
        k, p, N = self.k, self.p, self.N
        E = lambda w: float(w.sum()) / N
        cs = lambda idx, w, nc: np.bincount(idx, weights=w, minlength=nc) / N
        # HL deciles of eta_0 ~ N(-0.5, var_eta)
        self.hlcut = -0.5 + np.sqrt(self.var_eta) * ndtri(np.arange(1, 10) / 10)
        # accumulators (relevant block explicit; irrelevant block through u)
        Ev = Ev_u = Ev_u2 = 0.0
        Evx = np.zeros(2); Evux = np.zeros(2); Evxx = np.zeros((2, 2))
        Dv = np.zeros(K * K); Gx = np.zeros((K * K, 3)); Gu = np.zeros(K * K)
        Dhl = np.zeros(10); Ghlx = np.zeros((10, 3)); Ghlu = np.zeros(10); Phl = np.zeros(10); Ppi = np.zeros(10)
        Ezz = np.zeros((2, 2)); Ezx = np.zeros((2, 3)); Ezu = np.zeros(2)
        self.nodes = []
        for u, wu in zip(self.U, self.WU):
            eta = self._eta(u); pi0 = expit(eta); v = pi0 * (1 - pi0)
            Xt = np.column_stack([np.ones(N), self.X1, self.X2])
            hl = np.searchsorted(self.hlcut, eta)
            Zs = np.column_stack([0.5 * eta ** 2 * (eta >= 0), -0.5 * eta ** 2 * (eta < 0)])
            self.nodes.append((u, wu, pi0, v, hl, Zs))
            Ev += wu * E(v); Ev_u += wu * E(v) * u; Ev_u2 += wu * E(v) * u * u
            Evx += wu * np.array([E(v * self.X1), E(v * self.X2)])
            Evux += wu * u * np.array([E(v * self.X1), E(v * self.X2)])
            Evxx += wu * (Xt[:, 1:] * v[:, None]).T @ Xt[:, 1:] / N
            Dv += wu * cs(self.CELL, v, K * K)
            Gx += wu * np.column_stack([cs(self.CELL, v * Xt[:, j], K * K) for j in range(3)])
            Gu += wu * u * cs(self.CELL, v, K * K)
            Dhl += wu * cs(hl, v, 10); Phl += wu * cs(hl, np.ones(N), 10); Ppi += wu * cs(hl, pi0, 10)
            Ghlx += wu * np.column_stack([cs(hl, v * Xt[:, j], 10) for j in range(3)])
            Ghlu += wu * u * cs(hl, v, 10)
            Ezz += wu * (Zs * v[:, None]).T @ Zs / N
            Ezx += wu * (Zs * v[:, None]).T @ Xt / N
            Ezu += wu * u * (Zs * v[:, None]).T @ np.ones(N) / N
        # information matrix, (p+1) x (p+1): order (1, x1, x2, x3..xp)
        I = np.zeros((p + 1, p + 1))
        I[0, 0] = Ev; I[0, 1:3] = Evx; I[1:3, 0] = Evx; I[1:3, 1:3] = Evxx
        if k > 0:
            I[0, 3:] = Ev_u / np.sqrt(k); I[3:, 0] = I[0, 3:]
            I[1:3, 3:] = (Evux / np.sqrt(k))[:, None]; I[3:, 1:3] = I[1:3, 3:].T
            I[3:, 3:] = Ev_u2 / k + Ev * (np.eye(k) - 1.0 / k)
        self.I, self.Iinv = I, np.linalg.inv(I)
        self.Dv = Dv
        self.Gm = np.column_stack([Gx] + ([np.tile(Gu / np.sqrt(k), (k, 1)).T] if k > 0 else []))
        self.Dm12 = 1 / np.sqrt(Dv)
        self.Sigma = np.eye(K * K) - (self.Dm12[:, None] * self.Gm) @ self.Iinv @ (self.Gm.T * self.Dm12[None, :])
        ev, evec = np.linalg.eigh(self.Sigma)
        self.Sroot = evec * np.sqrt(np.clip(ev, 0, None)); self.sig_ev = ev
        # HL
        self.Ghl = np.column_stack([Ghlx] + ([np.tile(Ghlu / np.sqrt(k), (k, 1)).T] if k > 0 else []))
        pbar = Ppi / Phl; self.whl = Phl * pbar * (1 - pbar)
        e2, v2 = np.linalg.eigh(np.diag(Dhl) - self.Ghl @ self.Iinv @ self.Ghl.T)
        self.Rhl = v2 * np.sqrt(np.clip(e2, 0, None))
        # Stukel
        self.Ezx = np.column_stack([Ezx] + ([np.tile(Ezu / np.sqrt(k), (k, 1)).T] if k > 0 else []))
        self.Oinv = np.linalg.inv(Ezz - self.Ezx @ self.Iinv @ self.Ezx.T)
        self.Ev = Ev

    def drift(self, hfun):
        """hfun(x1, x2, pi0, v, u) -> h on the grid (unit tau). Returns a (cells), a_hl, E[z h], b (p+1)."""
        k, N = self.k, self.N
        cs = lambda idx, w, nc: np.bincount(idx, weights=w, minlength=nc) / N
        a = np.zeros(K * K); ahl = np.zeros(10); zh = np.zeros(2); bx = np.zeros(3); bu = 0.0
        for (u, wu, pi0, v, hl, Zs) in self.nodes:
            h = hfun(self.X1, self.X2, pi0, v, u)
            a += wu * cs(self.CELL, h, K * K); ahl += wu * cs(hl, h, 10)
            zh += wu * Zs.T @ h / N
            bx += wu * np.array([h.sum(), (h * self.X1).sum(), (h * self.X2).sum()]) / N
            bu += wu * u * h.sum() / N
        b = np.concatenate([bx, np.full(k, bu / np.sqrt(k))]) if k > 0 else bx
        return a, ahl, zh, b

    def delta_dg(self, a, b):
        return self.Dm12 * (a - self.Gm @ self.Iinv @ b)

    def drift_hl(self, ahl, b):
        return ahl - self.Ghl @ self.Iinv @ b

    def ncp_stukel(self, zh, b):
        s = zh - self.Ezx @ self.Iinv @ b
        return float(s @ self.Oinv @ s)

def h_local(name):
    """Theorem-L drift at unit tau: h = v g."""
    return lambda x1, x2, pi0, v, u: v * gfun(name, x1, x2)

def h_finite(name, C, n, k, cirr):
    """The study's exact finite-n drift, sqrt(n) {pi(eta_0 + C g) - pi_0}, as a function on the grid
    (unit 'tau': multiply by nothing; pass tau = 1 to the power functions)."""
    def f(x1, x2, pi0, v, u):
        eta = BASE[0] + BASE[1] * x1 + BASE[2] * x2 + cirr * np.sqrt(k) * u
        return np.sqrt(n) * (expit(eta + C * gfun(name, x1, x2)) - expit(eta))
    return f

# ----- limit powers -----
def score_chunked(m, chunk=10000):
    """The shipped score in chunks (the whole 100k batch needs 4 GB of activations)."""
    return np.concatenate([score(m[s:s + chunk]) for s in range(0, len(m), chunk)])

class DGPower:
    """Limit power of the deployed DeepGOF-1 (B = 199, k = 10) at drift delta, given Sigma."""
    def __init__(self, design, NF=100000, seed=90210):
        rng = np.random.default_rng(seed)
        self.F0 = np.sort(score_chunked(rng.normal(size=(NF, K * K)) @ design.Sroot.T))
        self.ZF = rng.normal(size=(NF, K * K)) @ design.Sroot.T
        self.NF = NF

    def power(self, delta):
        W_ = 1 - np.searchsorted(self.F0, score_chunked(self.ZF + delta[None, :]), side="right") / self.NF
        return float(binom.cdf(KREJ - 1, B, W_).mean())

class HLPower:
    def __init__(self, design, NF=200000, seed=777):
        rng = np.random.default_rng(seed)
        self.d = design
        ZH = rng.normal(size=(NF, 10)) @ design.Rhl.T
        self.A = (ZH ** 2 / design.whl).sum(1); self.ZH = ZH
        self.crit = chi2.ppf(1 - ALPHA, 8)

    def power(self, tau, drift):
        Bq = (self.ZH * (drift / self.d.whl)[None, :]).sum(1)
        Cq = float((drift ** 2 / self.d.whl).sum())
        return float((self.A + 2 * tau * Bq + tau ** 2 * Cq > self.crit).mean())

def power_stukel(tau, ncp1):
    return float(ncx2.sf(chi2.ppf(1 - ALPHA, 2), 2, tau ** 2 * ncp1)) if tau > 0 else ALPHA

# =====================================================================================================
# 2. Quadratic-form tests at a reference sample: the projection test
# =====================================================================================================
def proj_kernel(x):
    """A of proj_verbatim.R (Liu et al. 2024 / Escanciano): A_ij = sum_r |x_i-x_r||x_j-x_r| (pi - angle_ijr),
    on standardised columns."""
    X = (x - x.mean(0)) / x.std(0, ddof=1); n = len(X)
    A = np.zeros((n, n))
    for r in range(n):
        D = X - X[r]; nm = np.sqrt((D * D).sum(1)); Gm = D @ D.T; NN = np.outer(nm, nm)
        cs = np.where(NN > 1e-12, Gm / np.where(NN > 1e-12, NN, 1), 0.0)
        A += NN * (np.pi - np.arccos(np.clip(cs, -1, 1)))
    return A

def project_kernel(R, X, mu2):
    """R.cor = (I-H)' R (I-H), H = V X (X'VX)^-1 X'; plus E.Q, Var.Q of the
    moment-matched reference (Bernoulli fourth moment included)."""
    A = mu2[:, None] * X
    Bm = np.linalg.inv(X.T @ A)
    RA = R @ A; XB, XBt = X @ Bm, X @ Bm.T
    Rc = R - XBt @ RA.T - RA @ XB.T + XBt @ (RA.T @ A) @ XB.T
    Rc = (Rc + Rc.T) / 2
    dR = np.diag(Rc)
    EQ = float((dR * mu2).sum())
    mu4 = mu2 * (1 - 3 * mu2)
    VarQ = float((dR ** 2 * (mu4 - 3 * mu2 ** 2)).sum() + 2 * mu2 @ (Rc * Rc) @ mu2)
    return Rc, EQ, VarQ

class QuadForm:
    """Limit law of a residual quadratic form r' R r at a reference sample: sum_j lam_j (Z_j + tau c_j)^2.
    The Gaussian draws Z are kept in float32 until drop(); prep(c) caches the two per-draw numbers that
    a drift direction c needs, so that the power at every tau is then free."""
    def __init__(self, R, X, v, NZ=20000, seed=1234, chunk=4000, offdiag=False):
        """offdiag=True: U-statistic discretisation -- the diagonal of the projected M x M matrix is set to zero.
        Needed for the projection test, whose weight function at x = x' (about 2 pi p) exceeds its typical
        off-diagonal value ((4 pi / 3) p) by O(p): with the diagonal kept, the diagonal's share of the drift and
        of sum lambda^2 is O(p/M) and O(p^2/M), a finite-M bias (it shifts only the trace in the limit)."""
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
        self.Z = np.concatenate([rng.normal(size=(min(chunk, NZ - s), len(self.lam))).astype(np.float32)
                                 for s in range(0, NZ, chunk)])
        self.A0 = (self.Z.astype(np.float64) ** 2 @ self.lam)            # null draws of the form
        self.F0 = np.sort(self.A0)
        self.df = 2 * self.EQ ** 2 / self.VarQ
        self.crit_mm = chi2.ppf(1 - ALPHA, self.df) * self.VarQ / (2 * self.EQ)
        self._cache = {}

    def drop(self):
        self.Z = None

    def coef(self, hvec):
        """c per unit tau from a drift h at the sample points: d = h / sqrt(M), c = U' (d / sqrt(v))."""
        return self.U.T @ (hvec / np.sqrt(self.v) / np.sqrt(self.M))

    def prep(self, hvec, key):
        """Cache (c, Bq, signal) for the drift direction hvec under a name; needs Z."""
        c = self.coef(hvec)
        Bq = self.Z @ (self.lam * c).astype(np.float32)
        self._cache[key] = (c, Bq.astype(np.float64), float((self.lam * c ** 2).sum()))
        return key

    def signal(self, key):
        return self._cache[key][2]                                      # drift per tau^2

    def cross_sd(self, key):
        c = self._cache[key][0]
        return float(2 * np.sqrt((self.lam ** 2 * c ** 2).sum()))

    def noise(self):
        return float(np.sqrt(2 * self.sum_lam2))

    def draws(self, tau, key):
        c, Bq, sig = self._cache[key]
        return self.A0 + 2 * tau * Bq + tau ** 2 * sig

    def power_momentmatched(self, tau, key):
        """Moment-matched rule: Q * 2EQ/VarQ referred to chi2 with 2EQ^2/VarQ df (no arm of this archive uses it)."""
        return float((self.draws(tau, key) > self.crit_mm).mean())

    def power_bootstrap(self, tau, key, Bboot=B, kk=KREJ):
        """Deployed proj rule: parametric bootstrap with B = 199, p = (1 + #{T* >= T})/(B + 1), at its limit
        (Dufour's formula applied to the null law of the quadratic form)."""
        W_ = 1 - np.searchsorted(self.F0, self.draws(tau, key), side="right") / len(self.F0)
        return float(binom.cdf(kk - 1, Bboot, W_).mean())

def ref_sample(k, M, seed, cirr=C_IRR):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(M, 2 + k)); X = np.column_stack([np.ones(M), x])
    pi0 = expit(eta0_of(x, cirr)); v = pi0 * (1 - pi0)
    return x, X, pi0, v

def tau_for(fn, target=0.5, iters=18):
    hi = 1.0
    while fn(hi) < target:
        hi *= 2
        if hi > 65536: return float("inf")
    lo = 0.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        if fn(mid) < target: lo = mid
        else: hi = mid
    return hi
