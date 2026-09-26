"""sst_core.py -- the local limit of Stukel's test in the SstBoth form of Liu et al. (2024).

SstBoth adds the two Stukel covariates z_a = eta^2/2 1{eta >= 0}, z_b = -eta^2/2 1{eta < 0} ONE AT A TIME:
each has a score z-statistic Z_j = u_j / sqrt(Omega_jj), with u = n^{-1/2} z'(y - pi_hat) and Omega the
efficient score covariance after the fit. The statistic is T = Z_a^2 + Z_b^2, referred to chi-square(2)
(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal")), statmod::glm.scoretest).

Under pi_n = pi_0 + tau h / sqrt(n), u -> N(tau s(h), Omega) with s(h) = E[z h] - E[v z x'] I^-1 b(h)
(the same Omega and s(h) as the joint form). So W = D^-1/2 u, D = diag(Omega), has law
N(tau D^-1/2 s(h), R), R the correlation matrix of Omega, and T -> ||W||^2: a generalised chi-square,
sum_i lam_i (U_i + c_i)^2 with lam = 1 -/+ |rho| the eigenvalues of R and c = Lam^-1/2 Q' mu. Unless rho = 0
this is NOT chi-square(2) under the null, so the limit size of the deployed rule differs from .05.

The rejection probability P(T > q) is computed by one-dimensional quadrature over U_1 (conditional on U_1,
the second term is a scaled noncentral chi-square(1), whose cdf is a difference of two normal cdfs).
"""
import numpy as np
from scipy.stats import norm, chi2
from scipy.integrate import quad

ALPHA = 0.05
Q95 = float(chi2.ppf(1 - ALPHA, 2))


def gen_chisq_sf(lam, c, q=Q95):
    """P(lam1 (U1 + c1)^2 + lam2 (U2 + c2)^2 > q), U ~ N(0, I_2), lam > 0."""
    l1, l2 = float(lam[0]), float(lam[1]); c1, c2 = float(c[0]), float(c[1])
    if l1 < l2:
        l1, l2, c1, c2 = l2, l1, c2, c1                          # integrate over the larger weight
    r = np.sqrt(q / l1)
    def inner(u):                                               # P(term2 <= q - l1 (u + c1)^2)
        a = q - l1 * (u + c1) ** 2
        if a <= 0:
            return 0.0
        t = np.sqrt(a / l2)
        return (norm.cdf(t - c2) - norm.cdf(-t - c2)) * norm.pdf(u)
    lo, hi = -c1 - r, -c1 + r
    acc, _ = quad(inner, lo, hi, limit=200, epsabs=1e-12, epsrel=1e-10)
    return float(1.0 - acc)


class SstBoth:
    """Limit law of SstBoth given Omega (2 x 2 efficient score covariance)."""
    def __init__(self, Omega):
        self.Omega = np.asarray(Omega, float)
        self.dm12 = 1 / np.sqrt(np.diag(self.Omega))
        self.R = self.dm12[:, None] * self.Omega * self.dm12[None, :]
        self.rho = float(self.R[0, 1])
        lam, Q = np.linalg.eigh(self.R)
        self.lam, self.Q = lam, Q

    def power(self, tau, s):
        mu = tau * self.dm12 * np.asarray(s, float)
        c = (self.Q.T @ mu) / np.sqrt(self.lam)
        return gen_chisq_sf(self.lam, c)

    def size(self):
        return self.power(0.0, np.zeros(2))


def tau_for(fn, target=0.5, iters=40):
    hi = 1.0
    while fn(hi) < target:
        hi *= 2
        if hi > 65536:
            return float("inf")
    lo = 0.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        if fn(mid) < target:
            lo = mid
        else:
            hi = mid
    return hi


def sstboth_stat(X, y, b):
    """Finite-sample SstBoth from a logistic fit (X with intercept, beta b), as statmod::glm.scoretest:
    Z_j = u_j / sqrt(I_jj.x), u = z'(y - p), I_jj.x = z_j' V z_j - z_j' V X (X'VX)^-1 X' V z_j.
    Returns (T, p); NaN if one column is identically zero (no fitted risk on that side of .5)."""
    eh = X @ b; ph = 1 / (1 + np.exp(-eh)); vh = ph * (1 - ph); r = y - ph
    Z = np.column_stack([0.5 * eh ** 2 * (eh >= 0), -0.5 * eh ** 2 * (eh < 0)])
    if not ((Z[:, 0] != 0).any() and (Z[:, 1] != 0).any()):
        return float("nan"), float("nan")
    u = Z.T @ r
    A = Z.T @ (vh[:, None] * X)
    Cv = Z.T @ (vh[:, None] * Z) - A @ np.linalg.solve(X.T @ (vh[:, None] * X), A.T)
    T = float((u ** 2 / np.diag(Cv)).sum())
    return T, float(chi2.sf(T, 2))
