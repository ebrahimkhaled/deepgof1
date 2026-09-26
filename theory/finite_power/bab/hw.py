"""hw.py -- Hunter-Worsley (Kounias) bound on P(E^c) for the probabilistic CROWN band event.

Each one-sided band statement fails on a halfspace A_i = {w_i'e > z ||w_i||}, P(A_i) = Phi(-z).
For ANY spanning tree T on the events (Hunter 1976; Worsley 1982),
    P(union A_i) <= sum_i P(A_i) - sum_{(i,j) in T} P(A_i cap A_j),
and P(A_i cap A_j) = P(X > z, Y > z), corr(X, Y) = rho_ij, is an exact bivariate normal orthant
probability:  Phi2(-z, -z; rho) = Phi(-z) - 2 T(-z, sqrt((1 - rho)/(1 + rho)))  (Owen's T).
We take the maximum-weight spanning tree (Prim, dense).  The bound is sound for whatever tree is used;
the pairwise probabilities are computed to ~1e-16 and then rounded DOWN by 1e-12 each for safety."""
import numpy as np
from scipy.stats import norm
from scipy.special import owens_t


def orthant_pair(z, rho):
    rho = np.clip(rho, -1 + 1e-12, 1 - 1e-12)
    a = np.sqrt((1 - rho) / (1 + rho))
    return np.maximum(norm.cdf(-z) - 2 * owens_t(-z, a), 0.0)


def hw_bound(F, z, block=800):
    """F: (n x 36) matrix of band-form directions (rows need not be normalised).  Returns (HW, union)."""
    nrm = np.linalg.norm(F, axis=1)
    F = F[nrm > 0] / nrm[nrm > 0, None]
    n = len(F)
    p = norm.sf(z)
    union = n * p
    if n <= 1:
        return union, union
    # Prim's algorithm for the MAXIMUM spanning tree; weights computed row-block on demand
    inT = np.zeros(n, bool); inT[0] = True
    best = orthant_pair(z, F @ F[0])
    best[0] = -1
    total = 0.0
    for _ in range(n - 1):
        cand = np.where(inT, -np.inf, best)
        j = int(np.argmax(cand))
        total += max(cand[j] - 1e-12, 0.0)
        inT[j] = True
        w = orthant_pair(z, F @ F[j])
        best = np.maximum(best, w)
    return max(union - total, p), union
