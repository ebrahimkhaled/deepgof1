"""bab_core.py -- probabilistic CROWN on a STRATIFIED Gaussian cell.

Model.  e ~ N(0, I_36).  Fix an orthonormal U (36 x k) and write s = U'e ~ N(0, I_k) and
e_perp = (I - UU') e, which is independent of s.  A CELL is a box  C = {s : lo <= s <= hi}.
For any linear form f(e) = lam'e + a we have  f(e) = a + g's + lam_perp'e_perp  with g = U'lam and
lam_perp'e_perp ~ N(0, ||lam||^2 - ||g||^2) INDEPENDENT of s.  Hence, on C,
    f(e) >= a + sum_i min(g_i lo_i, g_i hi_i) - z * sd_perp
fails with conditional probability  Phi(-z)  given s in C (exactly, one-sided).

cell_pcrown runs CROWN back-substitution layer by layer, setting every ReLU unit's band from the two
one-sided forms above at level z.  Let E_C be the event that all 2 * N_RELU one-sided band statements
hold.  Conditionally on s in C, P(E_C^c) <= 2 * N_RELU * Phi(-z) (union bound), and on E_C every
relaxation is valid, so S(c0 + e) >= a_fin + g_fin's + lam_perp'e_perp.  Therefore
    P(S <= t_bar | s in C) <= 2 N_RELU Phi(-z) + P(a_fin + min_box g_fin's + N(0, sd^2) <= t_bar).
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from crown_core import NET, fwd, K, S_ref

T_BAR = 2.4865397005924965
N_RELU = int(sum(l["relu"].sum() for l in NET))

stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": np.repeat(stripe, K), "step_x1": np.repeat(np.asarray(step, float), K),
        "saddle": np.outer(step, step).ravel().astype(float), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}
DIRS = {k: v / np.linalg.norm(v) for k, v in DIRS.items()}


def box_min(g, lo, hi):
    """min over the box of g's, row-wise; g is (rows x k). Infinite ends are allowed."""
    if g.shape[1] == 0:
        return np.zeros(g.shape[0])
    with np.errstate(invalid="ignore"):
        a = np.where(g > 0, g * lo[None, :], np.where(g < 0, g * hi[None, :], 0.0))
    return a.sum(1)


def cell_pcrown(c0, U, lo, hi, z, want_bnds=False, want_forms=False):
    """Returns (a_min, sd_perp) of the final lower function on the cell (a_min includes the box
    minimum of the pinned part), and optionally the ReLU bands."""
    bnds = []; forms = []
    def backsub(t, sign):
        Lam = sign * NET[t]["W"]; const = sign * NET[t]["b"].copy()
        for s in range(t - 1, -1, -1):
            l, u = bnds[s]; rl = NET[s]["relu"]
            dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
            off = rl & (u <= 0); dl[off] = 0; du[off] = 0
            un = rl & (l < 0) & (u > 0)
            du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
            dl[un] = (u[un] >= -l[un]).astype(float)
            Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
            const = const + Ln @ cu
            Lam = Lp * dl[None, :] + Ln * du[None, :]
            const = const + Lam @ NET[s]["b"]
            Lam = Lam @ NET[s]["W"]
        g = Lam @ U
        sdp = np.sqrt(np.maximum((Lam ** 2).sum(1) - (g ** 2).sum(1), 0.0))
        if want_forms and t < len(NET) - 1:
            rlm = NET[t]["relu"]
            forms.append(-(Lam[rlm] - g[rlm] @ U.T))     # Gaussian (perp) part; violation = this > z sd
        return Lam @ c0 + const + box_min(g, lo, hi), sdp
    for t in range(len(NET) - 1):
        a_lo, s_lo = backsub(t, 1.0); a_up, s_up = backsub(t, -1.0)
        l = a_lo - z * s_lo; u = -a_up + z * s_up
        # identity (non-ReLU) units: their bands are never used by the relaxation
        rl = NET[t]["relu"]
        l = np.where(rl, l, -np.inf); u = np.where(rl, u, np.inf)
        bnds.append((l, u))
    a, sd = backsub(len(NET) - 1, 1.0)
    if want_forms:
        return float(a[0]), float(sd[0]), bnds, np.vstack(forms)
    return (float(a[0]), float(sd[0]), bnds) if want_bnds else (float(a[0]), float(sd[0]))


def cell_fail_bound(a, sd, z):
    """Conditional bound on P(S <= t_bar | cell) for band level z."""
    band = 2 * N_RELU * norm.sf(z)
    if sd == 0:
        fin = 1.0 if a <= T_BAR else 0.0
    else:
        fin = norm.cdf((T_BAR - a) / sd)
    return min(1.0, band + fin)


def cell_pcrown_ua(c0, U, lo, hi, z, details=False):
    """USAGE-AWARE probabilistic CROWN on a cell.  Same bounds as cell_pcrown, but the band budget is
    charged only for the statements the relaxations actually rely on:
      * stable-active unit (l >= 0): relaxation relu(x) = x needs only x >= 0, implied by its lower
        linear form being >= 0: charge P(N(0, s_lo^2) < -a_lo) = Phi(-a_lo / s_lo)  (<= Phi(-z));
      * stable-inactive unit (u <= 0): needs x <= 0: charge Phi(-(-a_up) / s_up) likewise;
      * unstable unit whose triangle UPPER line is used with a nonzero weight in ANY back-substitution
        (a negative coefficient at that unit): needs l <= x <= u: charge 2 Phi(-z);
      * unstable unit only ever met with nonnegative coefficients: its lower relaxation (0 or x) is valid
        for every x, so no statement is needed: charge 0.
    Everything is a deterministic function of (c0, U, box, z).  Returns (a, sd, charge)."""
    bnds, lo_parts, used = [], [], []
    def backsub(t, sign):
        Lam = sign * NET[t]["W"]; const = sign * NET[t]["b"].copy()
        for s in range(t - 1, -1, -1):
            l, u = bnds[s]; rl = NET[s]["relu"]
            dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
            off = rl & (u <= 0); dl[off] = 0; du[off] = 0
            un = rl & (l < 0) & (u > 0)
            du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
            dl[un] = (u[un] >= -l[un]).astype(float)
            Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
            used[s] |= un & (Ln != 0).any(0)
            const = const + Ln @ cu
            Lam = Lp * dl[None, :] + Ln * du[None, :]
            const = const + Lam @ NET[s]["b"]
            Lam = Lam @ NET[s]["W"]
        g = Lam @ U
        sdp = np.sqrt(np.maximum((Lam ** 2).sum(1) - (g ** 2).sum(1), 0.0))
        fin[:] = [Lam, const]
        return Lam @ c0 + const + box_min(g, lo, hi), sdp
    fin = [None, None]
    for t in range(len(NET) - 1):
        a_lo, s_lo = backsub(t, 1.0); a_up, s_up = backsub(t, -1.0)
        l = a_lo - z * s_lo; u = -a_up + z * s_up
        rl = NET[t]["relu"]
        l = np.where(rl, l, -np.inf); u = np.where(rl, u, np.inf)
        bnds.append((l, u)); lo_parts.append((a_lo, s_lo, a_up, s_up)); used.append(np.zeros(len(l), bool))
    a, sd = backsub(len(NET) - 1, 1.0)
    charge = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        for t in range(len(NET) - 1):
            l, u = bnds[t]; rl = NET[t]["relu"]; a_lo, s_lo, a_up, s_up = lo_parts[t]
            act = rl & (l >= 0); ina = rl & (u <= 0) & ~act; un = rl & (l < 0) & (u > 0)
            charge += np.where(s_lo[act] > 0, norm.cdf(-a_lo[act] / s_lo[act]), (a_lo[act] < 0) * 1.0).sum()
            # upper form = -(G + a_up), G ~ N(0, s_up^2); x <= 0 is implied by G >= -a_up
            charge += np.where(s_up[ina] > 0, norm.cdf(-a_up[ina] / s_up[ina]), (a_up[ina] < 0) * 1.0).sum()
            charge += 2 * norm.sf(z) * (un & used[t]).sum()
            if not np.all(np.isfinite(l[rl])) or not np.all(np.isfinite(u[rl])):
                charge = np.inf
    if details:
        return float(a[0]), float(sd[0]), float(charge), dict(bnds=bnds, used=used, lam=fin[0][0], const=float(fin[1][0]))
    return float(a[0]), float(sd[0]), float(charge)
