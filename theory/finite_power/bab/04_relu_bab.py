"""04_relu_bab.py -- probabilistic branch-and-bound by PRE-ACTIVATION (ReLU) splitting, BaB-prob style
(Wang & Tsiotras 2025, arXiv 2509.25647), in its simplest sound form.

ROOT.  Probabilistic CROWN at band level z (bab_core.cell_pcrown with k = 0): every ReLU unit j gets
linear forms  Lf_j(e) <= z_j <= Uf_j(e)  and bands [l_j, u_j]; E = {all 2 N_RELU one-sided band
statements hold},  P(E^c) <= 2 N_RELU Phi(-z) =: delta_int.  On E every relaxation below is valid.

LEAF.  A leaf R is a set of split decisions  z_j >= 0 (active) or z_j <= 0 (inactive).  On R cap E the
split units are LINEAR (slope 1 or 0, no offset) instead of relaxed, so the final back-substitution
gives a leaf lower function  S >= F_R(e) = lam_R'e + a_R  on R cap E  (intermediate bands are the root's,
still valid on E).  Also on E:  z_j >= 0  implies  Uf_j(e) >= 0,  z_j <= 0 implies  Lf_j(e) <= 0.
Hence R cap E is contained in the polyhedron  H_R = {Uf_j >= 0 (active), Lf_j <= 0 (inactive)}  and
    P(S <= t_bar, R, E) <= P(F_R <= t_bar, H_R) <= min( P(F_R <= t_bar), min_j P(F_R <= t_bar, H_j) ),
the last terms being EXACT bivariate normal probabilities (Owen's T).  Leaves partition R^36, so
    P(S <= t_bar) <= delta_int + sum_R  leafbound(R).
Power >= .95 is certified when this is <= .05.

Branching: split the leaf with the largest leafbound, on the not-yet-split unstable unit with the largest
final-backsubstitution looseness score  |coef_j| * (u_j - l_j) (upper-relaxation offset for negative
coefficients, lower-relaxation gap for positive ones).

Usage: python -u 04_relu_bab.py <direction> <t> <max_leaves>
"""
import os, sys, json, time, heapq
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.stats import norm
from scipy.special import owens_t
from bab_core import NET, DIRS, T_BAR, N_RELU

ZG = (4.3, 4.4, 4.5, 4.65, 4.8)
MODE = 'any'
DELTA = 0.05


def bvn_cdf(h, k, rho):
    """P(X <= h, Y <= k), standard bivariate normal with correlation rho; exact via Owen's T."""
    if rho > 1 - 1e-12: return norm.cdf(min(h, k))
    if rho < -1 + 1e-12: return max(0.0, norm.cdf(h) + norm.cdf(k) - 1)
    if h == 0 and k == 0: return 0.25 + np.arcsin(rho) / (2 * np.pi)
    s = np.sqrt(1 - rho * rho)
    if h == 0: h = 1e-300
    if k == 0: k = 1e-300
    ah = (k - rho * h) / (h * s); ak = (h - rho * k) / (k * s)
    beta = 0.0 if h * k > 0 or (h * k == 0 and h + k >= 0) else 0.5
    return 0.5 * norm.cdf(h) + 0.5 * norm.cdf(k) - owens_t(h, ah) - owens_t(k, ak) - beta


def root(c0, z):
    """Probabilistic CROWN at level z, keeping every ReLU unit's lower/upper linear forms."""
    bnds, forms = [], []
    def backsub(t, sign, fin_split=None):
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
        return Lam, Lam @ c0 + const
    for t in range(len(NET) - 1):
        Ll, al = backsub(t, 1.0); Lu, au = backsub(t, -1.0)
        l = al - z * np.linalg.norm(Ll, axis=1); u = -au + z * np.linalg.norm(Lu, axis=1)
        rl = NET[t]["relu"]
        l = np.where(rl, l, -np.inf); u = np.where(rl, u, np.inf)
        bnds.append((l, u)); forms.append((Ll, al, -Lu, -au))    # Lf = Ll e + al ; Uf = (-Lu) e + (-au)
    return bnds, forms


def leaf_final(c0, bnds, split):
    """Final lower function with split units made exact; returns lam, a and per-unit looseness scores."""
    T = len(NET) - 1
    Lam = NET[T]["W"].copy(); const = NET[T]["b"].copy(); scores = {}
    for s in range(T - 1, -1, -1):
        l, u = bnds[s]; rl = NET[s]["relu"]
        dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
        off = rl & (u <= 0); dl[off] = 0; du[off] = 0
        un = rl & (l < 0) & (u > 0)
        du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
        dl[un] = (u[un] >= -l[un]).astype(float)
        for (ls, j), sg in split.items():
            if ls == s:
                un[j] = False; cu[j] = 0.0
                dl[j] = du[j] = (1.0 if sg > 0 else 0.0)
        c = Lam[0]
        idx = np.where(un)[0]
        sc = np.abs(c[idx]) * (u[idx] - l[idx])
        for j, v in zip(idx, sc): scores[(s, int(j))] = float(v)
        Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
        const = const + Ln @ cu
        Lam = Lp * dl[None, :] + Ln * du[None, :]
        const = const + Lam @ NET[s]["b"]
        Lam = Lam @ NET[s]["W"]
    return Lam[0], float((Lam @ c0 + const)[0]), scores


def leaf_bound_chain(funcs, split, forms):
    """Every ancestor's lower function is valid on this leaf (its region contains the leaf): take the min."""
    return min(leaf_bound(l_, a_, split, forms) for l_, a_ in funcs)


def leaf_bound(lam, a, split, forms):
    sd = np.linalg.norm(lam)
    h = (T_BAR - a) / sd
    best = norm.cdf(h)
    for (s, j), sg in split.items():
        Ll, al, Ul, au = forms[s]
        if sg > 0:     # region inside {Uf >= 0} = {-Uf <= 0}
            w, c = -Ul[j], -au[j]
        else:          # region inside {Lf <= 0}
            w, c = Ll[j], al[j]
        sw = np.linalg.norm(w)
        if sw == 0:
            if c > 0: return 0.0
            continue
        rho = float(lam @ w / (sd * sw))
        best = min(best, bvn_cdf(h, -c / sw, rho))
    return max(best, 0.0)


def run(direction, t, max_leaves, log=print):
    v = DIRS[direction]; c0 = t * v
    t0 = time.time()
    # choose the root band level (deterministic, before any branching)
    best = None
    for z in ZG:
        bnds, forms = root(c0, z)
        lam, a, _ = leaf_final(c0, bnds, {})
        tot = 2 * N_RELU * norm.sf(z) + leaf_bound(lam, a, {}, forms)
        if best is None or tot < best[0]: best = (tot, z, bnds, forms)
    _, z, bnds, forms = best
    d_int = 2 * N_RELU * norm.sf(z)
    lam, a, sc = leaf_final(c0, bnds, {})
    b0 = leaf_bound(lam, a, {}, forms)
    heap = [(-b0, 0, {}, sc, b0, a, [(lam, a)])]; uid = 1; total = d_int + b0
    hist = [(1, total)]
    log("  t=%.2f root z=%.2f delta_int=%.4f root leafbound=%.4f a=%.2f sd=%.2f  unstable=%d"
        % (t, z, d_int, b0, a, np.linalg.norm(lam), len(sc)))
    while total > DELTA and len(heap) < max_leaves:
        negb, _, split, sc, b, a, funcs = heapq.heappop(heap)
        cand = [kk for kk in sorted(sc, key=sc.get, reverse=True) if kk not in split and (MODE == 'any' or kk[0] == 0)]
        if not cand or b <= 0:
            heapq.heappush(heap, (negb, _, split, sc, b, a, funcs)); break
        unit = cand[0]
        total -= b
        for sg in (1, -1):
            ns = dict(split); ns[unit] = sg
            lam2, a2, sc2 = leaf_final(c0, bnds, ns)
            f2 = funcs + [(lam2, a2)]
            b2 = leaf_bound_chain(f2, ns, forms)
            total += b2
            heapq.heappush(heap, (-b2, uid, ns, sc2, b2, a2, f2)); uid += 1
        hist.append((len(heap), total))
        if len(heap) % 50 == 0:
            log("   leaves %5d  bound %.4f  worst leaf %.4f depth %d  (%.0fs)"
                % (len(heap), total, -heap[0][0], max(len(h[2]) for h in heap), time.time() - t0))
    return dict(direction=direction, t=t, z=z, delta_int=d_int, bound=float(total),
                certified=bool(total <= DELTA), n_leaves=len(heap), wall_s=time.time() - t0,
                history=hist[:: max(1, len(hist) // 200)],
                leaves=[dict(split=[[k[0], k[1], sg] for k, sg in h[2].items()], b=h[4], a=h[5]) for h in heap]), bnds, forms


if __name__ == "__main__":
    d, t, ml = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    MODE = sys.argv[4] if len(sys.argv) > 4 else "any"
    r, _, _ = run(d, t, ml)
    print("RESULT %s t=%.2f bound=%.4f certified=%s leaves=%d wall=%.0fs"
          % (d, t, r["bound"], r["certified"], r["n_leaves"], r["wall_s"]))
    json.dump(r, open("relu_%s_%s_t%g.json" % (MODE, d, t), "w"), indent=1)
