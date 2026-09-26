"""04_proj_rate.py -- numerical check of Proposition P (the projection test of Liu et al. 2024 under k
irrelevant covariates), pure design (c_irr = 0), as proved in PROPOSITION.md, section P.

The deployed statistic is T = e'Ae/n^2 with A_ij = sum_r F(|x_i-x_r|^2, |x_j-x_r|^2, (x_i-x_r).(x_j-x_r)),
F(al, be, ga) = sqrt(al be) (pi - arccos(ga / sqrt(al be))), on standardised covariates (proj_verbatim.R).
Its population kernel is K(x, x') = E_z F(|x-z|^2, |x'-z|^2, (x-z).(x'-z)), z an independent copy of x.

Claims checked (p = 2 + k):
  signal  sigma_k = E[h~(x) h~(x') K(x, x')] = Q_P / p (1 + O(p^-1/2)),
          Q_P = F_ab (tr M)^2 + F_gg/2 ||M||_F^2,  F_ab = (pi - 1/sqrt3)/12,  F_gg = 2/(3 sqrt3)
  noise   2 (E v)^2 C_P (1 + o(1)) <= nu_k^2 <= C_P / 8 (1 + o(1)),   C_P = 4 F_ab^2 + F_gg^2 / 2
  and the standardized .95 null quantile q_k used by the corollary (condition: q_k > 1).

Signal: two independent samples of the relevant covariates (x_R and x_R'), each with h~ projected in-sample
on (1, x_R) so that every kernel term that is constant or linear in one argument cancels exactly; the
irrelevant block and z enter only through the Gram entries of (x_I - z_I, x_I' - z_I) and z_R, drawn NW times.
Noise: the M x M deployed kernel at a reference sample (M = 1000), hat-projected as in ir_core.QuadForm.

Resumable: every finished (part, p, seed) row is appended to PROJ_RATE.jsonl; PROJ_RATE.json is assembled
at the end.   Run:  python -u 04_proj_rate.py
"""
import os, sys, json, time
import numpy as np
from ir_core import GridDesign, gfun, BASE, DEPS, proj_kernel, project_kernel, ref_sample, expit

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "PROJ_RATE.jsonl")
PS_SIG = [2, 3, 4, 6, 10, 20, 50, 100, 200, 400]
PS_NOISE = [2, 3, 4, 6, 10, 20, 50]
MS, NW, MN, NZ = 1500, 64, 1000, 20000

def Fgram(al, be, ga):
    S = np.sqrt(al * be)
    rho = np.clip(np.where(S > 1e-12, ga / np.where(S > 1e-12, S, 1.0), 0.0), -1.0, 1.0)
    return S * (np.pi - np.arccos(rho))

# ---------------------------------------------------------------- constants (analytic, checked by differences)
s3 = np.sqrt(3.0)
F_ab, F_gg, F_ag = (np.pi - 1 / s3) / 12, 2 / (3 * s3), -1 / (6 * s3)
eps = 1e-3
fd_ab = (Fgram(2 + eps, 2 + eps, 1) - Fgram(2 + eps, 2 - eps, 1) - Fgram(2 - eps, 2 + eps, 1) + Fgram(2 - eps, 2 - eps, 1)) / (4 * eps * eps)
fd_gg = (Fgram(2, 2, 1 + eps) - 2 * Fgram(2, 2, 1) + Fgram(2, 2, 1 - eps)) / (eps * eps)
fd_ag = (Fgram(2 + eps, 2, 1 + eps) - Fgram(2 + eps, 2, 1 - eps) - Fgram(2 - eps, 2, 1 + eps) + Fgram(2 - eps, 2, 1 - eps)) / (4 * eps * eps)
C_P = 4 * F_ab ** 2 + 0.5 * F_gg ** 2
print("F_ab %.6f (fd %.6f)  F_gg %.6f (fd %.6f)  F_ag %.6f (fd %.6f)  C_P %.6f" % (F_ab, fd_ab, F_gg, fd_gg, F_ag, fd_ag, C_P))

des0 = GridDesign(0, cirr=0.0)
(u0, wu0, pi0g, vg, hl0, Zs0) = des0.nodes[0]
Xtg = np.column_stack([np.ones(des0.N), des0.X1, des0.X2])
Ev = float(des0.Ev)
MM, QP, BB = {}, {}, {}
for dep in DEPS:
    a, ahl, zh, b = des0.drift(lambda x1, x2, pi0, v, u, dep=dep: v * gfun(dep, x1, x2))
    BB[dep] = des0.Iinv @ b
    ht = vg * gfun(dep, des0.X1, des0.X2) - vg * (Xtg @ BB[dep])
    Mm = (Xtg[:, 1:] * ht[:, None]).T @ Xtg[:, 1:] / des0.N
    MM[dep] = Mm.tolist()
    QP[dep] = float(F_ab * np.trace(Mm) ** 2 + 0.5 * F_gg * (Mm ** 2).sum())
    print("%s: M = %s  tr M = %.5f  ||M||_F^2 = %.5f  Q_P = %.5f" % (dep, np.round(Mm, 4).tolist(), np.trace(Mm), (Mm ** 2).sum(), QP[dep]))
lo_nu2, hi_nu2 = 2 * Ev ** 2 * C_P, C_P / 8
print("E v = %.5f | asymptotic nu^2 bounds [%.5f, %.5f] -> nu in [%.4f, %.4f]" % (Ev, lo_nu2, hi_nu2, np.sqrt(lo_nu2), np.sqrt(hi_nu2)))
for dep in DEPS:
    print("  %s: p e_k in [%.4f, %.4f] asymptotically" % (dep, QP[dep] / np.sqrt(hi_nu2), QP[dep] / np.sqrt(lo_nu2)))
sys.stdout.flush()

def load():
    rows = []
    if os.path.exists(RAW):
        for line in open(RAW):
            try:
                rows.append(json.loads(line))
            except ValueError:
                pass
    return rows

def log(rec):
    with open(RAW, "a") as fh:
        fh.write(json.dumps(rec) + "\n")

done = {(r["part"], r["p"], r["seed"]) for r in load()}

# ---------------------------------------------------------------- signal
def htilde(x):
    X = np.column_stack([np.ones(len(x)), x])
    v = expit(X @ BASE); v = v * (1 - v)
    P = X @ np.linalg.solve(X.T @ X, X.T)
    out = {}
    for dep in DEPS:
        h = v * gfun(dep, x[:, 0], x[:, 1]) - v * (X @ BB[dep])
        out[dep] = h - P @ h
    return out

def signal(p, seed):
    k = p - 2
    rng = np.random.default_rng(seed)
    xa, xb = rng.normal(size=(MS, 2)), rng.normal(size=(MS, 2))
    Ha, Hb = htilde(xa), htilde(xb)
    per = {dep: [] for dep in DEPS}
    for w in range(NW):
        zR = rng.normal(size=2)
        if k > 0:
            xI, xIp, zI = rng.normal(size=(3, k))
            aI, bI = xI - zI, xIp - zI
            alI, beI, gaI = aI @ aI, bI @ bI, aI @ bI
        else:
            alI = beI = gaI = 0.0
        da, db = xa - zR, xb - zR
        F = Fgram(alI + (da * da).sum(1)[:, None], beI + (db * db).sum(1)[None, :], gaI + da @ db.T)
        for dep in DEPS:
            per[dep].append(float(Ha[dep] @ F @ Hb[dep]) / MS ** 2)
    return {dep: (float(np.mean(per[dep])), float(np.std(per[dep], ddof=1) / np.sqrt(NW))) for dep in DEPS}

t0 = time.time()
print("\nsignal: p sigma_k against Q_P (mean over %d draws of the irrelevant block, its SE)" % NW)
for p in PS_SIG:
    for seed in (21, 22):
        if ("signal", p, seed) in done:
            continue
        s = signal(p, seed)
        rec = dict(part="signal", p=p, seed=seed, **{"sigma_" + d: s[d][0] for d in DEPS}, **{"se_" + d: s[d][1] for d in DEPS})
        log(rec)
        print("  p %4d seed %d | inter p*sigma %.5f (SE %.5f; Q_P %.5f) | quad %.5f (SE %.5f; Q_P %.5f)  (%.0f s)" % (
            p, seed, p * s["inter"][0], p * s["inter"][1], QP["inter"], p * s["quad"][0], p * s["quad"][1], QP["quad"],
            time.time() - t0)); sys.stdout.flush()

# The two relevant-covariate samples above differ by about 25% (the U-statistic in x_R has sampling error of
# that size at MS = 1500), so the constant is pinned down with 16 further samples at p = 50 and 200
# (NW = 16 draws of the irrelevant block each suffice there: the per-draw spread is small at large p).
print("\nsignal, replicated samples of the relevant covariates (p = 50, 200; 16 samples, NW = 16)")
NW_MAIN = NW
NW = 16
for p in (50, 200):
    for seed in range(100, 116):
        if ("signal_rep", p, seed) in done:
            continue
        s = signal(p, seed)
        log(dict(part="signal_rep", p=p, seed=seed, **{"sigma_" + d: s[d][0] for d in DEPS}, **{"se_" + d: s[d][1] for d in DEPS}))
    rs = [r for r in load() if r["part"] == "signal_rep" and r["p"] == p]
    for d in DEPS:
        vals = np.array([p * r["sigma_" + d] for r in rs])
        print("  p %d %s: mean p*sigma %.5f, SE %.5f over %d samples; Q_P %.5f" % (p, d, vals.mean(), vals.std(ddof=1) / np.sqrt(len(vals)), len(vals), QP[d]))
    sys.stdout.flush()
NW = NW_MAIN

# ---------------------------------------------------------------- noise, largest share, standardized quantile
print("\nnoise: nu_k (reference sample M = %d), bounds, largest eigenvalue share, standardized .95 quantile q_k" % MN)
for p, seed in [(p, 31) for p in PS_NOISE] + [(p, 32) for p in (10, 20, 50)]:
    if ("noise_u", p, seed) in done:
        continue
    k = p - 2
    x, X, pi0, v = ref_sample(k, MN, 5000 + p + 1000 * (seed - 31), cirr=0.0)
    Rc, EQ, VarQ = project_kernel(proj_kernel(x), X, v)
    sv = np.sqrt(v)
    Bm = sv[:, None] * Rc * sv[None, :]
    del Rc
    # U-statistic (off-diagonal) discretisation: the diagonal of the M x M matrix carries the kernel's value at
    # x = x' (about 2 pi p against a typical (4 pi / 3) p off the diagonal); its share of sum lambda^2 is
    # O(p^2 / M) and vanishes only as M -> infinity, so it is removed (it shifts the trace only).
    nu_with_diag = float(np.sqrt(2 * (Bm ** 2).sum()) / MN ** 2)
    Bm[np.diag_indices(MN)] = 0.0
    lam = np.linalg.eigvalsh(Bm)
    del Bm
    lam = lam[np.abs(lam) > 1e-12 * np.abs(lam).max()]
    s2 = float((lam ** 2).sum()); nu = np.sqrt(2 * s2) / MN ** 2
    rng = np.random.default_rng(seed)
    draws = np.concatenate([(rng.normal(size=(2000, len(lam))) ** 2) @ lam for _ in range(NZ // 2000)])
    q95 = float((np.quantile(draws, 0.95) - lam.sum()) / np.sqrt(2 * s2))
    JI2 = (4 * k * k * F_ab ** 2 + 0.5 * F_gg ** 2 * k * (k + 1) + 4 * k * F_ab * F_gg) / p ** 2
    rec = dict(part="noise_u", p=p, seed=seed, nu=float(nu), nu2=float(nu ** 2), lower_finite=float(2 * Ev ** 2 * JI2),
               lower_asym=float(lo_nu2), upper_asym=float(hi_nu2), max_share=float((lam ** 2).max() / s2),
               n_neg=int((lam < -1e-9 * np.abs(lam).max()).sum()), n_keep=int(len(lam)), q95_std=q95,
               nu_with_diag=nu_with_diag, df_mm=float(2 * EQ ** 2 / VarQ))
    log(rec)
    print("  p %3d | nu %.4f  nu^2 %.5f  (lower %.5f finite-k / %.5f asym, upper %.5f) | max share %.3f | neg %d/%d | q95 %.3f  (%.0f s)" % (
        p, nu, nu ** 2, rec["lower_finite"], lo_nu2, hi_nu2, rec["max_share"], rec["n_neg"], rec["n_keep"], q95, time.time() - t0))
    sys.stdout.flush()

rows = load()
out = {"constants": dict(F_ab=F_ab, F_gg=F_gg, F_ag=F_ag, fd_ab=float(fd_ab), fd_gg=float(fd_gg), fd_ag=float(fd_ag), C_P=C_P,
                         Ev=Ev, M=MM, Q_P=QP, nu2_lower_asym=lo_nu2, nu2_upper_asym=hi_nu2),
       "settings": dict(MS=MS, NW=NW, M_noise=MN, NZ=NZ, design="pure, c_irr = 0"),
       "signal": sorted([r for r in rows if r["part"] == "signal"], key=lambda r: (r["p"], r["seed"])),
       "signal_rep": sorted([r for r in rows if r["part"] == "signal_rep"], key=lambda r: (r["p"], r["seed"])),
       "noise": sorted([r for r in rows if r["part"] == "noise_u"], key=lambda r: (r["p"], r["seed"]))}
json.dump(out, open(os.path.join(HERE, "PROJ_RATE.json"), "w"), indent=1)
print("wrote PROJ_RATE.json")
