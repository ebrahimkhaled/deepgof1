"""touch_sim.py -- W5(ii): how often does the covariate-only deletion rule fire on CLEAN data
when there are d covariates, and how often does it still catch a x10 decimal slip?

The rule (PREDECLARATION_TRIM.md): o_i = max_j |x_ij - med_j| / (1.4826 MAD_j), sample median
and MAD per covariate; if max_i o_i > c (c = 4.5 as declared) delete the record attaining the
max. By Corollary S.18 the touch and miss probabilities depend on the covariate law and the slip
mechanism only, so covariates alone are simulated.

Per replicate four numbers are kept, so ANY threshold can be evaluated afterwards:
  omax_clean : max_i o_i on the clean design
  o0_slip    : o of record 0 after x_01 <- 10 x_01 (record 0 w.l.o.g., the law is exchangeable)
  oth_slip   : max over records != 0 of o after the slip (the slip moves med and MAD of column 1)
  x01        : the original value of the slipped entry (gross / mild split)
The rule catches the slip at threshold t iff o0_slip > t and o0_slip > oth_slip.

Laws (all coordinates i.i.d.): normal N(0,1) [small-n studies]; t4 [Liu Setting 2];
unif(-3,3) [Liu Settings 1,3,4, corruption benchmark]; lognormal exp(N(0,1)) [skewed sensitivity,
the application's ratio covariates]; bern(.5) [Liu Setting 3].  Plus liu2 = (unif, N(0,1.5), t4),
d = 3 fixed, the paper's own three-covariate design.

One .npz per cell, written atomically; rerunning skips finished cells. 20,000 replicates per cell
in blocks of 2,000.  Run:  python touch_sim.py
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cells")
os.makedirs(OUT, exist_ok=True)
R, BLOCK = 20000, 2000
NS = (30, 50, 100, 200, 500)
DS = (2, 5, 10)
LAWS = ("normal", "t4", "unif", "lognormal", "bern")
SEED0 = 20260924


def draw(rng, law, size):
    if law == "normal":
        return rng.standard_normal(size)
    if law == "t4":
        return rng.standard_t(4, size)
    if law == "unif":
        return rng.uniform(-3, 3, size)
    if law == "lognormal":
        return np.exp(rng.standard_normal(size))
    if law == "bern":
        return rng.binomial(1, 0.5, size).astype(float)
    if law == "liu2":                     # size = (b, n, 3)
        b, n, d = size
        assert d == 3
        return np.stack([rng.uniform(-3, 3, (b, n)), rng.normal(0, 1.5, (b, n)),
                         rng.standard_t(4, (b, n))], axis=2)
    raise ValueError(law)


def o_of(X):
    """robust outlyingness per record, shape (b, n); MAD = 0 gives inf (as R's x / 0 does)."""
    med = np.median(X, axis=1, keepdims=True)
    dev = np.abs(X - med)
    mad = 1.4826 * np.median(dev, axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = dev / mad
        z = np.where(mad > 0, z, np.where(dev > 0, np.inf, 0.0))
    return z.max(2)


def one_cell(law, d, n, cell_id):
    rng = np.random.default_rng([SEED0, cell_id])
    omax_clean = np.empty(R); o0 = np.empty(R); oth = np.empty(R); x01 = np.empty(R)
    for s in range(0, R, BLOCK):
        X = draw(rng, law, (BLOCK, n, d))
        omax_clean[s:s + BLOCK] = o_of(X).max(1)
        Xs = X.copy(); Xs[:, 0, 0] *= 10.0
        o = o_of(Xs)
        o0[s:s + BLOCK] = o[:, 0]
        oth[s:s + BLOCK] = o[:, 1:].max(1)
        x01[s:s + BLOCK] = X[:, 0, 0]
    return dict(omax_clean=omax_clean, o0_slip=o0, oth_slip=oth, x01=x01)


def cells():
    cid = 0
    for law in LAWS:
        for d in DS:
            for n in NS:
                cid += 1; yield cid, law, d, n
    for n in NS:
        cid += 1; yield cid, "liu2", 3, n


if __name__ == "__main__":
    for cid, law, d, n in cells():
        f = os.path.join(OUT, "%s_d%02d_n%03d.npz" % (law, d, n))
        if os.path.exists(f):
            continue
        res = one_cell(law, d, n, cid)
        tmp = f + ".part.npz"
        np.savez(tmp, law=law, d=d, n=n, R=R, seed=[SEED0, cid], **res)
        os.replace(tmp, f)
        touch = float((res["omax_clean"] > 4.5).mean())
        caught = float(((res["o0_slip"] > 4.5) & (res["o0_slip"] > res["oth_slip"])).mean())
        print("%-9s d=%2d n=%3d  clean touched %.4f  slip caught %.4f" % (law, d, n, touch, caught))
        sys.stdout.flush()
    print("done")
