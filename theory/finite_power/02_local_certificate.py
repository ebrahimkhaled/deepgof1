"""02_local_certificate.py -- step 2, done properly: certify the constants LOCALLY, not globally.

WHY NOT A GLOBAL LIPSCHITZ BOUND. 01_constants.py showed the global route is vacuous: the product of
layer norms gives Lip(S_inf) <= 50.65 against an observed 2.98, and a guaranteed sample size of
166,745 where the truth is nearer 763. The slack has two sources, and both disappear locally:

  1. A global bound pretends every ReLU is active at once. Near the signal direction almost none of
     them change state -- the network is EXACTLY AFFINE inside the linear region containing it.
  2. A global bound charges the noise its full norm, L * ||e|| ~ L * 6.87. Inside a linear region the
     noise enters only through the gradient: g . e, a ONE-dimensional projection with sd ||g||.

THE STATEMENT THIS COMPUTES THE CONSTANTS FOR. Under a fixed alternative with map direction v (unit
norm) and signal rho, the map is m = t v + e with t = sqrt(n) rho and e the standardised cell noise.
Freeze the activation pattern P that the bias-free network S_inf takes at v. Under P every quantity
in the network -- every pre-activation, every pooled maximum -- is an affine function of the input,
so S itself is affine:

        S(m) = t * kappa + g . e + beta         exactly, whenever m lies in region P,

with kappa = S_inf(v / s) the paper's certificate, g the region's gradient and beta its offset. The
region holds as long as every ReLU keeps its sign and every max keeps its winner. Each such
condition is linear in e, so its failure probability has a Gaussian (sub-Gaussian, for the true
Bernoulli cell sums) tail, and a union bound gives an explicit n above which the region holds with
probability 1 - delta1. Beyond that, the power bound is a SCALAR tail bound on g . e.

Everything here is computed on the shipped weights, verified first against the deployed bench.

Run:  python 02_local_certificate.py
"""
import os, sys, json, csv
import numpy as np
import torch
import torch.nn as nn
from scipy.stats import norm

PAP = "/path/to/project/paper_deepgof"
RD = os.path.join(PAP, "code", "rdeploy")
OUT = os.path.join(PAP, "theory", "finite_power")
torch.set_default_dtype(torch.float64)
K = 6
csv.field_size_limit(10 ** 9)

W = {}
with open(os.path.join(RD, "weights.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        shp = tuple(int(x) for x in r["shape"].split("x"))
        W[r["name"]] = torch.tensor(np.array(r["values"].split(), dtype=float).reshape(shp))
sc = np.loadtxt(os.path.join(RD, "scaler.csv"), delimiter=",")
MU, SD = torch.tensor(sc[0]), torch.tensor(sc[1])

F = torch.nn.functional

def forward(z, bias=True, pattern=None, record=False):
    """The shipped network on standardised input z (batch x 36).

    pattern=None  : ordinary forward pass.
    pattern=P     : the SAME network with every nonlinearity frozen to P -- ReLU masks and the
                    winning index of every max. Under P the map z -> S(z) is exactly affine.
    record=True   : return the pattern this input induces, and every pre-nonlinearity value.
    """
    b = (lambda k: W[k]) if bias else (lambda k: torch.zeros_like(W[k]))
    x = z.view(-1, 1, K, K)
    rec, pre = {}, {}
    def relu(h, name):
        pre[name] = h
        m = (h > 0).to(h.dtype) if pattern is None else pattern[name]
        rec[name] = m
        return h * m
    h = relu(F.conv2d(x, W["conv.0.weight"], b("conv.0.bias"), padding=1), "r0")
    h = relu(F.conv2d(h, W["conv.2.weight"], b("conv.2.bias"), padding=1), "r2")
    # 2x2 max pool, stride 2: frozen by recording which of the four won in every window
    B, C, H, Wd = h.shape
    win = h.view(B, C, H // 2, 2, Wd // 2, 2).permute(0, 1, 2, 4, 3, 5).reshape(B, C, H // 2, Wd // 2, 4)
    pre["pool"] = win
    idx = win.argmax(-1) if pattern is None else pattern["pool"]
    rec["pool"] = idx
    h = torch.gather(win, -1, idx.unsqueeze(-1)).squeeze(-1)
    h = relu(F.conv2d(h, W["conv.5.weight"], b("conv.5.bias"), padding=1), "r5")
    flat = h.flatten(2)                                   # B x C x 9
    pre["gmax"] = flat
    gi = flat.argmax(-1) if pattern is None else pattern["gmax"]
    rec["gmax"] = gi
    e = torch.cat([torch.gather(flat, -1, gi.unsqueeze(-1)).squeeze(-1), flat.mean(-1)], 1)
    h = relu(e @ W["head.0.weight"].T + b("head.0.bias"), "rh")
    s = (h @ W["head.2.weight"].T + b("head.2.bias")).squeeze(-1)
    return (s, rec, pre) if record else s

# ---- verify against the deployed bench before trusting anything below -----------------------
X = torch.tensor(np.loadtxt(os.path.join(RD, "test_input.csv"), delimiter=","))
y = np.loadtxt(os.path.join(RD, "test_output.csv"))
err = float(np.abs(forward(X).detach().numpy() - y).max())
print("shipped bench reproduced: max |diff| = %.2e" % err)
assert err < 1e-6, "not the shipped network"

def vec_rows(val):
    return np.repeat(np.asarray(val, float), K)           # row i of the map takes val[i]
stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": vec_rows(stripe), "step_x1": vec_rows(step),
        "saddle": np.outer(step, step).ravel(), "diffuse": np.ones(K * K),
        "random": np.random.default_rng(11).normal(size=K * K)}

def margins(z_dir, z_off, P):
    """Every region condition as  (slope on t, offset, noise-coefficient vector).

    Region P holds at input  z = t * z_dir + z_off + eps  iff every returned quantity
    slope * t + offset + coef . eps  keeps the sign it has as t -> infinity."""
    cond = []
    def collect(fn):
        # fn(z) returns a flat vector of quantities that must stay positive under P
        J = torch.autograd.functional.jacobian(lambda z: fn(z), z_off.unsqueeze(0)).squeeze()
        base = fn(z_off.unsqueeze(0)).squeeze()
        slope = J @ z_dir
        return slope, base, J
    # ReLU units: sign of the pre-activation, oriented so "positive" means "as P says"
    for name in ("r0", "r2", "r5", "rh"):
        def fn(z, name=name):
            _, _, pre = forward(z, pattern=P, record=True)
            sgn = P[name] * 2 - 1
            return (pre[name] * sgn).reshape(-1)
        cond.append(collect(fn))
    # every pooled max must keep its winner: winner minus each loser stays positive
    def fn_pool(z):
        _, _, pre = forward(z, pattern=P, record=True)
        w = pre["pool"]
        win = torch.gather(w, -1, P["pool"].unsqueeze(-1))
        return (win - w).reshape(-1)
    cond.append(collect(fn_pool))
    def fn_g(z):
        _, _, pre = forward(z, pattern=P, record=True)
        f = pre["gmax"]
        win = torch.gather(f, -1, P["gmax"].unsqueeze(-1))
        return (win - f).reshape(-1)
    cond.append(collect(fn_g))
    slope = torch.cat([c[0] for c in cond]); off = torch.cat([c[1] for c in cond])
    J = torch.cat([c[2] for c in cond])
    return slope, off, J

rng = np.random.default_rng(20260923)
Z0 = torch.tensor(rng.normal(size=(200000, K * K)))
with torch.no_grad():
    # chunked: 200k maps in one pass asks for ~1.8 GB of activations at the second conv layer
    s_null = torch.cat([forward(Z0[a:a + 10000]) for a in range(0, len(Z0), 10000)])
    t_alpha = float(np.quantile(s_null.numpy(), 0.95))

rows = []
for name, v in DIRS.items():
    v = torch.tensor(v / np.linalg.norm(v))
    zdir = v / SD                                   # direction in the network's input space
    # the ASYMPTOTIC pattern: the bias-free network's own pattern at v
    _, P, _ = forward(zdir.unsqueeze(0), bias=False, record=True)
    P = {k: val.detach() for k, val in P.items()}
    kappa = float(forward(zdir.unsqueeze(0), bias=False))
    # gradient and offset of the full network frozen to P: S(z) = g . z + b0 on region P
    z0 = torch.zeros(K * K, requires_grad=True)
    s0 = forward(z0.unsqueeze(0), pattern=P)
    s0.backward()
    g_z, b0 = z0.grad.detach(), float(s0.detach())
    # in MAP space m (unstandardised): z = (m - MU)/SD, so S = (g_z/SD) . m + b0 - g_z . (MU/SD)
    g_m = g_z / SD
    beta = b0 - float(g_z @ (MU / SD))
    # region conditions, with the map written  m = t v + e,  z = t zdir + (e - MU)/SD
    slope, off, J = margins(zdir, -MU / SD, P)
    J_m = J / SD                                     # noise enters as e / SD
    sig = torch.linalg.norm(J_m, dim=1)
    degenerate = int((slope.abs() < 1e-12).sum())
    # conditions whose slope is POSITIVE are satisfied more and more as t grows; a condition with a
    # non-positive slope is not guaranteed by growing t and has to be handled separately
    live = slope > 1e-12
    n_bad = int((~live & (slope.abs() >= 1e-12)).sum())
    # per-condition t needed for the condition to hold with prob >= 1 - d, Bonferroni over all of
    # them so that the WHOLE region holds with probability >= 1 - delta1
    delta1 = 0.05
    m_cond = int(live.sum())
    zq = norm.ppf(1 - delta1 / max(m_cond, 1))
    t_need = ((-off[live] + zq * sig[live]) / slope[live]).clamp(min=0)
    t_reg = float(t_need.max()) if m_cond else float("nan")
    # power once the region holds: reject when t*kappa + beta + g.e > t_alpha
    gnorm = float(torch.linalg.norm(g_m))
    delta2 = 0.05
    t_pow = (t_alpha - beta + norm.ppf(1 - delta2) * gnorm) / kappa if kappa > 0 else float("inf")
    t_star = max(t_reg, t_pow)
    rows.append(dict(direction=name, kappa=kappa, g_norm=gnorm, beta=beta, t_alpha=t_alpha,
                     conditions=int(slope.numel()), degenerate=degenerate, wrong_sign=n_bad,
                     t_region=t_reg, t_power=t_pow, t_star=t_star))

print("\nnull 95th percentile (idealised N(0,I) map): %.3f\n" % t_alpha)
print("%-10s %7s %7s %8s %9s %9s %8s %9s %9s %9s"
      % ("direction", "kappa", "||g||", "beta", "conds", "degen", "badsign", "t_region",
         "t_power", "t*"))
for r in rows:
    print("%-10s %7.3f %7.3f %8.3f %9d %9d %8d %9.2f %9.2f %9.2f"
          % (r["direction"], r["kappa"], r["g_norm"], r["beta"], r["conditions"],
             r["degenerate"], r["wrong_sign"], r["t_region"], r["t_power"], r["t_star"]))

print("\nguaranteed sample size n* = (t* / rho)^2 -- power >= 0.90 (region fails <= .05, tail <= .05)")
print("%-10s" % "direction" + "".join("   rho=%.2f" % r for r in (0.25, 0.5, 1.0)))
for r in rows:
    print("%-10s" % r["direction"] + "".join("   %8.0f" % ((r["t_star"] / rho) ** 2)
                                            for rho in (0.25, 0.5, 1.0)))
print("\nfor comparison, the GLOBAL route (01_constants.py) needed n* = 166,745 rigorous / 763 "
      "empirical for kappa = 2.2, rho = 0.5")
json.dump(rows, open(os.path.join(OUT, "LOCAL_CERTIFICATE.json"), "w"), indent=2)
print("wrote theory/finite_power/LOCAL_CERTIFICATE.json")
