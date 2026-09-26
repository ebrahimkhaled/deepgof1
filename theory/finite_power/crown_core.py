"""crown_core.py -- generated from 03_certify.py by make_core.py: the shipped network in map
units, its dense linear pieces, the exact pure-ReLU rewrite (NET, NET0), the forward pass fwd,
CROWN (crown), ball sampling and a PGD attack. Importing it has no side effects beyond
loading and verifying the weights."""
import os, sys, json, csv, time
import numpy as np
import torch
import torch.nn.functional as F
from scipy.stats import chi2, binom

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
        W[r["name"]] = np.array(r["values"].split(), dtype=float).reshape(shp)
sc = np.loadtxt(os.path.join(RD, "scaler.csv"), delimiter=",")
MU, SD = sc[0], sc[1]

# ---- reference forward pass (torch, the architecture as trained) ----------------------------
Wt = {k: torch.tensor(v) for k, v in W.items()}
def S_ref(m, bias=True):
    """Score of raw maps m (batch x 36). bias=False gives S_inf: every bias AND the scaler shift
    set to zero, which is what the recession function of the composed map is."""
    b = (lambda k: Wt[k]) if bias else (lambda k: torch.zeros_like(Wt[k]))
    shift = torch.tensor(MU) if bias else torch.zeros(36)
    z = ((m - shift) / torch.tensor(SD)).view(-1, 1, K, K)
    h = F.relu(F.conv2d(z, Wt["conv.0.weight"], b("conv.0.bias"), padding=1))
    h = F.relu(F.conv2d(h, Wt["conv.2.weight"], b("conv.2.bias"), padding=1))
    h = F.max_pool2d(h, 2)
    h = F.relu(F.conv2d(h, Wt["conv.5.weight"], b("conv.5.bias"), padding=1))
    e = torch.cat([h.amax(dim=(2, 3)), h.mean(dim=(2, 3))], 1)
    h = F.relu(e @ Wt["head.0.weight"].T + b("head.0.bias"))
    return (h @ Wt["head.2.weight"].T + b("head.2.bias")).squeeze(-1)

X = np.loadtxt(os.path.join(RD, "test_input.csv"), delimiter=",")
y = np.loadtxt(os.path.join(RD, "test_output.csv"))
with torch.no_grad():
    err = float(np.abs(S_ref(torch.tensor(X * SD + MU)).numpy() - y).max())
print("shipped bench reproduced in map units: max |diff| = %.2e" % err)
assert err < 1e-8, "not the shipped network"

# ---- dense matrices for every linear piece ---------------------------------------------------
def conv_dense(w, cin, hw):
    n = cin * hw * hw
    basis = torch.eye(n).view(n, cin, hw, hw)
    return F.conv2d(basis, torch.tensor(w), padding=1).reshape(n, -1).T.numpy()
A1 = conv_dense(W["conv.0.weight"], 1, 6) / SD[None, :]        # 576 x 36, scaler folded in
A2 = conv_dense(W["conv.2.weight"], 16, 6)                     # 1152 x 576
A3 = conv_dense(W["conv.5.weight"], 32, 3)                     # 288 x 288
b1 = np.repeat(W["conv.0.bias"], 36) - A1 @ MU                 # the scaler shift is a bias
b2 = np.repeat(W["conv.2.bias"], 36)
b3 = np.repeat(W["conv.5.bias"], 9)
H0, hb0 = W["head.0.weight"], W["head.0.bias"]
H2, hb2 = W["head.2.weight"].reshape(-1), float(W["head.2.bias"].reshape(-1)[0])
POOLS = [[c * 36 + (2 * a + i) * 6 + (2 * bb + j) for i in (0, 1) for j in (0, 1)]
         for c in range(32) for a in range(3) for bb in range(3)]     # pooled index c*9 + a*3 + b

rng = np.random.default_rng(20260924)

def layer(Wm, b, relu):
    return dict(W=np.asarray(Wm, float), b=np.asarray(b, float), relu=np.asarray(relu, bool))

def build(bias=True):
    z = (lambda v: v) if bias else (lambda v: np.zeros_like(v))
    Ls = [layer(A1, z(b1), np.ones(576)), layer(A2, z(b2), np.ones(1152))]
    # pool level 1: per pool [relu(q1-q2), q2, relu(q3-q4), q4]
    M = np.zeros((1152, 1152)); rl = np.zeros(1152)
    for k, p in enumerate(POOLS):
        M[4 * k, p[0]], M[4 * k, p[1]] = 1, -1; rl[4 * k] = 1
        M[4 * k + 1, p[1]] = 1
        M[4 * k + 2, p[2]], M[4 * k + 2, p[3]] = 1, -1; rl[4 * k + 2] = 1
        M[4 * k + 3, p[3]] = 1
    Ls.append(layer(M, np.zeros(1152), rl))
    # pool level 2: m1 = x1 + x0, m2 = x3 + x2 ; [relu(m1 - m2), m2]
    M = np.zeros((576, 1152)); rl = np.zeros(576)
    for k in range(288):
        M[2 * k, [4 * k, 4 * k + 1]] = 1; M[2 * k, [4 * k + 2, 4 * k + 3]] = -1; rl[2 * k] = 1
        M[2 * k + 1, [4 * k + 2, 4 * k + 3]] = 1
    Ls.append(layer(M, np.zeros(576), rl))
    P = np.zeros((288, 576))
    for k in range(288):
        P[k, [2 * k, 2 * k + 1]] = 1                              # pooled = m2 + relu(m1 - m2)
    Ls.append(layer(A3 @ P, z(b3), np.ones(288)))
    # global max over 9 per channel by a tree, the mean carried alongside as an identity unit
    # level A: [r(v1-v2), v2, r(v3-v4), v4, r(v5-v6), v6, r(v7-v8), v8, v9, mean]
    M = np.zeros((320, 288)); rl = np.zeros(320)
    for c in range(32):
        o, v = 10 * c, 9 * c
        for q in range(4):
            M[o + 2 * q, v + 2 * q], M[o + 2 * q, v + 2 * q + 1] = 1, -1; rl[o + 2 * q] = 1
            M[o + 2 * q + 1, v + 2 * q + 1] = 1
        M[o + 8, v + 8] = 1
        M[o + 9, v:v + 9] = 1 / 9
    Ls.append(layer(M, np.zeros(320), rl))
    # level B: w_q = x[2q] + x[2q+1] (q < 4), w5 = x8 ; [r(w1-w2), w2, r(w3-w4), w4, w5, mean]
    M = np.zeros((192, 320)); rl = np.zeros(192)
    for c in range(32):
        o, v = 6 * c, 10 * c
        w = [[v + 0, v + 1], [v + 2, v + 3], [v + 4, v + 5], [v + 6, v + 7]]
        M[o, w[0]] = 1; M[o, w[1]] = -1; rl[o] = 1
        M[o + 1, w[1]] = 1
        M[o + 2, w[2]] = 1; M[o + 2, w[3]] = -1; rl[o + 2] = 1
        M[o + 3, w[3]] = 1
        M[o + 4, v + 8] = 1
        M[o + 5, v + 9] = 1
    Ls.append(layer(M, np.zeros(192), rl))
    # level C: y1 = x0 + x1, y2 = x2 + x3, y3 = x4 ; [r(y1-y2), y2, y3, mean]
    M = np.zeros((128, 192)); rl = np.zeros(128)
    for c in range(32):
        o, v = 4 * c, 6 * c
        M[o, [v, v + 1]] = 1; M[o, [v + 2, v + 3]] = -1; rl[o] = 1
        M[o + 1, [v + 2, v + 3]] = 1
        M[o + 2, v + 4] = 1
        M[o + 3, v + 5] = 1
    Ls.append(layer(M, np.zeros(128), rl))
    # level D: u1 = x0 + x1, u2 = x2 ; [r(u1-u2), u2, mean]
    M = np.zeros((96, 128)); rl = np.zeros(96)
    for c in range(32):
        o, v = 3 * c, 4 * c
        M[o, [v, v + 1]] = 1; M[o, v + 2] = -1; rl[o] = 1
        M[o + 1, v + 2] = 1
        M[o + 2, v + 3] = 1
    Ls.append(layer(M, np.zeros(96), rl))
    E = np.zeros((64, 96))
    for c in range(32):
        E[c, [3 * c, 3 * c + 1]] = 1                              # max = u2 + relu(u1 - u2)
        E[32 + c, 3 * c + 2] = 1
    Ls.append(layer(H0 @ E, z(hb0), np.ones(64)))
    Ls.append(layer(H2[None, :], np.array([z(np.array(hb2))]).reshape(1), np.zeros(1)))
    return Ls

def fwd(Ls, m, keep=False):
    x, pres = m.T, []
    for l in Ls:
        zz = l["W"] @ x + l["b"][:, None]
        pres.append(zz)
        x = np.where(l["relu"][:, None], np.maximum(zz, 0), zz)
    return (x[0], pres) if keep else x[0]

NET, NET0 = build(True), build(False)
mt = rng.normal(size=(2000, 36)) * 3
with torch.no_grad():
    d1 = np.abs(fwd(NET, mt) - S_ref(torch.tensor(mt)).numpy()).max()
    d0 = np.abs(fwd(NET0, mt) - S_ref(torch.tensor(mt), bias=False).numpy()).max()
print("\nC. pure-ReLU rewrite reproduces S and S_inf: max |diff| %.1e, %.1e" % (d1, d0))
assert max(d1, d0) < 1e-9

def crown(Ls, c0, r, p=2, want_pre=False):
    """Bounds on the output over {c0 + e : ||e||_p <= r}; CROWN back-substitution with adaptive
    lower slopes. Returns (lower, upper) and, if asked, the per-layer pre-activation bounds."""
    q = {2: 2, np.inf: 1}[p]
    bnds = []
    def backsub(t, sign):
        # bound sign * z_t from BELOW, for every unit of layer t at once
        Lam = sign * Ls[t]["W"]; const = sign * Ls[t]["b"].copy()
        for s in range(t - 1, -1, -1):
            l, u = bnds[s]; rl = Ls[s]["relu"]
            dl = np.ones_like(l); du = np.ones_like(l); cu = np.zeros_like(l)
            off = rl & (u <= 0); dl[off] = 0; du[off] = 0
            un = rl & (l < 0) & (u > 0)
            du[un] = u[un] / (u[un] - l[un]); cu[un] = -u[un] * l[un] / (u[un] - l[un])
            dl[un] = (u[un] >= -l[un]).astype(float)
            Lp, Ln = np.clip(Lam, 0, None), np.clip(Lam, None, 0)
            const = const + Ln @ cu
            Lam = Lp * dl[None, :] + Ln * du[None, :]
            const = const + Lam @ Ls[s]["b"]
            Lam = Lam @ Ls[s]["W"]
        return Lam @ c0 + const - r * np.linalg.norm(Lam, ord=q, axis=1)
    for t in range(len(Ls)):
        lo_t = backsub(t, 1.0); up_t = -backsub(t, -1.0)
        bnds.append((lo_t, up_t))
    lo_T, up_T = bnds[-1]
    return (float(lo_T[0]), float(up_T[0]), bnds) if want_pre else (float(lo_T[0]), float(up_T[0]))

def sample_ball(c0, r, n, p=2):
    if p == 2:
        d = rng.normal(size=(n, 36)); d /= np.linalg.norm(d, axis=1, keepdims=True)
        rad = r * rng.uniform(size=(n, 1)) ** (1 / 36)
        rad[: n // 2] = r                                       # half on the sphere itself
        return c0 + d * rad
    return c0 + rng.uniform(-r, r, size=(n, 36))

def attack_min(c0, r, iters=200, starts=20):
    """Projected gradient descent for min S over the l2 ball -- an UPPER bound on the true min."""
    best = np.inf
    for s in range(starts):
        e = rng.normal(size=36); e *= r / np.linalg.norm(e) * rng.uniform()
        for it in range(iters):
            z = torch.tensor(c0 + e, requires_grad=True)
            v = S_ref(z[None, :])[0]; v.backward()
            best = min(best, float(v))
            e = e - 0.2 * r * z.grad.numpy() / (np.linalg.norm(z.grad.numpy()) + 1e-12)
            nr = np.linalg.norm(e)
            if nr > r: e *= r / nr
    return best

