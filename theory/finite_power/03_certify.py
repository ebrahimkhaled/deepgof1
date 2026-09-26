"""03_certify.py -- step 2: CERTIFIED constants and a certified detection sample size.

01_constants.py showed that norm-product constants are vacuous (L 50.65 against 2.98 observed, c 98.5
against 7.2). 02_local_certificate.py showed that freezing one linear region does not help either,
because hundreds of units sit exactly on their kink along the signal direction. This file replaces
both with bounds that are RIGOROUS for the shipped weights and tight enough to use:

  A. c = sup_m |S(m) - S_inf(m)|, over ALL of R^36, by signed-interval propagation of the offset.
     The difference between the biased and the bias-free network, layer by layer, lies in an
     interval that does not depend on the input: relu(a + x) - relu(a) lies between min(0, x) and
     max(0, x) for every a, and a max moves by no more than its largest argument moves.

  B. L = sup_m ||grad S_inf(m)||_2 (and of S), over ALL of R^36, by interval back-propagation of the
     gradient: every ReLU multiplies by 0 or 1, every max routes to one argument, so each gradient
     coordinate lies in an interval computed once from the weights.

  C. A certified detection threshold WITHOUT any global constant. The max-pools and global maxima
     are rewritten exactly as ReLUs (max(a, b) = b + relu(a - b)), which turns the network into a
     plain ReLU network, and CROWN linear bound propagation (Zhang et al., 2018) then bounds
          min_{||e||_2 <= r}  S(t v + e)
     from below. If that certified minimum exceeds the null critical value, EVERY map within noise
     radius r of the signal is rejected -- so power >= P(||e|| <= r) -- and t = sqrt(n) rho turns
     it into a guaranteed sample size n*.

Every bound is validated against sampling and a gradient attack before it is reported: a
"certified" bound that a sample violates is a bug, and the script stops.

Everything is in MAP units m (the unstandardised 6x6 residual map); the shipped scaler is folded
into the first layer, so S(m) here is exactly what the package computes.

Run:  python 03_certify.py
"""
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

# ============================================================================================
# A. certified c = sup |S - S_inf|, over all of R^36
# ============================================================================================
def iaff(M, lo, hi, b=0.0):
    Mp, Mn = np.clip(M, 0, None), np.clip(M, None, 0)
    return Mp @ lo + Mn @ hi + b, Mp @ hi + Mn @ lo + b
def rclip(lo, hi):                    # relu(a + x) - relu(a), x in [lo, hi], any a
    return np.minimum(lo, 0), np.maximum(hi, 0)

lo, hi = rclip(b1, b1)                                        # layer 1: offset is exactly b1
lo, hi = rclip(*iaff(A2, lo, hi, b2))
lo = np.array([lo[p].min() for p in POOLS]); hi = np.array([hi[p].max() for p in POOLS])
lo, hi = rclip(*iaff(A3, lo, hi, b3))
L9, H9 = lo.reshape(32, 9), hi.reshape(32, 9)
lo = np.concatenate([L9.min(1), L9.mean(1)]); hi = np.concatenate([H9.max(1), H9.mean(1)])
lo, hi = rclip(*iaff(H0, lo, hi, hb0))
clo, chi_ = iaff(H2[None, :], lo, hi, hb2)
c_cert = float(max(abs(clo[0]), abs(chi_[0])))
print("\nA. certified  %.3f <= S - S_inf <= %.3f   everywhere   ->  c_cert = %.3f"
      % (clo[0], chi_[0], c_cert))

# ============================================================================================
# B. certified L = sup ||grad||_2, for S and S_inf alike (any activation pattern)
# ============================================================================================
def gclip(lo, hi):                    # multiply by 0 or 1
    return np.minimum(lo, 0), np.maximum(hi, 0)
g = H2.copy()
lo, hi = gclip(g, g)
lo, hi = iaff(H0.T, lo, hi)                                   # 64: [gmax_c, gmean_c]
gl, gh = gclip(lo[:32], hi[:32])
lo = np.repeat(gl, 9) + np.repeat(lo[32:], 9) / 9
hi = np.repeat(gh, 9) + np.repeat(hi[32:], 9) / 9
lo, hi = gclip(lo, hi)                                        # relu of conv5
lo, hi = iaff(A3.T, lo, hi)                                   # gradient on the 288 pooled units
pl, ph = gclip(lo, hi)                                        # each pool routes to one of four
lo = np.zeros(1152); hi = np.zeros(1152)
for k, p in enumerate(POOLS):
    lo[p] = pl[k]; hi[p] = ph[k]
lo, hi = gclip(lo, hi)
lo, hi = iaff(A2.T, lo, hi)
lo, hi = gclip(lo, hi)
glo, ghi = iaff(A1.T, lo, hi)
Gbox = np.maximum(np.abs(glo), np.abs(ghi))
L_cert = float(np.linalg.norm(Gbox))
print("B. certified  ||grad S||_2 <= %.3f  everywhere (per-cell sup |dS/dm_j| in [%.3f, %.3f])"
      % (L_cert, Gbox.min(), Gbox.max()))

# ---- validation of A and B against sampling and an attack ---------------------------------
rng = np.random.default_rng(20260924)
Zs = np.vstack([rng.normal(size=(40000, 36)) * s for s in (0.3, 1, 3, 10, 30, 100)])
with torch.no_grad():
    gap = np.concatenate([(S_ref(torch.tensor(Zs[a:a + 20000])) -
                           S_ref(torch.tensor(Zs[a:a + 20000]), bias=False)).numpy()
                          for a in range(0, len(Zs), 20000)])
def gnorm(Z, bias):
    out = []
    for a in range(0, len(Z), 20000):
        z = torch.tensor(Z[a:a + 20000], requires_grad=True)
        S_ref(z, bias).sum().backward()
        out.append(z.grad.norm(dim=1).numpy())
    return np.concatenate(out)
gI, gS = gnorm(Zs, False), gnorm(Zs, True)
# attack: push the gradient norm up by random local search over directions and scales
best = Zs[np.argmax(gI)].copy(); bestv = gI.max()
for it in range(3000):
    cand = best + rng.normal(size=36) * rng.choice([0.01, 0.1, 1.0])
    v = gnorm(cand[None, :], False)[0]
    if v > bestv:
        best, bestv = cand, v
c_emp, L_emp = float(np.abs(gap).max()), float(max(gI.max(), gS.max(), bestv))
print("   sampled/attacked:  max|S - S_inf| = %.3f (range %.3f..%.3f)   max ||grad|| = %.3f"
      % (c_emp, gap.min(), gap.max(), L_emp))
assert clo[0] - 1e-9 <= gap.min() and gap.max() <= chi_[0] + 1e-9, "c bound violated -- bug"
assert L_emp <= L_cert + 1e-9, "L bound violated -- bug"
print("   gap certified/observed:  c %.1fx   L %.1fx   (norm products gave 13.6x and 17.0x)"
      % (c_cert / c_emp, L_cert / L_emp))

# ============================================================================================
# C. CROWN on the network rewritten as a pure ReLU network
# ============================================================================================
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

T0 = time.time()
crown(NET, np.zeros(36), 1.0)
print("   one CROWN pass: %.2f s" % (time.time() - T0))

# ---- validate CROWN: every sampled value, and every sampled pre-activation, inside its bounds --
for c0, r in [(np.zeros(36), 1.0), (np.zeros(36), 7.14), (rng.normal(size=36) * 5, 3.0)]:
    lb, ub, bnds = crown(NET, c0, r, want_pre=True)
    pts = sample_ball(c0, r, 20000)
    val, pres = fwd(NET, pts, keep=True)
    ok = all(((bnds[t][0][:, None] - 1e-7 <= pres[t]) & (pres[t] <= bnds[t][1][:, None] + 1e-7)).all()
             for t in range(len(NET)))
    att = attack_min(c0, r, iters=60, starts=8)
    print("   ball r=%.2f at |c0|=%.1f: CROWN [%.3f, %.3f]  sampled [%.3f, %.3f]  attack min %.3f  %s"
          % (r, np.linalg.norm(c0), lb, ub, val.min(), val.max(), att,
             "all layers contained" if ok else "VIOLATION"))
    assert ok and lb <= min(val.min(), att) + 1e-7 and ub >= val.max() - 1e-7, "CROWN unsound -- bug"

# ============================================================================================
# D. the certified detection threshold
# ============================================================================================
# idealised null: m ~ N(0, I_36) in map units. Critical value = the 95% quantile of S, estimated by
# Monte Carlo and replaced by an exact 99.9% upper confidence bound from the order statistics.
N0 = 400000
M0 = rng.normal(size=(N0, 36))
with torch.no_grad():
    s0 = np.sort(np.concatenate([S_ref(torch.tensor(M0[a:a + 20000])).numpy()
                                 for a in range(0, N0, 20000)]))
kk = int(binom.ppf(0.999, N0, 0.95)) + 1
t_hat, t_bar = float(np.quantile(s0, .95)), float(s0[kk - 1])
print("\nD. null critical value: MC %.4f, certified (99.9%% UCB) %.4f" % (t_hat, t_bar))

R2 = float(np.sqrt(chi2.ppf(0.95, 36)))               # P(||e||_2 <= R2) = .95 for e ~ N(0, I)

def stripe_rows(val):
    return np.repeat(np.asarray(val, float), K)
stripe = [1.2, -.4, -1.1, -1.1, -.4, 1.2]; step = [-1, -1, -1, 1, 1, 1]
DIRS = {"stripe_x1": stripe_rows(stripe), "step_x1": stripe_rows(step),
        "saddle": np.outer(step, step).ravel(), "diffuse": np.ones(36),
        "random": np.random.default_rng(11).normal(size=36)}

def certified_t(v, r):
    """Smallest t on a fine grid with CROWN-min over the ball > t_bar, holding at every larger grid
    point up to 64x (checked, not assumed)."""
    grid = np.exp(np.linspace(np.log(1), np.log(20000), 60))
    ok = np.array([crown(NET, t * v, r)[0] > t_bar for t in grid])
    if not ok.any():
        return np.inf, grid, ok
    last_bad = np.where(~ok)[0]
    i0 = last_bad.max() + 1 if len(last_bad) else 0
    if i0 >= len(grid):
        return np.inf, grid, ok
    lo_t, hi_t = (grid[i0 - 1] if i0 > 0 else 0.0), grid[i0]
    for _ in range(30):                                   # refine between the last failure and success
        mid = 0.5 * (lo_t + hi_t)
        if crown(NET, mid * v, r)[0] > t_bar: hi_t = mid
        else: lo_t = mid
    return hi_t, grid, ok

def empirical_t(v, target=0.95, reps=20000):
    """The truth for comparison: smallest t with Monte Carlo power >= target at the MC critical value."""
    E = rng.normal(size=(reps, 36))
    def pw(t):
        with torch.no_grad():
            return float((S_ref(torch.tensor(t * v + E)).numpy() > t_hat).mean())
    lo_t, hi_t = 0.0, 1.0
    while pw(hi_t) < target:
        hi_t *= 2
        if hi_t > 1e5:
            return np.inf                                 # never reaches the target: say so
    for _ in range(40):
        mid = 0.5 * (lo_t + hi_t)
        if pw(mid) >= target: hi_t = mid
        else: lo_t = mid
    return hi_t

C_OLD = json.load(open(os.path.join(OUT, "CONSTANTS.json")))
rows = []
print("\n  direction    kappa   t_cert(ball r=%.2f)  t_emp(power .95)   ratio   n*_cert  n*_emp   global-constant n*"
      % R2)
for name, v in DIRS.items():
    v = v / np.linalg.norm(v)
    with torch.no_grad():
        kappa = float(S_ref(torch.tensor(v[None, :]), bias=False)[0])
    tc, grid, ok = certified_t(v, R2)
    te = empirical_t(v)
    # the old route with the NEW certified constants: t >= (t_bar + c + L R2) / kappa
    tg = (t_bar + c_cert + L_cert * R2) / kappa if kappa > 0 else np.inf
    rho = 0.5                                                 # ||mu||, as in 01_constants.py
    rows.append(dict(direction=name, kappa=kappa, t_cert=tc, t_emp=te, t_global=tg,
                     n_cert=(tc / rho) ** 2, n_emp=(te / rho) ** 2, n_global=(tg / rho) ** 2,
                     monotone_after=bool(ok[np.searchsorted(grid, tc):].all()) if np.isfinite(tc) else False))
    print("  %-10s %6.3f   %12.2f        %10.2f      %6.2f  %8.0f %7.0f   %10.0f"
          % (name, kappa, tc, te, tc / te, (tc / rho) ** 2, (te / rho) ** 2, (tg / rho) ** 2))

json.dump(dict(bench_err=err, c_cert=c_cert, c_range=[float(clo[0]), float(chi_[0])], c_emp=c_emp,
               L_cert=L_cert, L_emp=L_emp, grad_box=Gbox.tolist(),
               t_hat=t_hat, t_bar=t_bar, R2=R2, rho=0.5, rows=rows,
               old=dict(L_rig=C_OLD["L_rig"], c_rig=C_OLD["c_rig"])),
          open(os.path.join(OUT, "CERTIFIED.json"), "w"), indent=2)
print("\nwrote theory/finite_power/CERTIFIED.json")
