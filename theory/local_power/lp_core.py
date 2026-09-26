"""lp_core.py -- the light pieces shared by 01_local_power.py (formulas) and 02_simulate.py
(simulation workers): the shipped network, the departures, a batched logistic fit, and the three
tests (DeepGOF-1, Hosmer-Lemeshow, Stukel) exactly as the benchmark deploys them."""
import os, csv
import numpy as np
from scipy.special import expit
from scipy.stats import chi2
import torch
import torch.nn.functional as F

RD = os.environ.get("DEEPGOF_RDEPLOY", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "_support", "rdeploy"))
torch.set_default_dtype(torch.float64)
torch.set_num_threads(1)
csv.field_size_limit(10 ** 9)
K, B, KREJ, ALPHA = 6, 199, 10, 0.05
BETA0 = np.array([-0.5, 1.0, 0.6])
DEPS = ["quadratic", "interaction", "cubic", "bump", "threshold"]

W = {}
with open(os.path.join(RD, "weights.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        W[r["name"]] = torch.tensor(np.array(r["values"].split(), dtype=float).reshape(
            tuple(int(x) for x in r["shape"].split("x"))))
SC = np.loadtxt(os.path.join(RD, "scaler.csv"), delimiter=",")
MU, SD = torch.tensor(SC[0]), torch.tensor(SC[1])

def score(m):
    """S~(m) = S((m - a)/s) for raw maps m (batch x 36), row-major with x1 bins as rows."""
    with torch.no_grad():
        z = ((torch.as_tensor(np.asarray(m, float)) - MU) / SD).view(-1, 1, K, K)
        h = F.relu(F.conv2d(z, W["conv.0.weight"], W["conv.0.bias"], padding=1))
        h = F.relu(F.conv2d(h, W["conv.2.weight"], W["conv.2.bias"], padding=1))
        h = F.max_pool2d(h, 2)
        h = F.relu(F.conv2d(h, W["conv.5.weight"], W["conv.5.bias"], padding=1))
        e = torch.cat([h.amax(dim=(2, 3)), h.mean(dim=(2, 3))], 1)
        h = F.relu(e @ W["head.0.weight"].T + W["head.0.bias"])
        return (h @ W["head.2.weight"].T + W["head.2.bias"]).squeeze(-1).numpy()

def check_bench():
    X_ = np.loadtxt(os.path.join(RD, "test_input.csv"), delimiter=",")
    y_ = np.loadtxt(os.path.join(RD, "test_output.csv"))
    err = float(np.abs(score(X_ * SC[1] + SC[0]) - y_).max())
    assert err < 1e-8, "not the shipped network"
    return err

def gfun(name, x1, x2):
    """Bounded departures on the logit scale; h = tau * v * g."""
    g = {"quadratic": x1 ** 2 - 1,
         "interaction": x1 * x2,
         "cubic": x1 ** 3 - 3 * x1,
         "bump": 3 * np.exp(-((x1 - 1) ** 2 + (x2 + 0.5) ** 2) / 0.5),
         "threshold": 2.0 * (x2 > 1)}[name]
    return np.clip(g, -6, 6)

def newton(X, Y, beta, iters=12):
    """Batched logistic MLE: Y (b x n), X (n x p), beta (b x p)."""
    p = X.shape[1]
    XX = (X[:, :, None] * X[:, None, :]).reshape(len(X), p * p)   # n x p^2, built once
    for _ in range(iters):
        P = expit(beta @ X.T)
        grad = (Y - P) @ X
        H = ((P * (1 - P)) @ XX).reshape(-1, p, p)                  # no b x n x p x p intermediate
        beta = beta + np.linalg.solve(H, grad[:, :, None])[:, :, 0]
    return beta

def one_rep(args):
    """One dataset under pi_n = pi_0 + tau v g / sqrt(n); returns
    (name, n, p_DeepGOF-1, p_HL, p_Stukel, number of clipped records, seed)."""
    name, tau, n, seed = args[:4]
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 2)); X = np.column_stack([np.ones(n), x])
    pi0 = expit(X @ BETA0)
    h = 0.0 if name == "null" else tau * pi0 * (1 - pi0) * gfun(name, x[:, 0], x[:, 1])
    raw = pi0 + h / np.sqrt(n)
    nclip = int(((raw < 1e-9) | (raw > 1 - 1e-9)).sum())     # records where pi_n left (0,1)
    pin = np.clip(raw, 1e-9, 1 - 1e-9)
    y = (rng.uniform(size=n) < pin).astype(float)
    b = newton(X, y[None, :], np.zeros((1, 3)), 25)[0]
    ph = expit(X @ b); vh = ph * (1 - ph); r = y - ph
    # DeepGOF-1: rank cells (ties 'first'), x1 bins as rows, 199 parametric bootstrap refits
    r1 = np.empty(n, int); r1[np.argsort(x[:, 0], kind="stable")] = np.arange(n)
    r2 = np.empty(n, int); r2[np.argsort(x[:, 1], kind="stable")] = np.arange(n)
    cell = np.minimum(K * r1 // n, K - 1) * K + np.minimum(K * r2 // n, K - 1)
    O = np.zeros((n, K * K)); O[np.arange(n), cell] = 1
    m = (r @ O) / np.sqrt(np.maximum(vh @ O, 1e-8))
    T = score(m[None, :])[0]
    Yb = (rng.uniform(size=(B, n)) < ph[None, :]).astype(float)
    bb = newton(X, Yb, np.tile(b, (B, 1)), 10)
    Pb = expit(bb @ X.T)
    mb = ((Yb - Pb) @ O) / np.sqrt(np.maximum((Pb * (1 - Pb)) @ O, 1e-8))
    Tb = score(mb)
    Tb[~np.isfinite(Tb)] = np.inf                  # a failed refit counts against rejection
    p_dg = (1 + (Tb >= T).sum()) / (B + 1)
    # Hosmer-Lemeshow, 10 equal-size groups by fitted risk, chi-square 8
    grp = np.array_split(np.argsort(ph, kind="stable"), 10)
    HL = sum((y[g].sum() - ph[g].sum()) ** 2 / (ph[g].sum() * (1 - ph[g].sum() / len(g))) for g in grp)
    p_hl = chi2.sf(HL, 8)
    # Stukel, joint score test on both added covariates, chi-square 2
    eh = X @ b
    Z = np.column_stack([0.5 * eh ** 2 * (eh >= 0), -0.5 * eh ** 2 * (eh < 0)])
    s = Z.T @ r
    A = Z.T @ (vh[:, None] * X)
    Cv = Z.T @ (vh[:, None] * Z) - A @ np.linalg.solve(X.T @ (vh[:, None] * X), A.T)
    p_st = chi2.sf(float(s @ np.linalg.solve(Cv, s)), 2)
    return (name, n, float(p_dg), float(p_hl), float(p_st), nclip, seed)
