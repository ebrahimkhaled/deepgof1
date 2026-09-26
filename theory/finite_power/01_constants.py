"""01_constants.py -- would a finite-sample power guarantee be USEFUL, or vacuous?

THE CANDIDATE THEOREM. The consistency certificate today says only that power tends to one when
S_inf(mu/s) > 0. Both editors read that as thin ("consistency limited to visible residual map with
positive certificate"). The natural upgrade is non-asymptotic: because S is piecewise affine,

    S(z) >= S_inf(z) - c                  (Lemma 2: sup |S - S_inf| = c < inf)
    S_inf(a + b) >= S_inf(a) - L * ||b||  (S_inf is Lipschitz with constant L)

and the map under a fixed alternative is z = sqrt(n) mu/s + e with e a vector of standardised cell
sums of bounded residuals (sub-Gaussian by Hoeffding -- no normal approximation needed). So

    S(z) >= sqrt(n) * S_inf(mu/s) - L * ||e|| - c,

and the test rejects for certain once sqrt(n) * S_inf(mu/s) exceeds  t_alpha + c + L * q(||e||).
Every constant on the right is a property of the FROZEN WEIGHTS, so the result would be a sample
size n* above which detection of a named departure is GUARANTEED -- computed before any data exist.

WHETHER THAT IS WORTH PROVING depends on the size of L and c. A rigorous Lipschitz bound for a
network is usually a product of layer norms, which can be orders of magnitude above the truth; if
so, n* comes out astronomically large and the theorem is true but useless. This script measures:

    L_rig   rigorous upper bound on Lip(S_inf): product of exact layer operator norms
    L_emp   empirical lower bound: largest gradient norm of S_inf seen over many inputs
    c_rig   rigorous upper bound on sup|S - S_inf| by propagating the offsets
    c_emp   empirical lower bound on the same
    t_alpha 95th percentile of S under an idealised null map z ~ N(0, I)

and reports the required signal  G = t_alpha + c + L * q90(||e||)  under both the rigorous and the
empirical constants. If G_rig and G_emp differ by orders of magnitude, the theorem needs a tight
Lipschitz certificate (exact over the linear regions) before it is worth writing.

Run:  python 01_constants.py
"""
import os, sys, json
import numpy as np
import torch
import torch.nn as nn

PAP = "/path/to/project/paper_deepgof"
RD = os.path.join(PAP, "code", "rdeploy")
OUT = os.path.join(PAP, "theory", "finite_power")
torch.set_default_dtype(torch.float64)
K = 6

# ---- the shipped weights, row-major as PyTorch wrote them ------------------------------------
import csv
csv.field_size_limit(10 ** 9)          # one row holds a whole 32x32x3x3 kernel as a single field
W = {}
with open(os.path.join(RD, "weights.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        shp = tuple(int(x) for x in r["shape"].split("x"))
        W[r["name"]] = torch.tensor(np.array(r["values"].split(), dtype=float).reshape(shp))

class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(1, 16, 3, padding=1), nn.ReLU(),
            nn.Conv2d(16, 32, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 32, 3, padding=1), nn.ReLU())
        self.head = nn.Sequential(nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 1))
    def forward(self, x):
        h = self.conv(x.view(-1, 1, K, K))
        e = torch.cat([h.amax(dim=(2, 3)), h.mean(dim=(2, 3))], 1)
        return self.head(e).squeeze(1)

def build(zero_bias=False):
    net = CNN()
    sd = {k: v.clone() for k, v in W.items()}
    if zero_bias:
        for k in sd:
            if k.endswith(".bias"):
                sd[k].zero_()
    net.load_state_dict(sd); net.eval()
    return net

S, Sinf = build(False), build(True)

# ---- first, prove this is the shipped network: reproduce the deployed verification bench -------
X = np.loadtxt(os.path.join(RD, "test_input.csv"), delimiter=",")
y = np.loadtxt(os.path.join(RD, "test_output.csv"))
sc = np.loadtxt(os.path.join(RD, "scaler.csv"), delimiter=",")
MU, SD = sc[0], sc[1]
with torch.no_grad():
    # the bench inputs are ALREADY standardised: deepgof_forward.R scores them as given
    pred = S(torch.tensor(X)).numpy()
err = float(np.abs(pred - y).max())
print("reproduces the shipped bench: max |diff| = %.2e over %d maps" % (err, len(y)))
if err > 1e-6:
    sys.exit("NOT the shipped network -- refusing to report constants for something else")

# ---- rigorous Lipschitz bound: exact operator norm of each linear layer -----------------------
def conv_opnorm(w, hw, iters=300):
    """Largest singular value of a padded 3x3 convolution acting on a hw x hw input."""
    cin = w.shape[1]
    x = torch.randn(1, cin, hw, hw)
    for _ in range(iters):
        x = x / x.norm()
        y = torch.nn.functional.conv2d(x, w, padding=1)
        x = torch.nn.functional.conv_transpose2d(y, w, padding=1)
    x = x / x.norm()
    return float(torch.nn.functional.conv2d(x, w, padding=1).norm())

n0 = conv_opnorm(W["conv.0.weight"], 6)
n2 = conv_opnorm(W["conv.2.weight"], 6)
n5 = conv_opnorm(W["conv.5.weight"], 3)
h0 = float(torch.linalg.matrix_norm(W["head.0.weight"], 2))
h2 = float(torch.linalg.matrix_norm(W["head.2.weight"], 2))
# ReLU and 2x2 max-pooling are 1-Lipschitz in l2. Per channel, the global max is 1-Lipschitz and the
# mean over the 9 pooled positions is 1/3-Lipschitz, so the 64-vector readout is sqrt(1+1/9).
RD_FAC = np.sqrt(1 + 1 / 9)
L_rig = n0 * n2 * n5 * RD_FAC * h0 * h2
print("\nlayer operator norms: conv0 %.3f  conv2 %.3f  conv5 %.3f  head0 %.3f  head2 %.3f"
      % (n0, n2, n5, h0, h2))

# ---- empirical Lipschitz: largest gradient seen, over nulls, large-signal points and directions --
rng = np.random.default_rng(20260923)
def grad_norms(net, Z, bs=4096):
    out = []
    for a in range(0, len(Z), bs):
        z = torch.tensor(Z[a:a + bs], requires_grad=True)
        net(z).sum().backward()
        out.append(z.grad.norm(dim=1).detach().numpy())
    return np.concatenate(out)
Zs = []
for scale in (1.0, 3.0, 10.0, 30.0):
    Zs.append(rng.normal(size=(50000, 36)) * scale)
Zall = np.vstack(Zs)
gS = grad_norms(S, Zall)
gI = grad_norms(Sinf, Zall)
L_emp = float(gI.max())

# ---- the offset gap c = sup|S - S_inf| -------------------------------------------------------
with torch.no_grad():
    t = torch.tensor(Zall)
    gap = (S(t) - Sinf(t)).abs().numpy()
c_emp = float(gap.max())
# rigorous: propagate the offsets through the Lipschitz layers
b = {k: float(W[k].norm()) for k in W if k.endswith(".bias")}
# a conv bias adds the same value at every position: its l2 contribution is |b| * sqrt(positions)
d = b["conv.0.bias"] * 6
d = n2 * d + b["conv.2.bias"] * 6
d = n5 * d + b["conv.5.bias"] * 3
d = RD_FAC * d
d = h0 * d + b["head.0.bias"]
c_rig = h2 * d + b["head.2.bias"]

# ---- the idealised null critical value -------------------------------------------------------
with torch.no_grad():
    s0 = S(torch.tensor(rng.normal(size=(200000, 36)))).numpy()
t_alpha = float(np.quantile(s0, 0.95))

# ---- the noise term: ||e|| for e ~ N(0, I_36), 90th percentile --------------------------------
q90 = float(np.quantile(np.linalg.norm(rng.normal(size=(200000, 36)), axis=1), 0.90))

G_rig = t_alpha + c_rig + L_rig * q90
G_emp = t_alpha + c_emp + L_emp * q90

print("\n                     rigorous (upper)     empirical (lower)")
print("  Lip(S_inf)         %14.3f      %14.3f     ratio %.0fx" % (L_rig, L_emp, L_rig / L_emp))
print("  sup|S - S_inf|     %14.3f      %14.3f     ratio %.0fx" % (c_rig, c_emp, c_rig / max(c_emp, 1e-12)))
print("  t_0.95 (null)      %14.3f" % t_alpha)
print("  q90 ||e||          %14.3f" % q90)
print("  REQUIRED SIGNAL G  %14.1f      %14.1f" % (G_rig, G_emp))
print("\n  Lip(S) itself, empirical max gradient: %.3f (99.9th pct %.3f)"
      % (gS.max(), np.quantile(gS, .999)))

# what G means in sample size, for the paper's own certificate values (Table 1, per unit map norm)
print("\n  sample size to GUARANTEE detection, n* = (G / (kappa * ||mu||))^2, for ||mu|| = 0.5:")
for kappa in (2.200, 1.729, 1.075):
    print("    kappa %.3f :  rigorous n* = %12.0f     empirical n* = %8.0f"
          % (kappa, (G_rig / (kappa * .5)) ** 2, (G_emp / (kappa * .5)) ** 2))

json.dump(dict(bench_err=err, L_rig=L_rig, L_emp=L_emp, c_rig=c_rig, c_emp=c_emp,
               t_alpha=t_alpha, q90=q90, G_rig=G_rig, G_emp=G_emp,
               layer_norms=dict(conv0=n0, conv2=n2, conv5=n5, head0=h0, head2=h2)),
          open(os.path.join(OUT, "CONSTANTS.json"), "w"), indent=2)
print("\nwrote theory/finite_power/CONSTANTS.json")
