"""map_train.py -- GATE M1: can a 2-D residual map beat the 10-group summary on POWER? (2026-08-20)

Three arms, all residuals-only at inference (no z-scores ever reach the score):
  SUMMARY   the incumbent 22 grouped-residual features, MLP
  CNN       12x12 residual map over covariate ranks, small conv net
  CNN_KD    same + a distillation head predicting what each of the 11 classical tests would say.
            The battery teaches the representation during TRAINING; it is never an input.

Size is set by ranking against held-out nulls (exact rank p-value), so all three arms are compared
at the same realised level and the power comparison is size-matched by construction.

⚠️ LIMITATION OF THIS GATE, stated rather than hidden: speckit/10 declared "power >= .25 AND
size-spread <= .025", but the M1 corpus contains ONE design (n=400, p=2). Spread across designs is
therefore NOT measurable here and is not reported. M1 tests the POWER half only. That is the
decisive half -- the residual channel is already known to be well calibrated (spread .016-.019 from
decompose_features.py); what is unknown is whether it can be made powerful. If M1 passes on power,
the multi-design version must then re-test spread before any claim is made.
"""
import os, sys, time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

PAP = r"C:\Users\ebrah\.gemini\Projects\PDFs\Paper_Ebrahim_Frangiton\paper_deepgof"
MAP = os.path.join(PAP, "results", os.environ.get("DEEPGOF_MAPDIR", "map"))
DEV = "cuda" if torch.cuda.is_available() else "cpu"
K = 12
SUMCOLS = [f"g{i}" for i in range(1, 11)] + [f"s{i}" for i in range(1, 11)] + ["hl", "mabs"]
MAPCOLS = [f"m{i}" for i in range(1, K * K + 1)]
TESTS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
         "DEF.poly2", "DEF.poly3", "DEF.stukel", "A.poly", "A.spline"]

def zof(p):
    from scipy.stats import norm
    return norm.ppf(1 - np.clip(np.asarray(p, float), 1e-6, 1 - 1e-6))

class MLP(nn.Module):
    def __init__(self, d_in, hidden=64):
        super().__init__()
        self.f = nn.Sequential(nn.Linear(d_in, hidden), nn.ReLU(),
                               nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, 1))
    def forward(self, x, **kw): return self.f(x).squeeze(1)

class CNN(nn.Module):
    """Convolution finds a local residual bump ANYWHERE in covariate space -- the property the
    risk-ordered 1-D summary cannot have. Global max AND mean pooling: max for a localized bump,
    mean for a diffuse trend."""
    def __init__(self, n_aux=0, grid=K):
        # `grid` added 2026-08-21 for gate S4's resolution sweep. Default is the module-level K, so
        # every existing caller is bit-for-bit unchanged. Global max/mean pooling makes the head
        # size grid-independent (always 64), so only the reshape below needs to know the grid.
        super().__init__()
        self.grid = grid
        self.conv = nn.Sequential(
            nn.Conv2d(1, 16, 3, padding=1), nn.ReLU(),
            nn.Conv2d(16, 32, 3, padding=1), nn.ReLU(),
            nn.MaxPool2d(2),
            nn.Conv2d(32, 32, 3, padding=1), nn.ReLU())
        self.head = nn.Sequential(nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 1))
        self.aux = nn.Linear(64, n_aux) if n_aux else None
    def embed(self, x):
        h = self.conv(x.view(-1, 1, self.grid, self.grid))
        return torch.cat([h.amax(dim=(2, 3)), h.mean(dim=(2, 3))], 1)   # 64
    def forward(self, x, want_aux=False):
        e = self.embed(x)
        s = self.head(e).squeeze(1)
        return (s, self.aux(e)) if (want_aux and self.aux is not None) else s

def train(model, X, y, aux=None, epochs=25, bs=1024, lr=2e-3, kd=0.3):
    X = torch.from_numpy(X).to(DEV); y = torch.from_numpy(y.astype(np.float32)).to(DEV)
    A = torch.from_numpy(aux.astype(np.float32)).to(DEV) if aux is not None else None
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    bce, mse = nn.BCEWithLogitsLoss(), nn.MSELoss()
    n = len(y)
    for ep in range(epochs):
        perm = torch.randperm(n, device=DEV)
        for k in range(0, n, bs):
            i = perm[k:k + bs]
            opt.zero_grad()
            if A is not None:
                s, a = model(X[i], want_aux=True)
                loss = bce(s, y[i]) + kd * mse(a, A[i])
            else:
                loss = bce(model(X[i]), y[i])
            loss.backward(); opt.step()
    model.eval(); return model

@torch.no_grad()
def score(model, X, bs=8192):
    out = np.empty(len(X))
    for k in range(0, len(X), bs):
        out[k:k+bs] = model(torch.from_numpy(X[k:k+bs]).to(DEV)).cpu().numpy()
    return out

def main():
    t0 = time.time()
    tn = pd.read_csv(os.path.join(MAP, "train_null.csv"))
    ta = pd.read_csv(os.path.join(MAP, "train_alt.csv"))
    hn = pd.read_csv(os.path.join(MAP, "held_null.csv"))
    ha = pd.read_csv(os.path.join(MAP, "held_alt.csv"))
    fin = lambda d, c: np.isfinite(d[c].to_numpy(float)).all(1)
    cols = SUMCOLS + MAPCOLS + TESTS
    tn, ta = tn[fin(tn, cols)], ta[fin(ta, cols)]
    hn, ha = hn[fin(hn, cols)], ha[fin(ha, cols)]
    print(f"device {DEV} | train {len(tn)}+{len(ta)} | held {len(hn)}+{len(ha)}", flush=True)

    def mat(d, c): return d[c].to_numpy(dtype=np.float32)
    ytr = np.concatenate([np.zeros(len(tn)), np.ones(len(ta))])
    AUX = np.vstack([zof(tn[TESTS].to_numpy()), zof(ta[TESTS].to_numpy())]).astype(np.float32)
    AUX = np.clip(AUX, -6, 6)

    results = {}
    for name in ("SUMMARY", "CNN", "CNN_KD"):
        cols_i = SUMCOLS if name == "SUMMARY" else MAPCOLS
        Xtr = np.vstack([mat(tn, cols_i), mat(ta, cols_i)])
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-8
        Xtr = (Xtr - mu) / sd
        Xhn = (mat(hn, cols_i) - mu) / sd
        Xha = (mat(ha, cols_i) - mu) / sd
        torch.manual_seed(7001)
        if name == "SUMMARY":
            m = train(MLP(Xtr.shape[1]).to(DEV), Xtr, ytr)
        elif name == "CNN":
            m = train(CNN().to(DEV), Xtr, ytr)
        else:
            m = train(CNN(n_aux=len(TESTS)).to(DEV), Xtr, ytr, aux=AUX)
        Sn, Sa = score(m, Xhn), score(m, Xha)
        # exact rank p-value against the held-out nulls -> all arms at the same realised level
        cut = np.quantile(Sn, 0.95)
        size = float((Sn > cut).mean()); power = float((Sa > cut).mean())
        results[name] = (size, power)
        print(f"  {name:8s} size {size:.4f} | POWER {power:.4f}   ({time.time()-t0:.0f}s)", flush=True)

    print("\n" + "=" * 66)
    print("GATE M1: does a residual MAP beat the 10-group SUMMARY on power?")
    print("  reference points: residual-only channel .08-.16 | z-score channel .41")
    print("  M1 criterion (power half): best map arm >= .25")
    print("=" * 66)
    for k, (s, p) in results.items():
        print(f"  {k:8s} size {s:.4f}  power {p:.4f}")
    best = max(results['CNN'][1], results['CNN_KD'][1])
    base = results['SUMMARY'][1]
    print(f"\n  map gain over summary: {best - base:+.4f}")
    print(f"  M1 power criterion (>= .25): {'PASS' if best >= .25 else 'FAIL'}  (best map arm {best:.4f})")
    print("  NOTE: size-spread across designs NOT measurable here (single design). If the power")
    print("        half passes, the multi-design rerun must test spread before any claim.")

if __name__ == "__main__":
    main()
