"""map_train.py -- the network and the training loop of DeepGOF-1 (PyTorch).

CNN     1x6x6 -> Conv3x3(16) -> Conv3x3(32) -> MaxPool2 -> Conv3x3(32) -> GlobalMax + GlobalMean (64)
        -> FC64 -> FC1; the architecture of the shipped network (18,273 parameters at K = 6).
train   Adam (lr 2e-3, weight decay 1e-4), binary cross-entropy, 25 epochs, batch 1024.
score   batched forward pass.

The class and the two functions are unchanged from the module that trained the shipped network and
the per-block networks of the benchmark (benchmark/liu/liu_score.py, benchmark/grid/s3_score.py);
only a gate driver and its unused feature list were removed from this copy. training/02_train.py is
the same file.
"""
import numpy as np
import torch
import torch.nn as nn

DEV = "cuda" if torch.cuda.is_available() else "cpu"
K = 12


class CNN(nn.Module):
    """Convolution finds a local residual bump ANYWHERE in covariate space -- the property the
    risk-ordered 1-D summary cannot have. Global max AND mean pooling: max for a localized bump,
    mean for a diffuse trend."""
    def __init__(self, n_aux=0, grid=K):
        # Global max/mean pooling makes the head size grid-independent (always 64), so only the
        # reshape below needs to know the grid.
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
