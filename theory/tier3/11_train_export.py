"""11_train_export.py -- Tier 3 item B: train the held-out-family networks and export them for R.

For every nets/<variant>_s<s>/corpus.csv written by 10_corpus.R:
  * normalization constants from the corpus: mu = column mean, sd = column std + 1e-8 (float32), the
    procedure of code/liu_score.py and tier3_timing/09_train_time.py
  * the deposited CNN and train() (release/deepgof1_v2_final/training/02_train.py: Adam lr 2e-3,
    weight decay 1e-4, BCE, 25 epochs, batch 1024), on CPU, torch.manual_seed(9000 + s)
  * export with the deposited export_to_R() (03_export_to_R.py) -> weights.csv, scaler.csv,
    verify_input.csv / verify_output.csv (20 corpus maps, z-space, PyTorch float32 scores)
  * a second bench: 200 corpus maps scored in float64 by the same network (verify64_*.csv), so the
    R port can be checked to ~1e-12 as well as to the deposited 1e-6 gate
Writes nets/<net>/train_info.csv. Skips a network whose weights.csv exists.
"""
import glob
import importlib.util
import os
import time

import numpy as np
import pandas as pd
import torch

TR = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/training"
HERE = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


m = load("dg_train", os.path.join(TR, "02_train.py"))
ex = load("dg_export", os.path.join(TR, "03_export_to_R.py"))
m.DEV = "cpu"
torch.use_deterministic_algorithms(True)

for d in sorted(glob.glob(os.path.join(HERE, "nets", "*_s[123]"))):
    net_name = os.path.basename(d)
    if os.path.exists(os.path.join(d, "weights.csv")):
        print("exists", net_name)
        continue
    s = int(net_name[-1])
    D = pd.read_csv(os.path.join(d, "corpus.csv")).to_numpy(np.float32)
    y, X = D[:, 0], D[:, 1:]
    mu, sd = X.mean(0), X.std(0) + 1e-8
    Xs = ((X - mu) / sd).astype(np.float32)
    torch.manual_seed(9000 + s)
    t0 = time.time()
    net = m.train(m.CNN(grid=6).to("cpu"), Xs, y)
    secs = time.time() - t0
    sc = m.score(net, Xs)
    r = np.argsort(np.argsort(sc)) + 1
    n1, n0 = (y == 1).sum(), (y == 0).sum()
    auc = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    # bench maps: 10 alternatives and 10 nulls (raw map space; export_to_R standardizes them)
    rng = np.random.default_rng(100 + s)
    bench = np.vstack([X[rng.choice(np.where(y == 1)[0], 10, replace=False)],
                       X[rng.choice(np.where(y == 0)[0], 10, replace=False)]])
    ex.export_to_R(net, mu, sd, bench, out_dir=d, prefix="", n_test=20)
    # float64 bench
    idx = rng.choice(len(X), 200, replace=False)
    Z64 = (X[idx].astype(np.float64) - mu.astype(np.float64)) / sd.astype(np.float64)
    net64 = m.CNN(grid=6)
    net64.load_state_dict(net.state_dict())
    net64 = net64.double().eval()
    with torch.no_grad():
        S64 = net64(torch.from_numpy(Z64)).numpy().reshape(-1)
    np.savetxt(os.path.join(d, "verify64_input.csv"), Z64, delimiter=",", fmt="%.17g")
    np.savetxt(os.path.join(d, "verify64_output.csv"), S64, fmt="%.17g")
    torch.save({"state": net.state_dict(), "mu": mu, "sd": sd}, os.path.join(d, "net.pt"))
    pd.DataFrame([dict(net=net_name, rows=len(y), alt=int(n1), null=int(n0), torch_seed=9000 + s,
                       train_secs=round(secs, 2), insample_auc=round(float(auc), 4))]).to_csv(
        os.path.join(d, "train_info.csv"), index=False)
    print(f"{net_name}: {len(y)} maps, {secs:.1f} s, in-sample AUC {auc:.4f}", flush=True)
