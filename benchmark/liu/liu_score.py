"""liu_score.py -- Simulation II, stage 2: the DeepGOF-1 p-values of one Liu block, and its summary.

The block's network is trained on that block's own training corpus (half of its nulls, as bootstrap
replicates, against its random training departures, kind 1), with the architecture and training loop
of benchmark/common/map_train.py; the p-value of every dataset is then the rank of its observed score
among its B bootstrap scores. The two corpus seeds of the benchmark (101000, 102000) are the two
independent training corpora of the paper.

Two readings are written to summary.csv:
  raw      each arm at its own realised level (rejection at p <= .05; NA = no rejection);
  matched  every arm randomized to exactly .05 on the block's held-out nulls, so height is power alone.
Table 2 of the paper is recomputed from the per-replicate p-values by tables/make_tables.py.

Run:  python benchmark/liu/liu_score.py --dir results/liu/S2_n100_s101000 [--keep-bin]
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from map_train import CNN, train, DEV                       # noqa: E402
from matched_power import randomized_matched_power, ALPHA   # noqa: E402

K = 6
ARMS = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Stukel", "Osius-Rojek"]


@torch.no_grad()
def sc(net, X, mu, sd, bs=16384):
    o = np.empty(len(X))
    for k in range(0, len(X), bs):
        o[k:k+bs] = net(torch.from_numpy(((X[k:k+bs]-mu)/sd).astype(np.float32)).to(DEV)).cpu().numpy()
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--keep-bin", action="store_true")
    a = ap.parse_args()
    d = a.dir
    m = pd.read_csv(os.path.join(d, "blk_meta.csv")).iloc[0]
    B, n, seed, S, nd = int(m.B), int(m.n), int(m.seed), int(m.setting), int(m.ndatasets)
    idx = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    M = np.fromfile(os.path.join(d, "blk_maps.bin"), dtype=np.float32).reshape(nd, B+1, K*K)
    kind, rep = idx.kind.to_numpy(), idx.rep.to_numpy()
    half = int(m.nnull)//2
    tr = ((kind == 0) & (rep <= half)) | (kind == 1)
    ho = (kind == 0) & (rep > half)

    boot = M[tr][:, 1:, :].reshape(-1, K*K)
    Xalt = M[kind == 1][:, 0, :]
    rng = np.random.default_rng(17)
    Xnull = boot[rng.choice(len(boot), len(Xalt), replace=False)]
    X = np.vstack([Xnull, Xalt]).astype(np.float32)
    y = np.concatenate([np.zeros(len(Xnull)), np.ones(len(Xalt))])
    mu, sd = X.mean(0), X.std(0)+1e-8
    torch.manual_seed(7001 + seed % 1000)
    net = train(CNN(grid=K).to(DEV), (X-mu)/sd, y)
    Sc = sc(net, M.reshape(-1, K*K), mu, sd).reshape(nd, B+1)
    p_deep = (1 + (Sc[:, 1:] >= Sc[:, [0]]).sum(1)) / (B+1)

    arms = {"DEEP": p_deep}
    arms.update({t: idx[t].to_numpy() for t in ARMS})
    CS = sorted(idx.C[kind == 2].unique())
    rows = []
    for arm, p in arms.items():
        pc = np.nan_to_num(p, nan=1.0)
        rows.append(dict(setting=S, n=n, seed=seed, arm=arm, C=0.0, raw=float((pc[ho] <= ALPHA).mean()),
                         matched=np.nan))
        for C in CS:
            msk = (kind == 2) & (idx.C.to_numpy() == C)
            mp = randomized_matched_power(p[ho], p[msk])[0]
            rows.append(dict(setting=S, n=n, seed=seed, arm=arm, C=float(C),
                             raw=float((pc[msk] <= ALPHA).mean()), matched=mp))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(d, "summary.csv"), index=False)
    dep = idx[["kind", "rep", "C"]].copy(); dep["p_deep"] = p_deep
    dep.to_csv(os.path.join(d, "deep_pvalues.csv"), index=False)
    piv = T.pivot_table(index="arm", columns="C", values="raw")
    print(f"\n[LIU Setting {S}, n={n}, seed={seed}]  RAW rejection at nominal .05")
    print(piv.round(3).to_string())
    if not a.keep_bin:
        os.remove(os.path.join(d, "blk_maps.bin")); print("  block .bin deleted")


if __name__ == "__main__":
    main()
