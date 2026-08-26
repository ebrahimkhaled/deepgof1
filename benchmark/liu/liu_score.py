"""liu_score.py -- GATE LIU stage 2: our test on Liu et al.'s benchmark, against published tests.

Pre-declared in speckit/26_LIU_BENCHMARK.md. Read it before changing anything here.

Comparators are the EIGHT EXTERNAL published tests and their NA-robust Cauchy combination.
The author's own EF / DEF / EDGE family is deposited and printed for reference but is NEVER the
comparator (26 §2). `proj` and Stute-Zhu are absent; that is a stated limitation, and no
"best available" claim is permitted while they are.

Two readings are reported, because they answer different questions:
  RAW at nominal .05  -- what Liu et al. plot: each test at its own realised level, so an
                         over-sized test's curve is honestly inflated and visible.
  MATCHED level       -- every arm randomized to exactly .05 on the held-out nulls, so height is
                         power alone.
LIU2, the decisive criterion, is judged on the RAW curve because that is the quantity Liu et al.
report and the one their "power functions are not even increasing with C" statement refers to.
"""
import os, sys, argparse, glob
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_train import CNN, train, DEV
from probit_common import randomized_matched_power, narobust_cct, ALPHA

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
K = 6
EXT = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie",
       "Stukel", "Osius-Rojek", "le-Cessie"]
OWN = ["EF", "DEF.poly2", "DEF.poly3", "DEF.stukel", "A.poly", "A.spline"]


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
    d = a.dir if os.path.isabs(a.dir) else os.path.join(PAP, a.dir)
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

    arms = {"DEEP": p_deep, "CCT-8": narobust_cct(idx[EXT].to_numpy(float))}
    arms.update({t: idx[t].to_numpy() for t in EXT})
    arms.update({t: idx[t].to_numpy() for t in OWN})

    CS = sorted(idx.C[kind == 2].unique())
    rows = []
    for arm, p in arms.items():
        pc = np.nan_to_num(p, nan=1.0)
        rows.append(dict(setting=S, n=n, seed=seed, arm=arm, C=0.0, raw=float((pc[ho] <= ALPHA).mean()),
                         matched=np.nan, external=arm in EXT, own=arm in OWN))
        for C in CS:
            msk = (kind == 2) & (idx.C.to_numpy() == C)
            mp, F0, t0, g, na = randomized_matched_power(p[ho], p[msk])
            rows.append(dict(setting=S, n=n, seed=seed, arm=arm, C=float(C),
                             raw=float((pc[msk] <= ALPHA).mean()), matched=mp,
                             external=arm in EXT, own=arm in OWN))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(d, "summary.csv"), index=False)
    dep = idx[["kind", "rep", "C"]].copy(); dep["p_deep"] = p_deep
    dep.to_csv(os.path.join(d, "deep_pvalues.csv"), index=False)

    piv = T.pivot_table(index="arm", columns="C", values="raw")
    order = ["DEEP", "CCT-8"] + EXT
    print(f"\n[LIU Setting {S}, n={n}, seed={seed}]  RAW rejection at nominal .05")
    print(f"{'arm':16s}" + "".join(f"{('size' if c==0 else f'C={c:g}'):>9s}" for c in piv.columns)
          + f"{'rise':>8s}")
    for arm in order:
        if arm not in piv.index: continue
        r = piv.loc[arm]
        star = " *" if arm == "DEEP" else ""
        print(f"{arm:16s}" + "".join(f"{v:>9.3f}" for v in r) + f"{r.iloc[-1]-r.iloc[0]:>+8.3f}{star}")
    print(f"{'-- own family --':16s}")
    for arm in OWN:
        if arm not in piv.index: continue
        r = piv.loc[arm]
        print(f"{arm:16s}" + "".join(f"{v:>9.3f}" for v in r) + f"{r.iloc[-1]-r.iloc[0]:>+8.3f}")
    rise = piv.loc["DEEP"].iloc[-1] - piv.loc["DEEP"].iloc[0]
    if S == 2:
        print(f"\n  LIU2 (Setting 2 only): DEEP rise = {rise:+.3f}  "
              f"[{'PASS' if rise > .15 else 'FAIL'} vs the frozen +.15 bar]")
        ext_rise = {t: piv.loc[t].iloc[-1]-piv.loc[t].iloc[0] for t in EXT if t in piv.index}
        print("  external rises: " + "  ".join(f"{k} {v:+.3f}" for k, v in ext_rise.items()))
    if not a.keep_bin:
        os.remove(os.path.join(d, "blk_maps.bin")); print("  block .bin deleted")


if __name__ == "__main__":
    main()
