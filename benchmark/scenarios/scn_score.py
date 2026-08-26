"""scn_score.py -- GATE SCN stage 2: matched-level power per (scenario, severity, arm).

Pre-declared in speckit/25_SCENARIO_STUDY.md. Read it before changing anything here.

ONE deployment-matched rule per (n, corpus), trained on that block's nulls + random-departure
MIXTURE (the shipped prior) -- never on the single-kind scenarios it is then evaluated against, so
every scenario is a held-out KIND.

COMPARATORS (25 §2): the NINE classical tests and their Cauchy combination CCT-9.
`A.poly` / `A.spline` are computed and deposited but EXCLUDED from every arm used in a criterion or
figure -- they are the author's own unpublished EDGE/EDGES family.

Matched level: every arm placed at exactly alpha = .05 on the held-out nulls by the shared
randomized device, so a curve's height is power alone. Severity 0 is reported as each arm's RAW
size (matched level would force it to .05 by construction and destroy the information).
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
NINE = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie",
        "EF", "DEF.poly2", "DEF.poly3", "DEF.stukel"]
EXCLUDED = ["A.poly", "A.spline"]                       # 25 §2 -- deposited, never a comparator
PARTITION = ["HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie"]


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
    B, n, seed, nd = int(m.B), int(m.n), int(m.seed), int(m.ndatasets)
    idx = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    M = np.fromfile(os.path.join(d, "blk_maps.bin"), dtype=np.float32).reshape(nd, B+1, K*K)
    kind, rep = idx.kind.to_numpy(), idx.rep.to_numpy()
    half = int(m.nnull)//2
    tr = ((kind == 0) & (rep <= half)) | (kind == 1)
    ho = (kind == 0) & (rep > half)
    print(f"[SCN n={n} seed={seed}] train {tr.sum()} | held-out null {ho.sum()} | eval {(kind==2).sum()}",
          flush=True)

    # deployment-matched rule, mix pinned at 1.0
    boot = M[tr][:, 1:, :].reshape(-1, K*K)
    Xalt = M[kind == 1][:, 0, :]
    rng = np.random.default_rng(17)
    Xnull = boot[rng.choice(len(boot), len(Xalt), replace=False)]
    X = np.vstack([Xnull, Xalt]).astype(np.float32)
    y = np.concatenate([np.zeros(len(Xnull)), np.ones(len(Xalt))])
    mu, sd = X.mean(0), X.std(0)+1e-8
    torch.manual_seed(7001 + seed % 1000)
    net = train(CNN(grid=K).to(DEV), (X-mu)/sd, y)
    S = sc(net, M.reshape(-1, K*K), mu, sd).reshape(nd, B+1)
    p_deep = (1 + (S[:, 1:] >= S[:, [0]]).sum(1)) / (B+1)

    arms = {"DEEP": p_deep, "CCT-9": narobust_cct(idx[NINE].to_numpy(float))}
    arms.update({t: idx[t].to_numpy() for t in NINE})
    arms.update({t: idx[t].to_numpy() for t in EXCLUDED})    # deposited, flagged below

    rows = []
    # severity 0 = the null block: report RAW size per arm
    for arm, p in arms.items():
        rows.append(dict(n=n, seed=seed, scn="(null)", sev=0.0, arm=arm,
                         value=float((np.nan_to_num(p[ho], nan=1.0) <= ALPHA).mean()),
                         kindcol="size", excluded=arm in EXCLUDED, n_eval=int(ho.sum())))
    ev = idx[kind == 2]
    for scn in sorted(ev.scn.dropna().unique()):
        for sev in sorted(ev.sev[ev.scn == scn].unique()):
            msk = (kind == 2) & (idx.scn.to_numpy() == scn) & (idx.sev.to_numpy() == sev)
            for arm, p in arms.items():
                mp, F0, t0, g, na = randomized_matched_power(p[ho], p[msk])
                rows.append(dict(n=n, seed=seed, scn=scn, sev=float(sev), arm=arm,
                                 value=mp, kindcol="power", excluded=arm in EXCLUDED,
                                 n_eval=int(msk.sum())))
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(d, "summary.csv"), index=False)
    dep = idx[["kind", "rep", "scn", "sev"]].copy(); dep["p_deep"] = p_deep
    dep.to_csv(os.path.join(d, "deep_pvalues.csv"), index=False)

    P = T[(T.kindcol == "power") & (~T.excluded)]
    dp = P[P.arm == "DEEP"].set_index(["scn", "sev"]).value
    c9 = P[P.arm == "CCT-9"].set_index(["scn", "sev"]).value
    bp = P[P.arm.isin(PARTITION)].groupby(["scn", "sev"]).value.max()
    print(f"   size: DEEP {T[(T.arm=='DEEP')&(T.kindcol=='size')].value.iloc[0]:.4f}  "
          f"CCT-9 {T[(T.arm=='CCT-9')&(T.kindcol=='size')].value.iloc[0]:.4f}")
    for blk, names in [("BLOCK A (localised)", ["bump", "step", "inter", "abs"]),
                       ("BLOCK B (index/smooth)", ["indexq", "indexs", "quad", "cubic", "sine"])]:
        i = [k for k in dp.index if k[0] in names]
        print(f"   {blk:24s} vs CCT-9 {(dp[i]-c9[i]).mean():+.4f}   "
              f"vs best-partition {(dp[i]-bp[i]).mean():+.4f}", flush=True)
    if not a.keep_bin:
        os.remove(os.path.join(d, "blk_maps.bin")); print("   block .bin deleted", flush=True)


if __name__ == "__main__":
    main()
