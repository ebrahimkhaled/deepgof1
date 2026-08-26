"""s3_score.py -- GATE S3 stage 2: train ONE rule per (n, seed), apply it to all four designs.

Pre-declared in speckit/21_S3_PREDECLARATION.md. Read it before changing anything here.

Section 5 of that document fixes the rule-fitting choice: **one deployment-matched rule per
(n, seed), trained on the D1 corpus and applied UNCHANGED to D1-D4.** That is the honest test of
what the method actually offers -- a pretrained test is one an analyst downloads and applies to
their own design without retraining. Fitting a fresh rule per design would measure a capability the
method does not ship.

`--per-design` additionally trains a rule on each design's own corpus and reports it alongside, so
the COST of pooling is measured rather than assumed. Those numbers are a diagnostic and may never
be the headline (21 section 8).

Everything else is inherited from closed gates and not re-opened here: K=6, mix pinned at 1.0,
randomized matched level via probit_common, NA = a failed test that cannot reject.
"""
import os, sys, argparse, glob
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_train import CNN, train, DEV
from probit_common import randomized_matched_power, narobust_cct, ALPHA

PAP = "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
K = 6
DESIGNS = ["D1", "D2", "D3", "D4"]


@torch.no_grad()
def score(net, X, mu, sd, bs=16384):
    out = np.empty(len(X), dtype=np.float64)
    for k in range(0, len(X), bs):
        out[k:k + bs] = net(torch.from_numpy(((X[k:k + bs] - mu) / sd).astype(np.float32)).to(DEV)).cpu().numpy()
    return out


def load(n, seed, des):
    d = os.path.join(PAP, "results", "s3", f"{des}_n{n}_s{seed}")
    if not os.path.exists(os.path.join(d, "blk_meta.csv")):
        return None
    m = pd.read_csv(os.path.join(d, "blk_meta.csv")).iloc[0]
    idx = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    B, nd = int(m.B), int(m.ndatasets)
    M = np.fromfile(os.path.join(d, "blk_maps.bin"), dtype=np.float32).reshape(nd, B + 1, K * K)
    return dict(dir=d, meta=m, idx=idx, maps=M, B=B, half=int(m.nnull) // 2)


def fit_rule(blk, seed):
    """deployment-matched: null class = bootstrap replicates of the TRAINING datasets (mix = 1.0)."""
    idx, M, half = blk["idx"], blk["maps"], blk["half"]
    kind, rep = idx.kind.to_numpy(), idx.rep.to_numpy()
    tr = ((kind == 0) & (rep <= half)) | (kind == 1)
    boot = M[tr][:, 1:, :].reshape(-1, K * K)
    Xalt = M[kind == 1][:, 0, :]
    rng = np.random.default_rng(17)
    Xnull = boot[rng.choice(len(boot), len(Xalt), replace=False)]
    X = np.vstack([Xnull, Xalt]).astype(np.float32)
    y = np.concatenate([np.zeros(len(Xnull)), np.ones(len(Xalt))])
    mu, sd = X.mean(0), X.std(0) + 1e-8
    torch.manual_seed(7001 + seed % 1000)
    return train(CNN(grid=K).to(DEV), (X - mu) / sd, y), mu, sd


def evaluate(blk, net, mu, sd, tag, n, seed, des):
    idx, M, B, half = blk["idx"], blk["maps"], blk["B"], blk["half"]
    kind, rep = idx.kind.to_numpy(), idx.rep.to_numpy()
    ho = (kind == 0) & (rep > half)
    ev = kind == 2
    TESTS = [c for c in idx.columns if c not in ("kind", "rep", "sep", "cct", "axes")]
    S = score(net, M.reshape(-1, K * K), mu, sd).reshape(len(idx), B + 1)
    p = (1 + (S[:, 1:] >= S[:, [0]]).sum(1)) / (B + 1)
    s_t = score(net, M[(kind == 0) & (rep <= half)][:20000, 0, :], mu, sd)
    s_b = score(net, M[ho][:, 1:, :].reshape(-1, K * K)[:20000], mu, sd)
    auc = float(max(roc_auc_score(np.r_[np.zeros(len(s_t)), np.ones(len(s_b))], np.r_[s_t, s_b]),
                    1 - roc_auc_score(np.r_[np.zeros(len(s_t)), np.ones(len(s_b))], np.r_[s_t, s_b])))
    arms = {"DEEP": p, "CCT-11": narobust_cct(idx[TESTS].to_numpy(float))}
    arms.update({t: idx[t].to_numpy() for t in TESTS})
    rows = []
    for arm, q in arms.items():
        mp, F0, t0, g, na = randomized_matched_power(q[ho], q[ev])
        rows.append(dict(n=n, seed=seed, design=des, rule=tag, arm=arm, matched_power=mp,
                         raw_size=float((np.nan_to_num(q[ho], nan=1.0) <= ALPHA).mean()),
                         raw_power=float((np.nan_to_num(q[ev], nan=1.0) <= ALPHA).mean()),
                         score_auc=auc if arm == "DEEP" else np.nan, na_rate=na,
                         firth_rate=float(idx.sep.mean()), n_null=int(ho.sum()), n_alt=int(ev.sum())))
    return rows, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--per-design", action="store_true", help="diagnostic: also fit a rule per design")
    ap.add_argument("--keep-bin", action="store_true")
    a = ap.parse_args()

    blks = {d: load(a.n, a.seed, d) for d in DESIGNS}
    missing = [d for d, b in blks.items() if b is None]
    if missing:
        sys.exit(f"missing blocks for n={a.n} seed={a.seed}: {missing}")
    print(f"[S3 n={a.n} seed={a.seed}] blocks " +
          " ".join(f"{d}:{len(blks[d]['idx'])}" for d in DESIGNS), flush=True)

    # THE pre-declared rule: trained on D1, applied unchanged everywhere
    net, mu, sd = fit_rule(blks["D1"], a.seed)
    rows, deps = [], {}
    for d in DESIGNS:
        r, p = evaluate(blks[d], net, mu, sd, "pooled_D1", a.n, a.seed, d)
        rows += r; deps[d] = {"p_deep_pooled": p}
        s = {x["arm"]: x for x in r}
        print(f"   {d}  size {s['DEEP']['raw_size']:.4f}  DEEP {s['DEEP']['matched_power']:.4f}  "
              f"CCT {s['CCT-11']['matched_power']:.4f}  gain {s['DEEP']['matched_power']-s['CCT-11']['matched_power']:+.4f}"
              f"  AUC {s['DEEP']['score_auc']:.4f}  Firth {100*s['DEEP']['firth_rate']:.2f}%", flush=True)

    if a.per_design:                       # diagnostic only -- never the headline (21 section 8)
        for d in DESIGNS:
            netd, mud, sdd = fit_rule(blks[d], a.seed)
            r, p = evaluate(blks[d], netd, mud, sdd, "per_design", a.n, a.seed, d)
            rows += r; deps[d]["p_deep_perdesign"] = p
            s = {x["arm"]: x for x in r}
            print(f"   {d}* per-design rule: size {s['DEEP']['raw_size']:.4f}  "
                  f"DEEP {s['DEEP']['matched_power']:.4f}", flush=True)

    for d in DESIGNS:
        out = blks[d]["idx"][["kind", "rep"]].copy()
        for k, v in deps[d].items(): out[k] = v
        out.to_csv(os.path.join(blks[d]["dir"], "deep_pvalues.csv"), index=False)
        pd.DataFrame([r for r in rows if r["design"] == d]).to_csv(
            os.path.join(blks[d]["dir"], "summary.csv"), index=False)
        if not a.keep_bin:
            b = os.path.join(blks[d]["dir"], "blk_maps.bin")
            if os.path.exists(b): os.remove(b)
    print("   scored; block .bin files deleted (p-values + summaries retained)", flush=True)


if __name__ == "__main__":
    main()
