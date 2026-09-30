"""Raw rejection at p <= .05 against size-adjusted (matched-level) power, for the link and covariate-misfit
studies. Size adjustment: each test's threshold is set on its own null datasets at the same n so its size is
exactly .05 (randomized at the boundary p-value), the rule used throughout the paper."""
import glob, numpy as np, pandas as pd
H = r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/smalln_link/"
rd = lambda d: pd.concat([pd.read_csv(f) for f in glob.glob(H + d + "/*.csv")]).drop_duplicates(["id", "rep"])
N = rd("rows"); L = N[N.link != "logit"]; N = N[N.link == "logit"]; C = rd("rows_power")
T = {"DG": "p_DG", "proj": "p_proj", "StukelJ": "p_StukelJ", "GiViTI": "p_GiViTI", "HL": "p_HL"}


def thr(p0, a=.05):
    p0 = np.sort(p0[~np.isnan(p0)]); u = np.unique(p0); F = np.array([np.mean(p0 <= c) for c in u])
    k = np.argmax(F >= a); c = u[k]; below = F[k - 1] if k > 0 else 0
    return c, (a - below) / np.mean(p0 == c)


def adj(p, t):
    p = p[~np.isnan(p)]; c, g = t
    return np.mean(np.where(p < c, 1, np.where(p == c, g, 0)))


def raw(p):
    p = p[~np.isnan(p)]; return np.mean(p <= .05)


rows = []
for n in sorted(N.n.unique()):
    th = {k: thr(N[N.n == n][v].to_numpy(float)) for k, v in T.items()}
    size = {k: raw(N[N.n == n][v].to_numpy(float)) for k, v in T.items()}
    rows.append(dict(study="null (size)", cell="logit", n=n, **{f"{k} raw": round(size[k], 3) for k in T}))
    for lab, z in [(lk, L[(L.link == lk) & (L.n == n)]) for lk in ["probit", "cloglog", "stukel_asym", "stukel_heavy"]] + \
                  [(f"{f} C={c}", C[(C.fam == f) & (C.C == c) & (C.n == n)]) for f in ["quad", "inter"] for c in [0.5, 1.0]]:
        if len(z) == 0:
            continue
        r = dict(study="link" if lab[0] in "pcs" else "covariate", cell=lab, n=n)
        for k, v in T.items():
            p = z[v].to_numpy(float)
            r[f"{k} raw"] = round(raw(p), 3); r[f"{k} adj"] = round(adj(p, th[k]), 3)
        rows.append(r)
D = pd.DataFrame(rows)
pd.set_option("display.width", 250)
cols = ["study", "cell", "n"] + [f"{k} {w}" for k in T for w in ("raw", "adj")]
D = D.reindex(columns=cols)
print(D.to_string(index=False, na_rep=""))
D.to_csv(H + "RAW_VS_ADJUSTED.csv", index=False)
