"""Would adding the small-n and link-misfit cells change the ranking? Five tests common to both sources:
DeepGOF-1 (shipped network, default reading), projection, Stukel (SstBoth, as in Table 3), GiViTI, HL.
  A = the benchmark of Liu et al., Settings 2-4, n 50/100/500, same datasets (RANK_SHIPPED_LIU.csv): 9 blocks
  B = the new study: 12 covariate-misfit cells (n 20/30/50) and 16 link cells (n 20/30/50/100)
Scores, all at matched level .05: mean power over cells (each cell weight 1); mean rank within cell (1 = best);
and the share of cells in which the test's own size at .05 lies in the band [.035, .065]."""
import pandas as pd, numpy as np
H = r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/"
TESTS = ["DeepGOF-1", "projection", "Stukel", "GiViTI", "HL"]

a = pd.read_csv(H + "shipped_level/RANK_SHIPPED_LIU.csv")
mapA = {"DeepGOF-1 shipped": "DeepGOF-1", "projection": "projection", "Stukel": "Stukel", "GiViTI": "GiViTI", "HL": "HL",
        "DeepGOF-1 all-pairs": "DeepGOF-1 all-pairs"}
a = a[a.test.isin(mapA)].assign(test=lambda d: d.test.map(mapA))
A = a.pivot_table(index=["setting", "n"], columns="test", values="power").reset_index()
A_size = a.pivot_table(index=["setting", "n"], columns="test", values="size").reset_index()
A["source"] = "A benchmark"; A["cell"] = "S" + A.setting.astype(str) + " n" + A.n.astype(str)
A_size["cell"] = "S" + A_size.setting.astype(str) + " n" + A_size.n.astype(str)

lv = pd.read_csv(H + "smalln_link/LEVEL.csv")
def num(s): return float(str(s).split()[0])
size_new = {int(r.n): {"DeepGOF-1": num(r.DG), "projection": num(r.proj), "Stukel": num(r.StukelS),
                       "GiViTI": num(r.GiViTI), "HL": num(r.HL)} for r in lv.itertuples()}
ren = {"DG": "DeepGOF-1", "proj": "projection", "StukelS": "Stukel", "GiViTI": "GiViTI", "HL": "HL"}
lk = pd.read_csv(H + "smalln_link/POWER.csv").rename(columns=ren)
lk["source"] = "B link"; lk["cell"] = lk.link + " n" + lk.n.astype(str)
cv = pd.read_csv(H + "smalln_link/POWER_COVARIATE_MISFIT.csv").rename(columns=ren)
cv["source"] = "B covariate"; cv["cell"] = cv.fam + " C" + cv.C.astype(str) + " n" + cv.n.astype(str)
B = pd.concat([lk, cv], ignore_index=True)
B["DeepGOF-1 all-pairs"] = B["DeepGOF-1"]          # two covariates: the two readings are the same test


def score(D, S_of, label, dg="DeepGOF-1"):
    cols = [dg] + TESTS[1:]
    P = D[cols].to_numpy(float)
    ranks = pd.DataFrame(P).rank(axis=1, ascending=False).to_numpy()
    inband = np.array([[.035 <= S_of(i, t) <= .065 for t in cols] for i in range(len(D))])
    out = pd.DataFrame({"test": ["DeepGOF-1"] + TESTS[1:], "mean power": P.mean(0).round(3),
                        "mean rank": ranks.mean(0).round(2), "size in band": inband.mean(0).round(2)})
    out = out.sort_values("mean power", ascending=False)
    print(f"\n== {label}  ({len(D)} cells)")
    print(out.to_string(index=False))


sA = lambda D: (lambda i, t: float(A_size.set_index("cell").loc[D.cell.iloc[i], t]))
sB = lambda D: (lambda i, t: size_new[int(D.n.iloc[i])]["DeepGOF-1" if t.startswith("DeepGOF") else t])
score(A, sA(A), "A: benchmark only, default reading")
score(A, sA(A), "A: benchmark only, all-pairs reading", dg="DeepGOF-1 all-pairs")
score(B[B.source == "B covariate"], sB(B[B.source == "B covariate"]), "B: small-n covariate misfit only")
score(B[B.source == "B link"], sB(B[B.source == "B link"]), "B: link misfit only")
score(B, sB(B), "B: all new cells")
AB = pd.concat([A.assign(src="A"), B.assign(src="B")], ignore_index=True)
def sAB(i, t):
    r = AB.iloc[i]
    if r.src == "A":
        return float(A_size.set_index("cell").loc[r.cell, t])
    return size_new[int(r.n)]["DeepGOF-1" if t.startswith("DeepGOF") else t]
score(AB, sAB, "A + B pooled, default reading")
score(AB, sAB, "A + B pooled, all-pairs reading", dg="DeepGOF-1 all-pairs")
