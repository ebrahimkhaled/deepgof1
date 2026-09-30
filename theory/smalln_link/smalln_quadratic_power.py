"""Matched-level power at n = 30 and 50 against an omitted quadratic, from the clean cells of the
corrupted-record Study 1 (theory/small_n/rows; eta = -0.5 + x1 + 0.6 x2 + C (x1^2 - 1), no corruption).
Each test is randomized to exactly .05 on its own clean nulls at the same n (the paper's rule)."""
import glob, numpy as np, pandas as pd
R = pd.concat([pd.read_csv(f) for f in glob.glob(r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/small_n/rows/*.csv")])
R = R[R.type == "none"]
T = {"DeepGOF-1": "p_DG", "HL": "p_HL", "Stukel (joint)": "p_Stukel"}


def thr(p0, a=.05):
    p0 = np.sort(p0[~np.isnan(p0)]); u = np.unique(p0); F = np.array([np.mean(p0 <= c) for c in u])
    k = np.argmax(F >= a); c = u[k]; below = F[k - 1] if k > 0 else 0
    return c, (a - below) / np.mean(p0 == c)


rows = []
for n in (30, 50):
    nul = R[(R.n == n) & (R.C == 0)]
    th = {k: thr(nul[v].to_numpy(float)) for k, v in T.items()}
    for C in (0.5, 1.0, 1.5):
        alt = R[(R.n == n) & (R.C == C)]
        row = {"n": n, "C": C, "reps": len(alt)}
        for k, v in T.items():
            p = alt[v].to_numpy(float); p = p[~np.isnan(p)]; c, g = th[k]
            row[k] = round(float(np.mean(np.where(p < c, 1, np.where(p == c, g, 0)))), 3)
        rows.append(row)
print(pd.DataFrame(rows).to_string(index=False))
