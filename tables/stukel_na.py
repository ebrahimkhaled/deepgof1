"""stukel_na.py -- Online Resource 1, Section 18, Table S12: the share of datasets on which Stukel's test
(SstBoth) returns no p-value, in the extended Liu benchmark (pooled over the two corpora), and the same
share in the simulation with two active covariates among ten.
Run from the archive root:  python tables/stukel_na.py
"""
import glob, os
import pandas as pd

ARC = os.environ.get("DEEPGOF_ARCHIVE", ".")
rows = []
for d in sorted(glob.glob(os.path.join(ARC, "results", "liu", "S*_n*_s*"))):
    tag = os.path.basename(d)
    b = pd.read_csv(os.path.join(d, "blk_pvalues.csv"))
    rows.append(dict(setting=int(tag[1]), n=int(tag.split("_")[1][1:]),
                     null_na=int(b[b.kind == 0].Stukel.isna().sum()), null_N=int((b.kind == 0).sum()),
                     alt_na=int(b[b.kind == 2].Stukel.isna().sum()), alt_N=int((b.kind == 2).sum())))
R = pd.DataFrame(rows).groupby(["setting", "n"]).sum()
half_up = lambda x: (x * 10 + 0.5 + 1e-9) // 1 / 10          # the table rounds halves up
R["nulls_pct"] = half_up(100 * R.null_na / R.null_N)
R["alternatives_pct"] = half_up(100 * R.alt_na / R.alt_N)
print(R[["null_N", "nulls_pct", "alt_N", "alternatives_pct"]].to_string())
tot = R.sum()
print("all: nulls %.1f%%  alternatives %.1f%%" % (100 * tot.null_na / tot.null_N, 100 * tot.alt_na / tot.alt_N))
N = pd.read_csv(os.path.join(ARC, "theory", "niche_paper", "NICHE_ROWS.csv"), keep_default_na=False, na_values=["NA", ""])
print("two active covariates among ten: Stukel (SstBoth) missing in %d of %d datasets" % (int(N.Stukel_marginal.isna().sum()), len(N)))
