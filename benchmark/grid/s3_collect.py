"""s3_collect.py -- Simulation I, stage 3: collect the 60 block summaries into results/s3/VERDICT.csv.

Keeps the pooled rule (one network per (n, seed), trained on D1 and applied unchanged to D1-D4), which
is the rule the paper reports. Figure 6 (code in figures/fig6_s3grid_sc_v2.py) and the grid mean of
the deployed level (tables/make_tables.py) read VERDICT.csv.

Run from the archive root:  python benchmark/grid/s3_collect.py
"""
import glob, os
import pandas as pd

ARC = os.environ.get("DEEPGOF_ARCHIVE", ".")
f = sorted(glob.glob(os.path.join(ARC, "results", "s3", "*", "summary.csv")))
V = pd.concat([pd.read_csv(x) for x in f], ignore_index=True)
V = V[V.rule == "pooled_D1"]
out = os.path.join(ARC, "results", "s3", os.environ.get("S3_VERDICT_NAME", "VERDICT.csv"))
V.to_csv(out, index=False)
print("blocks:", V.groupby(["n", "design", "seed"]).ngroups, "| rows:", len(V), "->", out)
