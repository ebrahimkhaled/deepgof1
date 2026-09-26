#!/bin/sh
# smoke_test.sh -- the benchmark pipeline end to end on tiny blocks (about two minutes, two processes).
# Run from the archive root:  sh benchmark/smoke_test.sh
# Needs R with ebrahim.gof (>= 2.8.0) and brglm2, and Python with numpy, pandas, torch, scikit-learn.
# The classical p-values of the regenerated datasets must equal the deposited rows exactly; the
# comparison is printed at the end (max |difference| 0 in every column).
set -e
ROOT=$(pwd)
OUT=$ROOT/smoke_out
mkdir -p "$OUT/results/s3"
cd "$OUT"
# Simulation II (Liu et al. settings): Settings 1 and 3, n = 50, first corpus, a handful of datasets
for S in 1 3; do
  LU_SET=$S LU_N=50 LU_NULL=3 LU_RAND=1 LU_EVAL=1 LU_B=9 LU_W=1 LU_OUT=liu_S$S Rscript "$ROOT/benchmark/liu/liu_grid.R"
done
AU_SET=3 AU_N=50 AU_NULL=3 AU_RAND=1 AU_EVAL=1 AU_W=1 AU_DIR=liu_S3 Rscript "$ROOT/benchmark/liu/liu_augment.R"
python "$ROOT/benchmark/liu/liu_score.py" --dir liu_S1 --keep-bin
# Simulation I (the design grid): all four designs at n = 50, then the pooled rule and the collection
for D in D1 D2 D3 D4; do
  S3_DESIGN=$D S3_N=50 S3_NULL=4 S3_RAND=2 S3_EVAL=2 S3_B=9 S3_W=1 S3_OUT=results/s3/${D}_n50_s101000 Rscript "$ROOT/benchmark/grid/s3_grid.R"
done
DEEPGOF_ARCHIVE=. python "$ROOT/benchmark/grid/s3_score.py" --n 50 --seed 101000
DEEPGOF_ARCHIVE=. python "$ROOT/benchmark/grid/s3_collect.py"
# compare with the deposited rows
python - "$ROOT" <<'PY'
import sys, numpy as np, pandas as pd
root = sys.argv[1]
def cmp(new, dep, keys, cols, tag):
    n, d = pd.read_csv(new), pd.read_csv(dep)
    for f in (n, d):
        if "C" in f: f["_C"] = f["C"].fillna(-1)
    k = keys + (["_C"] if "_C" in n else [])
    m = n.merge(d, on=k, suffixes=("_n", "_d"))
    print(tag, len(m), "rows |", ", ".join("%s %.1e" % (c, float(np.nanmax(np.abs(m[c + "_n"] - m[c + "_d"])))) for c in cols))
P = ["HL", "HL-equalwidth", "Pigeon-Heyse"]
for s in (1, 3):
    cmp("liu_S%d/blk_pvalues.csv" % s, root + "/results/liu/S%d_n50_s101000/blk_pvalues.csv" % s, ["kind", "rep"], P + ["Stukel", "Osius-Rojek"], "Liu S%d" % s)
cmp("liu_S3/augment_sz.csv", root + "/results/liu/S3_n50_s101000/augment_pvalues.csv", ["kind", "rep"], ["SZ"], "Stute-Zhu S3")
for D in ("D2", "D4"):
    cmp("results/s3/%s_n50_s101000/blk_pvalues.csv" % D, root + "/results/s3/%s_n50_s101000/blk_pvalues.csv" % D, ["kind", "rep"], P, "grid " + D)
PY
