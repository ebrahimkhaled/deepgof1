# fdic_roworder.R -- the paper's FDIC adequacy calls (Table tab:fdic, column DG-1 observed p) re-run
# in the stored row order and in 30 random row orders, with the old (row-order ties) and the new
# (random ties) statistic. Same fits, same seed as code/fdic_deepgof.R (771000 + cell index), B = 199.
# Resumable: one row per (cell, order, version) in fdic_roworder.csv.
source("/path/to/project/paper_deepgof/theory/tie_fix/common.R")
d  <- read.csv(file.path(PAP, "data", "fdic_2008_banks.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(PAP, "data", "fdic_state_counts.csv"), stringsAsFactors = FALSE)
CELLS <- c("GA", setdiff(st$state[st$failures >= 15], "GA"))
FRM <- failed ~ eq_ta + npa_ta + roa + bro_dep + cre_ta
SEED0 <- 771000L
OUT <- file.path(OUTD, "fdic_roworder.csv")
dk <- done_keys(OUT)
for (S in CELLS) {
  ds0 <- d[d$state == S, ]
  set.seed(20260925L + match(S, CELLS))
  perms <- replicate(30, sample(nrow(ds0)), simplify = FALSE)
  for (k in 0:30) {
    ds <- if (k == 0) ds0 else ds0[perms[[k]], ]
    fit <- suppressWarnings(glm(FRM, data = ds, family = binomial()))
    for (v in c("old", "new")) {
      key <- paste(S, k, v, sep = "_"); if (key %in% dk) next
      set.seed(SEED0 + match(S, CELLS))
      r <- if (v == "old") dg_old(fit) else dg_new(fit)
      append_row(OUT, data.frame(key = key, state = S, order = k, version = v,
                                 axes = paste(r$axes, collapse = "|"), S = r$statistic, p = r$p.value))
    }
  }
  cat(S, "done\n"); flush.console()
}
cat("FDIC ROWORDER COMPLETE\n")
