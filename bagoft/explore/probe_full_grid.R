# =============================================================================
# probe_full_grid.R -- smoke test before extending the EXPLORATORY trial to the
# whole Liu grid (settings 1-4, n in {50,100,500}, C in {0,.25,.5,.75,1}).
#
# *** EXPLORATORY, NOT CONFIRMATORY. ***
#
# The first trial only ran settings 2 and 3 at n = 50 and 100. This checks the
# code paths it never touched -- setting 1 (ONE covariate), setting 4 (no
# intercept in the truth), and n = 500 -- with a deliberately cheap BAGofT
# (nsplits = 4, nsim = 4, ntree = 60). These p-values mean nothing; the point is
# whether the call errors, and how long it takes. A full-default dataset does
# nsplits * (nsim + 1) = 100 * 101 = 10,100 partition fits against 4 * 5 = 20
# here, so the timing is scaled by 505 to estimate one real dataset.
# =============================================================================
source("/path/to/project/paper_deepgof/explore_bagoft/common.R")

SEED0 <- 101000L
SCALE <- (100 * 101) / (4 * 5)
probe <- expand.grid(SET = 1:4, n = c(50L, 500L), KEEP.OUT.ATTRS = FALSE)

out <- do.call(rbind, lapply(seq_len(nrow(probe)), function(k) {
  SET <- probe$SET[k]; n <- probe$n[k]
  dd <- liu_dataset(SEED0, SET, n, 1.0, 1L)          # rep 1 of the C = 1 cell
  if (is.null(dd$data)) return(data.frame(SET, n, status = "degenerate_y", p = NA, secs = NA, est_full_s = NA))
  set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)
  r <- bagoft_p(dd$data, FRM(dd$data), nsplits = 4L, nsim = 4L, ntree = 60L)
  data.frame(SET, n, covariates = ncol(dd$data) - 1L,
             status = if (is.na(r$p)) "ERROR_OR_NA" else "ok",
             p = round(r$p, 3), stat_min = signif(r$stat_min, 3),
             secs = round(r$time_s, 2), est_full_s = round(r$time_s * SCALE))
}))
print(out, row.names = FALSE)
cat("\nmean estimated full-default seconds per dataset: n=50",
    round(mean(out$est_full_s[out$n == 50], na.rm = TRUE)), "| n=500",
    round(mean(out$est_full_s[out$n == 500], na.rm = TRUE)), "\n")
