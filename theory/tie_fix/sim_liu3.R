# sim_liu3.R -- old vs new statistic, paired, on Liu et al. (2024) Setting 3 (x2 ~ Bin(.5), so the
# second axis is tied), n = 100, B = 199. Each replicate: one data set, both versions, same seed.
#   COND = null    C = 0,   rows in generation (random) order   -> size, must agree
#   COND = power   C = 0.5, rows in generation (random) order   -> power, must agree
#   COND = power2  C = 1.0, rows in generation (random) order   -> power at a larger departure
#   COND = sorted  C = 0,   rows SORTED by y                    -> old should be liberal, new not
# Usage: Rscript sim_liu3.R COND REPS. Resumable: one row per replicate in liu3_<COND>.csv.
source("/path/to/project/paper_deepgof/theory/tie_fix/common.R")
a <- commandArgs(TRUE); COND <- a[1]; REPS <- as.integer(a[2])
C <- switch(COND, power = 0.5, power2 = 1.0, 0)
n <- 100L
OUT <- file.path(OUTD, paste0("liu3_", COND, ".csv")); dk <- done_keys(OUT)
base <- c(null = 510000L, power = 520000L, sorted = 530000L, power2 = 540000L)[[COND]]
for (i in seq_len(REPS)) {
  key <- paste(COND, i, sep = "_"); if (key %in% dk) next
  set.seed(base + i)
  x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5)
  y <- rbinom(n, 1, plogis(-1.972 + 0.135 * x1 + 1.117 * x2 + C * x1 * x2))
  dd <- data.frame(y, x1, x2)
  if (COND == "sorted") dd <- dd[order(dd$y), ]
  fit <- suppressWarnings(glm(y ~ x1 + x2, data = dd, family = binomial()))
  set.seed(base + 100000L + i); ro <- dg_old(fit)
  set.seed(base + 100000L + i); rn <- dg_new(fit)
  append_row(OUT, data.frame(key = key, rep = i, p_old = ro$p.value, p_new = rn$p.value,
                             S_old = ro$statistic, S_new = rn$statistic,
                             ax_old = paste(ro$axes, collapse = "|"), ax_new = paste(rn$axes, collapse = "|")))
  if (i %% 50 == 0) { cat(COND, i, "\n"); flush.console() }
}
cat("LIU3", COND, "COMPLETE\n")
