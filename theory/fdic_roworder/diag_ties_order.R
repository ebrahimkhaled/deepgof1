# diag_ties_order.R -- are the FDIC model covariates tied, and is the row order related to failure?
# Reads the deposited analysis file only. No test is run here.
PAP <- "/path/to/project/paper_deepgof"
d  <- read.csv(file.path(PAP, "data", "fdic_2008_banks.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(PAP, "data", "fdic_state_counts.csv"), stringsAsFactors = FALSE)
CELLS <- c("GA", setdiff(st$state[st$failures >= 15], "GA"))
V <- c("eq_ta", "npa_ta", "roa", "bro_dep", "cre_ta")

cat("whole file: n =", nrow(d), " sorted by cert (numeric)?", !is.unsorted(d$cert),
    " sorted by cert (character)?", !is.unsorted(as.character(d$cert)), "\n\n")

for (S in CELLS) {
  ds <- d[d$state == S, ]
  n <- nrow(ds); idx <- seq_len(n)
  cat(sprintf("== %s  n = %d  failures = %d\n", S, n, sum(ds$failed)))
  for (v in V) {
    x <- ds[[v]]
    tb <- table(x); tied <- sum(tb[tb > 1])
    top <- if (length(tb)) names(tb)[which.max(tb)] else NA
    cat(sprintf("   %-8s distinct %4d | rows in a tie group %4d (%.1f%%) | largest tie: value %s x %d\n",
                v, length(tb), tied, 100 * tied / n, substr(top, 1, 10), max(tb)))
  }
  # is the row order related to the outcome?
  w  <- suppressWarnings(wilcox.test(idx ~ ds$failed))
  g  <- glm(ds$failed ~ idx, family = binomial())
  cat(sprintf("   row index vs failed: mean index failed %.1f vs survived %.1f | Wilcoxon p %.3g | logit slope z %.2f\n",
              mean(idx[ds$failed == 1]), mean(idx[ds$failed == 0]), w$p.value,
              summary(g)$coefficients[2, 3]))
  # within the largest tie group of each axis, does row order carry y? (what the old map could see)
  for (v in V) {
    x <- ds[[v]]; tb <- table(x); if (max(tb) < 10) next
    g0 <- which(x == as.numeric(names(tb)[which.max(tb)]))
    yy <- ds$failed[g0]
    if (length(unique(yy)) < 2) { cat(sprintf("   %-8s largest tie group (%d rows): all failed = %d\n", v, length(g0), yy[1])); next }
    wg <- suppressWarnings(wilcox.test(seq_along(g0) ~ yy))
    cat(sprintf("   %-8s largest tie group (%d rows, %d failed): order-vs-failed Wilcoxon p %.3g\n",
                v, length(g0), sum(yy), wg$p.value))
  }
}
