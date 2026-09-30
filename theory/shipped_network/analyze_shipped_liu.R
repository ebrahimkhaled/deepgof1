# analyze_shipped_liu.R -- shipped network vs the per-block network of the benchmark, on the same datasets
# (Settings 2-4, n in {50,100,500}, corpus 101000). Level: rejection at p <= .05 on the 1,000 held-out nulls.
# Power: exact matched level .05 by a randomized threshold on those nulls (as in the paper), read on the
# first 500 alternatives at each C; paired over datasets, so the difference has a paired standard error.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level"
DEP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/release/deepgof1_v2_final/results/liu"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "rows_liu"), "\\.csv$", full.names = TRUE), utils::read.csv))

## matched-level rejection probability for each dataset: 1 below the critical p, gamma at it, 0 above
matched <- function(p0, p1, a = .05) {
  u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0
  g <- (a - below) / (mean(p0 == c))
  ifelse(p1 < c, 1, ifelse(p1 == c, g, 0))
}
out <- list()
for (S in 2:4) for (n in c(50, 100, 500)) {
  blk <- R[R$setting == S & R$n == n, ]
  dp <- utils::read.csv(file.path(DEP, sprintf("S%d_n%d_s101000", S, n), "deep_pvalues.csv"))
  key <- function(d) paste(d$kind, d$rep, ifelse(is.na(d$C), "", format(d$C, nsmall = 2)))
  blk$p_block <- dp$p_deep[match(key(blk), key(dp))]
  N0 <- blk[blk$kind == 0, ]; A <- blk[blk$kind == 2, ]
  if (nrow(N0) < 200 || !nrow(A)) next                  # block not (yet) run
  rs <- matched(N0$p_shipped, A$p_shipped); rb <- matched(N0$p_block, A$p_block); d <- rs - rb
  out[[length(out) + 1]] <- data.frame(setting = S, n = n, nulls = nrow(N0), alts = nrow(A),
    size_shipped = mean(N0$p_shipped <= .05), size_block = mean(N0$p_block <= .05, na.rm = TRUE),
    power_shipped = mean(rs), power_block = mean(rb, na.rm = TRUE), diff = mean(d, na.rm = TRUE),
    paired_se = stats::sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d))), matched_missing = sum(is.na(blk$p_block)))
}
S <- do.call(rbind, out)
print(format(S, digits = 3), row.names = FALSE)
cat(sprintf("\noverall (9 blocks): power shipped %.3f vs per-block %.3f, difference %+.3f (SE %.3f); size shipped %.4f, per-block %.4f\n",
            mean(S$power_shipped), mean(S$power_block), mean(S$diff), sqrt(sum(S$paired_se^2)) / nrow(S),
            mean(S$size_shipped), mean(S$size_block)))
utils::write.csv(S, file.path(HERE, "SHIPPED_LIU.csv"), row.names = FALSE)
