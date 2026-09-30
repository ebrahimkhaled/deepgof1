# analyze_ushape.R -- summary of run_ushape.R: size at .05 on the nulls, matched-level power (each test placed
# at exactly .05 on its own nulls at the same n by the randomized threshold of the paper), how often each
# DeepGOF-1 reading's map contains x1, and the paired difference all-pairs minus axis rule.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/ushape"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), "\\.csv$", full.names = TRUE), utils::read.csv))
R <- R[!duplicated(R[, c("id", "rep")]), ]
rej <- function(p0, p1, a = .05) {
  p0 <- p0[!is.na(p0)]; u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0; g <- (a - below) / mean(p0 == c)
  ifelse(is.na(p1), 0, ifelse(p1 < c, 1, ifelse(p1 == c, g, 0)))
}
tests <- c(v1 = "p_v1", allpairs = "p_all", HL = "p_HL", Stukel = "p_Stukel", projection = "p_proj")
out <- list()
for (n in sort(unique(R$n))) {
  nul <- R[R$n == n & R$C == 0, ]
  for (C in sort(unique(R$C[R$C > 0]))) {
    alt <- R[R$n == n & R$C == C, ]
    rr <- sapply(tests, function(t) rej(nul[[t]], alt[[t]]))
    d <- rr[, "allpairs"] - rr[, "v1"]
    out[[length(out) + 1]] <- data.frame(n = n, C = C, reps = nrow(alt),
      t(round(colMeans(rr), 3)), diff_all_v1 = round(mean(d), 3), se_diff = round(stats::sd(d) / sqrt(length(d)), 3),
      v1_picks_x1 = round(mean(alt$v1_x1), 3), all_max_on_x1 = round(mean(alt$all_x1), 3))
  }
  out[[length(out) + 1]] <- data.frame(n = n, C = 0, reps = nrow(nul),
    t(round(sapply(tests, function(t) mean(nul[[t]] <= .05, na.rm = TRUE)), 3)), diff_all_v1 = NA, se_diff = NA,
    v1_picks_x1 = round(mean(nul$v1_x1), 3), all_max_on_x1 = round(mean(nul$all_x1), 3))
}
O <- do.call(rbind, out); O <- O[order(O$n, O$C), ]
cat("C = 0 rows: raw size at .05; C > 0 rows: matched-level power\n")
print(O, row.names = FALSE)
cat("\nNA p-values per test:", colSums(is.na(R[, tests])), "\n")
utils::write.csv(O, file.path(HERE, "USHAPE.csv"), row.names = FALSE)
