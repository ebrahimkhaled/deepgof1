# analyze_smalln_power.R -- matched-level power of run_smalln_power.R. The thresholds come from the logistic
# nulls of run_smalln_link.R at the same n (the same design with C = 0), by the paper's randomized rule; each
# test is evaluated on the datasets where it returns a p-value. Paired differences DeepGOF-1 minus each rival
# are taken on the datasets where both return a p-value.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/smalln_link"
rd <- function(d) { R <- do.call(rbind, lapply(list.files(file.path(HERE, d), "\\.csv$", full.names = TRUE), utils::read.csv))
  R[!duplicated(R[, c("id", "rep")]), ] }
N <- rd("rows"); N <- N[N$link == "logit", ]
A <- rd("rows_power")
tests <- c(DG = "p_DG", HL = "p_HL", StukelJ = "p_StukelJ", StukelS = "p_StukelS", GiViTI = "p_GiViTI", proj = "p_proj")
thr <- function(p0, a = .05) {
  p0 <- p0[!is.na(p0)]; u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0
  list(c = c, g = (a - below) / mean(p0 == c))
}
rej <- function(t, p1) ifelse(p1 < t$c, 1, ifelse(p1 == t$c, t$g, 0))
out <- list()
for (n in sort(unique(A$n))) {
  th <- lapply(tests, function(t) thr(N[N$n == n, t]))
  for (fam in c("quad", "inter")) for (C in sort(unique(A$C))) {
    z <- A[A$n == n & A$fam == fam & A$C == C, ]
    pw <- sapply(names(tests), function(k) { p <- z[[tests[k]]]; mean(rej(th[[k]], p[!is.na(p)])) })
    dd <- sapply(setdiff(names(tests), "DG"), function(k) {
      ok <- !is.na(z$p_DG) & !is.na(z[[tests[k]]])
      d <- rej(th$DG, z$p_DG[ok]) - rej(th[[k]], z[[tests[k]]][ok])
      sprintf("%+.3f (%.3f)", mean(d), stats::sd(d) / sqrt(sum(ok))) })
    out[[length(out) + 1]] <- data.frame(fam = fam, C = C, n = n, reps = nrow(z), sep = round(mean(z$sep == 1, na.rm = TRUE), 3),
      t(round(pw, 3)), t(stats::setNames(dd, paste0("DG_minus_", names(dd)))), check.names = FALSE)
  }
}
O <- do.call(rbind, out); O <- O[order(O$fam, O$C, O$n), ]
print(O[, c("fam", "C", "n", "reps", "sep", names(tests))], row.names = FALSE)
cat("\npaired differences, DeepGOF-1 minus each test (SE):\n")
print(O[, c("fam", "C", "n", grep("^DG_minus", names(O), value = TRUE))], row.names = FALSE)
utils::write.csv(O, file.path(HERE, "POWER_COVARIATE_MISFIT.csv"), row.names = FALSE)
