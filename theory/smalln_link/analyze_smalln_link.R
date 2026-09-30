# analyze_smalln_link.R -- summary of run_smalln_link.R.
# Level: size at .05 on the logistic nulls, per n and test; the share of datasets with no p-value; separation.
# Power: matched-level power per link and n (each test placed at exactly .05 on its own nulls at the same n by
# the paper's randomized threshold, and evaluated on the datasets where it returns a p-value), and paired
# differences DeepGOF-1 minus Stukel (joint) and minus projection, on datasets where both return a p-value.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/smalln_link"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), "\\.csv$", full.names = TRUE), utils::read.csv))
R <- R[!duplicated(R[, c("id", "rep")]), ]
tests <- c(DG = "p_DG", HL = "p_HL", StukelJ = "p_StukelJ", StukelS = "p_StukelS", GiViTI = "p_GiViTI", proj = "p_proj")
thr <- function(p0, a = .05) {
  p0 <- p0[!is.na(p0)]; u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0
  list(c = c, g = (a - below) / mean(p0 == c))
}
rej <- function(t, p1) ifelse(p1 < t$c, 1, ifelse(p1 == t$c, t$g, 0))

cat("== LEVEL: size at .05 on logistic nulls (share with no p-value in brackets)\n")
L <- list()
for (n in sort(unique(R$n))) {
  z <- R[R$link == "logit" & R$n == n, ]
  L[[length(L) + 1]] <- data.frame(n = n, nulls = nrow(z), separated = round(mean(z$sep == 1, na.rm = TRUE), 3),
    no_fit = round(mean(is.na(z$sep)), 3),
    t(sapply(tests, function(t) sprintf("%.3f [%.2f]", mean(z[[t]] <= .05, na.rm = TRUE), mean(is.na(z[[t]]))))))
}
L <- do.call(rbind, L); print(L, row.names = FALSE)

cat("\n== POWER at matched level .05 (each test on the datasets where it returns a p-value)\n")
P <- list()
for (n in sort(unique(R$n))) {
  z0 <- R[R$link == "logit" & R$n == n, ]
  th <- lapply(tests, function(t) thr(z0[[t]]))
  for (lk in setdiff(unique(R$link), "logit")) {
    z <- R[R$link == lk & R$n == n, ]
    pw <- sapply(names(tests), function(k) { p <- z[[tests[k]]]; ok <- !is.na(p); mean(rej(th[[k]], p[ok])) })
    both <- function(k) { ok <- !is.na(z$p_DG) & !is.na(z[[tests[k]]])
      d <- rej(th$DG, z$p_DG[ok]) - rej(th[[k]], z[[tests[k]]][ok]); c(mean(d), stats::sd(d) / sqrt(sum(ok))) }
    dS <- both("StukelJ"); dP <- both("proj")
    P[[length(P) + 1]] <- data.frame(link = lk, n = n, reps = nrow(z), t(round(pw, 3)),
      DG_minus_StukelJ = sprintf("%+.3f (%.3f)", dS[1], dS[2]), DG_minus_proj = sprintf("%+.3f (%.3f)", dP[1], dP[2]))
  }
}
P <- do.call(rbind, P); P <- P[order(P$link, P$n), ]; print(P, row.names = FALSE)
utils::write.csv(L, file.path(HERE, "LEVEL.csv"), row.names = FALSE)
utils::write.csv(P, file.path(HERE, "POWER.csv"), row.names = FALSE)
