# analyze_ablation.R -- matched-level power (each statistic placed at exactly .05 on its own nulls of the same
# design and n, by the paper's randomized threshold) for the five statistics of run_ablation.R, per cell and
# averaged per design; paired differences with standard errors; raw size per design and n.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/ablation"
R <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), "\\.csv$", full.names = TRUE), utils::read.csv))
R$fam[is.na(R$fam)] <- "null"
R <- R[!duplicated(R[, c("design", "S", "fam", "C", "n", "rep")]), ]
ST <- c(axes = "p_axes", allpairs = "p_allpairs", combined = "p_combined", ss = "p_ss", max = "p_max")
thr <- function(p0, a = .05) {
  p0 <- p0[!is.na(p0)]; u <- sort(unique(p0)); F <- vapply(u, function(c) mean(p0 <= c), 0)
  k <- which(F >= a)[1]; c <- u[k]; below <- if (k > 1) F[k - 1] else 0
  list(c = c, g = (a - below) / mean(p0 == c))
}
rej <- function(t, p) ifelse(is.na(p), NA, ifelse(p < t$c, 1, ifelse(p == t$c, t$g, 0)))
cells <- list(); sizes <- list()
for (des in unique(R$design)) for (S in unique(R$S[R$design == des])) for (n in unique(R$n[R$design == des])) {
  z0 <- R[R$design == des & R$S == S & R$n == n & R$fam == "null", ]
  if (!nrow(z0)) next
  sizes[[length(sizes) + 1]] <- data.frame(design = des, S = S, n = n, nulls = nrow(z0),
    t(round(sapply(ST, function(k) mean(z0[[k]] <= .05, na.rm = TRUE)), 3)))
  th <- lapply(ST, function(k) thr(z0[[k]]))
  alts <- unique(R[R$design == des & R$S == S & R$n == n & R$fam != "null", c("fam", "C")])
  for (a in seq_len(nrow(alts))) {
    z <- R[R$design == des & R$S == S & R$n == n & R$fam == alts$fam[a] & R$C == alts$C[a], ]
    rr <- sapply(names(ST), function(k) rej(th[[k]], z[[ST[k]]]))
    d <- function(x, y) { v <- rr[, x] - rr[, y]; v <- v[!is.na(v)]; c(mean(v), stats::sd(v) / sqrt(length(v))) }
    dss <- d("axes", "ss"); dmx <- d("axes", "max"); dca <- d("combined", "axes"); dcp <- d("combined", "allpairs")
    cells[[length(cells) + 1]] <- data.frame(design = des, S = S, n = n, fam = alts$fam[a], C = alts$C[a], reps = nrow(z),
      t(round(colMeans(rr, na.rm = TRUE), 3)),
      net_minus_ss = round(dss[1], 3), se_ss = round(dss[2], 3), net_minus_max = round(dmx[1], 3), se_max = round(dmx[2], 3),
      comb_minus_axes = round(dca[1], 3), se_ca = round(dca[2], 3), comb_minus_all = round(dcp[1], 3), se_cp = round(dcp[2], 3),
      axes_hit = round(mean(z$axes_hit, na.rm = TRUE), 3), top_hit = round(mean(z$top_hit, na.rm = TRUE), 3))
  }
}
C <- do.call(rbind, cells); Z <- do.call(rbind, sizes)
cat("== size at .05 (raw) by design, setting and n\n"); print(Z, row.names = FALSE)
cat("\n== matched-level power per cell\n"); print(C, row.names = FALSE)
cat("\n== mean matched power per design (each cell weight 1)\n")
M <- aggregate(cbind(axes, allpairs, combined, ss, max) ~ design, C, mean); print(M, digits = 3, row.names = FALSE)
cat("\n== cells where the network (axis rule) is ahead of / behind the sum of squares by more than two SE\n")
print(table(design = C$design, vs_ss = ifelse(C$net_minus_ss > 2 * C$se_ss, "ahead", ifelse(C$net_minus_ss < -2 * C$se_ss, "behind", "level"))))
print(table(design = C$design, vs_max = ifelse(C$net_minus_max > 2 * C$se_max, "ahead", ifelse(C$net_minus_max < -2 * C$se_max, "behind", "level"))))
cat("\n== combined against the better of the two readings in each cell\n")
C$best <- pmax(C$axes, C$allpairs); C$comb_gap <- round(C$combined - C$best, 3)
print(aggregate(cbind(combined, best, comb_gap) ~ design, C, mean), digits = 3, row.names = FALSE)
utils::write.csv(C, file.path(HERE, "ABLATION_POWER.csv"), row.names = FALSE)
utils::write.csv(Z, file.path(HERE, "ABLATION_SIZE.csv"), row.names = FALSE)
