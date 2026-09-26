## 12_geom.R -- geometry of the Bernstein score confidence set at n = 100, D1.
## No pB evaluations; this is pure set geometry, so it is cheap.
##   * is C_gamma bounded? (radial bisection in 200 random directions, cap 1e4 Wald SEs)
##   * how big is it relative to the Wald ellipsoid?
##   * what fraction of an 11^3 grid over its bounding box lands inside it?

EX <- "/path/to/project/paper_deepgof/theory/exact_mmc"
source(file.path(EX, "10_core.R"))

ND <- 40L
res <- vector("list", ND)
set.seed(20260916)
for (i in seq_len(ND)) {
  g <- gen_dataset(0L, i)
  P <- prep_D1(g$data); y <- as.numeric(g$data$y)
  fo <- fast_fit(y, P$mm); if (is.null(fo) || is.null(fo$ph)) next
  th <- as.numeric(fo$coef)
  Iinv <- solve(crossprod(P$mm, fo$ph * (1 - fo$ph) * P$mm))
  se <- sqrt(diag(Iinv))
  ## radial extent in 200 random directions, measured in units of the Wald SE
  Dir <- matrix(rnorm(200 * 3), 200, 3); Dir <- Dir / sqrt(rowSums(Dir^2))
  Dir <- Dir * rep(se, each = 200)                 # isotropic in SE units
  rad <- numeric(200); cap <- 1e4
  for (k in 1:200) {
    lo <- 0; hi <- 1
    while (hi < cap && in_set_slack(matrix(th + hi * Dir[k, ], 1L), y, P$mm) <= 1) { lo <- hi; hi <- hi * 2 }
    if (hi >= cap) { rad[k] <- cap; next }
    for (it in 1:40) { mid <- (lo + hi) / 2
      if (in_set_slack(matrix(th + mid * Dir[k, ], 1L), y, P$mm) <= 1) lo <- mid else hi <- mid }
    rad[k] <- lo
  }
  bx <- set_box(y, P$mm, th, fo$ph, 7000L + i, nmc = 20000L)
  ax <- lapply(1:3, function(k) seq(bx$lo[k], bx$hi[k], length.out = NGRID))
  G <- as.matrix(expand.grid(ax[[1]], ax[[2]], ax[[3]]))
  nin <- sum(in_set_slack(G, y, P$mm) <= 1)
  res[[i]] <- c(i = i, rmin = min(rad), rmed = median(rad), rmax = max(rad),
                runbdd = mean(rad >= cap),
                h1 = bx$half[1]/se[1], h2 = bx$half[2]/se[2], h3 = bx$half[3]/se[3],
                F = bx$F, nmcin = bx$nin, ngridin = nin,
                theta0in = as.numeric(in_set_slack(matrix(g$theta0, 1L), y, P$mm) <= 1))
}
R <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(R, file.path(EX, "geom_probe.csv"), row.names = FALSE)
cat(sprintf("datasets: %d\n", nrow(R)))
cat(sprintf("radial extent (Wald-SE units)  min %.1f  median %.1f  max %.1f  (over 200 dirs/dataset)\n",
            median(R[, "rmin"]), median(R[, "rmed"]), median(R[, "rmax"])))
cat(sprintf("directions hitting the 1e4 cap: mean %.4f, datasets with any: %d\n",
            mean(R[, "runbdd"]), sum(R[, "runbdd"] > 0)))
cat(sprintf("box half-width / Wald SE   median  %.1f  %.1f  %.1f\n",
            median(R[, "h1"]), median(R[, "h2"]), median(R[, "h3"])))
cat(sprintf("11^3 grid points inside C_gamma: median %.0f, min %.0f, max %.0f\n",
            median(R[, "ngridin"]), min(R[, "ngridin"]), max(R[, "ngridin"])))
cat(sprintf("MC points inside / 20000: median %.0f\n", median(R[, "nmcin"])))
cat(sprintf("theta0 in C_gamma: %d / %d\n", sum(R[, "theta0in"]), nrow(R)))
