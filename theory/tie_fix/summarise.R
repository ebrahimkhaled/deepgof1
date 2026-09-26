# summarise.R -- tables for RESULT.md from the deposited per-unit CSVs.
D <- "/path/to/project/paper_deepgof/theory/tie_fix"
rd <- function(f) { p <- file.path(D, f); if (file.exists(p)) read.csv(p) else NULL }

f <- rd("fdic_roworder.csv")
if (!is.null(f)) {
  cat("== FDIC: stored order (0) and 30 random orders, B = 199, paper seed ==\n")
  for (S in unique(f$state)) for (v in c("old", "new")) {
    z <- f[f$state == S & f$version == v, ]; r <- z[z$order > 0, ]
    cat(sprintf("%s %s  stored: p %.3f S %.3f axes %s | random (%d): p median %.3f range %.3f-%.3f, p<=.05 %d | S range %.3f to %.3f | distinct axes %s\n",
                S, v, z$p[z$order == 0], z$S[z$order == 0], z$axes[z$order == 0], nrow(r),
                median(r$p), min(r$p), max(r$p), sum(r$p <= .05), min(r$S), max(r$S),
                paste(unique(r$axes), collapse = ",")))
  }
}

g <- rd("glow_roworder.csv")
if (!is.null(g)) {
  cat("\n== GLOW, new statistic ==\n")
  for (mn in c("M1", "M2")) {
    a <- g[g$part == "A" & g$model == mn, ]; r <- a[a$order > 0, ]
    cat(sprintf("%s part A: file order p %.3f S %.3f | 30 random orders p median %.3f range %.3f-%.3f, p<=.05 %d\n",
                mn, a$p[a$order == 0], a$S[a$order == 0], median(r$p), min(r$p), max(r$p), sum(r$p <= .05)))
    b1 <- g$p[g$part == "B_fileorder_seeds" & g$model == mn]
    b2 <- g$p[g$part == "B_randomorder_fixedseed" & g$model == mn]
    if (length(b1) && length(b2)) {
      ks <- suppressWarnings(ks.test(b1, b2))
      cat(sprintf("%s part B: file order x %d seeds: median %.3f, p<=.05 %.3f (SE %.3f) | %d random orders: median %.3f, p<=.05 %.3f (SE %.3f) | KS p %.3f\n",
                  mn, length(b1), median(b1), mean(b1 <= .05), sqrt(mean(b1 <= .05) * (1 - mean(b1 <= .05)) / length(b1)),
                  length(b2), median(b2), mean(b2 <= .05), sqrt(mean(b2 <= .05) * (1 - mean(b2 <= .05)) / length(b2)), ks$p.value))
    }
  }
}

cat("\n== Liu Setting 3, n = 100, B = 199, paired old vs new ==\n")
for (cn in c("null", "power", "power2", "sorted")) {
  z <- rd(paste0("liu3_", cn, ".csv")); if (is.null(z)) next
  R <- nrow(z); ro <- z$p_old <= .05; rn <- z$p_new <= .05
  se <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))
  dd <- rn - ro
  b <- sum(ro & !rn); cc <- sum(!ro & rn)
  mc <- if (b + cc > 0) binom.test(cc, b + cc)$p.value else 1
  cat(sprintf("%-6s R = %4d | old %.3f (SE %.3f) | new %.3f (SE %.3f) | new-old %+.3f (paired SE %.3f) | discordant old-only %d new-only %d, exact McNemar p %.3f | KS on p %.3f | cor(p_old,p_new) %.2f\n",
              cn, R, mean(ro), se(ro), mean(rn), se(rn), mean(dd), sd(dd) / sqrt(R), b, cc, mc,
              suppressWarnings(ks.test(z$p_old, z$p_new))$p.value, cor(z$p_old, z$p_new)))
}
