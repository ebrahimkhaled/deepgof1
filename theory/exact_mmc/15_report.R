## 15_report.R -- every number that goes into FEASIBILITY.md, printed.

EX <- "/path/to/project/paper_deepgof/theory/exact_mmc"
ALPHA <- 0.05; GAMMA <- 0.01

fs <- list.files(file.path(EX, "rows"), pattern = "^k\\d+_r\\d+\\.csv$", full.names = TRUE)
A <- do.call(rbind, lapply(fs, function(f) tryCatch(utils::read.csv(f), error = function(e) NULL)))
A <- A[!duplicated(A[, c("kind", "rep")]), ]
A <- A[order(A$kind, A$rep), ]
utils::write.csv(A, file.path(EX, "per_dataset.csv"), row.names = FALSE)

## use a BALANCED prefix of the two arms, so the arms are comparable
mx <- min(max(A$rep[A$kind == 0]), max(A$rep[A$kind == 2]))
nn <- sort(A$rep[A$kind == 0]); na <- sort(A$rep[A$kind == 2])
keep <- intersect(nn[nn <= mx], na[na <= mx])
N0 <- A[A$kind == 0 & A$rep %in% keep, ]; N1 <- A[A$kind == 2 & A$rep %in% keep, ]

cp <- function(x, n) { if (n == 0) return(c(NA, NA)); ci <- binom.test(x, n)$conf.int; as.numeric(ci) }
ln <- function(lab, x, n) cat(sprintf("  %-34s %3d/%3d = %.4f   [%.4f, %.4f]\n",
                                      lab, x, n, x / n, cp(x, n)[1], cp(x, n)[2]))
## paired difference of two correlated proportions (discordant-pair Wald CI)
pdiff <- function(r1, r2, lab) {
  n <- length(r1); b <- sum(r1 & !r2); c0 <- sum(!r1 & r2); d <- (b - c0) / n
  se <- sqrt((b + c0 - (b - c0)^2 / n)) / n
  cat(sprintf("  %-34s %+.4f   [%+.4f, %+.4f]   (discordant %d vs %d of %d)\n",
              lab, d, d - 1.96 * se, d + 1.96 * se, b, c0, n))
}
qs <- function(v) sprintf("median %.1f  q90 %.1f  max %.1f", median(v), quantile(v, .9), max(v))
sf <- function(f, v) if (!length(v)) NA_real_ else f(v)

cat("=================================================================================\n")
cat(sprintf("exact_mmc | D1, p=2, indep, 50%% events, n=%d | B=%d | gamma=%.2f | alpha=%.2f\n",
            A$n[1], A$B[1], GAMMA, ALPHA))
cat(sprintf("rows on disk: null %d, alt %d ; balanced paired prefix used: %d per arm\n",
            sum(A$kind == 0), sum(A$kind == 2), length(keep)))

cat("\n--- 1. the Bernstein confidence set -----------------------------------------------\n")
rel <- cbind(A$box_half_0 / A$se_0, A$box_half_1 / A$se_1, A$box_half_2 / A$se_2)
cat(sprintf("  bounding-box half-width / Wald SE, median over %d datasets: %.1f  %.1f  %.1f\n",
            nrow(A), median(rel[, 1]), median(rel[, 2]), median(rel[, 3])))
cat(sprintf("  ... q90:  %.1f  %.1f  %.1f      max: %.1f  %.1f  %.1f\n",
            quantile(rel[, 1], .9), quantile(rel[, 2], .9), quantile(rel[, 3], .9),
            max(rel[, 1]), max(rel[, 2]), max(rel[, 3])))
cat(sprintf("  truncation at %d Wald SE binds (C_gamma numerically unbounded): %d / %d = %.3f\n",
            200L, sum(A$box_unbdd), nrow(A), mean(A$box_unbdd)))
cat(sprintf("  theta0 inside C_gamma (null arm only): %d / %d\n",
            sum(A$theta0_in[A$kind == 0], na.rm = TRUE), sum(!is.na(A$theta0_in))))
cat(sprintf("  11^3 = 1331 grid points inside C_gamma: %s ; zero for %d datasets\n",
            qs(A$n_grid_in), sum(A$n_grid_in == 0)))
cat(sprintf("  uniform points inside / 20000 drawn from the box: %s\n", qs(A$box_nin)))
gp <- file.path(EX, "geom_probe.csv")
if (file.exists(gp)) {
  G <- read.csv(gp)
  cat(sprintf("  [12_geom.R, 40 null datasets, 200 random directions each]\n"))
  cat(sprintf("  radial extent of C_gamma in Wald-SE units: min-direction %.1f, median %.1f, max-direction %.1f\n",
              median(G$rmin), median(G$rmed), median(G$rmax[G$runbdd == 0])))
  cat(sprintf("  datasets with at least one unbounded direction: %d / %d\n",
              sum(G$runbdd > 0), nrow(G)))
}

cat("\n--- 2. size (null arm) ------------------------------------------------------------\n")
ln("p_plugin <= .05 (deployed)",       sum(N0$p_plugin <= ALPHA), nrow(N0))
ln("p_exact  <= .05",                  sum(N0$p_exact  <= ALPHA), nrow(N0))
ln("p_exact  <= .04 (= alpha - gamma)", sum(N0$p_exact <= ALPHA - GAMMA), nrow(N0))
pdiff(N0$p_plugin <= ALPHA, N0$p_exact <= ALPHA, "size: plugin.05 - exact.05")
cat(sprintf("  p_exact >= p_plugin on all null datasets: %s\n", all(N0$p_exact >= N0$p_plugin - 1e-12)))
cat(sprintf("  median p_plugin %.3f, median p_exact %.3f, median (p_exact - p_plugin) %.3f\n",
            median(N0$p_plugin), median(N0$p_exact), median(N0$p_exact - N0$p_plugin)))

cat("\n--- 3. power (alternative arm) ----------------------------------------------------\n")
ln("p_plugin <= .05 (deployed)",       sum(N1$p_plugin <= ALPHA), nrow(N1))
ln("p_exact  <= .05",                  sum(N1$p_exact  <= ALPHA), nrow(N1))
ln("p_exact  <= .04 (guaranteed .05)", sum(N1$p_exact <= ALPHA - GAMMA), nrow(N1))
pdiff(N1$p_plugin <= ALPHA, N1$p_exact <= ALPHA, "power: plugin.05 - exact.05")
pdiff(N1$p_plugin <= ALPHA, N1$p_exact <= ALPHA - GAMMA, "power: plugin.05 - exact.04")
cat(sprintf("  median p_plugin %.3f, median p_exact %.3f, median inflation %.3f\n",
            median(N1$p_plugin), median(N1$p_exact), median(N1$p_exact - N1$p_plugin)))

cat("\n--- 4. where the supremum sat -----------------------------------------------------\n")
for (k in c(0L, 2L)) {
  R <- A[A$kind == k, ]
  dhat <- sqrt(((R$sup_0 - R$that_0) / R$se_0)^2 + ((R$sup_1 - R$that_1) / R$se_1)^2 +
               ((R$sup_2 - R$that_2) / R$se_2)^2)
  cat(sprintf("  kind %d: sup at theta-hat %d/%d = %.3f ; when elsewhere, ||sup - theta-hat|| in Wald-SE units: median %.2f, max %.2f\n",
              k, sum(R$at_hat), nrow(R), mean(R$at_hat),
              sf(median, dhat[R$at_hat == 0]), sf(max, dhat[R$at_hat == 0])))
}
cat(sprintf("  Nelder-Mead improved on the grid: %d / %d = %.3f\n",
            sum(A$nm_improved), nrow(A), mean(A$nm_improved)))
cat(sprintf("  p_exact > p_grid: median gain when it improved %.4f\n",
            sf(median, (A$p_exact - A$p_grid)[A$nm_improved == 1])))
cat(sprintf("  search hit the %.0f s budget (p_exact is then a strict LOWER bound): %d / %d = %.3f\n",
            as.numeric(Sys.getenv("EX_BUDGET", "480")), sum(A$trunc), nrow(A), mean(A$trunc)))

cat("\n--- 5. compute --------------------------------------------------------------------\n")
cat(sprintf("  theta evaluations per dataset : %s\n", qs(A$n_eval)))
cat(sprintf("  infeasible NM probes (free)   : %s\n", qs(A$n_infeas)))
cat(sprintf("  refits per dataset            : %s   (B=199 per FULL evaluation)\n", qs(A$nrep_tot)))
cat(sprintf("  seconds per dataset           : %s   mean %.0f\n", qs(A$secs), mean(A$secs)))
cat(sprintf("  the deployed plug-in test is ONE evaluation = 199 refits ~ %.2f s\n",
            median(A$secs) / median(A$n_eval)))
cat(sprintf("  cost ratio exact/plug-in, by refits: median %.0f x, max %.0f x\n",
            median(A$nrep_tot) / 199, max(A$nrep_tot) / 199))

dm <- file.path(EX, "deepmax.csv")
if (file.exists(dm)) {
  D <- read.csv(dm)
  M <- merge(D, A[, c("kind", "rep", "p_plugin", "p_exact", "n_eval", "secs")],
             by = c("kind", "rep"), suffixes = c("_deep", "_prot"))
  cat("\n--- 6. deep-max audit: is the protocol supremum the real supremum? ---------------\n")
  cat(sprintf("  %d datasets re-maximised with 6 extremes + %.0f rejection points + 8 NM restarts\n",
              nrow(M), median(D$n_pts)))
  cat(sprintf("  p_deep > p_exact(protocol): %d / %d = %.3f ; median shortfall %.4f, max %.4f\n",
              sum(M$p_deep > M$p_exact + 1e-12), nrow(M), mean(M$p_deep > M$p_exact + 1e-12),
              median(M$p_deep - M$p_exact), max(M$p_deep - M$p_exact)))
  cat(sprintf("  p_deep = 1.0 on %d / %d datasets (the sup is attained at a theta where EVERY\n",
              sum(M$p_deep >= 1 - 1e-12), nrow(M)))
  cat("  replicate beats the observed statistic, so the exact test can never reject)\n")
  cat(sprintf("  contribution of the six coordinate extremes alone: p_ext > p_hat on %d / %d\n",
              sum(M$p_ext > M$p_hat + 1e-12), nrow(M)))
  cat(sprintf("  rejection-sample stage: acceptance rate median %.5f, points drawn median %.0f\n",
              median(D$accept), median(D$n_drawn)))
  cat(sprintf("  deep search: median %.0f evaluations, median %.0f s per dataset\n",
              median(D$n_eval), median(D$secs)))
  cat(sprintf("  size at .05 under the deep sup (null arm): %d / %d\n",
              sum(M$p_deep[M$kind == 0] <= ALPHA), sum(M$kind == 0)))
  cat(sprintf("  power at .05 under the deep sup (alt arm): %d / %d\n",
              sum(M$p_deep[M$kind == 2] <= ALPHA), sum(M$kind == 2)))
}
cat("\n=================================================================================\n")
