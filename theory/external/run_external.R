# run_external.R -- external validation: a published risk model with FIXED coefficients is checked on new
# patients. H0: each y_i ~ Bernoulli(p_i) with p_i given; nothing is refitted. Five N(0,1) covariates;
# the published model is eta0 = -0.5 + 0.5 x1 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5 (the localization design).
# The truth in the validation population is eta0 + d:
#   null    d = 0                          n = 250, 500 (500 reps), 1,000 (300 reps)
#   large   d = C                          C in {.2, .4}   calibration-in-the-large
#   slope   d = (C - 1) eta0               C in {.8, .6}   an overfitted model (calibration slope C)
#   ushape  d = C (x1^2 - 1)               C in {.4, .8}
#   thresh  d = C max(x2 - 1, 0)           C in {2, 4}
#   inter   d = C x1 x3                    C in {.6, 1}
# Alternatives at n = 500, 300 reps each.
# Tests. Monte Carlo (B = 199 draws y* ~ Bernoulli(p), shared by all MC tests; exact level at any n):
#   DeepGOF-1 axis rule (axes fixed by the published coefficients), all-pairs, combined (exact min-p);
#   the projection statistic of Liu et al. (2024); le Cessie-van Houwelingen's kernel statistic Q (their
#   kernel and bandwidth); Pearson's X^2.
# Standard external-validation references: Hosmer-Lemeshow deciles against chi-square(10); Stukel's two
# directions added to the offset eta0, likelihood-ratio chi-square(2); Cox's recalibration (intercept and
# slope), LR chi-square(2); Spiegelhalter's z; the GiViTI belt with devel = "external".
# Output rows/<cell>_r<first>.csv, 25 replicates per file (resumable).
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/external"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L

CELLS <- rbind(data.frame(fam = "null", C = 0, n = c(250L, 500L, 1000L), reps = c(500L, 500L, 300L)),
               data.frame(fam = rep(c("large", "slope", "ushape", "thresh", "inter"), each = 2),
                          C = c(.2, .4, .8, .6, .4, .8, 2, 4, .6, 1), n = 500L, reps = 300L))
CELLS$id <- seq_len(nrow(CELLS))
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (k in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[k]), ceiling(seq_len(CELLS$reps[k]) / 25))) {
  f <- file.path(OUT, sprintf("c%02d_%s_C%.1f_n%04d_r%04d.csv", CELLS$id[k], CELLS$fam[k], CELLS$C[k], CELLS$n[k], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[k, ], reps = ch, f = f)
}
if (Sys.getenv("EXT_REV") == "1") jobs <- rev(jobs)

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
  K <- 6L; b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  R <- do.call(rbind, lapply(j$reps, function(r) {
    n <- j$cc$n
    set.seed(7000000L + j$cc$id * 10000L + r)
    X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta0 <- as.numeric(cbind(1, X) %*% b0); p <- stats::plogis(eta0)
    d <- switch(j$cc$fam, null = 0, large = j$cc$C, slope = (j$cc$C - 1) * eta0,
                ushape = j$cc$C * (X[, 1]^2 - 1), thresh = j$cc$C * pmax(X[, 2] - 1, 0),
                inter = j$cc$C * X[, 1] * X[, 3])
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
    ## --- the Monte Carlo statistics, as functions of an outcome vector
    cells <- ebrahim.gof:::.dg_cells(X, K, NULL)
    sdc <- abs(b0[-1]) * apply(X, 2, stats::sd)                 # the axis rule, published coefficients
    ax <- sort(order(sdc, decreasing = TRUE)[1:2])
    kax <- which(cells$pairs[, 1] == colnames(X)[ax[1]] & cells$pairs[, 2] == colnames(X)[ax[2]])
    A <- ebrahim.gof:::.proj_weight(X)
    Dm <- as.matrix(stats::dist(scale(X)) * sqrt(0.5))           # le Cessie: sqrt of summed 0.5 d^2
    Rk <- pmax(1 - Dm / mean(Dm), 0)
    v <- p * (1 - p)
    stats_of <- function(yy) {
      e <- yy - p
      sc <- vapply(cells$idx, function(i) score(ebrahim.gof:::.dg_cellmap(i, e, p, K)), 0)
      c(axes = sc[kax], allpairs = max(sc), proj = ebrahim.gof:::.proj_stat(e, A),
        lecessie = sum(as.numeric(e %*% Rk) * e), pearson = sum(e^2 / v))
    }
    obs <- stats_of(y)
    bt <- t(vapply(seq_len(B), function(b) stats_of(stats::rbinom(n, 1L, p)), obs))
    pmc <- vapply(names(obs), function(k) (1 + sum(bt[, k] >= obs[[k]])) / (B + 1), 0)
    all_a <- c(obs[["axes"]], bt[, "axes"]); all_p <- c(obs[["allpairs"]], bt[, "allpairs"])
    within <- function(v) vapply(seq_along(v), function(i) sum(v >= v[i]) / length(v), 0)
    mins <- pmin(within(all_a), within(all_p)); pcomb <- sum(mins <= mins[1]) / (B + 1)
    ## --- standard external-validation tests
    g <- pmin(ceiling(rank(p, ties.method = "first") / (n / 10)), 10)
    O <- tapply(y, g, sum); E <- tapply(p, g, sum); ng <- tabulate(g, 10)
    phl <- stats::pchisq(sum((O - E)^2 / (E * (1 - E / ng))), 10, lower.tail = FALSE)
    dev0 <- -2 * sum(y * log(p) + (1 - y) * log(1 - p))
    za <- 0.5 * eta0^2 * (p >= .5); zb <- -0.5 * eta0^2 * (p < .5)
    gs <- suppressWarnings(stats::glm(y ~ 0 + za + zb, offset = eta0, family = stats::binomial()))
    pstuk <- stats::pchisq(dev0 - gs$deviance, 2, lower.tail = FALSE)
    gc <- suppressWarnings(stats::glm(y ~ eta0, family = stats::binomial()))
    pcox <- stats::pchisq(dev0 - gc$deviance, 2, lower.tail = FALSE)
    zsp <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))
    pspz <- 2 * stats::pnorm(-abs(zsp))
    pgiv <- tryCatch(suppressWarnings(givitiR::givitiCalibrationBelt(o = y, e = p, devel = "external")$p.value),
                     error = function(e) NA_real_)
    data.frame(fam = j$cc$fam, C = j$cc$C, n = n, rep = r, events = sum(y),
               DG_axes = pmc[["axes"]], DG_allpairs = pmc[["allpairs"]], DG_combined = pcomb,
               proj = pmc[["proj"]], lecessie = pmc[["lecessie"]], pearson = pmc[["pearson"]],
               HL = phl, stukel = pstuk, cox = pcox, spiegelhalter = pspz, giviti = pgiv)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("external%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("EXT_WORKERS", "8")))
  parallel::clusterExport(cl, c("LIB", "B", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
