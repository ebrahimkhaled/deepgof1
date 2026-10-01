# run_localization.R -- referee M4: operating characteristics of "where the model fails". Five N(0,1)
# covariates, n = 500, eta0 = -0.5 + 0.5 x1 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5, fitted linear in all five.
#   ushape  eta0 + C (x1^2 - 1),        C in {.4, .8}   truth: x1, cells at both ends of x1 (sixths 1, 6)
#   thresh  eta0 + C max(x2 - 1, 0),     C in {2, 4}     truth: x2, cells at the top of x2 (sixths 5, 6)
#   inter   eta0 + C x1 x3,              C in {.6, 1}    truth: the pair (x1, x3), its corner cells
#   null    eta0                                          (500 datasets; every flag is false)
# 300 datasets per alternative cell. One set of B = 199 parametric-bootstrap refits per dataset serves:
#   DeepGOF-1, axis rule and all-pairs reading: p-value, the pair displayed, and the cells beyond the
#     simultaneous 95% threshold of the displayed map (the 95% quantile of the largest |cell| of the same
#     pair's map over the bootstrap replicates)
#   Lin-Wei-Ying cumulative residuals: for each covariate the sup of |cumulative residual| over its ranks,
#     p-value from the same refits; covariates flagged at Bonferroni .05/5
#   the spline screen: for each covariate, the likelihood-ratio test of ns(x, 3) in place of x (2 df);
#     covariates flagged at Bonferroni .05/5
# ebrahim.gof 2.9.0 internals from the private library. Output rows/<cell>_r<first>.csv, 25 per file.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/localization"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L; N <- 500L

CELLS <- rbind(data.frame(fam = "null", C = 0, reps = 500L),
               data.frame(fam = c("ushape", "ushape", "thresh", "thresh", "inter", "inter"),
                          C = c(.4, .8, 2, 4, .6, 1), reps = 300L))
CELLS$id <- seq_len(nrow(CELLS))
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (k in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[k]), ceiling(seq_len(CELLS$reps[k]) / 25))) {
  f <- file.path(OUT, sprintf("c%02d_%s_C%.1f_r%04d.csv", CELLS$id[k], CELLS$fam[k], CELLS$C[k], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[k, ], reps = ch, f = f)
}
if (Sys.getenv("LO_REV") == "1") jobs <- rev(jobs)

one_chunk <- function(j) {
  suppressMessages({ library(ebrahim.gof, lib.loc = LIB); library(splines) })
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
  K <- 6L
  ## is cell (i, j) of the map over pair (a, b) inside the true region?
  in_region <- function(fam, pair, i, j) {
    bin_of <- function(v) if (pair[1] == v) i else if (pair[2] == v) j else NA
    switch(fam,
      ushape = "x1" %in% pair && bin_of("x1") %in% c(1, 6),
      thresh = "x2" %in% pair && bin_of("x2") %in% c(5, 6),
      inter  = all(c("x1", "x3") %in% pair) && bin_of("x1") %in% c(1, 2, 5, 6) && bin_of("x3") %in% c(1, 2, 5, 6),
      FALSE)
  }
  flags <- function(map, thr, fam, pair) {
    k <- which(abs(map) > thr)
    ins <- vapply(k, function(z) in_region(fam, pair, (z - 1L) %/% K + 1L, (z - 1L) %% K + 1L), TRUE)
    c(n = length(k), inside = sum(ins))
  }
  R <- do.call(rbind, lapply(j$reps, function(r) {
    s <- 9000000L + j$cc$id * 10000L + r
    set.seed(s)
    X0 <- matrix(stats::rnorm(N * 5), N, 5, dimnames = list(NULL, paste0("x", 1:5)))
    eta <- -0.5 + X0 %*% c(0.5, 0.8, 0.6, 0.4, 0.3) + switch(j$cc$fam, null = 0,
             ushape = j$cc$C * (X0[, 1]^2 - 1), thresh = j$cc$C * pmax(X0[, 2] - 1, 0), inter = j$cc$C * X0[, 1] * X0[, 3])
    d <- data.frame(X0, y = stats::rbinom(N, 1, stats::plogis(eta)))
    fit <- suppressWarnings(stats::glm(y ~ x1 + x2 + x3 + x4 + x5, data = d, family = stats::binomial()))
    mm <- stats::model.matrix(fit); X <- ebrahim.gof:::.dg_covariates(fit); vars <- colnames(X)
    tcols <- ebrahim.gof:::.dg_term_cols(fit, mm, vars)
    ph <- as.numeric(stats::fitted(fit)); y <- d$y
    bins <- apply(X, 2, function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / N)))
    cells <- ebrahim.gof:::.dg_cells(X, K, NULL)
    idx_of <- function(p) as.integer((bins[, p[1]] - 1L) * K + bins[, p[2]])
    cm <- function(p, rr, pp) ebrahim.gof:::.dg_cellmap(idx_of(p), rr, pp, K)
    ## observed
    axO <- ebrahim.gof:::.dg_pick_axes(mm, tcols, stats::coef(fit)); mA <- cm(axO, y - ph, ph); SA <- score(mA)
    so <- vapply(cells$idx, function(i) score(ebrahim.gof:::.dg_cellmap(i, y - ph, ph, K)), 0)
    topO <- cells$pairs[which.max(so), ]; mT <- cm(topO, y - ph, ph); ST <- max(so)
    ord <- lapply(vars, function(v) order(X[, v]))
    lwy <- function(rr) vapply(ord, function(o) max(abs(cumsum(rr[o]))) / sqrt(N), 0)
    LO <- lwy(y - ph)
    ## bootstrap
    itc <- TRUE; bA <- bT <- mxA <- mxT <- numeric(B); bL <- matrix(0, B, 5)
    for (b in seq_len(B)) {
      ys <- stats::rbinom(N, 1L, ph)
      fb <- suppressWarnings(stats::glm.fit(mm, ys, family = stats::binomial(), intercept = itc))
      pb <- as.numeric(fb$fitted.values); rb <- ys - pb
      axb <- ebrahim.gof:::.dg_pick_axes(mm, tcols, fb$coefficients)
      bA[b] <- score(cm(axb, rb, pb))
      bT[b] <- max(vapply(cells$idx, function(i) score(ebrahim.gof:::.dg_cellmap(i, rb, pb, K)), 0))
      mxA[b] <- max(abs(cm(axO, rb, pb))); mxT[b] <- max(abs(cm(topO, rb, pb)))
      bL[b, ] <- lwy(rb)
    }
    pA <- (1 + sum(bA >= SA)) / (B + 1); pT <- (1 + sum(bT >= ST)) / (B + 1)
    fA <- flags(mA, stats::quantile(mxA, .95, names = FALSE), j$cc$fam, axO)
    fT <- flags(mT, stats::quantile(mxT, .95, names = FALSE), j$cc$fam, topO)
    pL <- vapply(1:5, function(k) (1 + sum(bL[, k] >= LO[k])) / (B + 1), 0)
    ## spline screen, one covariate at a time
    pS <- vapply(vars, function(v) {
      f2 <- stats::as.formula(paste("y ~", paste(ifelse(vars == v, sprintf("ns(%s, 3)", v), vars), collapse = " + ")))
      g <- suppressWarnings(stats::glm(f2, data = d, family = stats::binomial()))
      stats::pchisq(fit$deviance - g$deviance, 2, lower.tail = FALSE)
    }, 0)
    data.frame(fam = j$cc$fam, C = j$cc$C, rep = r, events = sum(y),
               p_axes = pA, axes = paste(axO, collapse = "|"), flagA = fA[["n"]], flagA_in = fA[["inside"]],
               p_all = pT, top = paste(topO, collapse = "|"), flagT = fT[["n"]], flagT_in = fT[["inside"]],
               lwy_flags = paste(vars[pL <= .05 / 5], collapse = "|"), lwy_minp = min(pL),
               spl_flags = paste(vars[pS <= .05 / 5], collapse = "|"), spl_minp = min(pS))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("localization%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("LO_WORKERS", "6")))
  parallel::clusterExport(cl, c("LIB", "B", "N", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
