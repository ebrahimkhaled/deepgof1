# =============================================================================
# partA.R -- Part A of PREDECLARATION.md: every paper test on the REAL outcome of each included
# main-effects model, rows permuted with the pre-declared seed (set.seed(20260925 + k)).
#   Rscript partA.R main        -> partA_pvalues_main.csv, partA_axes.csv, partA_perm_sensitivity.csv, perm_*.csv
#   Rscript partA.R bagoft k    -> partA_bagoft_<data>.csv (BAGofT, fast engine, optional arm)
# Test seeds: 20260926 + 10 k + {1 DeepGOF-1, 2 Stute-Zhu, 3 projection, 4 BAGofT}.
# =============================================================================
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche/common.R")
args <- commandArgs(TRUE); mode <- args[1]
SEEDA <- 20260926L

if (mode == "main") {
  rows <- list(); axes <- list(); sens <- list()
  for (k in 1:3) {
    d <- prep_data(k, write_perm = TRUE); frm <- frm_of(d)
    fit <- glm(frm, data = d, family = binomial())
    stopifnot(fit$converged)
    cat(sprintf("\n==== %s: n = %d, events = %d, p = %d, logLik %.4f\n", DNAMES[k], nrow(d), sum(d$y), ncol(d) - 1, as.numeric(logLik(fit))))
    ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
    res <- list(); tm <- list()
    tt <- function(nm, expr) { t0 <- proc.time()[[3]]; v <- expr; res[[nm]] <<- v; tm[[nm]] <<- proc.time()[[3]] - t0 }
    set.seed(SEEDA + 10L * k + 1L)
    tt("DeepGOF1", { dg <- ebrahim.gof::deepgof1(fit, B = 199L, K = 6L); dg$p.value })
    tt("HL", pv(ebrahim.gof:::gof_hl(ctx)))
    tt("Stukel_marginal", pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))))
    tt("Stukel_joint", pv(ebrahim.gof:::gof_stukel(ctx, list(form = "joint"))))
    set.seed(SEEDA + 10L * k + 2L)
    tt("StuteZhu", pv(ebrahim.gof:::gof_stutezhu(ctx)))
    tt("GiViTI", pv(ebrahim.gof:::gof_giviti(ctx)))
    tt("mHL", pv(ebrahim.gof:::gof_ftest(ctx)))
    tt("OsiusRojek", pv(ebrahim.gof:::gof_osius(ctx)))
    X <- model.matrix(fit)[, -1, drop = FALSE]
    set.seed(SEEDA + 10L * k + 3L)
    tt("proj_raw", { AL <- list(raw = proj_A_liu_fast(X, std = FALSE), std = proj_A_liu_fast(X, std = TRUE))
                     pr <- proj_p_multi(d, frm, fit, B = 999L, AL = AL); unname(pr$p["raw"]) })
    res$proj_std <- unname(pr$p["std"]); tm$proj_std <- NA
    for (a in names(res)) rows[[length(rows) + 1]] <- data.frame(data = DNAMES[k], test = a, p = res[[a]], secs = tm[[a]])
    b <- coef(fit)[-1]; sdx <- apply(X, 2, sd)
    axes[[k]] <- data.frame(data = DNAMES[k], axis1 = dg$axes[1], axis2 = dg$axes[2], statistic = dg$statistic,
                            ranking_by_abs_b_sd = paste(sprintf("%s %.3f", names(sort(abs(b) * sdx, decreasing = TRUE)),
                                                                 sort(abs(b) * sdx, decreasing = TRUE)), collapse = "; "))
    cat("  DeepGOF-1 axes:", dg$axes, " S =", format(dg$statistic, digits = 5), "\n")
    for (a in names(res)) cat(sprintf("  %-16s p = %.4f\n", a, res[[a]]))
    ## declared secondary: 20 further random row orders of the ORIGINAL file, same bootstrap seed
    d0 <- raw_data(k); d0 <- d0[stats::complete.cases(d0), ]; rownames(d0) <- NULL
    for (j in 1:20) {
      set.seed(20261000L + 100L * k + j); o <- sample(nrow(d0)); dj <- d0[o, ]; rownames(dj) <- NULL
      fj <- glm(frm, data = dj, family = binomial())
      set.seed(SEEDA + 10L * k + 1L); gj <- ebrahim.gof::deepgof1(fj, B = 199L, K = 6L)
      sens[[length(sens) + 1]] <- data.frame(data = DNAMES[k], order = j, p = gj$p.value, statistic = gj$statistic,
                                             axis1 = gj$axes[1], axis2 = gj$axes[2])
    }
    S <- do.call(rbind, sens); S <- S[S$data == DNAMES[k], ]
    cat(sprintf("  20 other row orders: p median %.3f, range %.3f-%.3f, #p<=.05 = %d\n", median(S$p), min(S$p), max(S$p), sum(S$p <= .05)))
    flush.console()
  }
  write.csv(do.call(rbind, rows), file.path(RD, "partA_pvalues_main.csv"), row.names = FALSE)
  write.csv(do.call(rbind, axes), file.path(RD, "partA_axes.csv"), row.names = FALSE)
  write.csv(do.call(rbind, sens), file.path(RD, "partA_perm_sensitivity.csv"), row.names = FALSE)
  cat("\nDONE main\n")
} else if (mode == "bagoft") {
  k <- as.integer(args[2])
  source(file.path(THEORY, "addon_arms", "bagoft_fast_pre.R"))
  d <- prep_data(k); frm <- frm_of(d)
  set.seed(SEEDA + 10L * k + 4L)
  t0 <- proc.time()[[3]]; bg <- bagoft_fast_pre(d, frm); tb <- proc.time()[[3]] - t0
  out <- data.frame(data = DNAMES[k], test = c("BAGofT", "BAGofT_pvalue3"), p = c(bg$p.value, bg$p.value3), secs = c(tb, NA))
  write.csv(out, file.path(RD, sprintf("partA_bagoft_%s.csv", DNAMES[k])), row.names = FALSE)
  print(out)
}
