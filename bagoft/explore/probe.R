# =============================================================================
# probe.R -- STEP 1: feasibility + timing probe.   *** EXPLORATORY, NOT
# CONFIRMATORY. *** Its numbers may never be published as a head-to-head.
#
# (a) does BAGofT load, and which version
# (b) is ne = floor(5*sqrt(n)) well defined / degenerate at n = 50 and n = 100
# (c) seconds per dataset at n = 50 and n = 100 on Liu Setting 2
# =============================================================================
source("/path/to/project/paper_deepgof/explore_bagoft/common.R")
SEED0 <- 101000L
sink(file.path(EXPL, "PROBE.log"), split = TRUE)

cat("=== EXPLORATORY BAGofT small-n probe ===\n")
cat("BAGofT version:", as.character(packageVersion("BAGofT")), "\n")
cat("randomForest  :", as.character(packageVersion("randomForest")), "\n")
cat(R.version.string, "\n\n")

cat("--- (b1) split arithmetic at the published defaults ---\n")
geo <- do.call(rbind, lapply(c(50, 100, 500), bag_geometry))
print(geo, row.names = FALSE)
cat("\n  ne      = size of the VALIDATION set (where the chi-square is formed)\n")
cat("  n_train = n - ne = rows left to (i) REFIT the logistic model and\n")
cat("            (ii) grow the parRF partition forest.\n")
cat("  df of the per-split chi-square = ngp (the number of groups), BAGofT_sin():\n")
cat("            P = 1 - pchisq(chisq, ngp)\n\n")

## --- (b2) what actually happens inside one split at each n -------------------
## nsim = 0 => the statistic only, no bootstrap: ~1/101 of the cost. This is
## where we can see ngp, the refit, and any warnings.
diag_one <- function(n, SET, C, i) {
  dd <- liu_dataset(SEED0, SET, n, C, i)
  if (is.null(dd$data)) return(NULL)
  d <- as.data.frame(lapply(dd$data, as.vector))
  frm <- FRM(d); environment(frm) <- globalenv()
  ne <- floor(5*sqrt(n)); nt <- n - ne
  ## how bad is the refit on nt rows? replicate BAGofT_sin's split 100x
  set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)
  sep <- 0L; nev <- integer(0)
  for (k in 1:100) {
    tri <- sample(seq_len(n), nt); dt <- d[tri, ]
    nev <- c(nev, sum(dt$y))
    g <- suppressWarnings(glm(frm, data = dt, family = binomial()))
    if (!g$converged || any(abs(coef(g)[-1]) > 8) || length(unique(dt$y)) < 2L) sep <- sep + 1L
  }
  ## the real thing, statistic only
  set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)
  warns <- character(0)
  t0 <- Sys.time()
  r <- withCallingHandlers(
    tryCatch(invisible(capture.output(res <- suppressMessages(
      BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = frm, link = "logit"),
                     parFun = BAGofT::parRF(ntree = 60), data = d,
                     nsplits = 100, ne = ne, nsim = 0)))), error = function(e) {
        warns <<- c(warns, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!exists("res")) return(NULL)
  ngp <- vapply(res$singleSplit.results, function(z) z$ngp, numeric(1))
  chs <- vapply(res$singleSplit.results, function(z) z$chisq, numeric(1))
  pvs <- vapply(res$singleSplit.results, function(z) z$p.value, numeric(1))
  list(n = n, i = i, ne = ne, nt = nt,
       events_train_med = median(nev), events_train_min = min(nev),
       pct_refit_bad = sep,
       ngp_tab = table(ngp), chisq_med = median(chs), pv_mean = mean(pvs),
       n_nonfinite_chisq = sum(!is.finite(chs)),
       nwarn = length(warns), warn1 = if (length(warns)) warns[1] else "",
       secs_nsim0 = el)
}

for (n in c(50, 100)) {
  cat(sprintf("--- (b2) inside the splits, Liu S2, C=0, n=%d (5 datasets, nsim=0) ---\n", n))
  for (i in 1:5) {
    z <- diag_one(n, 2L, 0, i)
    if (is.null(z)) { cat("  dataset", i, ": DROPPED (degenerate y), as liu_grid.R drops it\n"); next }
    cat(sprintf("  ds %d: ne=%d n_train=%d | train events med=%d min=%d | refits flagged bad %d/100 | ngp: %s | median chisq=%.2f | mean split p=%.3f | nonfinite chisq=%d | warns=%d [%s] | %.1fs (nsim=0)\n",
                i, z$ne, z$nt, z$events_train_med, z$events_train_min, z$pct_refit_bad,
                paste(sprintf("%s:%d", names(z$ngp_tab), as.integer(z$ngp_tab)), collapse=" "),
                z$chisq_med, z$pv_mean, z$n_nonfinite_chisq, z$nwarn,
                substr(z$warn1, 1, 60), z$secs_nsim0))
    flush.console()
  }
  cat("\n")
}

## --- (c) FULL published-default timing, 5 datasets per n ---------------------
cat("--- (c) FULL published-default timing (nsplits=100, nsim=100, ntree=60) ---\n")
tim <- list()
for (n in c(50, 100)) {
  for (i in 1:5) {
    dd <- liu_dataset(SEED0, 2L, n, 0, i)
    if (is.null(dd$data)) next
    frm <- FRM(dd$data)
    set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)
    t0 <- Sys.time()
    r <- bagoft_p(dd$data, frm, nsplits = 100L, nsim = 100L, ntree = 60L)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  n=%3d ds %d: p=%s p_med=%s | %.1f s\n", n, i,
                format(r$p), format(r$p_med), el)); flush.console()
    tim[[length(tim)+1]] <- data.frame(n = n, i = i, seed = dd$seed,
                                       p = r$p, p_med = r$p_med, secs = el)
  }
}
TT <- do.call(rbind, tim)
write.csv(TT, file.path(EXPL, "probe_timing.csv"), row.names = FALSE)
cat("\nmean seconds/dataset (SERIAL, 1 worker):\n")
print(aggregate(secs ~ n, TT, mean), row.names = FALSE)
cat("\nWritten:", file.path(EXPL, "probe_timing.csv"), "\n")
sink()
