# =============================================================================
# verify_fast.R -- is bagoft_fast() the SAME test as BAGofT 1.0.0, and how much faster?
#
# Three checks, all required to pass before bagoft_fast() is used for anything:
#
#   A. AGAINST THE PACKAGE. 24 datasets (4 Liu settings x n = 50/100/500 x C = 0 and .5),
#      the same seed into both, at nsplits = 10 and nsim = 5 so the package run is short.
#      Every output -- the three p-values, the three statistics and all 15 simulated
#      statistics -- must be identical(), not merely close. Setting 1 is compared with
#      parRF_dropfix.R, because BAGofT 1.0.0 crashes on a single covariate.
#
#   B. AGAINST THE TRIAL'S DEPOSIT, AT THE PUBLISHED DEFAULTS. One dataset per setting x n
#      block (12), regenerated from the seed in bagoft_pvalues.csv, seeded exactly as
#      run_trial.R seeded it, run at nsplits = 100, nsim = 100. The deposited p, p_med and
#      pmin statistic -- produced by the package itself during the trial -- must be reproduced
#      exactly. This is the check that matters: full defaults, real trial data, and no
#      possibility that the reference was computed by the code under test.
#
#   C. SPEED. One observed-data pass (100 splits, 1/101 of a test) by each, same dataset,
#      one core, nothing else in parallel with it.
#
# Run:  Rscript verify_fast.R          (writes VERIFY_FAST.log beside itself)
# =============================================================================
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/explore_bagoft"
source(file.path(HERE, "common.R"))
source(file.path(HERE, "parRF_dropfix.R"))
source(file.path(HERE, "bagoft_fast.R"))
suppressPackageStartupMessages({ library(BAGofT); library(randomForest) })

LOG <- file.path(HERE, "VERIFY_FAST.log")
sink(LOG, split = TRUE)
cat("verify_fast.R --", format(Sys.time()), "\n",
    "R", R.version.string, "| BAGofT", as.character(packageVersion("BAGofT")),
    "| randomForest", as.character(packageVersion("randomForest")), "\n\n")

prep <- function(d) as.data.frame(lapply(d, as.vector))      # what bagoft_p() does first

orig <- function(dat, frm, set1, nsplits, nsim) {
  environment(frm) <- globalenv()
  out <- NULL
  invisible(capture.output(out <- suppressMessages(suppressWarnings(
    BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = frm, link = "logit"),
                   parFun = if (set1) parRF_dropfix(ntree = 60) else BAGofT::parRF(ntree = 60),
                   data = dat, nsplits = nsplits,
                   ne = floor(5 * sqrt(nrow(dat))), nsim = nsim)))))
  out
}
fast <- function(dat, frm, nsplits, nsim) {
  environment(frm) <- globalenv()
  suppressWarnings(bagoft_fast(dat, frm, nsplits = nsplits, nsim = nsim, ntree = 60L))
}
same <- function(a, b) {
  k <- c("p.value", "p.value2", "p.value3", "pmean", "pmedian", "pmin")
  all(vapply(k, function(x) identical(unname(a[[x]]), unname(b[[x]])), TRUE)) &&
    identical(unname(a$simRes$pmeanSim),   unname(b$simRes$pmeanSim)) &&
    identical(unname(a$simRes$pmediansim), unname(b$simRes$pmediansim)) &&
    identical(unname(a$simRes$pminsim),    unname(b$simRes$pminsim))
}

# ---------------------------------------------------------------- A. against the package
cat("== A. against BAGofT 1.0.0 itself (nsplits = 10, nsim = 5) ==\n")
okA <- 0L; nA <- 0L
for (SET in 1:4) for (N in c(50L, 100L, 500L)) for (C in c(0, 0.5)) {
  dd <- liu_dataset(101000L, SET, N, C, 7L)
  if (is.null(dd$data)) next
  dat <- prep(dd$data); frm <- FRM(dat)
  set.seed(dd$seed + 500000000L); a <- orig(dat, frm, SET == 1L, 10L, 5L)
  set.seed(dd$seed + 500000000L); b <- fast(dat, frm, 10L, 5L)
  ok <- !is.null(a) && same(a, b)
  nA <- nA + 1L; okA <- okA + ok
  cat(sprintf("  S%d n=%3d C=%.2f  p %.2f/%.2f  pmean %.12f/%.12f  -> %s\n", SET, N, C,
              a$p.value, b$p.value, a$pmean, b$pmean, if (ok) "IDENTICAL" else "*** DIFFERS ***"))
}
cat(sprintf("A: %d of %d identical\n\n", okA, nA))

# ------------------------------------------------- B. against the trial deposit, full defaults
cat("== B. against the trial's deposited results (nsplits = 100, nsim = 100) ==\n")
dep <- utils::read.csv(file.path(HERE, "bagoft_pvalues.csv"))
pick <- do.call(rbind, lapply(split(dep[dep$C == 0.5, ], list(dep$setting[dep$C == 0.5],
                                                              dep$n[dep$C == 0.5])),
                              function(g) g[which.min(g$rep), ]))
cl <- parallel::makeCluster(as.integer(Sys.getenv("VF_WORKERS", "4")))   # the machine may be shared
parallel::clusterExport(cl, c("HERE", "prep", "fast"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(HERE, "common.R")); source(file.path(HERE, "bagoft_fast.R")); NULL }))
resB <- parallel::parLapplyLB(cl, split(pick, seq_len(nrow(pick))), function(r) {
  dd  <- liu_dataset(101000L, r$setting, r$n, r$C, r$rep)
  dat <- prep(dd$data); frm <- FRM(dat)
  set.seed(dd$seed + 500000000L)
  t0 <- proc.time()[["elapsed"]]
  b  <- fast(dat, frm, 100L, 100L)
  data.frame(setting = r$setting, n = r$n, rep = r$rep, seed_ok = dd$seed == r$seed,
             p_dep = r$p, p_fast = b$p.value, pmed_dep = r$p_med, pmed_fast = b$p.value2,
             stat_dep = r$bag_pmin_stat, stat_fast = b$pmin,
             secs_dep = r$time_s, secs_fast = proc.time()[["elapsed"]] - t0)
})
parallel::stopCluster(cl)
B <- do.call(rbind, resB)
# The deposit is a CSV, and R writes a double to text at 15 significant digits. The pmin
# statistic is continuous, so the value read back can differ from the full-precision one in
# its last bit (a first run compared them raw and failed 7 of 12 on a relative difference of
# 1.4e-16). The right comparison is at the precision the deposit actually holds: the fast
# value, written to text the way the trial wrote it, must read back as the deposited number.
B$match <- B$seed_ok & B$p_dep == B$p_fast & B$pmed_dep == B$pmed_fast &
           as.numeric(as.character(B$stat_fast)) == B$stat_dep
print(B[order(B$setting, B$n), ], row.names = FALSE, digits = 6)
cat(sprintf("B: %d of %d deposited full-default results reproduced exactly\n\n",
            sum(B$match), nrow(B)))

# ---------------------------------------------------------------------------- C. speed
cat("== C. speed: one observed-data pass (100 splits = 1/101 of a test), one core ==\n")
for (N in c(50L, 100L, 500L)) {
  dd  <- liu_dataset(101000L, 4L, N, 0.5, 3L); dat <- prep(dd$data); frm <- FRM(dat)
  set.seed(1L); t0 <- proc.time()[["elapsed"]]; a <- orig(dat, frm, FALSE, 100L, 0L)
  to <- proc.time()[["elapsed"]] - t0
  set.seed(1L); t0 <- proc.time()[["elapsed"]]; b <- fast(dat, frm, 100L, 0L)
  tf <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  n=%3d  package %.2f s | fast %.2f s | x%.1f faster | same statistic: %s\n",
              N, to, tf, to / tf, identical(unname(a$pmean), unname(b$pmean))))
}
sink()
