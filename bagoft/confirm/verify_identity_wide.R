# =============================================================================
# verify_identity_wide.R -- widen the certificate that bagoft_fast() IS BAGofT 1.0.0.
#
# WHY. The confirmatory run used bagoft_fast(), which was verified on 36 datasets: 24 at reduced
# nsplits/nsim and 12 at full defaults. An adversarial review put the arithmetic plainly -- that is
# 0.5% of the 7,200 datasets and about 0.002% of the splits, no dataset was checked at full
# defaults UNDER THE NULL, and three of the four departure strengths were never checked at all.
# Thirty-six clean comparisons are consistent with a per-dataset divergence rate as high as 8%.
# Nothing suggests a divergence exists; the point is that the deposit could not rule one out.
#
# The exposure is data-dependent branches that few datasets reach: a constant forest prediction
# (length(urt) == 1), the is.nan carry-forward in the K search, the merge loop over undersized
# groups with its tied sort(), and a rank-deficient glm.fit taking the pivot path. Those fire most
# at n = 50 and at strong departures -- the region the original check never entered.
#
# WHAT THIS DOES. Draws a stratified random sample of the confirmatory run's own datasets, covering
# every setting, sample size and departure strength, and runs the RELEASED package at the published
# defaults on each, comparing its p-value with the deposited one exactly. Setting 1 uses
# parRF_dropfix.R, because BAGofT 1.0.0 crashes on a single covariate -- and that is itself worth
# stating in the paper rather than leaving inside an "identical to BAGofT 1.0.0" stamp.
#
# Run:  Rscript verify_identity_wide.R        Env: VI_PER_CELL (default 4), VI_WORKERS (20)
# =============================================================================
PAP  <- "/path/to/project/paper_deepgof"
EXPL <- file.path(PAP, "explore_bagoft")
HERE <- file.path(PAP, "confirm_bagoft")
source(file.path(EXPL, "common.R"))
source(file.path(EXPL, "parRF_dropfix.R"))
suppressPackageStartupMessages(library(BAGofT))

PER  <- as.integer(Sys.getenv("VI_PER_CELL", "4"))     # datasets per (setting, n, C) cell
NW   <- as.integer(Sys.getenv("VI_WORKERS", "20"))
OUT  <- file.path(HERE, "IDENTITY_WIDE.csv")
## NOT the file the launcher redirects into. On Windows that redirect holds an exclusive handle,
## so a second writer gets "Permission denied" -- which is exactly how the first attempt died, at
## startup, before computing anything. Progress goes to stdout, which the launcher already captures.
LOG  <- file.path(HERE, "IDENTITY_WIDE_progress.log")

dep <- utils::read.csv(file.path(HERE, "rows", "..", "..", "confirm_bagoft", "CONFIRM_ROW.csv"))
rows <- do.call(rbind, lapply(list.files(file.path(HERE, "rows"), full.names = TRUE),
                              utils::read.csv))
rows <- rows[rows$status == "ok", ]

set.seed(717171L)                                       # the sample itself is reproducible
pick <- do.call(rbind, lapply(split(rows, list(rows$setting, rows$n, rows$C), drop = TRUE),
                              function(g) g[sample(nrow(g), min(PER, nrow(g))), ]))
## VI_MAXROWS exists to smoke-test the plumbing on a couple of datasets before committing three
## hours to the real thing; leave it unset for the run that certifies anything.
if (nzchar(Sys.getenv("VI_MAXROWS", ""))) {
  pick <- pick[seq_len(min(nrow(pick), as.integer(Sys.getenv("VI_MAXROWS")))), ]
  OUT <- file.path(HERE, "IDENTITY_WIDE_smoke.csv")
}
cat(sprintf("%s  identity re-check: %d datasets over %d cells, %d workers\n",
            format(Sys.time()), nrow(pick), nrow(unique(pick[c("setting", "n", "C")])), NW),
    file = LOG, append = TRUE)

one <- function(r) {
  dd <- liu_dataset(101000L, r$setting, r$n, r$C, r$rep)
  dat <- as.data.frame(lapply(dd$data, as.vector))
  frm <- FRM(dat); environment(frm) <- globalenv()
  set.seed(dd$seed + 500000000L)                        # the stream the run itself used
  t0 <- proc.time()[["elapsed"]]
  b <- tryCatch({
    out <- NULL
    invisible(capture.output(out <- suppressMessages(suppressWarnings(
      BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = frm, link = "logit"),
                     parFun = if (r$setting == 1L) parRF_dropfix(ntree = 60L)
                              else BAGofT::parRF(ntree = 60L),
                     data = dat, nsplits = 100L, nsim = 100L,
                     ne = floor(5 * sqrt(nrow(dat))))))))
    out
  }, error = function(e) NULL)
  data.frame(setting = r$setting, n = r$n, C = r$C, rep = r$rep,
             p_deposit = r$p, p_package = if (is.null(b)) NA_real_ else b$p.value,
             pmed_deposit = r$p_med, pmed_package = if (is.null(b)) NA_real_ else b$p.value2,
             stat_deposit = r$bag_pmin_stat, stat_package = if (is.null(b)) NA_real_ else b$pmin,
             secs = round(proc.time()[["elapsed"]] - t0, 1),
             reference = if (r$setting == 1L) "parRF_dropfix" else "BAGofT 1.0.0")
}

cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("EXPL", "one"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(EXPL, "common.R")); source(file.path(EXPL, "parRF_dropfix.R"))
  suppressPackageStartupMessages(library(BAGofT)); NULL }))
res <- do.call(rbind, parallel::parLapplyLB(cl, split(pick, seq_len(nrow(pick))), one))
parallel::stopCluster(cl)

## the deposit stores 15 significant digits, so the statistic is compared at that precision
res$match <- res$p_deposit == res$p_package & res$pmed_deposit == res$pmed_package &
  as.numeric(as.character(res$stat_package)) == res$stat_deposit
utils::write.csv(res, OUT, row.names = FALSE)
bad <- sum(!res$match | is.na(res$match))
msg <- sprintf("%s  %d of %d reproduced exactly (%d mismatched or errored); median %.0f s\n",
               format(Sys.time()), sum(res$match, na.rm = TRUE), nrow(res), bad, median(res$secs))
cat(msg, file = LOG, append = TRUE)
cat(msg)
if (bad) print(res[!res$match | is.na(res$match), ], row.names = FALSE)
