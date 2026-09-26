# =============================================================================
# niche_bagoft_fast.R -- the BAGofT arm on the niche, n = 200 scenarios only, 100 reps each
# (7 scenarios, 700 datasets), with bagoft_fast_pre() in place of the package call.
#
# Same datasets as theory/niche/rows/: niche_addon.R's scenario table, gen() and frm() are
# parsed out of niche_addon.R (itself a verbatim copy of run_niche.R), seed
# 123000 + id * 100000 + rep. BAGofT's own RNG is seeded with dataset seed + 500000000, as in
# niche_addon.R's bagoft mode and confirm_bagoft/run_confirm.R. Settings: testGlmBi(logit),
# parRF(ntree = 60) with its dcPre() pre-selection (ten covariates > 5), nsplits = 100,
# nsim = 100, ne = floor(5 sqrt(n)).
#
# GATE. It refuses to run unless verify_fast_pre.csv shows at least 20 stage-B datasets
# (full settings) with p.value, p.value2, p.value3 and the statistics identical to the package.
#
# Output: one file per dataset, niche_rows_bagoft/<scen>_r<rep>.csv, the same columns as
# niche_addon.R's bagoft mode (keys scen, fam, C, n, rep, events, dseed; p = p.value, the
# headline; p_med = p.value2; p_min = p.value3; bag_pmin_stat = the statistic pmin, NOT a
# p-value) plus impl = "fast_pre". Rep-major order, so a cut-short run stays balanced.
# Resumable: existing files are skipped; each file is written .part and renamed.
#
# Run:  Rscript niche_bagoft_fast.R        Env: AD_W workers (8), AD_MAXREP (100)
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
NW   <- as.integer(Sys.getenv("AD_W", "8"))
MAXR <- as.integer(Sys.getenv("AD_MAXREP", "100"))
OUTD <- Sys.getenv("AD_OUT", file.path(HERE, "niche_rows_bagoft"))
SMOKE <- nzchar(Sys.getenv("AD_SMOKE"))   # smoke test only: nsplits = 2, nsim = 1, no gate; needs AD_OUT
if (SMOKE && !nzchar(Sys.getenv("AD_OUT"))) stop("AD_SMOKE needs AD_OUT")
NSPL <- if (SMOKE) 2L else 100L
NSIM <- if (SMOKE) 1L else 100L

## ---- gate: the verification must have passed ----
vf <- file.path(HERE, "verify_fast_pre.csv")
if (!SMOKE && !file.exists(vf)) stop("verify_fast_pre.csv missing: run the verification first")
if (!SMOKE) { V <- utils::read.csv(vf); VB <- V[V$stage == "B", ]
okB <- VB$p_identical & VB$p3_identical & VB$stats_identical &
       VB$p2_pkg == VB$p2_fast
if (nrow(VB) < 20 || !all(okB)) stop(sprintf("verification not passed: %d stage-B rows, %d identical", nrow(VB), sum(okB))) }
dir.create(OUTD, showWarnings = FALSE)

## ---- the niche generator, parsed out of niche_addon.R ----
ex <- parse(file.path(HERE, "niche_addon.R"))
keep <- vapply(ex, function(e) {
  s <- paste(deparse(e), collapse = " ")
  grepl("^(CH|SEEDB|NNULL|NALT|SC) <-", s) || grepl("^SC\\$", s) || grepl("^niche_defs <-", s)
}, TRUE)
for (e in ex[keep]) eval(e, globalenv())
SC <- SC[SC$n == 200L, ]
stopifnot(nrow(SC) == 7L)

jobs <- list()
for (r in seq_len(MAXR)) for (i in seq_len(nrow(SC))) {
  sc <- SC[i, ]
  f <- file.path(OUTD, sprintf("%s_r%04d.csv", sc$scen, r))
  if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(sc = sc, r = r, f = f)
}
cat(sprintf("niche_bagoft_fast | %d scenarios | %d reps | %d jobs left | %d workers | %s\n",
            nrow(SC), MAXR, length(jobs), NW, format(Sys.time())))
flush.console()
if (!length(jobs)) { cat("nothing to do\n"); quit(save = "no") }

cl <- parallel::makeCluster(min(NW, length(jobs)))
parallel::clusterExport(cl, c("HERE", "niche_defs", "SEEDB", "NSPL", "NSIM"))
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages({ library(randomForest); library(dcov) })
  source(file.path(HERE, "bagoft_fast_pre.R"))
  eval(niche_defs)
  one <- function(sc, r, f) {
    s <- SEEDB + sc$id * 100000L + r
    set.seed(s); dat <- gen(sc)
    dat <- as.data.frame(lapply(dat, as.vector))
    fo  <- frm(dat); environment(fo) <- globalenv()
    set.seed(s + 500000000L)
    t0 <- proc.time()[["elapsed"]]
    b  <- tryCatch(suppressWarnings(bagoft_fast_pre(dat, fo, nsplits = NSPL, nsim = NSIM, ntree = 60L)),
                   error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    g  <- function(x) if (!is.null(x) && is.numeric(x) && length(x) == 1) as.numeric(x) else NA_real_
    row <- if (inherits(b, "error"))
      cbind(key_of(sc, r, dat, s), p = NA, p_med = NA, p_min = NA, bag_pmin_stat = NA, time_s = el,
            status = "error", err = substr(conditionMessage(b), 1, 200), impl = "fast_pre")
    else
      cbind(key_of(sc, r, dat, s), p = g(b$p.value), p_med = g(b$p.value2), p_min = g(b$p.value3),
            bag_pmin_stat = g(b$pmin), time_s = el, status = "ok", err = "", impl = "fast_pre")
    tmp <- paste0(f, ".part"); utils::write.csv(row, tmp, row.names = FALSE)
    if (file.exists(f)) file.remove(f); file.rename(tmp, f)
    NULL
  }
  NULL }))
t0 <- Sys.time(); done <- 0L
for (bt in split(jobs, ceiling(seq_along(jobs) / (2L * NW)))) {
  parallel::parLapplyLB(cl, bt, function(jb) one(jb$sc, jb$r, jb$f))
  done <- done + length(bt)
  el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  cat(sprintf("%s  %d / %d datasets | %.2f h | ETA %s\n", format(Sys.time(), "%Y-%m-%d %H:%M"), done, length(jobs), el,
              format(Sys.time() + (length(jobs) - done) * el / done * 3600, "%Y-%m-%d %H:%M")))
  flush.console()
}
parallel::stopCluster(cl)
cat("niche_bagoft_fast DONE", format(Sys.time()), "\n")
