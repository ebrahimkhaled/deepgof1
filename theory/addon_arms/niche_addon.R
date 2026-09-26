# =============================================================================
# niche_addon.R -- mHL and GiViTI (and, optionally, BAGofT) on EXACTLY the datasets of the niche
# study, theory/niche/rows/ (14 scenarios, 16,000 datasets).
#
# Same data. theory/niche/run_niche.R is hashed (PREDECLARATION.sha256) and is only READ here,
# never sourced or edited. Its scenario table, gen(), frm() and fit_pol() are copied VERBATIM
# below, with the same per-dataset seed:
#   set.seed(SEEDB + sc$id * 100000 + r),   SEEDB = 123000
# gen() is the first thing drawn after the seed and nothing after the fit uses the RNG, so the
# arms cannot shift the data. HL_check (HL recomputed here) must equal the stored HL column
# (verify_pairing.R checks it).
#
# Two modes (AD_MODE):
#   main    HL_check, Stukel (joint and marginal), mHL, GiViTI. Same blocks as run_niche.R:
#           one file per (scenario, 20 reps), niche_rows/<scen>_r<first rep>.csv, keyed
#           (scen, fam, C, n, rep) like theory/niche/rows/.
#   bagoft  BAGofT 1.0.0 at the settings of the paper's BAGofT arm (jasa_study bagoft_p() /
#           confirm_bagoft/run_confirm.R): testGlmBi(logit), parRF(ntree = 60), nsplits = 100,
#           nsim = 100, ne = floor(5 sqrt(n)); its own RNG seeded with dataset seed + 500000000,
#           as in run_confirm.R. The PACKAGE function is called, not bagoft_fast(): with ten
#           covariates parRF pre-selects with dcPre(), a branch bagoft_fast() does not reproduce.
#           Headline p-value = p.value (splits combined by the mean); p.value2, p.value3 and the
#           statistic pmin are kept, named as in the confirmatory rows. One file per dataset,
#           niche_rows_bagoft/<scen>_r<rep>.csv.
# Resumable in both modes: a file that exists is skipped; files are written .part then renamed.
#
# Run:   Rscript niche_addon.R
# Env:   AD_MODE  main (default) / bagoft
#        AD_W     workers (4)
#        AD_SCEN  comma list of scenarios (default all 14)
#        AD_MAXREP cap on reps per scenario, for smoke tests
#        AD_OUT   output folder (default niche_rows or niche_rows_bagoft)
# =============================================================================
HERE  <- "/path/to/project/paper_deepgof/theory/addon_arms"
MODE  <- Sys.getenv("AD_MODE", "main"); stopifnot(MODE %in% c("main", "bagoft"))
NW    <- as.integer(Sys.getenv("AD_W", "4"))
MAXR  <- as.integer(Sys.getenv("AD_MAXREP", "100000"))
OUTD  <- file.path(HERE, Sys.getenv("AD_OUT", if (MODE == "main") "niche_rows" else "niche_rows_bagoft"))
dir.create(OUTD, showWarnings = FALSE)

## ---- VERBATIM from run_niche.R (lines 22-32, full-study values) ----
CH    <- 20L
SEEDB <- 123000L
NNULL <- 2000L
NALT  <- 1000L
SC <- expand.grid(fam = c("inter", "quad"), C = c(0.3, 0.5, 0.8), n = c(200L, 500L), stringsAsFactors = FALSE)
SC <- rbind(data.frame(fam = "null", C = 0, n = c(200L, 500L)), SC)
SC$id   <- seq_len(nrow(SC))
SC$reps <- ifelse(SC$fam == "null", NNULL, NALT)
SC$scen <- ifelse(SC$fam == "null", sprintf("null_n%d", SC$n), sprintf("%s_C%.1f_n%d", SC$fam, SC$C, SC$n))
if (nzchar(Sys.getenv("AD_SCEN"))) SC <- SC[SC$scen %in% strsplit(Sys.getenv("AD_SCEN"), ",")[[1]], ]

niche_defs <- quote({
  ## ---- VERBATIM from run_niche.R 'defs' (lines 39-58) ----
  P <- 10L
  gen <- function(sc) {
    n <- sc$n; X <- matrix(stats::rnorm(n * P), n, P); colnames(X) <- paste0("x", 1:P)
    f <- switch(sc$fam, null = 0, inter = X[, 1] * X[, 2], quad = X[, 1]^2 - 1)
    eta <- -0.5 + 1.0 * X[, 1] + 0.8 * X[, 2] + 0.15 * rowSums(X[, 3:P]) + sc$C * f
    data.frame(X, y = stats::rbinom(n, 1, stats::plogis(eta)))
  }
  frm <- function(dat) stats::as.formula(paste("y ~", paste(setdiff(names(dat), "y"), collapse = " + ")))
  fit_pol <- function(dat) {
    f <- frm(dat)
    g <- tryCatch(suppressWarnings(stats::glm(f, data = dat, family = stats::binomial())), error = function(e) NULL)
    if (is.null(g)) return(NULL)
    if (!g$converged || any(!is.finite(stats::coef(g)))) {
      h <- tryCatch(suppressWarnings(stats::glm(f, data = dat, family = stats::binomial(), method = "brglmFit", type = "AS_mean")),
                    error = function(e) NULL)
      if (!is.null(h)) return(h)
    }
    g
  }
  key_of <- function(sc, r, dat, s) data.frame(scen = sc$scen, fam = sc$fam, C = sc$C, n = sc$n, rep = r,
                                               events = sum(dat$y), dseed = s)
  niche_main <- function(sc, r) {
    s <- SEEDB + sc$id * 100000L + r
    set.seed(s); dat <- gen(sc); fit <- fit_pol(dat)
    t0 <- proc.time()[["elapsed"]]
    a <- if (is.null(fit)) c(HL_check = NA, Stukel_joint = NA, Stukel_marginal = NA, mHL = NA, GiViTI = NA,
                             giviti_status = NA, giviti_secs = NA) else arms_for_fit(fit, do_giviti = TRUE)
    firth <- if (is.null(fit)) NA else as.numeric(inherits(fit, "brglmFit") || identical(fit$method, "brglmFit"))
    cbind(key_of(sc, r, dat, s), firth = firth, t(a), secs = proc.time()[["elapsed"]] - t0)
  }
  niche_bagoft <- function(sc, r) {
    s <- SEEDB + sc$id * 100000L + r
    set.seed(s); dat <- gen(sc)
    dat <- as.data.frame(lapply(dat, as.vector))
    fo  <- frm(dat); environment(fo) <- globalenv()
    set.seed(s + 500000000L)                                   # the BAGofT arm's own RNG stream
    t0 <- proc.time()[["elapsed"]]
    b <- tryCatch({
      invisible(utils::capture.output(res <- suppressMessages(suppressWarnings(
        BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = fo, link = "logit"),
                       parFun = BAGofT::parRF(ntree = 60), data = dat, nsplits = 100L,
                       ne = floor(5 * sqrt(nrow(dat))), nsim = 100L)))))
      res }, error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    g  <- function(x) if (!is.null(x) && is.numeric(x) && length(x) == 1) as.numeric(x) else NA_real_
    if (inherits(b, "error"))
      return(cbind(key_of(sc, r, dat, s), p = NA, p_med = NA, p_min = NA, bag_pmin_stat = NA, time_s = el,
                   status = "error", err = substr(conditionMessage(b), 1, 200)))
    cbind(key_of(sc, r, dat, s), p = g(b$p.value), p_med = g(b$p.value2), p_min = g(b$p.value3),
          bag_pmin_stat = g(b$pmin), time_s = el, status = "ok", err = "")
  }
})

jobs <- list()
for (i in seq_len(nrow(SC))) {
  sc <- SC[i, ]; R <- min(sc$reps, MAXR)
  if (MODE == "main") {
    for (s in seq(1L, R, by = CH)) {
      f <- file.path(OUTD, sprintf("%s_r%04d.csv", sc$scen, s))
      if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(sc = sc, reps = s:min(R, s + CH - 1L), f = f)
    }
  } else {
    ## rep-major round robin, so a cut-short BAGofT run still covers every scenario evenly
    for (r in seq_len(R)) {
      f <- file.path(OUTD, sprintf("%s_r%04d.csv", sc$scen, r))
      if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(sc = sc, reps = r, f = f)
    }
  }
}
if (MODE == "bagoft") jobs <- jobs[order(vapply(jobs, function(j) j$reps[1], 1L))]
cat(sprintf("niche_addon %s | %d scenarios | %d jobs | %d workers | %s\n", MODE, nrow(SC), length(jobs), NW,
            format(Sys.time())))
flush.console()
if (!length(jobs)) { cat("nothing to do\n"); quit(save = "no") }

cl <- parallel::makeCluster(min(NW, length(jobs)))
parallel::clusterExport(cl, c("HERE", "niche_defs", "SEEDB", "MODE"))
invisible(parallel::clusterEvalQ(cl, {
  if (MODE == "main") source(file.path(HERE, "addon_common.R")) else {
    suppressMessages({ library(BAGofT); library(randomForest) })
    write_atomic <- function(df, f) { tmp <- paste0(f, ".part"); utils::write.csv(df, tmp, row.names = FALSE)
      if (file.exists(f)) file.remove(f); file.rename(tmp, f) }
  }
  eval(niche_defs); NULL }))
t0 <- Sys.time(); done <- 0L
for (bt in split(jobs, ceiling(seq_along(jobs) / (4L * NW)))) {
  parallel::parLapplyLB(cl, bt, function(jb) {
    R <- do.call(rbind, lapply(jb$reps, function(r) if (MODE == "main") niche_main(jb$sc, r) else niche_bagoft(jb$sc, r)))
    write_atomic(R, jb$f); NULL
  })
  done <- done + length(bt)
  el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  cat(sprintf("%s  %d / %d jobs | %.2f h | ETA %s\n", format(Sys.time(), "%Y-%m-%d %H:%M"), done, length(jobs), el,
              format(Sys.time() + (length(jobs) - done) * el / done * 3600, "%H:%M")))
  flush.console()
}
if (MODE == "main") invisible(parallel::clusterEvalQ(cl, giviti_close()))
parallel::stopCluster(cl)
cat("niche_addon DONE", format(Sys.time()), "\n")
