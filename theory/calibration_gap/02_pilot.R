## 02_pilot.R -- the calibration-gap pilot. Pre-declared protocol, design D1, n in {50, 200}.
##
##   L0 = law of T(X, y') with y' ~ Bern(pi_theta0(x_i))      (true null law, theta0 captured
##        from the generator)
##   L1 = law of T(X, y*) with y* ~ Bern(pi_thetahat(x_i))    (the deployed bootstrap reference)
##
## Both at the SAME observed X, both through the whole deployed pipeline (separation policy,
## rank bins, residual map, shipped scaler, shipped network). M scores each; two-sample
## empirical DFs; Delta_n, Delta_n^+, eps_tau (tau = .30); DKW inflation at gamma = .01.
##
## CHECKPOINTED AND RESUMABLE. The moment a dataset finishes, its row is (a) written atomically
## to rows/n<N>_rep<i>.csv and (b) appended to per_dataset.csv under a lock. On start, any
## dataset whose (n, rep, M) row is already on disk is skipped, so the run can be killed and
## restarted with the same command and it picks up where it stopped. cell_summary.csv and
## run_meta.rds are rebuilt from the rows at the end of every invocation.
##
## Run:  Rscript 02_pilot.R      env: CG_M, CG_D, CG_W, CG_N, CG_CHUNK, CG_OUTDIR

suppressMessages(library(brglm2))
RD  <- "/path/to/project/paper_deepgof/code/rdeploy"
PAP <- "/path/to/project/paper_deepgof"
SRC <- file.path(PAP, "theory", "calibration_gap")          # where the scripts live
OUT <- Sys.getenv("CG_OUTDIR", "")                          # where the outputs go
if (!nzchar(OUT)) OUT <- SRC

M      <- as.integer(Sys.getenv("CG_M", "100000"))     # scores per law per dataset
D      <- as.integer(Sys.getenv("CG_D", "40"))         # datasets per cell
NW     <- as.integer(Sys.getenv("CG_W", "12"))
CELLS  <- as.integer(strsplit(Sys.getenv("CG_N", "50,200"), ",")[[1]])
CHUNK  <- as.integer(Sys.getenv("CG_CHUNK", "2500"))   # maps scored per network batch
TAU    <- 0.30
GAMMA  <- 0.01
ALPHA  <- 0.05
SEED0  <- 901000L
DKW    <- 2 * sqrt(log(4 / GAMMA) / (2 * M))

dir.create(file.path(OUT, "rows"), recursive = TRUE, showWarnings = FALSE)
PDCSV <- file.path(OUT, "per_dataset.csv")
PROG  <- file.path(OUT, "progress.log")

source(file.path(RD, "deepgof_forward.R"))
source(file.path(SRC, "01_batched_net.R"))
W  <- load_weights(file.path(RD, "weights.csv"))
ms <- as.matrix(utils::read.csv(file.path(RD, "scaler.csv"), header = FALSE))
MU <- ms[1, ]; SDV <- ms[2, ]
BFAM <- binomial()

## ---- checkpointing -------------------------------------------------------------------------------
## A short-lived lock directory (dir.create is atomic on NTFS) serialises the appends to
## per_dataset.csv; the per-dataset row file is written to a temp name and renamed, so a row
## file is never seen half-written.
with_lock <- function(dir, expr) {
  lk <- file.path(dir, ".pdlock")
  got <- FALSE
  for (k in 1:6000) {
    if (isTRUE(suppressWarnings(dir.create(lk, showWarnings = FALSE)))) { got <- TRUE; break }
    Sys.sleep(0.02)
  }
  on.exit(if (got) unlink(lk, recursive = TRUE), add = TRUE)
  force(expr)
}

hb <- function(msg) {
  cat(sprintf("%s pid%-6d %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getpid(), msg),
      file = PROG, append = TRUE)
}

save_row <- function(row) {
  df <- as.data.frame(as.list(row), check.names = FALSE)
  names(df) <- sub("\\..*$", "", names(df))
  rf <- file.path(OUT, "rows", sprintf("n%d_rep%d.csv",
                                       as.integer(row[["n"]]), as.integer(row[["rep"]])))
  tf <- sprintf("%s.tmp%d", rf, Sys.getpid())
  utils::write.csv(df, tf, row.names = FALSE)
  if (file.exists(rf)) unlink(rf)
  file.rename(tf, rf)
  with_lock(OUT, {
    if (!file.exists(PDCSV)) utils::write.csv(df, PDCSV, row.names = FALSE)
    else utils::write.table(df, PDCSV, sep = ",", append = TRUE,
                            col.names = FALSE, row.names = FALSE)
  })
  invisible(TRUE)
}

## every (n, rep, M) already on disk, from the row files and from per_dataset.csv
done_keys <- function() {
  k <- character(0)
  fs <- list.files(file.path(OUT, "rows"), pattern = "^n\\d+_rep\\d+\\.csv$", full.names = TRUE)
  for (f in fs) {
    r <- tryCatch(utils::read.csv(f), error = function(e) NULL)
    if (!is.null(r) && nrow(r)) k <- c(k, sprintf("%d_%d_%d", r$n[1], r$rep[1], r$M[1]))
  }
  if (file.exists(PDCSV)) {
    r <- tryCatch(utils::read.csv(PDCSV), error = function(e) NULL)
    if (!is.null(r) && nrow(r)) k <- c(k, sprintf("%d_%d_%d", r$n, r$rep, r$M))
  }
  unique(k)
}

## ---- the deployed separation policy (rule verbatim from code/s3_grid.R) -------------------------
fast_fit <- function(y, mm) {
  g <- tryCatch(suppressWarnings(glm.fit(mm, y, family = BFAM)), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- g$coefficients
  firth <- 0L
  if (any(!is.finite(b)) || any(abs(b[-1]) > 8) || !g$converged) {
    f <- tryCatch(suppressWarnings(
           brglm2::brglmFit(x = mm, y = y, family = BFAM,
                            control = brglm2::brglmControl(type = "AS_mean"), intercept = TRUE)),
         error = function(e) NULL)
    if (!is.null(f) && all(is.finite(f$coefficients))) { g <- f; firth <- 1L }
    else return(list(ph = NULL, firth = 0L, coef = b))
  }
  list(ph = as.numeric(g$fitted.values), firth = firth, coef = g$coefficients)
}

## fixed model matrix + fixed 36-cell incidence (p = 2, so the shipped axis rule always returns
## c(x1, x2) and the rank bins depend only on X, which is held fixed within a dataset)
prep_D1 <- function(d, K = 6L) {
  mm <- model.matrix(~ x1 + x2, data = d)
  n  <- nrow(mm)
  bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
  cell <- (bn(mm[, "x1"]) - 1L) * K + bn(mm[, "x2"])
  A <- matrix(0, K * K, n)
  A[cbind(cell, seq_len(n))] <- 1
  list(mm = mm, A = A, n = n)
}
fast_map <- function(y, ph, P) {
  s  <- as.numeric(P$A %*% (y - ph))
  vv <- as.numeric(P$A %*% (ph * (1 - ph)))
  s / sqrt(pmax(vv, 1e-8))
}

## ---- generator: gen_design_null() of code/universal_design.R for D1, with theta0 captured -------
## (local copy; the original file is not touched. p = 2, law "indep", cat 0, rate .5.)
gen_D1_null <- function(n, rate = 0.5) {
  Z <- matrix(rnorm(n * 2L), n, 2L)
  Z[, 2] <- runif(n, -3, 3)
  d <- as.data.frame(Z); names(d) <- c("x1", "x2")
  MM  <- stats::model.matrix(~ ., data = d)
  bet <- c(0, stats::runif(ncol(MM) - 1L, -0.8, 0.8))
  lin <- as.numeric(MM %*% bet)
  cscl <- 1 / max(stats::sd(lin), 1e-8) * 1.2
  lin <- lin * cscl
  f  <- function(b0) mean(plogis(b0 + lin)) - rate
  b0 <- tryCatch(stats::uniroot(f, c(-25, 25))$root, error = function(e) 0)
  p0 <- plogis(b0 + lin)
  d$y <- stats::rbinom(n, 1, p0)
  list(data = d, p0 = p0, theta0 = c(b0, bet[-1] * cscl))
}

## ---- M scores from one law -----------------------------------------------------------------------
draw_law <- function(pvec, P, M, seed, tag = "") {
  set.seed(seed)
  n <- P$n
  S <- numeric(M); nfirth <- 0L; nfail <- 0L
  maps <- matrix(0, CHUNK, 36L)
  done <- 0L; nb <- 0L
  while (done < M) {
    cs <- min(CHUNK, M - done)
    bad <- logical(cs)
    for (j in seq_len(cs)) {
      y <- rbinom(n, 1L, pvec)
      f <- fast_fit(y, P$mm)
      if (is.null(f) || is.null(f$ph)) { bad[j] <- TRUE; nfail <- nfail + 1L; next }
      nfirth <- nfirth + f$firth
      maps[j, ] <- fast_map(y, f$ph, P)
    }
    sc <- deepgof_score_batch((maps[seq_len(cs), , drop = FALSE] -
                               rep(MU, each = cs)) / rep(SDV, each = cs), W, chunk = CHUNK)
    sc[bad] <- Inf                                   # shipped convention: a failed refit scores +Inf
    S[(done + 1L):(done + cs)] <- sc
    done <- done + cs
    nb <- nb + 1L
    if (nb %% 4L == 0L || done == M) hb(sprintf("%s %d/%d scores", tag, done, M))
  }
  list(S = S, firth = nfirth / M, fail = nfail / M)
}

## ---- the gap functionals -------------------------------------------------------------------------
gap_stats <- function(z0, z1, tau = TAU) {
  M <- length(z0)
  s0 <- sort(z0); s1 <- sort(z1)
  u  <- sort(unique(c(s0, s1)))
  F0 <- findInterval(u, s0) / M                       # #{z0 <= u}/M
  F1 <- findInterval(u, s1) / M
  d  <- F1 - F0
  elig <- F0 >= 1 - tau                               # 1 - F0(t) <= tau, t just above u
  c(Delta   = max(abs(d)),
    Deltap  = max(0, max(d)),
    epstau  = if (any(elig)) max(0, max(d[elig])) else 0)
}

## ---- one dataset ----------------------------------------------------------------------------------
one_dataset <- function(i, n) {
  t0 <- Sys.time()
  hb(sprintf("START n=%d rep=%d", n, i))
  set.seed(SEED0 + i)
  g <- gen_D1_null(n)
  P <- prep_D1(g$data)
  fo <- fast_fit(g$data$y, P$mm)
  if (is.null(fo) || is.null(fo$ph)) { hb(sprintf("SKIP  n=%d rep=%d observed fit failed", n, i)); return(NULL) }
  Tobs <- deepgof_score_batch(matrix((fast_map(g$data$y, fo$ph, P) - MU) / SDV, 1L, 36L), W)

  cellcode <- if (n == 50L) 0L else 1L
  L0 <- draw_law(g$p0, P, M, SEED0 + i + 3000000L + 10000000L * cellcode,
                 sprintf("n=%d rep=%d L0", n, i))
  L1 <- draw_law(fo$ph, P, M, SEED0 + i + 6000000L + 10000000L * cellcode,
                 sprintf("n=%d rep=%d L1", n, i))

  ## the M scores of each law for the first dataset of a cell, kept for plotting
  if (i == 1L) {
    gzf <- file.path(OUT, sprintf("scores_n%d_rep1.csv.gz", n))
    if (!file.exists(gzf)) {
      gz <- gzfile(gzf, "w")
      utils::write.csv(data.frame(law = rep(c("L0", "L1"), each = M),
                                  score = c(L0$S, L1$S)), gz, row.names = FALSE)
      close(gz)
    }
  }

  gs <- gap_stats(L0$S, L1$S)
  V  <- mean(L1$S >= Tobs)
  V0 <- mean(L0$S >= Tobs)
  row <- c(n = n, rep = i, M = M,
           Tobs = Tobs,
           Delta = gs["Delta"], Deltap = gs["Deltap"], epstau = gs["epstau"],
           Delta_dkw = gs["Delta"] + DKW, Deltap_dkw = gs["Deltap"] + DKW,
           epstau_dkw = gs["epstau"] + DKW,
           V = V, V0 = V0, VmV0 = V - V0,
           rej_boot = as.numeric(V <= ALPHA), rej_true = as.numeric(V0 <= ALPHA),
           firth_obs = fo$firth, firth_L0 = L0$firth, firth_L1 = L1$firth,
           fail_L0 = L0$fail, fail_L1 = L1$fail,
           event_rate = mean(g$data$y),
           b0 = g$theta0[1], beta1 = g$theta0[2], beta2 = g$theta0[3],
           b0hat = fo$coef[1], beta1hat = fo$coef[2], beta2hat = fo$coef[3])
  save_row(row)
  hb(sprintf("DONE  n=%d rep=%d  %.1f min", n, i,
             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  row
}

## ---- run -------------------------------------------------------------------------------------------
cat(sprintf("calibration-gap pilot | D1 (p=2, indep, 50%% events) | M=%d D=%d tau=%.2f gamma=%.2f\n",
            M, D, TAU, GAMMA))
cat(sprintf("DKW term 2*sqrt(log(4/gamma)/(2M)) = %.5f | workers = %d | chunk = %d\n", DKW, NW, CHUNK))
cat(sprintf("outputs -> %s\nstarted %s\n", OUT, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
HAVE <- done_keys()
cat(sprintf("resume: %d dataset rows already on disk\n", length(HAVE)))
TODO <- lapply(CELLS, function(n) {
  need <- Filter(function(i) !(sprintf("%d_%d_%d", n, i, M) %in% HAVE), seq_len(D))
  cat(sprintf("  cell n=%3d: %d of %d datasets to compute\n", n, length(need), D))
  need
})
names(TODO) <- as.character(CELLS)
flush.console()

if (sum(lengths(TODO)) > 0L) {
  cl <- parallel::makeCluster(NW)
  parallel::clusterExport(cl, c("RD", "SRC", "OUT", "PDCSV", "PROG", "M", "TAU", "GAMMA", "ALPHA",
                                "SEED0", "CHUNK", "DKW", "with_lock", "hb", "save_row",
                                "fast_fit", "prep_D1", "fast_map", "gen_D1_null", "draw_law",
                                "gap_stats", "one_dataset"), envir = environment())
  invisible(parallel::clusterEvalQ(cl, {
    suppressMessages(library(brglm2))
    source(file.path(RD, "deepgof_forward.R"))
    source(file.path(SRC, "01_batched_net.R"))
    W  <- load_weights(file.path(RD, "weights.csv"))
    ms <- as.matrix(utils::read.csv(file.path(RD, "scaler.csv"), header = FALSE))
    MU <- ms[1, ]; SDV <- ms[2, ]
    BFAM <- binomial()
    NULL }))
} else cl <- NULL

timing <- list()
for (n in CELLS) {
  need <- TODO[[as.character(n)]]
  if (!length(need)) { cat(sprintf("cell n=%3d already complete\n", n)); flush.console(); next }
  t0 <- Sys.time()
  invisible(parallel::parLapplyLB(cl, need, one_dataset, n = n))
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  timing[[as.character(n)]] <- el
  cat(sprintf("cell n=%3d done: %d datasets computed this invocation, %.1f min (%s)\n",
              n, length(need), el, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  flush.console()
}
if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)

## ---- rebuild per_dataset.csv, cell_summary.csv and run_meta.rds from the checkpointed rows ------
fs <- list.files(file.path(OUT, "rows"), pattern = "^n\\d+_rep\\d+\\.csv$", full.names = TRUE)
ALLROWS <- if (length(fs)) do.call(rbind, lapply(fs, function(f)
             tryCatch(utils::read.csv(f), error = function(e) NULL))) else NULL
if (is.null(ALLROWS) && file.exists(PDCSV))
  ALLROWS <- tryCatch(utils::read.csv(PDCSV), error = function(e) NULL)
if (!is.null(ALLROWS) && nrow(ALLROWS)) {
  ALLROWS <- ALLROWS[ALLROWS$M == M, , drop = FALSE]
  ALLROWS <- ALLROWS[!duplicated(ALLROWS[, c("n", "rep")]), , drop = FALSE]
  ALLROWS <- ALLROWS[order(ALLROWS$n, ALLROWS$rep), , drop = FALSE]
  with_lock(OUT, utils::write.csv(ALLROWS, PDCSV, row.names = FALSE))
}
saveRDS(list(M = M, D = D, NW = NW, CHUNK = CHUNK, TAU = TAU, GAMMA = GAMMA, DKW = DKW,
             CELLS = CELLS, timing = timing, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        file.path(OUT, "run_meta.rds"))

q <- function(v, p) as.numeric(stats::quantile(v, p, type = 7, names = FALSE))
if (!is.null(ALLROWS) && nrow(ALLROWS)) {
  SUM <- do.call(rbind, lapply(sort(unique(ALLROWS$n)), function(n) {
    R <- ALLROWS[ALLROWS$n == n, ]
    data.frame(n = n, D = nrow(R), M = M, DKW = DKW,
               eps_med = median(R$epstau), eps_q90 = q(R$epstau, .90), eps_q95 = q(R$epstau, .95),
               epsdkw_med = median(R$epstau_dkw), epsdkw_q90 = q(R$epstau_dkw, .90),
               epsdkw_q95 = q(R$epstau_dkw, .95),
               dp_med = median(R$Deltap), dp_q90 = q(R$Deltap, .90), dp_q95 = q(R$Deltap, .95),
               dpdkw_q90 = q(R$Deltap_dkw, .90), dpdkw_q95 = q(R$Deltap_dkw, .95),
               D2_med = median(R$Delta), D2_q90 = q(R$Delta, .90), D2_q95 = q(R$Delta, .95),
               mean_VmV0 = mean(R$VmV0),
               size_boot = mean(R$rej_boot), size_true = mean(R$rej_true),
               firth_obs = mean(R$firth_obs), firth_L0 = mean(R$firth_L0),
               firth_L1 = mean(R$firth_L1),
               fail_L0 = mean(R$fail_L0), fail_L1 = mean(R$fail_L1),
               cert_q90 = ALPHA + q(R$epstau, .90), cert_q90_dkw = ALPHA + q(R$epstau_dkw, .90),
               mins = if (!is.null(timing[[as.character(n)]])) timing[[as.character(n)]] else NA_real_)
  }))
  utils::write.csv(SUM, file.path(OUT, "cell_summary.csv"), row.names = FALSE)
  print(t(SUM))
}
cat(sprintf("\nPILOT INVOCATION COMPLETE -- %d rows on disk at M=%d (protocol D=%d per cell)\n",
            if (is.null(ALLROWS)) 0L else nrow(ALLROWS), M, D))
