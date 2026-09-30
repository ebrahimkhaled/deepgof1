# =============================================================================
# run_phd.R -- supervised map axes on Liu Settings 2, 3 and 4.
#
# Four readings of ONE frozen network, on the SAME datasets, from the SAME bootstrap draws:
#   v1        the deployed rule: the two covariates with the largest |beta_j| * sd_j
#   phd       the top two eigenvectors of M = (1/n) sum_i r_i x_i x_i'   (Li 1992)
#   phdw      the same, built on the Pearson residual
#   allpairs  the max over every coordinate pair -- Setting 2 only (at p = 2 it IS v1)
#
# v1 is recomputed here rather than read from results/liu/: the deposited column was scored with a
# network trained per block (code/liu_score.py:65), not with the shipped weights, so taking it would
#
# The data-chosen rule is recomputed inside every bootstrap replicate. That is the whole reason the
# rank p-value needs no correction for having chosen: observed and replicates are treated alike.
#
# Run:  Rscript run_phd.R      Env: PH_WORKERS (18), PH_SEEDS, PH_N, PH_SET, PH_SMOKE
# =============================================================================
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
HERE <- file.path(PAP, "phd")
RD   <- file.path(PAP, "code", "rdeploy")
ROWD <- file.path(HERE, "rows"); dir.create(ROWD, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(brglm2))
source(file.path(HERE, "run_phd_defs.R"))     # gen_liu, fit_pol, pick_axes, map_uv, readings, K

B <- as.integer(Sys.getenv("PH_B", "199"))
NW <- as.integer(Sys.getenv("PH_WORKERS", "18"))
SEEDS <- as.integer(strsplit(Sys.getenv("PH_SEEDS", "101000,102000"), ",")[[1]])
NS  <- as.integer(strsplit(Sys.getenv("PH_N", "50,100,500"), ",")[[1]])
SETS <- as.integer(strsplit(Sys.getenv("PH_SET", "2,3,4"), ",")[[1]])
SMOKE <- Sys.getenv("PH_SMOKE", "0") == "1"
NNULL <- if (SMOKE) 4L else 2000L
NEVAL <- if (SMOKE) 2L else 1000L
# liu_grid.R seeds with ci = the index of C in c(0, .25, .5, .75, 1), so the departure strengths are
# indices 2..5, NOT 1..4. The all-pairs run got this wrong first time; verify_same_datasets.R caught
# it. Keep the full vector and index into it.
CS_ALL <- c(0, 0.25, 0.5, 0.75, 1.0)

one_dataset <- function(SEED0, NOBS, SET, kind, ci, C, i, score) {
  set.seed(SEED0 + SET * 500000L + kind * 50000L + ci * 5000L + i)   # liu_grid.R's seeding
  d <- gen_liu(NOBS, SET, if (kind == 0L) 0 else C)
  if (length(unique(d$y)) < 2L) return(NULL)
  xn <- setdiff(names(d), "y")
  frm <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
  fp <- fit_pol(frm, d); if (is.null(fp)) return(NULL)
  prs <- utils::combn(xn, 2, simplify = FALSE)
  obs <- tryCatch(readings(fp$fit, d, xn, prs, score), error = function(e) NULL)
  if (is.null(obs)) return(NULL)
  ph <- as.numeric(stats::fitted(fp$fit)); d2 <- d
  Sb <- matrix(NA_real_, B, length(obs), dimnames = list(NULL, names(obs)))
  for (b in seq_len(B)) {
    d2$y <- stats::rbinom(NOBS, 1, ph)
    fb <- fit_pol(frm, d2)
    ## a failed refit scores 0 on every reading, as liu_grid.R does: it cannot push a maximum up
    Sb[b, ] <- if (is.null(fb)) rep(0, length(obs)) else
      tryCatch(readings(fb$fit, d2, xn, prs, score), error = function(e) rep(0, length(obs)))
  }
  p <- (1 + colSums(sweep(Sb, 2, obs, ">="))) / (B + 1)
  cbind(data.frame(setting = SET, n = NOBS, seed = SEED0, C = if (kind == 0L) 0 else C,
                   rep = i, kind = kind, firth = fp$firth,
                   v1_axes = paste(pick_axes(fp$fit, d), collapse = "+")),
        as.data.frame(as.list(p)), label = "EXPLORATORY_NOT_CONFIRMATORY")
}

jobs <- list()
for (SET in SETS) for (SEED0 in SEEDS) for (NOBS in NS) {
  f <- file.path(ROWD, sprintf("S%d_n%d_s%d.csv", SET, NOBS, SEED0))
  if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(SET = SET, SEED0 = SEED0, NOBS = NOBS, f = f)
}
cat(sprintf("phd%s | %d blocks to run | B=%d | %d workers | %s\n", if (SMOKE) " SMOKE" else "",
            length(jobs), B, NW, format(Sys.time()))); flush.console()
if (!length(jobs)) { cat("nothing to do\n"); quit(save = "no") }

run_block <- function(j) {
  setwd(RD); source("deepgof_forward.R")
  W <- load_weights("weights.csv")
  ms <- as.matrix(utils::read.csv("scaler.csv", header = FALSE))
  mu <- ms[1, ]; sdv <- ms[2, ]
  score <- function(m) deepgof_score((m - mu) / sdv, W, K)
  out <- list()
  for (i in seq_len(NNULL))
    out[[length(out) + 1L]] <- one_dataset(j$SEED0, j$NOBS, j$SET, 0L, 0L, 0, i, score)
  for (ci in seq_along(CS_ALL)) if (CS_ALL[ci] > 0) for (i in seq_len(NEVAL))
    out[[length(out) + 1L]] <- one_dataset(j$SEED0, j$NOBS, j$SET, 2L, ci, CS_ALL[ci], i, score)
  R <- do.call(rbind, Filter(Negate(is.null), out))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  sprintf("%s  S%d n=%d seed=%d: %d rows", format(Sys.time()), j$SET, j$NOBS, j$SEED0, nrow(R))
}

cl <- parallel::makeCluster(min(NW, length(jobs)))
parallel::clusterExport(cl, c("PAP", "RD", "HERE", "K", "B", "CS_ALL", "NNULL", "NEVAL",
                              "one_dataset"))
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(library(brglm2)); source(file.path(HERE, "run_phd_defs.R")); NULL }))
t0 <- Sys.time()
msg <- parallel::parLapplyLB(cl, jobs, run_block)
parallel::stopCluster(cl)
cat(paste(unlist(msg), collapse = "\n"), "\n")
cat(sprintf("DONE %s: %.2f h\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
