# =============================================================================
# run_s1.R -- Liu Setting 1 with the SHIPPED frozen weights.
#
# At p = 1 there are no covariate pairs and no two directions, so all-pairs, pHd and the deployed
# axis rule are the same function of the same data. One arm is computed, named v1, and the map is
# liu_grid.R's one-dimensional construction: 36 quantile cells on the ranks of x1, row-major -- NOT
# a 6x6 grid over two axes.
#
# This completes the 45-cell shipped-versus-deposited comparison of
#
# for the run that has already finished.
#
# Run:  Rscript run_s1.R      Env: S1_WORKERS (6), S1_SMOKE
# =============================================================================
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
HERE <- file.path(PAP, "phd")
RD   <- file.path(PAP, "code", "rdeploy")
ROWD <- file.path(HERE, "rows_s1"); dir.create(ROWD, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(brglm2))
source(file.path(HERE, "run_phd_defs.R"))       # gen_liu, fit_pol, K -- unchanged and hashed

B  <- as.integer(Sys.getenv("S1_B", "199"))
NW <- as.integer(Sys.getenv("S1_WORKERS", "6"))
SMOKE <- Sys.getenv("S1_SMOKE", "0") == "1"
NNULL <- if (SMOKE) 4L else 2000L
NEVAL <- if (SMOKE) 2L else 1000L
CS_ALL <- c(0, 0.25, 0.5, 0.75, 1.0)            # ci indexes THIS vector, as liu_grid.R does
SET <- 1L

## the p = 1 map, verbatim from code/liu_grid.R's resid_map() single-axis branch
map_1d <- function(fit, d, ax) {
  ph <- as.numeric(stats::fitted(fit)); r <- d$y - ph; n <- length(r)
  i <- pmin(K * K, 1L + floor(K * K * (rank(d[[ax]], ties.method = "first") - 1) / n))
  idx <- factor(i, levels = 1:(K * K))
  s <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8))
}

one_dataset <- function(SEED0, NOBS, kind, ci, C, i, score) {
  set.seed(SEED0 + SET * 500000L + kind * 50000L + ci * 5000L + i)   # liu_grid.R's seeding
  d <- gen_liu(NOBS, SET, if (kind == 0L) 0 else C)
  if (length(unique(d$y)) < 2L) return(NULL)
  frm <- y ~ x1
  fp <- fit_pol(frm, d); if (is.null(fp)) return(NULL)
  obs <- tryCatch(score(map_1d(fp$fit, d, "x1")), error = function(e) NULL)
  if (is.null(obs)) return(NULL)
  ph <- as.numeric(stats::fitted(fp$fit)); d2 <- d; sb <- numeric(B)
  for (b in seq_len(B)) {
    d2$y <- stats::rbinom(NOBS, 1, ph)
    fb <- fit_pol(frm, d2)
    sb[b] <- if (is.null(fb)) 0 else tryCatch(score(map_1d(fb$fit, d2, "x1")), error = function(e) 0)
  }
  data.frame(setting = SET, n = NOBS, seed = SEED0, C = if (kind == 0L) 0 else C,
             rep = i, kind = kind, firth = fp$firth, v1_axes = "x1",
             v1 = (1 + sum(sb >= obs)) / (B + 1), label = "SHIPPED_WEIGHTS")
}

jobs <- list()
for (SEED0 in c(101000L, 102000L)) for (NOBS in c(50L, 100L, 500L)) {
  f <- file.path(ROWD, sprintf("S1_n%d_s%d.csv", NOBS, SEED0))
  if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(SEED0 = SEED0, NOBS = NOBS, f = f)
}
cat(sprintf("S1%s | %d blocks | B=%d | %d workers | %s\n", if (SMOKE) " SMOKE" else "",
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
    out[[length(out) + 1L]] <- one_dataset(j$SEED0, j$NOBS, 0L, 0L, 0, i, score)
  for (ci in seq_along(CS_ALL)) if (CS_ALL[ci] > 0) for (i in seq_len(NEVAL))
    out[[length(out) + 1L]] <- one_dataset(j$SEED0, j$NOBS, 2L, ci, CS_ALL[ci], i, score)
  R <- do.call(rbind, Filter(Negate(is.null), out))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  sprintf("%s  S1 n=%d seed=%d: %d rows", format(Sys.time()), j$NOBS, j$SEED0, nrow(R))
}

cl <- parallel::makeCluster(min(NW, length(jobs)))
## CS_ALL is exported. Omitting it is what killed the all-pairs run on 2026-09-20, after every
parallel::clusterExport(cl, c("PAP", "RD", "HERE", "B", "SET", "CS_ALL", "NNULL", "NEVAL",
                              "map_1d", "one_dataset"))
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(library(brglm2)); source(file.path(HERE, "run_phd_defs.R")); NULL }))
t0 <- Sys.time()
msg <- parallel::parLapplyLB(cl, jobs, run_block)
parallel::stopCluster(cl)
cat(paste(unlist(msg), collapse = "\n"), "\n")
cat(sprintf("DONE %s: %.2f h\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
