## 00_verify_and_bench.R -- calibration-gap pilot, step 0.
##
## Three things, all of which must pass before a single score of the pilot is computed:
##   (1) the shipped base-R forward pass reproduces rdeploy/test_output.csv;
##   (2) a FAST re-expression of the deployed per-replicate score (fixed model matrix, fixed rank
##       bins, glm.fit / brglmFit.fit instead of glm(formula, data)) agrees BIT-FOR-BIT with the
##       shipped .deepgof_map + deepgof_score path on real D1 datasets;
##   (3) timings, so M can be chosen honestly.
##
## Why (2) is legitimate and not a re-implementation of the test: for design D1, p = 2, so
## .deepgof_map's axis rule ("two largest |beta_j| * sd(x_j), returned in model order") selects
## BOTH covariates and returns them as c(x1, x2) for every fit, always. The rank bins therefore
## depend only on X, which is held fixed across all M replicates of a dataset. The fast path
## precomputes them; everything downstream (residual map, scaler, network) is the shipped code.

suppressMessages(library(brglm2))
RD  <- "/path/to/project/paper_deepgof/code/rdeploy"
PAP <- "/path/to/project/paper_deepgof"
OUT <- file.path(PAP, "theory", "calibration_gap")

source(file.path(RD, "deepgof_forward.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
## deepgof.R sources its own forward pass off sys.frame(1)$ofile; supply the pieces it needs and
## take .deepgof_map out of it verbatim.
.deepgof_map <- local({
  txt <- readLines(file.path(RD, "deepgof.R"))
  i0  <- grep("^\\.deepgof_map <- function", txt)
  i1  <- grep("^\\}", txt); i1 <- min(i1[i1 > i0])
  e   <- new.env(parent = globalenv())
  eval(parse(text = txt[i0:i1]), envir = e)
  get(".deepgof_map", envir = e)
})
cat("pulled .deepgof_map() out of the shipped deepgof.R\n")

## ---- (1) forward pass vs PyTorch ---------------------------------------------------------------
W  <- load_weights(file.path(RD, "weights.csv"))
ms <- as.matrix(utils::read.csv(file.path(RD, "scaler.csv"), header = FALSE))
MU <- ms[1, ]; SD <- ms[2, ]
Xt <- as.matrix(utils::read.csv(file.path(RD, "test_input.csv"), header = FALSE))
Sp <- as.numeric(readLines(file.path(RD, "test_output.csv")))
Sr <- deepgof_scores(Xt, W)
v1 <- max(abs(Sr - Sp))
cat(sprintf("(1) forward pass vs test_output.csv: n=%d maps, max|diff| = %.3e\n", length(Sp), v1))

score_shipped <- function(m) deepgof_score((m - MU) / SD, W, 6L)

## ---- the deployed separation policy (verbatim from code/s3_grid.R) ------------------------------
fit_pol <- function(frm, d) {
  g <- tryCatch(suppressWarnings(glm(frm, data = d, family = binomial())), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1]) > 8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(frm, data = d, family = binomial(),
                                       method = "brglmFit", type = "AS_mean")), error = function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(list(fit = f, firth = 1L))
  }
  list(fit = g, firth = 0L)
}

## ---- the fast path -----------------------------------------------------------------------------
## precompute, once per dataset: model matrix, the 36-cell incidence matrix, the scaler.
prep_D1 <- function(d, K = 6L) {
  mm <- model.matrix(~ x1 + x2, data = d)
  n  <- nrow(mm)
  bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
  cell <- (bn(mm[, "x1"]) - 1L) * K + bn(mm[, "x2"])
  A <- matrix(0, K * K, n)
  A[cbind(cell, seq_len(n))] <- 1
  list(mm = mm, A = A, n = n, K = K)
}
## one score, fast: y -> fit (with the separation rule) -> map -> scaler -> network
fast_fit <- function(y, P) {
  g <- tryCatch(suppressWarnings(glm.fit(P$mm, y, family = binomial())), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- g$coefficients
  firth <- 0L
  if (any(!is.finite(b)) || any(abs(b[-1]) > 8) || !g$converged) {
    f <- tryCatch(suppressWarnings(
           brglm2::brglmFit(x = P$mm, y = y, family = binomial(),
                            control = brglm2::brglmControl(type = "AS_mean"), intercept = TRUE)),
         error = function(e) NULL)
    if (!is.null(f) && all(is.finite(f$coefficients))) { g <- f; firth <- 1L }
  }
  list(ph = as.numeric(g$fitted.values), firth = firth)
}
fast_map <- function(y, ph, P) {
  r  <- y - ph
  s  <- as.numeric(P$A %*% r)
  vv <- as.numeric(P$A %*% (ph * (1 - ph)))
  s / sqrt(pmax(vv, 1e-8))
}
fast_score <- function(y, P) {
  f <- fast_fit(y, P)
  if (is.null(f)) return(list(S = Inf, firth = 0L, fail = 1L))
  list(S = score_shipped(fast_map(y, f$ph, P)), firth = f$firth, fail = 0L)
}

## ---- the generator, with theta0 captured (local copy; the original is untouched) ----------------
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

## ---- (2) fast path vs shipped path --------------------------------------------------------------
set.seed(7)
dif <- numeric(0); nfirth <- 0L; ntot <- 0L
for (n in c(50L, 200L)) {
  for (rep in 1:8) {
    g <- gen_D1_null(n)
    P <- prep_D1(g$data)
    for (b in 1:30) {
      yb <- if (b == 1L) g$data$y else rbinom(n, 1, g$p0)
      if (length(unique(yb)) < 2L) next
      dd <- g$data; dd$y <- yb
      fp <- fit_pol(y ~ x1 + x2, dd)
      if (is.null(fp)) next
      S_ship <- score_shipped(.deepgof_map(fp$fit, 6L)$map)
      fs <- fast_score(yb, P)
      dif <- c(dif, abs(S_ship - fs$S))
      nfirth <- nfirth + (fp$firth != fs$firth)
      ntot <- ntot + 1L
    }
  }
}
cat(sprintf("(2) fast vs shipped: %d comparisons, max|diff| = %.3e, firth-flag mismatches = %d\n",
            ntot, max(dif), nfirth))

## ---- (3) timing ---------------------------------------------------------------------------------
tm <- list()
for (n in c(50L, 200L)) {
  g <- gen_D1_null(n); P <- prep_D1(g$data)
  R <- 400L
  ys <- lapply(seq_len(R), function(i) rbinom(n, 1, g$p0))
  t0 <- Sys.time(); for (i in seq_len(R)) fast_score(ys[[i]], P)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / R
  ## split: fit only, and network only
  t1 <- Sys.time(); for (i in seq_len(R)) fast_fit(ys[[i]], P)
  e1 <- as.numeric(difftime(Sys.time(), t1, units = "secs")) / R
  m0 <- fast_map(ys[[1]], fast_fit(ys[[1]], P)$ph, P)
  t2 <- Sys.time(); for (i in seq_len(R)) score_shipped(m0)
  e2 <- as.numeric(difftime(Sys.time(), t2, units = "secs")) / R
  tm[[as.character(n)]] <- c(total = el, fit = e1, net = e2)
  cat(sprintf("(3) n=%3d : %.3f ms / score   (fit %.3f ms, network %.3f ms)\n",
              n, 1000*el, 1000*e1, 1000*e2))
}
## shipped-path timing, for the record
g <- gen_D1_null(200L)
t3 <- Sys.time(); for (i in 1:100) { dd <- g$data; dd$y <- rbinom(200, 1, g$p0)
  fp <- fit_pol(y ~ x1 + x2, dd); score_shipped(.deepgof_map(fp$fit, 6L)$map) }
cat(sprintf("(3) shipped path n=200: %.3f ms / score\n",
            1000*as.numeric(difftime(Sys.time(), t3, units="secs"))/100))

saveRDS(list(v1 = v1, fast_max_diff = max(dif), ncmp = ntot, firth_mismatch = nfirth, timing = tm),
        file.path(OUT, "verify.rds"))
cat("\nverify.rds written\n")
