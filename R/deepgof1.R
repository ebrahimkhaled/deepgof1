## deepgof1.R -- DeepGOF-1: the practitioner-facing test. The whole user interface:
##
##     source("R/deepgof1.R")          # from the release root
##     fit <- glm(y ~ x1 + x2, data = d, family = binomial())
##     deepgof1(fit)
##
## No training, no Python, no GPU, no tuning. The network was trained ONCE, offline, on simulated
## data and ships as 18,273 frozen numbers (inst/extdata/deepgof1_weights.csv); the analyst only
## runs the forward pass. Runtime ~5 s per test (200 forward passes + 199 logistic refits).
##
## What happens per call:
##   1. build the residual map of the fitted model (6 x 6 grid over covariate RANKS)
##   2. score it with the frozen network                       -> S_obs
##   3. simulate B datasets from the fitted model, refit each, score each -> S_1 .. S_B
##   4. p = (1 + #{S_b >= S_obs}) / (B + 1)
##
## Step 3 makes the level a rank within the analyst's own bootstrap: it does not depend on the
## network being right about anything (Lemma S.2 / Proposition S.3 of the paper's supplement).
##
## TIES (version 2 of this archive). Ties among covariate values are broken at random, once per
## call: a random permutation of the rows is drawn only when some covariate column has a tied
## value, and the same permutation serves the observed map and every bootstrap map. This is the
## rule of ebrahim.gof 2.8.0 and of the paper's Algorithm 1. With continuous covariates nothing is
## drawn and the p-value is bit-identical to the earlier row-order rule. ties = "roworder" restores
## the earlier rule (rank(ties.method = "first")); it exists only to reproduce numbers computed
## before the fix, and must not be used on real data: a file sorted by the outcome then puts tied
## rows into cells by their outcome.
##
## VERSIONING. This file defines deepgof1() -- the frozen first model. deepgof() is an alias for
## the current default version. A DeepGOF-2 ships new weights + a deepgof2() with this interface.

## locate this file so the release works from any working directory: source()'s frame carries
## `ofile`; scripts get it from --file=; the last resort assumes the release root is the cwd.
.dg1_here <- local({
  of <- tryCatch(sys.function(-1), error = function(e) NULL)   # placate R CMD check
  h <- NULL
  for (fr in rev(sys.frames())) {
    o <- tryCatch(get("ofile", envir = fr), error = function(e) NULL)
    if (is.character(o)) { h <- dirname(normalizePath(o)); break }
  }
  if (is.null(h)) {
    fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    h <- if (length(fa)) dirname(normalizePath(fa[1])) else file.path(getwd(), "R")
  }
  h
})
source(file.path(.dg1_here, "deepgof1_forward.R"))
.deepgof1_root <- normalizePath(file.path(.dg1_here, ".."))

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- the residual map, identical to the training-time builder ----------------------------------
## the axis rule (reading = "v1"): the two model-matrix columns with the largest |beta_j| * sd(x_j)
.deepgof1_map_core <- function(mm, y, ph, coefs, K = 6L, tb = NULL) {
  vars <- colnames(mm)[-1]
  if (length(vars) < 2L) stop("deepgof1(): needs at least two covariates")
  b  <- coefs[vars]
  sc <- abs(b) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
  ax <- vars[sort(order(sc, decreasing = TRUE)[1:2])]
  r  <- y - ph
  n  <- length(r)
  ## tb: the tie-breaking permutation drawn by deepgof1(); NULL = ties by row order
  rk <- function(v) {
    if (is.null(tb)) return(rank(v, ties.method = "first"))
    o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr
  }
  bn <- function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n))
  idx <- factor((bn(mm[, ax[1]]) - 1L) * K + bn(mm[, ax[2]]), levels = 1:(K * K))
  s  <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  list(map = as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8)), axes = ax)
}
.deepgof1_map <- function(fit, K = 6L, tb = NULL)
  .deepgof1_map_core(stats::model.matrix(fit), fit$y, as.numeric(stats::fitted(fit)),
                     stats::coef(fit), K, tb)

## the covariates of the model, untransformed: ns(x, 3), poly(x, 2) and I(x^2) are ranked by x, and an
## interaction adds no axis. Transformed covariates come from the data the model was fitted to, on the
## rows the fit kept; factors enter by their level codes. NULL when a covariate cannot be recovered.
.deepgof1_covariates <- function(fit) {
  mf <- stats::model.frame(fit)
  tl <- attr(stats::terms(mf), "term.labels")
  vars <- unique(unlist(lapply(tl, function(l) all.vars(str2lang(l)))))
  if (!length(vars)) return(matrix(numeric(0), nrow(mf), 0L))
  src <- fit$data
  get1 <- function(v) {
    if (v %in% names(mf) && is.null(dim(mf[[v]]))) return(mf[[v]])
    if (is.data.frame(src) && v %in% names(src)) {
      if (nrow(src) == nrow(mf) && is.null(fit$na.action)) return(src[[v]])
      ids <- rownames(mf)
      if (!anyDuplicated(rownames(src)) && all(ids %in% rownames(src)))
        return(src[ids, v, drop = TRUE])
    }
    if (is.null(fit$na.action)) {
      e <- tryCatch(get(v, envir = environment(stats::formula(fit))), error = function(e) NULL)
      if (is.atomic(e) && is.null(dim(e)) && length(e) == nrow(mf)) return(e)
    }
    NULL
  }
  cols <- lapply(vars, get1)
  if (any(vapply(cols, is.null, TRUE))) return(NULL)
  X <- vapply(cols, function(v) {
    if (is.character(v)) v <- factor(v)
    as.numeric(if (is.factor(v)) as.integer(v) else v)
  }, numeric(nrow(mf)))
  X <- matrix(X, nrow(mf), dimnames = list(NULL, vars))
  X[, apply(X, 2, function(v) length(unique(v)) > 1L), drop = FALSE]
}

## the cells of every pair (or, with one covariate, 36 quantile cells along it), found once per call
.deepgof1_cells <- function(X, K = 6L, tb = NULL) {
  n <- nrow(X)
  rk <- function(v) {
    if (is.null(tb)) return(rank(v, ties.method = "first"))
    o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr
  }
  if (ncol(X) == 1L) {
    i <- pmin(K * K, 1L + floor(K * K * (rk(X[, 1]) - 1) / n))
    return(list(idx = list(as.integer(i)), pairs = matrix(colnames(X), 1L, 2L)))
  }
  bins <- apply(X, 2, function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n)))
  prs <- utils::combn(colnames(X), 2L)
  list(idx = lapply(seq_len(ncol(prs)), function(j)
         as.integer((bins[, prs[1L, j]] - 1L) * K + bins[, prs[2L, j]])),
       pairs = t(prs))
}
.deepgof1_cellmap <- function(i, r, ph, K = 6L) {
  f <- factor(i, levels = seq_len(K * K))
  s  <- vapply(split(r, f), sum, 0)
  vv <- vapply(split(ph * (1 - ph), f), sum, 0)
  s / sqrt(pmax(vv, 1e-8))
}

## ---- the test ----------------------------------------------------------------------------------
## reading = "v1" (default): one map over the two strongest columns. reading = "allpairs": the largest
## score over the maps of every pair of covariates, taken again in every bootstrap replicate. With one
## covariate both readings are the 36-cell map along it. The same code as deepgof1() in ebrahim.gof 2.9.0.
deepgof1 <- function(fit, B = 199L, K = 6L, weights_dir = NULL, ties = c("random", "roworder"),
                     reading = c("v1", "allpairs"), covariates = NULL) {
  stopifnot(inherits(fit, "glm"), fit$family$family == "binomial")
  if (any(fit$prior.weights != 1) || !all(fit$y %in% c(0, 1)))
    stop("deepgof1(): needs a 0/1 outcome, one row per observation, without prior weights")
  ties <- match.arg(ties); reading <- match.arg(reading)
  if (!is.null(covariates) && reading == "v1") stop("'covariates' applies to reading = \"allpairs\" only")
  wd <- weights_dir %||% file.path(.deepgof1_root, "inst", "extdata")
  W  <- load_weights(file.path(wd, "deepgof1_weights.csv"))
  ms <- as.matrix(utils::read.csv(file.path(wd, "deepgof1_scaler.csv"), header = FALSE))
  mu <- ms[1, ]; sdv <- ms[2, ]
  score <- function(m) deepgof_score((m - mu) / sdv, W, K)

  mm <- stats::model.matrix(fit)
  if (reading == "allpairs") {
    X <- .deepgof1_covariates(fit)
    if (is.null(X)) {
      warning("deepgof1(): could not recover the covariates; the maps use model-matrix columns")
      v <- colnames(mm)[colnames(mm) != "(Intercept)"]
      X <- mm[, v[!is.na(stats::coef(fit)[v])], drop = FALSE]
    }
    if (!is.null(covariates)) X <- X[, covariates, drop = FALSE]
    if (ncol(X) < 1L) stop("deepgof1(): needs at least one covariate")
  } else {
    v <- colnames(mm)[colnames(mm) != "(Intercept)"]
    if (length(v) < 1L) stop("deepgof1(): needs at least one covariate")
    X <- mm[, v, drop = FALSE]
  }
  pairs_path <- reading == "allpairs" || ncol(X) == 1L

  ## ties among covariate values: one random permutation per call, drawn only when some
  ## covariate is tied, shared by the observed map and every bootstrap map
  tied <- any(apply(X, 2, anyDuplicated) > 0L)
  tb   <- if (tied && ties == "random") sample.int(nrow(mm)) else NULL

  ph <- as.numeric(stats::fitted(fit)); y <- as.numeric(fit$y)
  if (pairs_path) {
    cells <- .deepgof1_cells(X, K, tb)
    stat_all <- function(r, p) vapply(cells$idx, function(i) score(.deepgof1_cellmap(i, r, p, K)), 0)
    so <- stat_all(y - ph, ph); Sobs <- max(so); top <- which.max(so)
    axes <- if (ncol(X) == 1L) colnames(X) else cells$pairs[top, ]
    map  <- .deepgof1_cellmap(cells$idx[[top]], y - ph, ph, K)
  } else {
    obs <- .deepgof1_map(fit, K, tb); Sobs <- score(obs$map); axes <- obs$axes; map <- obs$map
  }

  ## the analyst's own parametric bootstrap: the source of the level. Each refit uses the fitted
  ## model's design matrix (archive versions up to 2.0.0 refitted the formula on the model frame,
  ## which fails for every transformed term, so every replicate scored +Inf and p was 1).
  itc <- attr(stats::terms(fit), "intercept") > 0L
  Sb  <- numeric(B)
  for (b in seq_len(B)) {
    ys <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(mm, ys, weights = fit$prior.weights,
                                                   offset = fit$offset, family = stats::binomial(),
                                                   control = fit$control, intercept = itc)),
                   error = function(e) NULL)
    ## a FAILED refit scores +Inf: it then counts AGAINST rejection (conservative direction --
    ## -Inf would deflate the p-value). Measured activation on the paper's designs: 0.03%.
    Sb[b] <- if (is.null(fb)) Inf else if (pairs_path) {
      pb <- as.numeric(fb$fitted.values); max(stat_all(ys - pb, pb))
    } else score(.deepgof1_map_core(mm, ys, as.numeric(fb$fitted.values), fb$coefficients, K, tb)$map)
  }
  p <- (1 + sum(Sb >= Sobs)) / (B + 1)
  structure(list(statistic = Sobs, p.value = p, B = B, K = K, reading = reading, axes = axes,
                 map = matrix(map, K, K, byrow = TRUE), boot = Sb,
                 method = "DeepGOF-1: pretrained residual-map goodness-of-fit test",
                 data.name = deparse(substitute(fit))), class = "deepgof1")
}

## the current default version
deepgof <- deepgof1

print.deepgof1 <- function(x, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("S = %.4f, B = %d, p-value = %.4f\n", x$statistic, x$B, x$p.value))
  if (length(x$axes) == 1L) cat(sprintf("map: %d quantile cells over ranks of %s\n", x$K^2, x$axes))
  else cat(sprintf("map: %d x %d over ranks of %s and %s\n", x$K, x$K, x$axes[1], x$axes[2]))
  cat("alternative hypothesis: the logistic model is misspecified\n\n")
  invisible(x)
}
