## deepgof.R -- the practitioner-facing test. This is the whole user interface:
##
##     fit <- glm(y ~ x1 + x2, data = d, family = binomial())
##     deepgof(fit)
##
## No training, no Python, no GPU, no tuning. The network was trained ONCE by us on simulated data
## and ships as 18,273 frozen numbers inside the package; the analyst only ever runs the forward
## pass. This file is what would become R/deepgof.R in `ebrahim.gof`, with `weights.csv` and
## `scaler.csv` becoming a single `sysdata.rda` object.
##
## What happens per call:
##   1. build the residual map of the fitted model (K x K grid over covariate RANKS)
##   2. score it with the frozen network                       -> S_obs
##   3. simulate B datasets from the fitted model, refit each, score each -> S_1 .. S_B
##   4. p = (1 + #{S_b >= S_obs}) / (B + 1)
##
## Step 3 is what makes the level exact (Prop. 1): the p-value is a rank within the analyst's own
## bootstrap, so it does not depend on the network being right about anything.

source(file.path(dirname(sys.frame(1)$ofile %||% "."), "deepgof_forward.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- the residual map, identical to the training-time builder ----------------------------------
.deepgof_map <- function(fit, K = 6L) {
  d  <- stats::model.frame(fit)
  mm <- stats::model.matrix(fit)
  vars <- colnames(mm)[-1]
  if (length(vars) < 2L) stop("deepgof(): needs at least two covariates")
  ## the two grid axes: largest |beta_j| * sd(x_j), returned in model order
  b  <- stats::coef(fit)[vars]
  sc <- abs(b) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
  ax <- vars[sort(order(sc, decreasing = TRUE)[1:2])]
  ph <- as.numeric(stats::fitted(fit))
  r  <- stats::model.response(d) - ph
  n  <- length(r)
  bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
  idx <- factor((bn(mm[, ax[1]]) - 1L) * K + bn(mm[, ax[2]]), levels = 1:(K * K))
  s  <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  list(map = as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8)), axes = ax)
}

## ---- the test ----------------------------------------------------------------------------------
deepgof <- function(fit, B = 199L, K = 6L, weights_dir = NULL) {
  stopifnot(inherits(fit, "glm"), fit$family$family == "binomial")
  wd <- weights_dir %||% getOption("deepgof.weights",
        "/path/to/project/paper_deepgof/code/rdeploy")
  W  <- load_weights(file.path(wd, "weights.csv"))
  ms <- as.matrix(utils::read.csv(file.path(wd, "scaler.csv"), header = FALSE))
  mu <- ms[1, ]; sdv <- ms[2, ]
  score <- function(m) deepgof_score((m - mu) / sdv, W, K)

  obs  <- .deepgof_map(fit, K)
  Sobs <- score(obs$map)

  ## the analyst's own parametric bootstrap: the source of the exact level
  ph  <- as.numeric(stats::fitted(fit))
  dat <- stats::model.frame(fit)
  yn  <- names(dat)[1]
  Sb  <- numeric(B)
  for (b in seq_len(B)) {
    dat[[yn]] <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(stats::formula(fit), data = dat,
                                               family = stats::binomial())),
                   error = function(e) NULL)
    ## a FAILED refit scores +Inf: it then counts AGAINST rejection (conservative direction --
    ## -Inf would deflate the p-value). Measured activation on the paper's designs: 0.03%.
    Sb[b] <- if (is.null(fb)) Inf else score(.deepgof_map(fb, K)$map)
  }
  p <- (1 + sum(Sb >= Sobs)) / (B + 1)
  structure(list(statistic = Sobs, p.value = p, B = B, K = K, axes = obs$axes,
                 boot = Sb, method = "Deep-GOF: pretrained residual-map test",
                 data.name = deparse(substitute(fit))), class = "deepgof")
}

print.deepgof <- function(x, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("S = %.4f, B = %d, p-value = %.4f\n", x$statistic, x$B, x$p.value))
  cat(sprintf("grid: %d x %d over ranks of %s and %s\n", x$K, x$K, x$axes[1], x$axes[2]))
  cat("alternative hypothesis: the logistic model is misspecified\n\n")
  invisible(x)
}
