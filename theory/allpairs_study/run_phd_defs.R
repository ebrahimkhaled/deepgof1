# =============================================================================
# run_phd_defs.R -- the generator, the fit policy and the readings, in ONE place.
#
# Both run_phd.R and verify_same_datasets.R source this. That is the point: if the verifier carried
# its own copy of the generator it would reproduce whatever mistake the runner made, and agree with
# it. Sourcing one definition means the only thing the verifier supplies independently is the
# deposited p-value it checks against.
#
# Everything from `gen_liu` through `pick_axes` is verbatim from code/liu_grid.R.
# =============================================================================
K <- 6L

gen_liu <- function(n, setting, C) {
  if (setting == 1L) {
    x1 <- runif(n, -3, 3); d <- data.frame(x1)
    eta <- -1.963 + 0.982 * x1 + C * x1^2
  } else if (setting == 2L) {
    x1 <- runif(n, -3, 3); x2 <- rnorm(n, 0, 1.5); x3 <- rt(n, 4)
    d <- data.frame(x1, x2, x3)
    eta <- -2 + 0.3 * x1 + 0.3 * x2 + 0.3 * x3 + C * x1^2
  } else if (setting == 3L) {
    x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); d <- data.frame(x1, x2)
    eta <- -1.972 + 0.135 * x1 + 1.117 * x2 + C * x1 * x2
  } else {
    x1 <- runif(n, -3, 3); x2 <- runif(n, -3, 3); d <- data.frame(x1, x2)
    eta <- 0.3 * x1 + 0.3 * x2 + C * x1 * x2                  # no intercept in the truth
  }
  d$y <- rbinom(n, 1, plogis(eta)); d
}

## the Firth policy, so a separated fit takes the same path the deposited column took
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

pick_axes <- function(fit, d) {
  vars <- setdiff(names(d), "y")
  mm <- stats::model.matrix(fit); b <- stats::coef(fit)
  sc <- vapply(vars, function(v) {
    cols <- which(startsWith(colnames(mm), v)); if (!length(cols)) return(0)
    nm <- colnames(mm)[cols]
    val <- abs(b[nm]) * apply(mm[, cols, drop = FALSE], 2, stats::sd)
    val <- val[is.finite(val)]; if (!length(val)) 0 else max(val)
  }, numeric(1))
  sc[!is.finite(sc)] <- 0
  sel <- vars[order(sc, decreasing = TRUE)][1:2]
  vars[sort(match(sel, vars))]
}

bn <- function(v, n) pmin(K, 1L + floor(K * (rank(as.numeric(v), ties.method = "first") - 1) / n))

map_uv <- function(u1, u2, r, ph) {
  n <- length(r)
  idx <- factor((bn(u1, n) - 1L) * K + bn(u2, n), levels = 1:(K * K))
  s <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8))
}

## every reading of one fitted model -> a named vector of scores
readings <- function(fit, d, xn, prs, score) {
  ph <- as.numeric(stats::fitted(fit)); r <- d$y - ph
  X  <- as.matrix(d[, xn, drop = FALSE])
  ax <- pick_axes(fit, d)
  out <- c(v1 = score(map_uv(X[, ax[1]], X[, ax[2]], r, ph)))
  ## pHd axes. A rank-deficient or degenerate M would make eigen() fail on a handful of separated
  ## fits; fall back to v1's axes there rather than dropping the dataset, so every reading stays
  ## defined on exactly the same rows.
  phd_axes <- function(rr) tryCatch({
    M <- crossprod(X * rr, X) / nrow(X)
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
    X %*% ev$vectors[, order(abs(ev$values), decreasing = TRUE)[1:2], drop = FALSE]
  }, error = function(e) X[, ax, drop = FALSE])
  Z  <- phd_axes(r)
  Zw <- phd_axes(r / sqrt(pmax(ph * (1 - ph), 1e-8)))
  out <- c(out, phd = score(map_uv(Z[, 1], Z[, 2], r, ph)),
                phdw = score(map_uv(Zw[, 1], Zw[, 2], r, ph)))
  if (length(prs) > 1L)
    out <- c(out, allpairs = max(vapply(prs, function(q) score(map_uv(X[, q[1]], X[, q[2]], r, ph)), 0)))
  out
}
