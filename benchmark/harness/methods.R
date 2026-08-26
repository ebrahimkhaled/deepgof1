# Baselines (run.all.gof) + Method A (learned-basis directed) + Method C (meta-learner).
# Assumes the ebrahim.gof R sources are already source()d (run.all.gof available).

BASELINE_TESTS <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                    "DEF.poly2", "DEF.poly3", "DEF.stukel",
                    "Ensemble.Vote(3DEF)", "Ensemble.Univ(3DEF+EF)")

# p-values for the baseline battery (one run.all.gof call).
baseline_pvals <- function(fit) {
  out <- setNames(rep(NA_real_, length(BASELINE_TESTS)), BASELINE_TESTS)
  r <- tryCatch(suppressWarnings(run.all.gof(fit)), error = function(e) NULL)
  if (!is.null(r)) out[] <- r$p_value[match(BASELINE_TESTS, r$Test)]
  out
}

# Features for the meta-learner: partition + DEF-family + EF (NOT the ensembles),
# transformed to one-sided evidence z = qnorm(1 - p) so larger = stronger misfit.
META_FEATS <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                "DEF.poly2", "DEF.poly3", "DEF.stukel")
feats_from_pvals <- function(pv) {
  z <- qnorm(pmin(pmax(pv[META_FEATS], 1e-6), 1 - 1e-6), lower.tail = FALSE)
  z[is.na(z)] <- 0
  setNames(as.numeric(z), make.names(META_FEATS))
}

# ---- Method A: cross-fitted learned-direction directed test ----
# Learns the misfit direction in covariate space (squares/cubes + pairwise
# products), fitted on one fold and applied to the other (cross-fit => valid).
methodA <- function(fit, K = 2) {
  mm <- model.matrix(fit); X <- mm[, -1, drop = FALSE]
  y <- fit$y; n <- length(y)
  make_dict <- function(Xm) {
    cols <- list()
    for (j in seq_len(ncol(Xm))) { cols[[length(cols) + 1L]] <- Xm[, j]^2
                                   cols[[length(cols) + 1L]] <- Xm[, j]^3 }
    if (ncol(Xm) >= 2)
      for (a in 1:(ncol(Xm) - 1)) for (b in (a + 1):ncol(Xm))
        cols[[length(cols) + 1L]] <- Xm[, a] * Xm[, b]
    do.call(cbind, cols)
  }
  folds <- sample(rep(1:K, length.out = n))
  num <- 0; den <- 0
  for (k in 1:K) {
    tr <- which(folds != k); te <- which(folds == k)
    g <- tryCatch(suppressWarnings(glm.fit(mm[tr, , drop = FALSE], y[tr],
                                           family = binomial())),
                  error = function(e) NULL)
    if (is.null(g)) next
    eta_te  <- as.numeric(mm[te, , drop = FALSE] %*% g$coefficients)
    phat_te <- plogis(eta_te); r_te <- y[te] - phat_te
    r_tr    <- y[tr] - g$fitted.values
    Ztr <- make_dict(X[tr, , drop = FALSE]); Zte <- make_dict(X[te, , drop = FALSE])
    ctr <- colMeans(Ztr); sds <- apply(Ztr, 2, sd); sds[sds == 0] <- 1
    Ztr <- scale(Ztr, ctr, sds); Zte <- scale(Zte, ctr, sds)
    df <- tryCatch(lm.fit(cbind(1, Ztr), r_tr), error = function(e) NULL)
    if (is.null(df)) next
    co <- df$coefficients; co[is.na(co)] <- 0
    s_te <- as.numeric(cbind(1, Zte) %*% co)
    num <- num + sum(s_te * r_te)
    den <- den + sum(s_te^2 * phat_te * (1 - phat_te))
  }
  z <- if (den > 0) num / sqrt(den) else 0
  list(stat = z, p = pnorm(z, lower.tail = FALSE))   # one-sided: learned dir. aligns with misfit
}

# ---- Method A, FIXED: valid covariate-space directed (DEF-style) score test ----
# Projects the standardized residual onto a covariate-space basis Z and calibrates
# with the Farrington Omega-projection (accounts for beta-hat) + weighted-chi-square
# (Satterthwaite/Imhof, reusing .def_pvalue). Valid by construction; reuses the exact
# machinery that makes def.gof valid, but with a basis in covariate space instead of
# fitted-probability space.  basis = "poly" (squares/cubes/interactions) or "spline".
methodA_def <- function(fit, basis = c("poly", "spline", "combined"),
                        method = "satterthwaite") {
  basis <- match.arg(basis)
  y  <- as.numeric(fit$y)
  ph <- pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6)
  X  <- stats::model.matrix(fit)
  V  <- ph * (1 - ph)
  r  <- (y - ph) / sqrt(V)
  Xc <- X[, -1, drop = FALSE]
  # covariate-space polynomial basis (squares, cubes, pairwise products)
  poly_cov <- function() {
    cols <- list()
    for (j in seq_len(ncol(Xc))) { cols[[length(cols)+1L]] <- Xc[,j]^2
                                   cols[[length(cols)+1L]] <- Xc[,j]^3 }
    if (ncol(Xc) >= 2) for (a in 1:(ncol(Xc)-1)) for (b in (a+1):ncol(Xc))
      cols[[length(cols)+1L]] <- Xc[,a] * Xc[,b]
    do.call(cbind, cols)
  }
  # fitted-probability-space bends (DEF-style calibration-curve shape, individual level)
  phat_shape <- function() { e <- stats::qlogis(ph); cbind(e^2, e^3) }
  if (basis == "poly") {
    Z <- poly_cov()
  } else if (basis == "spline") {
    if (!requireNamespace("splines", quietly = TRUE)) return(list(stat = NA, df = NA, p = NA))
    parts <- lapply(seq_len(ncol(Xc)), function(j) {
      xj <- Xc[, j]; if (length(unique(xj)) > 5) splines::ns(xj, df = 4) else NULL })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    inter <- if (ncol(Xc) >= 2) Xc[, 1] * Xc[, 2] else NULL
    Z <- do.call(cbind, c(parts, list(inter)))
  } else {                       # combined: covariate-space + fitted-probability-space
    Z <- cbind(poly_cov(), phat_shape())
  }
  Z <- scale(Z)
  Z <- Z[, apply(Z, 2, function(z) all(is.finite(z)) && stats::sd(z) > 1e-8), drop = FALSE]
  if (is.null(dim(Z)) || ncol(Z) < 1) return(list(stat = NA, df = NA, p = NA))
  qz <- qr(Z)                                   # drop collinear/rank-deficient columns
  if (qz$rank < ncol(Z)) Z <- Z[, qz$pivot[seq_len(qz$rank)], drop = FALSE]
  if (ncol(Z) < 1) return(list(stat = NA, df = NA, p = NA))
  ZtZ  <- crossprod(Z)
  Zr   <- crossprod(Z, r)
  S    <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
  ZtU  <- crossprod(Z, sqrt(V) * X)            # Z' W^{1/2} X
  XtWX <- crossprod(X, V * X)
  ZtOZ <- ZtZ - ZtU %*% solve(XtWX) %*% t(ZtU) # Z' Omega Z
  lam  <- Re(eigen(solve(ZtZ) %*% ZtOZ, only.values = TRUE)$values)
  lam  <- lam[lam > 1e-9]
  if (!length(lam)) return(list(stat = S, df = NA, p = NA))
  list(stat = S, df = sum(lam)^2 / sum(lam^2), p = .def_pvalue(S, lam, method))
}

# ---- Method C: meta-learner over the battery ----
train_meta <- function(FM, labels) {
  d <- as.data.frame(FM); d$.lab <- labels
  suppressWarnings(glm(.lab ~ ., data = d, family = binomial()))
}
score_meta <- function(model, FM) {
  d <- as.data.frame(FM, stringsAsFactors = FALSE)
  as.numeric(predict(model, newdata = d, type = "link"))
}
