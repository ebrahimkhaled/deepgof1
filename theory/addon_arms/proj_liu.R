# =============================================================================
# proj_liu.R -- the projection test of Liu et al. (2024, Stat Comput 34:175, Sec. 2.1) with the
# weight A_ijl as they define it (after Escanciano 2006, Econometric Theory 22:1030-1051).
#
# Definition.  T = n^-2 sum_i sum_j sum_l e_i e_j A_ijl,  e_i = y_i - p-hat_i,
#   A_ijl = integral over the unit sphere of I(X_i'w <= X_l'w) I(X_j'w <= X_l'w) dw.
# Write u = X_i - X_l, v = X_j - X_l. The indicators are I(u'w <= 0) and I(v'w <= 0), so A_ijl is the
# measure of the intersection of two closed half-spheres. With dw the uniform PROBABILITY measure:
#   u, v both nonzero:  A_ijl = (pi - angle(u, v)) / (2 pi),  angle = arccos(u'v / (|u| |v|)).
#     (Project w onto span(u, v); its direction is uniform on the circle. Each half-sphere is a
#     half-circle; two half-circles whose normals are at angle theta overlap on an arc of pi - theta.
#     Checks: u = v gives 1/2; u = -v gives 0.)
#   u = 0, v nonzero (X_i = X_l):  I(0 <= 0) = 1 always, so A_ijl = P(v'w <= 0) = 1/2.
#   u = v = 0:  A_ijl = 1.
# This is Escanciano's (2006, Appendix) A^(0)_ijr = pi - arccos(...) up to the constant 2 pi, with
# his stated special cases (X_i = X_j != X_r: pi; X_i = X_j = X_r: 2 pi; X_i = X_r != X_j: pi).
# The printed Liu et al. formula reads "C_p (1/pi) arccos(...)", without "pi -"; the integral
# above fixes the orientation. There is NO length factor |u| |v|: A_ijl depends on directions only.
#
# The constant (here 1/(2 pi) per term, i.e. the sphere taken with total mass 1) does not change
# a bootstrap p-value. The zero-vector convention and the orientation do.
#
# Two versions of the weight: on the RAW covariates (Liu et al.'s statistic exactly; the primary
# arm) and on covariates standardised with scale() (the deployed arm's choice; secondary). The angle
# weight is not invariant to rescaling a column, so the two are different tests.
# Everything else is kept exactly as in the deployed arm (code/liu_augment.R, proj_A / proj_T /
# the bootstrap block of one(); theory/niche/proj_verbatim.R): X = model matrix without the
# intercept; T = e'Ae / n^2; model-based bootstrap, B = 199, y* ~
# Bernoulli(fitted), refit with the same fit_pol (Firth fallback), a failed refit scores -Inf;
# p = (1 + #{T* >= T}) / (B + 1).
# =============================================================================

## ---- the deployed weight, VERBATIM from code/liu_augment.R (kept only for the pairing check) ----
proj_A_nn <- function(X) {
  X <- scale(X); n <- nrow(X)
  A <- matrix(0, n, n)
  for (r in seq_len(n)) {
    D  <- sweep(X, 2, X[r, ], "-")                 # n x p
    nm <- sqrt(rowSums(D*D))
    G  <- tcrossprod(D)                            # n x n inner products
    NN <- outer(nm, nm)
    cs <- ifelse(NN > 1e-12, G/NN, 0)
    A  <- A + NN * (pi - acos(pmin(pmax(cs, -1), 1)))
  }
  A
}

## ---- the corrected weight: A_ij = sum_l A_ijl, A_ijl = (pi - angle) / (2 pi) ----
## std = TRUE standardises the columns first, as the deployed arm does. tol: a difference vector
## shorter than tol (relative to the data scale) counts as the zero vector.
proj_A_liu <- function(X, std = TRUE, tol = 1e-10) {
  X <- as.matrix(X); if (std) X <- scale(X)
  n <- nrow(X); A <- matrix(0, n, n)
  eps <- tol * max(1, sqrt(max(rowSums(X * X))))
  for (l in seq_len(n)) {
    D  <- sweep(X, 2, X[l, ], "-")
    nm <- sqrt(rowSums(D * D)); z <- nm <= eps
    s  <- 1 / ifelse(z, 1, nm)
    cs <- tcrossprod(D * s)                        # cosines (rows of D scaled to unit length)
    K  <- (pi - acos(pmin(pmax(cs, -1), 1))) / (2 * pi)
    diag(K) <- 0.5                                 # i = j: angle 0 exactly (avoids acos rounding at 1)
    if (any(z)) { K[z, ] <- 0.5; K[, z] <- 0.5; K[z, z] <- 1 }
    A  <- A + K
  }
  A
}

## ---- both weights from ONE pass (shared differences and inner products) ----
## old is bit-identical to proj_A_nn (same operations in the same order on the same numbers).
proj_A_both <- function(X, tol = 1e-10) {
  X <- scale(X); n <- nrow(X)
  Aold <- matrix(0, n, n); Anew <- matrix(0, n, n)
  eps <- tol * max(1, sqrt(max(rowSums(X * X))))
  for (r in seq_len(n)) {
    D  <- sweep(X, 2, X[r, ], "-")
    nm <- sqrt(rowSums(D*D))
    G  <- tcrossprod(D)
    NN <- outer(nm, nm)
    cs <- ifelse(NN > 1e-12, G/NN, 0)
    Aold <- Aold + NN * (pi - acos(pmin(pmax(cs, -1), 1)))
    z  <- nm <= eps
    s  <- 1 / ifelse(z, 1, nm)
    cn <- tcrossprod(D * s)
    K  <- (pi - acos(pmin(pmax(cn, -1), 1))) / (2 * pi)
    diag(K) <- 0.5
    if (any(z)) { K[z, ] <- 0.5; K[, z] <- 0.5; K[z, z] <- 1 }
    Anew <- Anew + K
  }
  list(old = Aold, new = Anew)
}

## ---- compiled copy of proj_A_liu() (proj_liu_c.c; about 5x faster at n = 500) ----
## Same definition and zero rule; agrees with proj_A_liu() to rounding (check_proj_liu_c.R).
## Falls back to proj_A_liu() when the DLL is not built.
## std = TRUE: columns standardised first (the deployed arm's choice); std = FALSE: raw X, which is
## Liu et al.'s statistic exactly.
proj_A_liu_fast <- function(X, std = TRUE, tol = 1e-10) {
  dll <- file.path("/path/to/project/paper_deepgof/theory/addon_arms",
                   paste0("proj_liu_c", .Platform$dynlib.ext))
  if (!is.loaded("proj_A_liu_c")) { if (file.exists(dll)) dyn.load(dll) else return(proj_A_liu(X, std = std, tol = tol)) }
  X <- as.matrix(X); if (std) X <- scale(X)
  n <- nrow(X); p <- ncol(X)
  eps <- tol * max(1, sqrt(max(rowSums(X * X))))
  matrix(.C("proj_A_liu_c", as.double(X), as.integer(n), as.integer(p), as.double(eps),
            A = double(n * n))$A, n, n)
}

proj_T <- function(e, A) as.numeric(crossprod(e, A %*% e)) / (length(e)^2)

## ---- fit_pol of the deployed arm, VERBATIM (code/liu_augment.R; proj_verbatim.R proj_fit_pol) ----
proj_fit_pol <- function(frm, d) {
  g <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial())), error=function(e) NULL)
  if (is.null(g)) return(NULL); b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1])>8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(frm,data=d,family=binomial(),method="brglmFit",type="AS_mean")),
                  error=function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(f)
  }
  g
}

## ---- the bootstrap block of the deployed arm, run once for several weights ----
## Same draws, same refits, one quadratic form per weight. With AL = list(old = proj_A_nn(X))
## this returns exactly the deployed proj p-value (given the same RNG state on entry).
## Returns p-values, observed statistics and the number of failed refits.
proj_p_multi <- function(d, frm, fit, B, AL) {
  NOBS <- nrow(d)
  ph <- as.numeric(stats::fitted(fit))
  Tobs <- vapply(AL, function(A) proj_T(d$y - ph, A), 0)
  Tb <- matrix(0, B, length(AL)); d2 <- d; nfail <- 0L
  for (b in seq_len(B)) {
    d2$y <- rbinom(NOBS, 1, ph)
    fb <- proj_fit_pol(frm, d2)
    if (is.null(fb)) { Tb[b, ] <- -Inf; nfail <- nfail + 1L; next }
    eb <- d2$y - as.numeric(stats::fitted(fb))
    Tb[b, ] <- vapply(AL, function(A) proj_T(eb, A), 0)
  }
  p <- (1 + colSums(sweep(Tb, 2, Tobs, ">="))) / (B + 1)
  list(p = stats::setNames(p, names(AL)), T = stats::setNames(Tobs, names(AL)), nfail = nfail)
}
