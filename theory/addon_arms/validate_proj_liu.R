# =============================================================================
# validate_proj_liu.R -- checks proj_A_liu() against its integral definition by brute force.
#
#   T = integral over w of integral over u of M_n(u, w)^2 F_nw(du) dw,
#   M_n(u, w) = n^-1/2 sum_i e_i I(X_i'w <= u),  F_nw = empirical cdf of X_i'w,
# with w uniform on the unit sphere (probability measure). For one w the inner integral is
#   n^-1 sum_l M_n(X_l'w, w)^2 = n^-2 sum_l (sum_i e_i I(X_i'w <= X_l'w))^2,
# so averaging over random directions estimates T; the closed form is e'Ae / n^2 with A from
# proj_A_liu(). Also checked entrywise: A_ij = sum_l P_w(X_i'w <= X_l'w, X_j'w <= X_l'w).
# Designs: p = 3 continuous; p = 2 with a binary covariate (as Liu Setting 3); p = 1; and a
# design with duplicated rows, which exercises the zero-vector convention for i != l.
# Also: proj_A_both() reproduces proj_A_nn() and proj_A_liu() bit for bit.
# Output: validate_proj_liu.log
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
source(file.path(HERE, "proj_liu.R"))
M <- 200000L

mc_check <- function(label, X, e, seed, std = TRUE) {
  Xs <- if (std) scale(X) else as.matrix(X); n <- nrow(Xs); p <- ncol(Xs)
  A  <- proj_A_liu(X, std = std)                   # same scaling as the arm checked
  Tcf <- proj_T(e, A)
  set.seed(seed)
  Tw <- numeric(M); Amc <- matrix(0, n, n); Amc2 <- matrix(0, n, n)
  CHK <- 20000L
  for (s0 in seq(1L, M, by = CHK)) {
    m <- min(CHK, M - s0 + 1L)
    W <- matrix(stats::rnorm(p * m), p, m); W <- sweep(W, 2, sqrt(colSums(W * W)), "/")
    S <- Xs %*% W                                  # n x m projections
    for (k in seq_len(m)) {
      sk <- S[, k]
      I  <- outer(sk, sk, "<=")                    # I[i, l] = 1(X_i'w <= X_l'w), ties included
      cl <- colSums(I * e)                         # sum_i e_i I(i, l)
      Tw[s0 + k - 1L] <- sum(cl^2) / n^2
      Ak  <- tcrossprod(I)
      Amc <- Amc + Ak; Amc2 <- Amc2 + Ak * Ak
    }
  }
  Amc <- Amc / M; seA <- sqrt(pmax(Amc2 / M - Amc^2, 1e-300) / M)
  Tmc <- mean(Tw); se <- stats::sd(Tw) / sqrt(M)
  zA  <- ifelse(seA > 1e-12, (Amc - A) / seA, ifelse(abs(Amc - A) < 1e-6, 0, Inf))   # zero-variance entries (p = 1): acos rounding near +-1 is ~1e-8
  dA  <- max(abs(Amc - A))
  ## the printed orientation (arccos without "pi -"), same zero-vector rule, for contrast
  Aalt <- matrix(0, n, n)
  for (l in seq_len(n)) {
    D <- sweep(Xs, 2, Xs[l, ], "-"); nm <- sqrt(rowSums(D * D)); z <- nm <= 1e-10
    cs <- tcrossprod(D / ifelse(z, 1, nm)); K <- acos(pmin(pmax(cs, -1), 1)) / (2 * pi)
    if (any(z)) { K[z, ] <- 0.5; K[, z] <- 0.5; K[z, z] <- 1 }
    Aalt <- Aalt + K
  }
  cat(sprintf("%-34s n=%d p=%d | T closed %.6f | T MC %.6f (SE %.6f) | z = %+.2f\n",
              label, n, p, Tcf, Tmc, se, (Tmc - Tcf) / se))
  cat(sprintf("%-34s entrywise: max|A_MC - A| = %.4f, max |z| over the %d distinct entries = %.2f\n", "",
              dA, n * (n + 1) / 2, max(abs(zA[upper.tri(zA, diag = TRUE)]))))
  cat(sprintf("%-34s printed orientation (arccos, no 'pi -'): T = %.6f, z = %+.1f; max|A_MC - A_alt| = %.3f\n", "",
              proj_T(e, Aalt), (Tmc - proj_T(e, Aalt)) / se, max(abs(Amc - Aalt))))
  if (std) {                                       # deployed weight and one-pass version: standardised only
  Ao <- proj_A_nn(X)
  cat(sprintf("%-34s deployed NN weight: T_nn/T_def = %.3f, cor(A_nn, A_def) off-diagonal = %.3f\n", "",
              proj_T(e, Ao) / Tcf, stats::cor(Ao[upper.tri(Ao)], A[upper.tri(A)])))
  ## the one-pass version used by the add-on runner
  AB <- proj_A_both(X)
  cat(sprintf("%-34s proj_A_both: old identical to proj_A_nn: %s | new identical to proj_A_liu: %s\n", "",
              identical(AB$old, Ao), identical(AB$new, A)))
  }
  invisible(c(Tcf = Tcf, Tmc = Tmc, se = se, dA = dA))
}

set.seed(20260924)
n <- 25
X1 <- cbind(runif(n, -3, 3), rnorm(n, 0, 1.5), rt(n, 4))
e1 <- rbinom(n, 1, 0.3) - 0.3
X2 <- cbind(runif(n, -3, 3), rbinom(n, 1, 0.5))
e2 <- rbinom(n, 1, 0.4) - 0.4
X3 <- cbind(runif(n, -3, 3))
e3 <- rbinom(n, 1, 0.5) - 0.5
X4 <- rbind(X1[1:20, ], X1[c(3, 3, 7, 12, 12), ])   # duplicated rows: X_i = X_l with i != l
e4 <- rbinom(n, 1, 0.3) - 0.3

sink(file.path(HERE, "validate_proj_liu.log"), split = TRUE)
cat("validate_proj_liu |", format(Sys.time()), "| directions M =", M, "\n")
mc_check("p=3 continuous (Setting 2 design)", X1, e1, 1)
mc_check("p=2 with binary (Setting 3 design)", X2, e2, 2)
mc_check("p=1 (Setting 1 design)", X3, e3, 3)
mc_check("p=3 with duplicated rows", X4, e4, 4)
## raw covariates (Liu et al.'s statistic exactly, the primary arm)
mc_check("RAW p=3 continuous", X1, e1, 5, std = FALSE)
mc_check("RAW p=2 with binary", X2, e2, 6, std = FALSE)
mc_check("RAW p=3 with duplicated rows", X4, e4, 7, std = FALSE)
cat(sprintf("raw vs standardised, p=2 with binary: T_raw / T_std = %.3f\n",
            proj_T(e2, proj_A_liu(X2, std = FALSE)) / proj_T(e2, proj_A_liu(X2))))
## p = 1 exactly: the sphere is {-1, +1}, so A is the average over the two directions
Xs <- scale(X3); A3 <- proj_A_liu(X3)
ex <- (tcrossprod(outer(Xs[, 1], Xs[, 1], "<=")) + tcrossprod(outer(-Xs[, 1], -Xs[, 1], "<="))) / 2
cat(sprintf("p=1 EXACT (two directions) vs closed form: max |diff| = %.2e\n", max(abs(ex - A3))))
sink()
