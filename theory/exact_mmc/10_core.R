## 10_core.R -- the CMMC machinery. Sourced by 11_bench.R and 12_run.R; defines nothing at load
## time except functions and constants.
##
## Nothing here re-implements the deployed test. The residual map, the scaler and the network are
## taken from theory/calibration_gap/01_batched_net.R (verified bit-for-bit against the shipped
## deepgof_score() and against PyTorch's test_output.csv) and from code/rdeploy. The refit policy
## is the deployed one: glm, Firth (brglmFit AS_mean) on separation, exactly as code/s3_grid.R.
##
## What IS new here:
##   * pB(theta)  -- the common-random-number Monte Carlo p-value at an arbitrary theta,
##   * the Bernstein score confidence set C_gamma and its bounding box,
##   * the maximisation of pB over C_gamma (grid + Nelder-Mead).

suppressMessages(library(brglm2))

PAP <- "/path/to/project/paper_deepgof"
RD  <- file.path(PAP, "code", "rdeploy")
CG  <- file.path(PAP, "theory", "calibration_gap")
EX  <- file.path(PAP, "theory", "exact_mmc")

source(file.path(RD, "deepgof_forward.R"))
source(file.path(CG, "01_batched_net.R"))
source(file.path(PAP, "code", "universal_design.R"))

W    <- load_weights(file.path(RD, "weights.csv"))
.ms  <- as.matrix(utils::read.csv(file.path(RD, "scaler.csv"), header = FALSE))
MU   <- .ms[1, ]; SDV <- .ms[2, ]
BFAM <- binomial()

## ---- protocol constants -------------------------------------------------------------------------
BREP    <- 199L          # Monte Carlo replicates per theta
GAMMA   <- 0.01          # confidence-set miscoverage
ALPHA   <- 0.05          # nominal level
KBIN    <- 6L            # deployed grid resolution
NGRID   <- 11L           # points per axis of the maximisation grid  (11^3 = 1331)
NMSTART <- 3L            # Nelder-Mead restarts, from the best NMSTART grid points
NMMAXIT <- 60L
SEED_DATA <- 101000L     # s3_grid.R's SEED0
SEED_U    <- 770000L
SEED_BOX  <- 880000L
NOBS      <- 100L
D1        <- list(p = 2L, law = "indep", rho = 0.5, cat = 0L, rate = 0.50)
AMPMULT   <- 2 * sqrt(100 / NOBS)

## ---- deployed refit policy (rule verbatim from code/s3_grid.R fit_pol) ---------------------------
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
  list(ph = as.numeric(g$fitted.values), firth = firth, coef = as.numeric(g$coefficients))
}

## fixed model matrix + fixed 36-cell incidence. p = 2, so the shipped axis rule always returns
## c(x1, x2) and the rank bins depend only on X, which is fixed within a dataset.
prep_D1 <- function(d, K = KBIN) {
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
score_one <- function(map) deepgof_score_batch(matrix((map - MU) / SDV, 1L, 36L), W)

## ---- generators ----------------------------------------------------------------------------------
## D1 null, RNG-identical to gen_design_null(n, 2, "indep", 0.5, 0, 0.5) of code/universal_design.R,
## with theta0 captured. Verified against the original in 11_bench.R.
gen_D1_null_theta <- function(n, rate = 0.5) {
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
  d$y <- stats::rbinom(n, 1, plogis(b0 + lin))
  list(data = d, theta0 = c(b0, bet[-1] * cscl))
}

## one dataset of the confirmatory grid, seeded exactly as s3_grid.R's one(kind, i)
gen_dataset <- function(kind, i) {
  set.seed(SEED_DATA + kind * 300000L + i)
  if (kind == 0L) {
    g <- gen_D1_null_theta(NOBS, D1$rate)
    list(data = g$data, theta0 = g$theta0)
  } else {
    fs  <- if (kind == 1L) SEED_DATA * 17L + i else SEED_DATA * 17L + 700000L + i
    fam <- make_family_p(fs, D1$p); fam$amp <- fam$amp * AMPMULT
    g <- gen_design_alt(fam, NOBS, D1$p, D1$law, D1$rho, D1$cat, D1$rate)
    list(data = g$data, theta0 = rep(NA_real_, 3L))
  }
}

## ---- Bernstein score confidence set ---------------------------------------------------------------
## S_j(theta) = sum_i x_ij (y_i - pi_i(theta)) is a sum of n independent terms, mean zero at the
## true theta, with |term| <= |x_ij| and variance x_ij^2 pi_i (1 - pi_i). Bernstein's inequality:
##
##      P(|sum_i W_i| >= t) <= 2 exp( - t^2 / (2 (v + M t / 3)) ),   |W_i| <= M,  v = sum Var(W_i).
##
## Set the right side to gamma / (d+1) per coordinate and union-bound over the d+1 = 3 coordinates.
## With L = log(2 (d+1) / gamma), solving t^2 = 2 L (v + M t / 3) gives the threshold
##
##      t_j(theta) = L M_j / 3 + sqrt( (L M_j / 3)^2 + 2 L v_j(theta) ),   M_j = max_i |x_ij|.
##
## C_gamma = {theta : |S_j(theta)| <= t_j(theta) for all j} then has coverage >= 1 - gamma at every
## n, conditional on X. No asymptotics, no estimator condition.
bern_L <- function(dplus1 = 3L, gamma = GAMMA) log(2 * dplus1 / gamma)

## Theta: q x 3 matrix of candidate parameters. Returns q x 1 slack = max_j (|S_j| / t_j); a point
## is in C_gamma iff slack <= 1.
in_set_slack <- function(Theta, y, mm, L = bern_L(), Mj = NULL) {
  if (is.null(dim(Theta))) Theta <- matrix(Theta, 1L)
  if (is.null(Mj)) Mj <- apply(abs(mm), 2, max)
  Pi <- plogis(mm %*% t(Theta))                     # n x q
  R  <- y - Pi
  S  <- crossprod(mm, R)                            # 3 x q
  V  <- crossprod(mm^2, Pi * (1 - Pi))              # 3 x q
  a  <- L * Mj / 3
  Tt <- a + sqrt(a^2 + 2 * L * V)                   # 3 x q  (a recycles down columns)
  apply(abs(S) / Tt, 2, max)
}

## Bounding box of C_gamma: linear-approximation seed, Monte Carlo fill, then six penalised
## Nelder-Mead coordinate maximisations to catch extremes the Monte Carlo missed.
##
## TRUNCATION. C_gamma is NOT always bounded: 12_geom.R finds directions along which the score
## saturates inside the Bernstein thresholds on 3 of 40 D1 datasets at n = 100. Dufour's argument
## needs a compact parameter set, so the search is truncated at RCAP Wald standard errors from
## theta-hat and the truncation is flagged whenever it binds. That truncation is data-dependent,
## so on a flagged dataset the exactness claim is void, not merely approximate.
RCAP <- 200

set_box <- function(y, mm, thetahat, ph, seed, nmc = 20000L, rcap = RCAP) {
  L  <- bern_L(); Mj <- apply(abs(mm), 2, max)
  w  <- ph * (1 - ph)
  I  <- crossprod(mm, w * mm)
  Iinv <- tryCatch(solve(I), error = function(e) diag(3) * 1e3)
  Pi0 <- as.numeric(ph)
  v0  <- as.numeric(crossprod(mm^2, Pi0 * (1 - Pi0)))
  a0  <- L * Mj / 3
  t0  <- a0 + sqrt(a0^2 + 2 * L * v0)
  h   <- as.numeric(abs(Iinv) %*% t0)                # linear-approximation half-widths
  se  <- sqrt(pmax(diag(Iinv), 0))

  clo <- thetahat - rcap * se; chi <- thetahat + rcap * se     # truncation box

  set.seed(seed)
  F <- 3; ok <- integer(0); Th <- NULL
  for (it in 1:8) {
    lo <- pmax(thetahat - F * h, clo); hi <- pmin(thetahat + F * h, chi)
    Th <- cbind(runif(nmc, lo[1], hi[1]), runif(nmc, lo[2], hi[2]), runif(nmc, lo[3], hi[3]))
    sl <- in_set_slack(Th, y, mm, L, Mj)
    ok <- which(sl <= 1)
    if (!length(ok)) { F <- F * 1.6; next }
    blo <- apply(Th[ok, , drop = FALSE], 2, min)
    bhi <- apply(Th[ok, , drop = FALSE], 2, max)
    touch <- any((blo - lo) / (hi - lo) < 0.02) || any((hi - bhi) / (hi - lo) < 0.02)
    if (!touch || all(lo <= clo + 1e-12) && all(hi >= chi - 1e-12)) break
    F <- F * 1.6
  }
  if (!length(ok)) { blo <- thetahat - h; bhi <- thetahat + h }

  ## penalised coordinate maximisation over C_gamma intersected with the truncation box,
  ## six runs, started from the Monte Carlo extremes
  pen <- function(th, k, s) {
    ex <- max(c(0, (clo - th) / se, (th - chi) / se))
    if (ex > 0) return(1e6 + 1e3 * ex)
    sl <- in_set_slack(matrix(th, 1L), y, mm, L, Mj)
    if (sl > 1) return(1e3 + 1e3 * (sl - 1))
    -s * th[k]
  }
  for (k in 1:3) for (s in c(-1, 1)) {
    st <- if (length(ok)) Th[ok[which.max(s * Th[ok, k])], ] else thetahat
    o <- tryCatch(optim(st, pen, k = k, s = s, method = "Nelder-Mead",
                        control = list(maxit = 400, reltol = 1e-10,
                                       parscale = pmax(se, 1e-8))), error = function(e) NULL)
    if (!is.null(o) && o$value < 0) {
      if (s > 0) bhi[k] <- max(bhi[k], o$par[k]) else blo[k] <- min(blo[k], o$par[k])
    }
  }
  blo <- pmax(blo, clo); bhi <- pmin(bhi, chi)
  ## does the truncation bind? (i.e. is C_gamma numerically unbounded for this dataset?)
  unbdd <- as.integer(any((thetahat - blo) / se > 0.95 * rcap) ||
                      any((bhi - thetahat) / se > 0.95 * rcap))
  wdt <- bhi - blo
  blo <- pmax(pmin(blo, thetahat) - 0.02 * wdt, clo)
  bhi <- pmin(pmax(bhi, thetahat) + 0.02 * wdt, chi)
  list(lo = blo, hi = bhi, half = as.numeric((bhi - blo) / 2), se = as.numeric(se),
       lin_half = as.numeric(h), F = F, nin = length(ok), nmc = nmc, unbdd = unbdd,
       mc_in = if (length(ok)) Th[ok, , drop = FALSE] else NULL)
}

## ---- the CMMC p-value at one theta ------------------------------------------------------------------
## Common random numbers: U is n x B, drawn ONCE per dataset. y*_b(theta) = 1{U[,b] <= pi(theta)}.
## Each replicate is refit by the deployed policy, remapped and rescored; a replicate whose refit
## fails scores +Inf (the shipped convention: a failed refit cannot help the observed row).
##
## floor_count: exact early abandonment for the maximisation. If count + (B - k) <= floor_count the
## remaining replicates cannot make this theta the running maximum, so we stop. The returned p is
## then a LOWER bound, flagged by full = FALSE. This changes no supremum.
pB <- function(theta, Tobs, U, P, floor_count = -1L, chunk = 40L) {
  n <- P$n; B <- ncol(U)
  pi <- as.numeric(plogis(P$mm %*% theta))
  cnt <- 0L; k <- 0L; nfirth <- 0L; nfail <- 0L
  maps <- matrix(0, chunk, 36L)
  while (k < B) {
    cs <- min(chunk, B - k)
    bad <- logical(cs)
    for (j in seq_len(cs)) {
      yb <- as.numeric(U[, k + j] <= pi)
      f <- fast_fit(yb, P$mm)
      if (is.null(f) || is.null(f$ph)) { bad[j] <- TRUE; nfail <- nfail + 1L; next }
      nfirth <- nfirth + f$firth
      maps[j, ] <- fast_map(yb, f$ph, P)
    }
    sc <- deepgof_score_batch((maps[seq_len(cs), , drop = FALSE] -
                               rep(MU, each = cs)) / rep(SDV, each = cs), W, chunk = chunk)
    sc[bad] <- Inf
    cnt <- cnt + sum(sc >= Tobs)
    k <- k + cs
    if (cnt + (B - k) <= floor_count) break
  }
  list(p = (1 + cnt) / (B + 1), cnt = cnt, nrep = k, full = (k == B),
       firth = nfirth, fail = nfail)
}

## ---- the maximisation --------------------------------------------------------------------------
## sup over C_gamma of pB: 11^3 grid on the bounding box (points outside C_gamma dropped), then
## NMSTART Nelder-Mead runs from the best grid points. theta-hat is always evaluated first, so
## p_exact >= p_plugin by construction.
## SEARCH BUDGET. A theta evaluation costs 199 refits, and its cost is NOT constant over
## C_gamma: at parameter values far from theta-hat the simulated responses are near-separated,
## the deployed refit policy falls back to Firth (brglmFit, ~14 ms against glm.fit's ~2 ms), and
## one evaluation goes from ~0.6 s to ~3 s. Without a cap a single dataset can run past an hour.
## The search is therefore stopped after BUDGET seconds and the dataset flagged; on a flagged
## dataset p_exact is a strict lower bound on the supremum, so the exactness claim is void there.
BUDGET <- as.numeric(Sys.getenv("EX_BUDGET", "480"))

cmmc <- function(y, mm, P, U, Tobs, thetahat, ph, seed_box, budget = BUDGET) {
  L  <- bern_L(); Mj <- apply(abs(mm), 2, max)
  bx <- set_box(y, mm, thetahat, ph, seed_box)
  tstop <- Sys.time() + budget

  state <- new.env(parent = emptyenv())
  state$best_cnt <- -1L; state$best <- rep(NA_real_, 3L); state$n_eval <- 0L
  state$n_infeas <- 0L; state$nrep_tot <- 0L; state$trunc <- 0L
  evalp <- function(th) {
    if (state$trunc || Sys.time() > tstop) {
      state$trunc <- 1L
      return((1 + max(state$best_cnt, 0L)) / (BREP + 1))    # frozen: ends NM quickly
    }
    r <- pB(th, Tobs, U, P, floor_count = state$best_cnt)
    state$n_eval <- state$n_eval + 1L
    state$nrep_tot <- state$nrep_tot + r$nrep
    if (r$full && r$cnt > state$best_cnt) { state$best_cnt <- r$cnt; state$best <- th }
    r$p
  }

  ## theta-hat (the deployed test's parameter), same U. Always evaluated, budget or no budget:
  ## it is the plug-in comparator and it guarantees p_exact >= p_plugin.
  r0 <- pB(thetahat, Tobs, U, P, floor_count = -1L)
  state$n_eval <- 1L; state$nrep_tot <- r0$nrep
  state$best_cnt <- r0$cnt; state$best <- thetahat
  cnt_hat <- r0$cnt

  ## grid
  ax <- lapply(1:3, function(k) seq(bx$lo[k], bx$hi[k], length.out = NGRID))
  G  <- as.matrix(expand.grid(ax[[1]], ax[[2]], ax[[3]]))
  sl <- in_set_slack(G, y, mm, L, Mj)
  keep <- which(sl <= 1)
  n_grid_in <- length(keep)
  pg <- numeric(n_grid_in)
  if (n_grid_in) {
    ord <- keep[order(sl[keep])]                    # interior-first
    for (q in seq_along(ord)) { pg[q] <- evalp(G[ord[q], ]); if (state$trunc) break }
  } else ord <- integer(0)
  grid_best_cnt <- state$best_cnt
  grid_best <- state$best

  ## Nelder-Mead refinement from the best NMSTART grid points
  negp <- function(th) {
    s <- in_set_slack(matrix(th, 1L), y, mm, L, Mj)
    if (s > 1) { state$n_infeas <- state$n_infeas + 1L; return(1 + (s - 1)) }
    -evalp(th)
  }
  psc <- pmax((bx$hi - bx$lo) / (NGRID - 1) / 2, 1e-8)
  ## starts: the best NMSTART grid points. The 11^3 grid over the bounding box can contain FEW or
  ## NO interior points (12_geom.R: median 6.5 of 1331, min 0), so the list is topped up with
  ## theta-hat, which is always in C_gamma because the score equations hold there exactly.
  starts <- list()
  if (n_grid_in) {
    top <- ord[order(pg, decreasing = TRUE)[seq_len(min(NMSTART, n_grid_in))]]
    starts <- lapply(seq_along(top), function(q) as.numeric(G[top[q], ]))
  }
  if (!length(starts)) starts <- list(as.numeric(thetahat))
  for (s in starts) {
    tryCatch(optim(s, negp, method = "Nelder-Mead",
                   control = list(maxit = NMMAXIT, reltol = 1e-10, parscale = psc)),
             error = function(e) NULL)
  }
  nm_improved <- as.integer(state$best_cnt > grid_best_cnt)

  list(p_plugin = (1 + cnt_hat) / (BREP + 1),
       p_exact  = (1 + state$best_cnt) / (BREP + 1),
       p_grid   = (1 + grid_best_cnt) / (BREP + 1),
       argmax   = state$best, grid_argmax = grid_best,
       at_hat   = as.integer(isTRUE(all.equal(as.numeric(state$best), as.numeric(thetahat)))),
       n_eval   = state$n_eval, n_infeas = state$n_infeas, nrep_tot = state$nrep_tot,
       n_grid_in = n_grid_in, nm_improved = nm_improved, trunc = state$trunc, box = bx)
}

## ---- one dataset, end to end ------------------------------------------------------------------------
run_dataset <- function(kind, i) {
  t0 <- Sys.time()
  g <- tryCatch(gen_dataset(kind, i), error = function(e) NULL)
  if (is.null(g) || length(unique(g$data$y)) < 2L) return(NULL)
  d <- g$data
  P <- prep_D1(d)
  y <- as.numeric(d$y)
  fo <- fast_fit(y, P$mm)
  if (is.null(fo) || is.null(fo$ph)) return(NULL)
  Tobs <- score_one(fast_map(y, fo$ph, P))

  set.seed(SEED_U + kind * 300000L + i)
  U <- matrix(runif(P$n * BREP), P$n, BREP)

  r <- cmmc(y, P$mm, P, U, Tobs, as.numeric(fo$coef), fo$ph, SEED_BOX + kind * 300000L + i)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  c(kind = kind, rep = i, n = P$n, B = BREP,
    Tobs = Tobs, event_rate = mean(y), firth_obs = fo$firth,
    p_plugin = r$p_plugin, p_exact = r$p_exact, p_grid = r$p_grid,
    theta0_0 = g$theta0[1], theta0_1 = g$theta0[2], theta0_2 = g$theta0[3],
    that_0 = fo$coef[1], that_1 = fo$coef[2], that_2 = fo$coef[3],
    sup_0 = r$argmax[1], sup_1 = r$argmax[2], sup_2 = r$argmax[3],
    at_hat = r$at_hat, nm_improved = r$nm_improved, trunc = r$trunc,
    n_eval = r$n_eval, n_infeas = r$n_infeas, nrep_tot = r$nrep_tot,
    n_grid_in = r$n_grid_in,
    box_half_0 = r$box$half[1], box_half_1 = r$box$half[2], box_half_2 = r$box$half[3],
    se_0 = r$box$se[1], se_1 = r$box$se[2], se_2 = r$box$se[3],
    box_F = r$box$F, box_nin = r$box$nin, box_unbdd = r$box$unbdd,
    theta0_in = if (is.na(g$theta0[1])) NA_real_ else
      as.numeric(in_set_slack(matrix(g$theta0, 1L), y, P$mm) <= 1),
    secs = secs)
}
