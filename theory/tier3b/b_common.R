# b_common.R -- shared definitions for Tier 3b (sourced, never run on its own).
# Tier 3b: non-polynomial departures. Reuses theory/tier3/tier3_common.R (RNG helpers, score_stat).
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
T3B <- file.path(PAP, "theory/tier3b")
TABB <- file.path(T3B, "tables")

DESIGNS <- c("D2", "D10")
SHAPES  <- c("threshold", "bump", "oscillation", "corner")
NS      <- c(200L, 500L)
N_PILOT <- 316L                       # geometric midpoint of 200 and 500 (the "middle n" of the plan)

## departure shapes, uncentred
h_raw <- function(shape, x1, x2) switch(shape,
  threshold   = as.numeric(x2 > 1),
  bump        = exp(-((x1 - 0.5)^2 + (x2 - 1)^2) / 0.5),
  oscillation = sin(2 * x1),
  corner      = as.numeric(x1 > 0.8 & x2 > 0.8),
  null        = 0 * x1)

## population means (Monte Carlo, computed once by 00_hmeans.R, stored in tables/H_MEANS.csv)
h_means <- function() {
  f <- file.path(TABB, "H_MEANS.csv"); if (!file.exists(f)) stop("run 00_hmeans.R first")
  z <- utils::read.csv(f); stats::setNames(z$mc_mean, z$shape)
}

## data-generating process: x ~ N(0, I_P); eta = linear part + C * (h - E h); y ~ Bernoulli(plogis(eta))
gen_b <- function(design, shape, n, C, hm) {
  P <- if (design == "D2") 2L else 10L
  X <- matrix(stats::rnorm(n * P), n, P); colnames(X) <- paste0("x", 1:P)
  lin <- if (design == "D2") -0.5 + X[, 1] + 0.6 * X[, 2] else -0.5 + X[, 1] + 0.8 * X[, 2] + 0.15 * rowSums(X[, 3:P])
  hc <- if (shape == "null") 0 else h_raw(shape, X[, 1], X[, 2]) - hm[[shape]]
  eta <- lin + C * hc
  data.frame(X, y = stats::rbinom(n, 1, stats::plogis(eta)))
}

## seeds: data set, bootstrap (DeepGOF-1 loop), GRP. Ranges never overlap (all < 2^31 - 1).
##   pilot ci = 10 + grid index; main ci = 1..3; nulls shape "null", ci = 0; ni: 1 = n200, 2 = n500, 3 = pilot n
shape_id <- function(shape) match(shape, c("null", SHAPES)) - 1L
n_id     <- function(n) if (n == 200L) 1L else if (n == 500L) 2L else if (n == N_PILOT) 3L else stop("n")
data_seed <- function(design, shape, n, ci, rep)
  300000000L + match(design, DESIGNS) * 10000000L + shape_id(shape) * 1000000L + n_id(n) * 100000L + ci * 5000L + rep
boot_seed <- function(dseed) dseed + 1000000000L
grp_seed  <- function(dseed) dseed + 500000000L
