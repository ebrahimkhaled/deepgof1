## rate_common.R -- null designs for the level-rate study (2026-09-24)
## The generators are copied verbatim from code/universal_design.R (make_X, .b0_for_rate,
## gen_design_null) so this study needs nothing else in scope. Designs D1, D3, D4 are the
## S3 grid designs of Simulation I (benchmark/grid/s3_designs.R); D2 (factor) is not used here.

suppressPackageStartupMessages(library(ebrahim.gof))

DESIGNS <- list(
  D1 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.50),   # textbook p=2, 50% events
  D3 = list(p=3L, law="ar1",   rho=0.7, cat=0L, rate=0.50),   # p=3, AR(1) rho=.7
  D4 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.15))   # p=2, 15% events

make_X <- function(n, p, law = "indep", rho = 0.5, cat_levels = 0L) {
  pc <- max(1L, p - (cat_levels > 0L))
  if (law == "ar1") {
    S <- rho^abs(outer(seq_len(pc), seq_len(pc), "-"))
    L <- chol(S)
    Z <- matrix(rnorm(n * pc), n, pc) %*% L
  } else {
    Z <- matrix(rnorm(n * pc), n, pc)
    if (pc >= 2) Z[, 2] <- runif(n, -3, 3)
  }
  if (law == "skew") Z[, 1] <- (exp(Z[, 1] * 0.6) - exp(0.18)) / 1.2
  d <- as.data.frame(Z); names(d) <- paste0("x", seq_len(pc))
  if (cat_levels > 0L)
    d[[paste0("x", pc + 1L)]] <- factor(sample.int(cat_levels, n, replace = TRUE))
  d
}

.b0_for_rate <- function(lin, target) {
  f <- function(b0) mean(plogis(b0 + lin)) - target
  tryCatch(stats::uniroot(f, c(-25, 25))$root, error = function(e) 0)
}

gen_design_null <- function(n, p, law = "indep", rho = 0.5, cat_levels = 0L, rate = 0.5) {
  d <- make_X(n, p, law, rho, cat_levels)
  MM  <- stats::model.matrix(~ ., data = d)
  bet <- c(0, stats::runif(ncol(MM) - 1L, -0.8, 0.8))
  lin <- as.numeric(MM %*% bet)
  lin <- lin / max(stats::sd(lin), 1e-8) * 1.2
  b0  <- .b0_for_rate(lin, rate)
  d$y <- stats::rbinom(n, 1, plogis(b0 + lin))
  frm <- stats::as.formula(paste("y ~", paste(names(d)[names(d) != "y"], collapse = " + ")))
  list(data = d, fit = suppressWarnings(stats::glm(frm, data = d, family = binomial())),
       rate = rate)
}

## One replicate: the analyst's own plain-glm fit handed to the released test, B = 199.
## Returns the p-value, the statistic, a separation flag (any |coef| > 8 or non-finite),
## and the number of bootstrap refits that errored (scored +Inf by the released code).
one_rep <- function(design, n, B = 199L) {
  cfg <- DESIGNS[[design]]
  g <- gen_design_null(n, cfg$p, cfg$law, cfg$rho, cfg$cat, cfg$rate)
  cf <- coef(g$fit)
  sep <- as.integer(any(!is.finite(cf)) || any(abs(cf) > 8))
  r <- tryCatch(deepgof1(g$fit, B = B, K = 6L), error = function(e) NULL)
  if (is.null(r)) return(c(p = NA_real_, stat = NA_real_, sep = sep, ninf = NA_real_))
  c(p = r$p.value, stat = r$statistic, sep = sep, ninf = sum(!is.finite(r$boot)))
}
