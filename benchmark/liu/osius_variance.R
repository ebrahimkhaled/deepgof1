# =============================================================================
# osius_variance.R -- why the Osius-Rojek level fails on Liu settings 3 and 4.
#
# The test refers (X^2 - (J - p)) to a normal whose variance, when every covariate pattern is
# unique, is the WEIGHTED RESIDUAL sum of squares of
#       c_i = (1 - 2 p_i) / {p_i (1 - p_i)}
# after regressing c on the model's design matrix with weights p_i(1 - p_i). That residual is the
# whole denominator. It collapses when c is nearly a linear function of the covariates, or nearly
# zero -- and then |z| explodes and the test rejects true models.
#
# This script measures how much of c's weighted variance SURVIVES that projection, per setting and
# sample size, on the same null datasets the benchmark uses (liu_grid.R's seeding, both corpora).
# It writes results/liu/OSIUS_VARIANCE.csv, which tables/make_tables.py turns into the macros the
# discussion cites. Nothing here re-runs a test: it is arithmetic on refits of the null datasets.
#
# Run from the archive root:  Rscript benchmark/liu/osius_variance.R
# =============================================================================
PAP  <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
REPS <- as.integer(Sys.getenv("OV_REPS", "150"))       # null datasets per (setting, n, corpus)
SEEDS <- c(101000L, 102000L)

## gen_liu() and the null-block seeding, verbatim from benchmark/liu/liu_grid.R
gen_liu <- function(n, setting, C) {
  if (setting == 1L) { x1 <- runif(n,-3,3); d <- data.frame(x1); eta <- -1.963 + 0.982*x1 + C*x1^2
  } else if (setting == 2L) { x1 <- runif(n,-3,3); x2 <- rnorm(n,0,1.5); x3 <- rt(n,4)
    d <- data.frame(x1,x2,x3); eta <- -2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2
  } else if (setting == 3L) { x1 <- runif(n,-3,3); x2 <- rbinom(n,1,.5); d <- data.frame(x1,x2)
    eta <- -1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2
  } else { x1 <- runif(n,-3,3); x2 <- runif(n,-3,3); d <- data.frame(x1,x2); eta <- 0.3*x1 + 0.3*x2 + C*x1*x2 }
  d$y <- rbinom(n, 1, plogis(eta)); d
}
null_seed <- function(SEED0, SET, i) as.integer(SEED0) + as.integer(SET)*500000L + as.integer(i)

rows <- list()
for (SET in 1:4) for (NOBS in c(50L, 100L, 500L)) {
  share <- numeric(0)
  for (SEED0 in SEEDS) for (i in seq_len(REPS)) {
    set.seed(null_seed(SEED0, SET, i))
    d <- gen_liu(NOBS, SET, 0)
    if (length(unique(d$y)) < 2L) next
    frm <- stats::as.formula(paste("y ~", paste(setdiff(names(d), "y"), collapse = " + ")))
    fit <- tryCatch(suppressWarnings(stats::glm(frm, data = d, family = binomial())),
                    error = function(e) NULL)
    if (is.null(fit) || !fit$converged) next
    X  <- stats::model.matrix(fit)
    ph <- as.numeric(stats::fitted(fit)); V <- ph * (1 - ph)
    cv <- (1 - 2 * ph) / V
    w  <- tryCatch(stats::lm.wfit(X, cv, V), error = function(e) NULL)
    if (is.null(w)) next
    rss <- sum(V * w$residuals^2)                                   # the test's denominator
    tot <- sum(V * (cv - stats::weighted.mean(cv, V))^2)            # before the projection
    if (is.finite(rss) && is.finite(tot) && tot > 0) share <- c(share, rss / tot)
  }
  rows[[length(rows) + 1L]] <- data.frame(setting = SET, n = NOBS, ndatasets = length(share),
                                          pct_var_left = 100 * median(share))
}
R <- do.call(rbind, rows)
utils::write.csv(R, file.path(PAP, "results", "liu", "OSIUS_VARIANCE.csv"), row.names = FALSE)
print(R, row.names = FALSE, digits = 3)
cat("\nsettings 1-2:", sprintf("%.0f-%.0f%%", min(R$pct_var_left[R$setting <= 2]), max(R$pct_var_left[R$setting <= 2])),
    "| setting 3:", sprintf("%.1f%%", median(R$pct_var_left[R$setting == 3])),
    "| setting 4:", sprintf("%.1f%%", median(R$pct_var_left[R$setting == 4])), "\n")
