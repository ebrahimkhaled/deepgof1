## check_against_package.R -- R/deepgof1.R must return the same p-value as deepgof1() in the R package
## ebrahim.gof (2.9.0 or later) for the same seed: plain, spline, one-covariate and tied models, both readings.
## Run from the release root:  Rscript R/check_against_package.R
source("R/deepgof1.R")
repo <- deepgof1
stopifnot(utils::packageVersion("ebrahim.gof") >= "2.9.0")
pkg <- ebrahim.gof::deepgof1
set.seed(1); n <- 150
d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n), x3 = rnorm(n), g = rbinom(n, 1, .4))
d$y <- rbinom(n, 1, plogis(0.3 + 0.8 * d$x1 - 0.5 * d$x2 + 0.6 * (d$x1^2 - 3) + 0.4 * d$g))
models <- list(plain  = y ~ x1 + x2 + x3,
               spline = y ~ splines::ns(x1, 3) + x2,
               one    = y ~ x1,
               tied   = y ~ x1 + g)
ok <- TRUE
for (nm in names(models)) for (rd in c("v1", "allpairs")) {
  f <- glm(models[[nm]], data = d, family = binomial())
  set.seed(7); a <- repo(f, B = 39, reading = rd)$p.value
  set.seed(7); b <- pkg(f, B = 39, reading = rd)$p.value
  cat(sprintf("%-7s %-9s archive %.3f  package %.3f  %s\n", nm, rd, a, b, if (identical(a, b)) "same" else "DIFFERENT"))
  ok <- ok && identical(a, b)
}
stopifnot(ok)
cat("all identical\n")
