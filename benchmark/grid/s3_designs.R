## s3_designs.R -- the covariate designs and the random departure family of Simulation I (the
## 60-cell grid): make_X, gen_design_null, make_family_p, gen_design_alt. Copied unchanged from the
## generator module of the deposited run; sourced by s3_grid.R.

## ---- covariate designs ------------------------------------------------------------------
## law: "indep" (mixed normal/uniform), "ar1" (AR(1) correlation rho), "skew" (one log-normal)
## cat_levels: 0 = none, else one factor with that many levels
make_X <- function(n, p, law = "indep", rho = 0.5, cat_levels = 0L) {
  pc <- max(1L, p - (cat_levels > 0L))          # continuous columns
  if (law == "ar1") {
    S <- rho^abs(outer(seq_len(pc), seq_len(pc), "-"))
    L <- chol(S)
    Z <- matrix(rnorm(n * pc), n, pc) %*% L
  } else {
    Z <- matrix(rnorm(n * pc), n, pc)
    if (pc >= 2) Z[, 2] <- runif(n, -3, 3)      # keep the mixed-scale flavour of the old design
  }
  if (law == "skew") Z[, 1] <- (exp(Z[, 1] * 0.6) - exp(0.18)) / 1.2
  d <- as.data.frame(Z); names(d) <- paste0("x", seq_len(pc))
  if (cat_levels > 0L)
    d[[paste0("x", pc + 1L)]] <- factor(sample.int(cat_levels, n, replace = TRUE))
  d
}

## Solve the intercept so the marginal event rate hits the target.
.b0_for_rate <- function(lin, target) {
  f <- function(b0) mean(plogis(b0 + lin)) - target
  tryCatch(stats::uniroot(f, c(-25, 25))$root, error = function(e) 0)
}

## A null (correctly specified) dataset with a controlled event rate.
gen_design_null <- function(n, p, law = "indep", rho = 0.5, cat_levels = 0L, rate = 0.5) {
  d <- make_X(n, p, law, rho, cat_levels)
  MM  <- stats::model.matrix(~ ., data = d)
  bet <- c(0, stats::runif(ncol(MM) - 1L, -0.8, 0.8))
  lin <- as.numeric(MM %*% bet)
  lin <- lin / max(stats::sd(lin), 1e-8) * 1.2                 # fix signal scale
  b0  <- .b0_for_rate(lin, rate)
  d$y <- stats::rbinom(n, 1, plogis(b0 + lin))
  frm <- stats::as.formula(paste("y ~", paste(names(d)[names(d) != "y"], collapse = " + ")))
  list(data = d, fit = suppressWarnings(stats::glm(frm, data = d, family = binomial())),
       rate = rate)
}

## The p-generalized random departure family (same nine kinds as make_family()).
make_family_p <- function(fseed, p) {
  set.seed(fseed)
  K <- sample(1:3, 1)
  terms <- replicate(K, {
    kind <- sample(c("quad","cubic","inter","sine","bump","step","abs","indexq","indexs"), 1,
                   prob = c(.15,.10,.15,.12,.12,.10,.06,.10,.10))
    i <- sample(seq_len(p), 1); j <- sample(seq_len(p), 1)
    list(kind=kind, i=i, j=j, w=runif(1,.5,3), ph=runif(1,0,2*pi), mu=runif(1,-2,2),
         s=runif(1,.15,.8), t=runif(1,-1,1), c=sample(c(-1,1),1)*runif(1,.5,1.5))
  }, simplify = FALSE)
  list(fseed = fseed, K = K, terms = terms, amp = runif(1, 0.25, 0.75))
}

gen_design_alt <- function(fam, n, p, law = "indep", rho = 0.5, cat_levels = 0L, rate = 0.5) {
  d <- make_X(n, p, law, rho, cat_levels)
  MM  <- stats::model.matrix(~ ., data = d)
  bet <- c(0, stats::runif(ncol(MM) - 1L, -0.8, 0.8))
  lin <- as.numeric(MM %*% bet); lin <- lin / max(stats::sd(lin), 1e-8) * 1.2
  ## numeric columns only carry the departure shapes
  num <- vapply(d, is.numeric, logical(1))
  Xn  <- as.matrix(d[, num, drop = FALSE]); q <- ncol(Xn)
  delta <- rep(0, n)
  for (tm in fam$terms) {
    xi <- Xn[, ((tm$i - 1L) %% q) + 1L]; xj <- Xn[, ((tm$j - 1L) %% q) + 1L]
    delta <- delta + tm$c * switch(tm$kind,
      quad = xi^2, cubic = xi^3/3, inter = xi*xj, sine = sin(tm$w*xi + tm$ph),
      bump = exp(-(xi - tm$mu)^2/tm$s), step = as.numeric(xi > tm$t), abs = abs(xi),
      indexq = (lin - mean(lin))^2/2, indexs = sin(tm$w*lin))
  }
  sdd <- stats::sd(delta); if (sdd > 1e-8) delta <- delta * (fam$amp / sdd)
  b0 <- .b0_for_rate(lin + delta, rate)
  d$y <- stats::rbinom(n, 1, plogis(b0 + lin + delta))
  frm <- stats::as.formula(paste("y ~", paste(names(d)[names(d) != "y"], collapse = " + ")))
  list(data = d, fit = suppressWarnings(stats::glm(frm, data = d, family = binomial())))
}
