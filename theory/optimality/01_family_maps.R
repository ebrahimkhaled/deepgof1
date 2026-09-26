# =============================================================================
# 01_family_maps.R -- W2 data: held-out maps from the SHIPPED departure prior, per family.
#
# The prior and the map builder are copied verbatim from release/deepgof1/training/01_corpus.R
# (make_family, gen_random, map_of). Fresh dataset seeds throughout: none of these maps is in the
# training corpus. Two family sets:
#   train families  f = 1..100   (fseed = 60101 + f, the exact families the network was trained on)
#   fresh families  f = 101..200 (fseed = 70101 + f, new draws from the same prior)
# R_PER maps per family, plus NNULL bootstrap-null maps built exactly as the corpus builds them.
#
# Output: FAMILY_MAPS.csv (set, fam, y, c1..c36)
# Run:  Rscript 01_family_maps.R    Env: W2_RPER (300), W2_NNULL (30000), W2_WORKERS (8)
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/optimality"
RPER  <- as.integer(Sys.getenv("W2_RPER", "300"))
NNULL <- as.integer(Sys.getenv("W2_NNULL", "30000"))
NW    <- as.integer(Sys.getenv("W2_WORKERS", "8"))
K <- 6L; NOBS <- 200L

defs <- quote({
  K <- 6L; NOBS <- 200L
  map_of <- function(fit, K) {
    d  <- stats::model.frame(fit); mm <- stats::model.matrix(fit)
    vars <- colnames(mm)[-1]
    b  <- stats::coef(fit)[vars]
    sc <- abs(b) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
    ax <- vars[sort(order(sc, decreasing = TRUE)[1:2])]
    ph <- as.numeric(stats::fitted(fit)); r <- stats::model.response(d) - ph
    n  <- length(r)
    bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
    idx <- factor((bn(mm[, ax[1]]) - 1L) * K + bn(mm[, ax[2]]), levels = 1:(K * K))
    s <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
    s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
    as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8))
  }
  make_family <- function(fseed) {
    set.seed(fseed)
    four  <- runif(1) < 0.3
    vars  <- if (four) 4L else 2L
    Kt    <- sample(1:3, 1)
    KINDS <- c("quad","cubic","inter","sine","bump","step","abs","indexq","indexs")
    PR    <- c(.15,.10,.15,.12,.12,.10,.06,.10,.10)
    terms <- replicate(Kt, {
      kind <- sample(KINDS, 1, prob = PR)
      i <- sample(seq_len(vars), 1); j <- sample(setdiff(seq_len(vars), i), 1)
      list(kind = kind, i = i, j = j, w = runif(1, .5, 3), ph = runif(1, 0, 2*pi),
           mu = runif(1, -2, 2), s = runif(1, .15, .8), t = runif(1, -1, 1),
           c = sample(c(-1, 1), 1) * runif(1, .5, 1.5))
    }, simplify = FALSE)
    amp <- runif(1, 0.25, 0.75)
    list(fseed = fseed, four = four, K = Kt, terms = terms, amp = amp)
  }
  gen_random <- function(fam, n) {
    if (fam$four) { x1<-runif(n,-3,3); x2<-rnorm(n); x3<-rnorm(n); x4<-runif(n,-2,2)
      eta0 <- 0.2+0.6*x1-0.4*x2+0.3*x3-0.2*x4
      X <- cbind(x1,x2,x3,x4); frm <- y ~ x1+x2+x3+x4
    } else { x1<-runif(n,-3,3); x2<-rnorm(n)
      eta0 <- 0.3+0.8*x1-0.5*x2
      X <- cbind(x1,x2); frm <- y ~ x1+x2 }
    delta <- rep(0, n)
    for (tm in fam$terms) {
      xi <- X[, tm$i]; xj <- X[, tm$j %% ncol(X) + 1]
      delta <- delta + tm$c * switch(tm$kind,
        quad  = xi^2, cubic = xi^3/3, inter = xi*xj,
        sine  = sin(tm$w*xi + tm$ph),
        bump  = exp(-(xi - tm$mu)^2 / tm$s),
        step  = as.numeric(xi > tm$t), abs = abs(xi),
        indexq = (eta0 - mean(eta0))^2 / 2,
        indexs = sin(tm$w * eta0))
    }
    sdd <- stats::sd(delta); if (sdd > 1e-8) delta <- delta * (fam$amp / sdd)
    y <- rbinom(n, 1, plogis(eta0 + delta))
    d <- as.data.frame(X); d$y <- y
    list(data = d, fit = suppressWarnings(glm(frm, data = d, family = binomial())))
  }
  one_family <- function(job) {
    fam <- make_family(job$fseed); out <- matrix(NA_real_, 0, 36)
    for (i in seq_len(job$R)) {
      set.seed(job$dseed + i)
      g <- gen_random(fam, NOBS)
      if (is.null(g$fit) || !g$fit$converged) next
      m <- tryCatch(map_of(g$fit, K), error = function(e) NULL)
      if (!is.null(m)) out <- rbind(out, m)
    }
    data.frame(set = job$set, fam = job$f, y = 1L, out)
  }
  one_null_block <- function(job) {
    set.seed(job$dseed); out <- matrix(NA_real_, 0, 36)
    while (nrow(out) < job$R) {
      x1 <- runif(NOBS, -3, 3); x2 <- rnorm(NOBS)
      d0 <- data.frame(x1, x2, y = rbinom(NOBS, 1, plogis(0.3 + 0.8*x1 - 0.5*x2)))
      f0 <- suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial()))
      if (!f0$converged) next
      d0$y <- rbinom(NOBS, 1, fitted(f0))
      fb <- tryCatch(suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial())), error = function(e) NULL)
      if (is.null(fb) || !fb$converged) next
      m <- tryCatch(map_of(fb, K), error = function(e) NULL)
      if (!is.null(m)) out <- rbind(out, m)
    }
    data.frame(set = "null", fam = 0L, y = 0L, out)
  }
})

jobs <- c(lapply(1:100, function(f) list(set = "train", f = f, fseed = 60101L + f, dseed = 3100000L + f * 1000L, R = RPER, kind = "fam")),
          lapply(101:200, function(f) list(set = "fresh", f = f, fseed = 70101L + f, dseed = 3100000L + f * 1000L, R = RPER, kind = "fam")),
          lapply(1:30, function(b) list(dseed = 4200000L + b * 7919L, R = ceiling(NNULL / 30), kind = "null")))
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, "defs")
invisible(parallel::clusterEvalQ(cl, { eval(defs); NULL }))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, jobs, function(j) if (j$kind == "fam") one_family(j) else one_null_block(j))
parallel::stopCluster(cl)
D <- do.call(rbind, res); colnames(D) <- c("set", "fam", "y", paste0("c", 1:36))
utils::write.csv(D, file.path(HERE, "FAMILY_MAPS.csv"), row.names = FALSE)
cat(sprintf("wrote FAMILY_MAPS.csv: %d rows (%d null) in %.1f min\n", nrow(D), sum(D$y == 0),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
