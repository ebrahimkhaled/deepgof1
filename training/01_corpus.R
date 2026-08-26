## 01_corpus.R -- generate a DeepGOF training corpus. THE ONE FILE A DeepGOF-2 AUTHOR EDITS.
##
## The corpus defines the test. DeepGOF-1's measured behaviour -- correct size everywhere, power
## concentrated at small n, the specific blind spots -- is a property of THIS distribution over
## (null maps, alternative maps), not of the architecture (paper, Future Work direction 5: the
## same architecture scored +.090 under one departure prior and -.071 under another).
##
## Two design decisions are load-bearing. Change them knowingly:
##
##   1. DEPLOYMENT-MATCHED NULLS. The null class (y = 0) is NOT "data from a correct logistic
##      model". It is BOOTSTRAP REPLICATES: simulate from a FITTED model, refit, map. That is
##      exactly the object the deployed test compares against, and training against it pushes the
##      score toward pivotality -- the condition under which the bootstrap rank p-value is exact
##      (Lemma S.2). Training against plain correct-model nulls collapsed the deployed level in
##      our sweeps (score-AUC diagnostic, paper Fig. 4): the score learns to see the bootstrap
##      itself and calibration breaks.
##
##   2. THE DEPARTURE PRIOR. Alternatives are ~random FAMILIES (recipe below): 1-3 terms drawn
##      from {quad, cubic, inter, sine, bump, step, abs, indexq, indexs} on a 2- or 4-covariate
##      design, amplitude scaled to a target sd on the linear predictor. The prior is a
##      first-class, versioned object: publish yours with your weights.
##
## Output: ship-corpus CSV (y, cell_1..cell_36) consumed by 02_train.py.
## The EXACT corpus used for DeepGOF-1 is deposited alongside this script as ship_corpus.csv;
## this generator documents and reproduces its recipe (fresh seeds give a fresh corpus draw).

N_FAM   <- as.integer(Sys.getenv("DG_FAMILIES", "100"))   # departure families
PER_FAM <- as.integer(Sys.getenv("DG_PERFAM",   "42"))    # alternative maps per family
N_NULL  <- as.integer(Sys.getenv("DG_NULLS",    "4200"))  # bootstrap-replicate null maps
NOBS    <- as.integer(Sys.getenv("DG_N",        "200"))   # sample size per dataset
K       <- as.integer(Sys.getenv("DG_K",        "6"))     # map resolution (gate S4 selected 6)
SEED0   <- as.integer(Sys.getenv("DG_SEED",     "60101"))
OUT     <- Sys.getenv("DG_OUT", "corpus_new.csv")

## ---- the residual map: MUST stay identical to R/deepgof1.R's builder ---------------------------
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

## ---- the departure prior (verbatim from the DeepGOF-1 build; edit HERE for DeepGOF-2) ----------
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

## ---- build --------------------------------------------------------------------------------------
rows <- vector("list", N_FAM * PER_FAM + N_NULL); ri <- 0L
cat(sprintf("corpus: %d families x %d + %d bootstrap nulls, n=%d, K=%d\n",
            N_FAM, PER_FAM, N_NULL, NOBS, K))

## alternatives: label 1
for (f in seq_len(N_FAM)) {
  fam <- make_family(SEED0 + f)
  for (i in seq_len(PER_FAM)) {
    set.seed(SEED0 + f * 10000L + i)
    g <- gen_random(fam, NOBS)
    if (is.null(g$fit) || !g$fit$converged) next
    m <- tryCatch(map_of(g$fit, K), error = function(e) NULL)
    if (is.null(m)) next
    ri <- ri + 1L; rows[[ri]] <- c(1, m)
  }
  if (f %% 20 == 0) cat("  family", f, "of", N_FAM, "\n")
}

## nulls: label 0 -- BOOTSTRAP REPLICATES of correctly specified fits (deployment-matched)
set.seed(SEED0 * 7L)
b <- 0L
while (b < N_NULL) {
  ## a correct logistic dataset, fitted...
  x1 <- runif(NOBS, -3, 3); x2 <- rnorm(NOBS)
  eta0 <- 0.3 + 0.8*x1 - 0.5*x2
  d0 <- data.frame(x1, x2, y = rbinom(NOBS, 1, plogis(eta0)))
  f0 <- suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial()))
  if (!f0$converged) next
  ph <- fitted(f0)
  ## ...then ONE bootstrap replicate of that fit: simulate from pi-hat, refit, map
  d0$y <- rbinom(NOBS, 1, ph)
  fb <- tryCatch(suppressWarnings(glm(y ~ x1 + x2, data = d0, family = binomial())),
                 error = function(e) NULL)
  if (is.null(fb) || !fb$converged) next
  m <- tryCatch(map_of(fb, K), error = function(e) NULL)
  if (is.null(m)) next
  b <- b + 1L; ri <- ri + 1L; rows[[ri]] <- c(0, m)
}

M <- do.call(rbind, rows[seq_len(ri)])
colnames(M) <- c("y", paste0("c", seq_len(K * K)))
utils::write.csv(M, OUT, row.names = FALSE)
cat(sprintf("wrote %s: %d rows (%d alt, %d null)\n", OUT, nrow(M), sum(M[,1] == 1), sum(M[,1] == 0)))
