# =============================================================================
# bagoft_fast.R -- BAGofT 1.0.0, the same test, computed without the R overhead.
#
# WHAT THIS IS. A re-expression of BAGofT::BAGofT() exactly as the trial calls it:
#   testModel = testGlmBi(formula, link = "logit"), parFun = parRF(ntree = 60),
#   nsplits = 100, nsim = 100, ne = floor(5 * sqrt(n)).
# Every random number is drawn by the same call, in the same order, as in the package:
# sample() for each split, randomForest() for each forest, rbinom() for each simulated
# response. Everything deterministic is the same arithmetic on the same numbers, reached
# through cheaper entry points:
#   glm(formula, data)             -> glm.fit(X, y)            glm() calls glm.fit itself
#   predict(glm, se.fit = TRUE)    -> linkinv(X %*% beta)       predict.lm's own product,
#                                                               without the summary.glm()
#                                                               that se.fit = TRUE runs
#   resid(glm, "pearson")          -> residuals.glm's formula
#   xtabs(v ~ f)                   -> tapply(v, f, sum, default = 0L)   what xtabs calls
#   cut(x, breaks, TRUE)           -> .bincode(x, breaks, TRUE, TRUE)   cut's own codes;
#                                     cut's labels are always distinct ("Range_i" is its
#                                     fallback), so the grouping is the codes alone
#   randomForest(res ~ ., datRf)   -> randomForest(x, y)       the formula method's target
#   mtry = floor(1/3) = 0          -> mtry = 1                 parRF computes mtry from
#                                     length(".") = 1, so it is always 0; randomForest warns
#                                     and resets it to max(1, min(p, round(0))) = 1. Passing
#                                     1 grows the same forest without 10,100 warnings a test.
# NOTHING about the statistic, its defaults or its random draws is changed. verify_fast.R
# checks the outputs are IDENTICAL to BAGofT's -- identical(), not all.equal().
#
# WHY. A profile (2026-09-19) showed the package spends most of each test in formula
# parsing, model.frame(), xtabs(), label formatting inside cut() and the predict()
# wrappers, 10,100 times over. The forest-growing C code is about 4% of it at n = 50.
#
# SCOPE. Numeric covariates, at most five. Above five, parRF's dcPre() pre-selection runs,
# and that branch is not reproduced here, so the function stops rather than guess. One
# covariate works, and matches parRF_dropfix.R -- BAGofT 1.0.0 itself crashes there.
# =============================================================================

bagoft_fast <- function(dat, formula, nsplits = 100L, nsim = 100L, ntree = 60L,
                        ne = floor(5 * sqrt(nrow(dat)))) {
  Rsp  <- as.character(formula)[2]
  fam  <- stats::binomial(link = "logit")
  # parRF grows its forest on EVERY column except the response, whatever the formula says
  xcol <- setdiff(names(dat), Rsp)
  if (length(xcol) > 5L)
    stop("more than five covariates: parRF would pre-select with dcPre(), not reproduced here")

  mt <- stats::terms(formula, data = dat)
  X  <- stats::model.matrix(mt, stats::model.frame(mt, dat))    # glm's design, built once
  Z  <- dat[, xcol, drop = FALSE]                               # what the forest sees
  nr <- nrow(dat)
  nt <- nr - ne
  nmin     <- ceiling(sqrt(ne))                                 # parRF defaults
  Kmax     <- floor(ne / nmin)
  maxnodes <- min(nt, ceiling(5 * (ncol(dat) - 1)))

  # predict.lm's product for a fitted glm, pivoting exactly as it does when rank-deficient
  glm_pred <- function(fit, Xn) {
    p   <- fit$rank
    piv <- if (p) fit$qr$pivot[seq_len(p)]
    fam$linkinv(drop(Xn[, piv, drop = FALSE] %*% fit$coefficients[piv]))
  }

  # cut(x, breaks, include.lowest = TRUE) as a factor of bin codes
  cut_codes <- function(x, breaks) {
    nb <- length(breaks <- sort.int(as.double(breaks)))
    if (anyDuplicated(breaks)) stop("'breaks' are not unique")
    factor(.bincode(x, breaks, TRUE, TRUE), levels = seq_len(nb - 1L))
  }

  # BAGofT_sin + testGlmBi + parRF, for one split of the response y
  one_split <- function(y) {
    trainIn <- sample(c(1:nr), nt)
    testIdx <- c(1:nr)[-trainIn]              # datset[-trainIn, ] keeps the original order
    yT <- y[trainIn]
    yE <- y[testIdx]
    XT <- X[trainIn, , drop = FALSE]
    XE <- X[testIdx, , drop = FALSE]

    fit   <- stats::glm.fit(XT, yT, family = fam)
    predT <- glm_pred(fit, XT)
    predE <- glm_pred(fit, XE)
    mu    <- fit$fitted.values
    res   <- (fit$y - mu) * sqrt(fit$prior.weights) / sqrt(fam$variance(mu))

    ZT <- Z[trainIn, , drop = FALSE]
    rf <- randomForest::randomForest(x = ZT, y = res, ntree = ntree, maxnodes = maxnodes,
                                     mtry = 1, importance = TRUE)
    trainsetPred <- stats::predict(rf, newdata = ZT)
    rt  <- round(trainsetPred, digits = 10)
    urt <- unique(rt)

    if (length(urt) == 1) {
      gup <- as.factor(rep(1, ne))
    } else {
      Kmax_adj    <- min(Kmax, length(urt))
      chitrainVec <- numeric(Kmax_adj)
      dT <- predT - yT
      vT <- predT * (1 - predT)
      for (gt in c(1:Kmax_adj)) {
        gupt    <- cut_codes(rt, stats::quantile(urt, probs = seq(0, 1, 1/gt), names = FALSE))
        dift    <- abs(tapply(dT, gupt, sum, default = 0L))
        dent    <- tapply(vT, gupt, sum, default = 0L)
        contrit <- (dift)^2/dent
        chitrainVec[gt] <- sum(contrit)
        if (is.nan(sum(contrit))) chitrainVec[gt] <- chitrainVec[(gt - 1)]
      }
      chitrainVec4 <- chitrainVec[-1] - chitrainVec[-length(chitrainVec)]
      Ksel <- which.max(chitrainVec4) + 1
      testsetPred <- stats::predict(rf, newdata = Z[testIdx, , drop = FALSE])
      gup <- cut_codes(testsetPred, c(Inf, -Inf, stats::quantile(urt, probs = seq(0, 1, 1/Ksel),
                                                                 names = FALSE)))
      gup <- droplevels(gup)
      freqTab <- table(gup)
      while (min(freqTab) < nmin) {             # parRF's merge of undersized groups, verbatim
        levels(gup)[which(levels(gup) == names(sort(freqTab))[1])] <- names(sort(freqTab))[2]
        freqTab <- table(gup)
      }
    }

    ngp     <- length(levels(gup))
    dif     <- abs(tapply(predE - yE, gup, sum, default = 0L))
    den     <- tapply(predE * (1 - predE), gup, sum, default = 0L)
    contri  <- (dif)^2/den
    chisq   <- sum(contri)
    1 - stats::pchisq(chisq, ngp)
  }

  # BAGofT_multi
  multi <- function(y) {
    pvdat <- numeric(nsplits)
    for (j in c(1:nsplits)) pvdat[j] <- one_split(y)
    c(mean(pvdat), stats::median(pvdat), min(pvdat))
  }

  # BAGofT: the observed data first, then nsim responses simulated from the full-data fit
  y0  <- dat[[Rsp]]
  obs <- multi(y0)
  if (nsim < 1)
    return(list(pmean = obs[1], pmedian = obs[2], pmin = obs[3]))
  pdat2 <- glm_pred(stats::glm.fit(X, y0, family = fam), X)
  sims  <- matrix(NA_real_, nsim, 3L)
  for (i in c(1:nsim)) {
    ydat2 <- sapply(pdat2, function(x) stats::rbinom(1, 1, x))   # the package's own draw
    sims[i, ] <- multi(ydat2)
  }
  list(p.value  = mean(obs[1] > sims[, 1]),
       p.value2 = mean(obs[2] > sims[, 2]),
       p.value3 = mean(obs[3] > sims[, 3]),
       pmean = obs[1], pmedian = obs[2], pmin = obs[3],
       simRes = list(pmeanSim = sims[, 1], pmediansim = sims[, 2], pminsim = sims[, 3]))
}
