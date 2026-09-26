# =============================================================================
# bagoft_fast_pre.R -- BAGofT 1.0.0 without the R overhead, INCLUDING the dcPre() variable
# pre-selection that parRF() runs when there are more than five covariates.
#
# This is explore_bagoft/bagoft_fast.R (verified identical to the package at p <= 5, left
# untouched) with the one branch it refused to reproduce added. The call it reproduces is
#   BAGofT(testGlmBi(formula, "logit"), parRF(ntree = 60), data, nsplits = 100,
#          ne = floor(5 sqrt(n)), nsim = 100)
#
# THE PRE-SELECTION (BAGofT 1.0.0, parRF -> dcPre(npreSel = 5, type = "V")). In every split,
# after the working model is fitted on the training part:
#   datRf       <- training covariates (every column except the response, data order) + res
#   nParV       <- ncol(datRf) - 1                       (parVar = "." by default)
#   if nParV > 5:
#     dcRes     <- t(dcov::mdcor(res, covariates))       one distance correlation per column,
#                                                        V-statistic; compiled, no random numbers
#     parVarNew <- names[order(-dcRes)[1:5]]             the five largest, in decreasing order
#   mtry        <- floor(length(parVarNew) / 3) = 1
#   forest      <- randomForest(res ~ v1 + ... + v5, datRf, ntree, maxnodes, mtry, importance = TRUE)
# maxnodes stays min(n_train, ceiling(5 * (ncol(data) - 1))), computed from ALL covariates.
# The forest therefore sees the five selected columns IN THE ORDER OF THE FORMULA (decreasing
# distance correlation). Column order matters: with mtry = 1 the variable tried at each node is
# drawn by position. Here the forest gets Z[, parVarNew] in that same order, through the x/y
# entry point that the formula method calls. dcov::mdcor draws no random numbers, so the random
# stream is the same as the package's: sample() per split, the forest, rbinom() per simulated
# response.
#
# With five or fewer covariates this function does exactly what bagoft_fast() does.
# verify_fast_pre.R checks it against the package.
# =============================================================================

bagoft_fast_pre <- function(dat, formula, nsplits = 100L, nsim = 100L, ntree = 60L,
                            ne = floor(5 * sqrt(nrow(dat))), npreSel = 5L) {
  Rsp  <- as.character(formula)[2]
  fam  <- stats::binomial(link = "logit")
  xcol <- setdiff(names(dat), Rsp)          # parRF uses every column except the response
  presel <- length(xcol) > npreSel          # dcPre(): nParV > npreSel

  mt <- stats::terms(formula, data = dat)
  X  <- stats::model.matrix(mt, stats::model.frame(mt, dat))
  Z  <- dat[, xcol, drop = FALSE]
  nr <- nrow(dat)
  nt <- nr - ne
  nmin     <- ceiling(sqrt(ne))
  Kmax     <- floor(ne / nmin)
  maxnodes <- min(nt, ceiling(5 * (ncol(dat) - 1)))

  glm_pred <- function(fit, Xn) {
    p   <- fit$rank
    piv <- if (p) fit$qr$pivot[seq_len(p)]
    fam$linkinv(drop(Xn[, piv, drop = FALSE] %*% fit$coefficients[piv]))
  }
  cut_codes <- function(x, breaks) {
    nb <- length(breaks <- sort.int(as.double(breaks)))
    if (anyDuplicated(breaks)) stop("'breaks' are not unique")
    factor(.bincode(x, breaks, TRUE, TRUE), levels = seq_len(nb - 1L))
  }

  one_split <- function(y) {
    trainIn <- sample(c(1:nr), nt)
    testIdx <- c(1:nr)[-trainIn]
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
    if (presel) {
      dc   <- t(dcov::mdcor(res, ZT))                       # dcPre's call, same arguments
      vsel <- xcol[order(-as.numeric(dc))[c(1:npreSel)]]
      mtry <- floor(length(vsel) / 3)
    } else {
      vsel <- xcol
      mtry <- 1                                             # floor(1/3) = 0, reset to 1
    }
    ZTs <- ZT[, vsel, drop = FALSE]
    rf  <- randomForest::randomForest(x = ZTs, y = res, ntree = ntree, maxnodes = maxnodes,
                                      mtry = mtry, importance = TRUE)
    trainsetPred <- stats::predict(rf, newdata = ZTs)
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
      testsetPred <- stats::predict(rf, newdata = Z[testIdx, vsel, drop = FALSE])
      gup <- cut_codes(testsetPred, c(Inf, -Inf, stats::quantile(urt, probs = seq(0, 1, 1/Ksel),
                                                                 names = FALSE)))
      gup <- droplevels(gup)
      freqTab <- table(gup)
      while (min(freqTab) < nmin) {
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

  multi <- function(y) {
    pvdat <- numeric(nsplits)
    for (j in c(1:nsplits)) pvdat[j] <- one_split(y)
    c(mean(pvdat), stats::median(pvdat), min(pvdat))
  }

  y0  <- dat[[Rsp]]
  obs <- multi(y0)
  if (nsim < 1)
    return(list(pmean = obs[1], pmedian = obs[2], pmin = obs[3]))
  pdat2 <- glm_pred(stats::glm.fit(X, y0, family = fam), X)
  sims  <- matrix(NA_real_, nsim, 3L)
  for (i in c(1:nsim)) {
    ydat2 <- sapply(pdat2, function(x) stats::rbinom(1, 1, x))
    sims[i, ] <- multi(ydat2)
  }
  list(p.value  = mean(obs[1] > sims[, 1]),
       p.value2 = mean(obs[2] > sims[, 2]),
       p.value3 = mean(obs[3] > sims[, 3]),
       pmean = obs[1], pmedian = obs[2], pmin = obs[3],
       simRes = list(pmeanSim = sims[, 1], pmediansim = sims[, 2], pminsim = sims[, 3]))
}
