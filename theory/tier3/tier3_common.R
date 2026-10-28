# tier3_common.R -- shared definitions for the Tier 3 scripts (sourced, never run on its own).
PAP   <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
T3    <- file.path(PAP, "theory/tier3")
RDEP  <- file.path(PAP, "code/rdeploy")
NETS  <- c(paste0("Control_s", 1:3), paste0("NoInter_s", 1:3), paste0("NoQuad_s", 1:3), paste0("NoQuadStrict_s", 1:3))

## the default RNG, reset before every set.seed (the L'Ecuyer rule)
rng_default <- function() RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seed_at <- function(s) { rng_default(); set.seed(s) }
## run a third-party call without letting it move the caller's RNG stream
keep_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  kind <- RNGkind()
  on.exit({ do.call(RNGkind, as.list(kind)); if (had) assign(".Random.seed", old, envir = globalenv()) })
  expr
}

## ---- data-generating processes, verbatim ----
gen_liu <- function(n, setting, C) {                   # code/liu_grid.R
  if (setting == 1L) {
    x1 <- runif(n,-3,3); d <- data.frame(x1)
    eta <- -1.963 + 0.982*x1 + C*x1^2
  } else if (setting == 2L) {
    x1 <- runif(n,-3,3); x2 <- rnorm(n,0,1.5); x3 <- rt(n,4)
    d <- data.frame(x1,x2,x3)
    eta <- -2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2
  } else if (setting == 3L) {
    x1 <- runif(n,-3,3); x2 <- rbinom(n,1,.5); d <- data.frame(x1,x2)
    eta <- -1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2
  } else {
    x1 <- runif(n,-3,3); x2 <- runif(n,-3,3); d <- data.frame(x1,x2)
    eta <- 0.3*x1 + 0.3*x2 + C*x1*x2
  }
  d$y <- rbinom(n, 1, plogis(eta)); d
}
liu_seed  <- function(seed0, S, k, ci, i) seed0 + S * 500000L + k * 50000L + ci * 5000L + i
liu_dgseed <- function(S, n, k, ci, i) 880000L + S * 100000L + n * 10L + k * 7L + ci + i * 13L   # run_shipped_liu.R

## ---- the score statistic for adding columns Z to a fitted logistic model (02_pilot_A_C.R) ----
score_stat <- function(X, ph, r, Z) {
  w <- ph * (1 - ph)
  keep <- apply(Z, 2, function(z) stats::sd(z) > 1e-10)
  Z <- Z[, keep, drop = FALSE]; if (!ncol(Z)) return(c(T = NA, df = 0))
  XtW <- t(X * w); A <- solve(XtW %*% X, XtW %*% Z)
  Zr <- Z - X %*% A
  q <- qr(Zr * sqrt(w)); Zr <- Zr[, q$pivot[seq_len(q$rank)], drop = FALSE]
  if (!ncol(Zr)) return(c(T = NA, df = 0))
  U <- crossprod(Zr, r); V <- crossprod(Zr * w, Zr)
  c(T = as.numeric(crossprod(U, solve(V, U))), df = ncol(Zr))
}

## ---- batched forward pass (theory/niche/run_niche.R load_score_batch, generalized to any export) ----
make_batch_scorer <- function(wfile, sfile, K = 6L, fast = FALSE) {
  e <- new.env(); old <- setwd(RDEP); on.exit(setwd(old)); sys.source("deepgof_forward.R", envir = e)
  W <- e$load_weights(wfile); ms <- as.matrix(utils::read.csv(sfile, header = FALSE))
  flat <- function(Wt) { O <- dim(Wt)[1]; C <- dim(Wt)[2]; Wm <- matrix(0, O, C * 9L); k <- 0L
    for (di in 1:3) for (dj in 1:3) { k <- k + 1L; Wm[, ((k - 1L) * C + 1L):(k * C)] <- Wt[, , di, dj, drop = FALSE][, , 1, 1] }
    Wm }
  W0 <- flat(W[["conv.0.weight"]]); W2 <- flat(W[["conv.2.weight"]]); W5 <- flat(W[["conv.5.weight"]])
  convb <- function(x, Wm, b) {
    C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]; N <- dim(x)[4]
    xp <- array(0, c(C, H + 2L, Wd + 2L, N)); xp[, 2:(H + 1L), 2:(Wd + 1L), ] <- x
    cols <- matrix(0, C * 9L, H * Wd * N); k <- 0L
    for (di in 0:2) for (dj in 0:2) { k <- k + 1L
      cols[((k - 1L) * C + 1L):(k * C), ] <- matrix(xp[, (1L + di):(H + di), (1L + dj):(Wd + dj), , drop = FALSE], C, H * Wd * N) }
    array(pmax(Wm %*% cols + b, 0), c(nrow(Wm), H, Wd, N))
  }
  zscore <- function(Z) {                               # rows of Z: maps ALREADY standardized
    N <- nrow(Z)
    x <- array(aperm(array(t(Z), c(K, K, N)), c(2, 1, 3)), c(1L, K, K, N))
    x <- convb(x, W0, W[["conv.0.bias"]]); x <- convb(x, W2, W[["conv.2.bias"]])
    o <- seq(1L, K, 2L); ee <- o + 1L
    x <- pmax(x[, o, o, , drop = FALSE], x[, ee, o, , drop = FALSE], x[, o, ee, , drop = FALSE], x[, ee, ee, , drop = FALSE])
    x <- convb(x, W5, W[["conv.5.bias"]])
    if (fast) {                                          # global max/mean without apply(): same values, the
      O <- dim(x)[1]; N <- dim(x)[4]                     # mean differing from apply(mean) by <= 1 ulp
      xm <- matrix(aperm(x, c(1, 4, 2, 3)), O * N, 9L)
      mx <- xm[, 1]; for (j in 2:9) mx <- pmax(mx, xm[, j])
      E <- rbind(matrix(mx, O, N), matrix(rowMeans(xm), O, N))
    } else
    E <- rbind(apply(x, c(1, 4), max), apply(x, c(1, 4), mean))
    h <- pmax(W[["head.0.weight"]] %*% E + as.numeric(W[["head.0.bias"]]), 0)
    as.numeric(W[["head.2.weight"]] %*% h + as.numeric(W[["head.2.bias"]]))
  }
  list(score = function(M) zscore(sweep(sweep(M, 2, ms[1, ], "-"), 2, ms[2, ], "/")),
       zscore = zscore, W = W, ms = ms, unbatched = function(z) e$deepgof_score(z, W, K))
}
load_extra_nets <- function() {
  L <- lapply(NETS, function(nm) make_batch_scorer(file.path(T3, "nets", nm, "weights.csv"), file.path(T3, "nets", nm, "scaler.csv"), fast = TRUE)$score)
  names(L) <- NETS; L
}
## p-value of a score against its bootstrap (failed refits are rows of NA in M -> +Inf, as deepgof1 does)
boot_p <- function(s) { s[is.na(s)] <- Inf; (1 + sum(s[-1] >= s[1])) / length(s) }
