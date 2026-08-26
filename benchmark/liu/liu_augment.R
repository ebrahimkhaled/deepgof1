suppressMessages(library(brglm2))
## liu_augment.R -- add Stute-Zhu and proj (Escanciano RMEP) to an existing GATE LIU block.
##
## Pre-declared in speckit/26_LIU_BENCHMARK.md, AMENDMENT (2026-08-23).
##
## WHY THIS IS CHEAP. Every dataset in a LIU block is generated from a DETERMINISTIC seed
## (`set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i)`), so the exact same data can be
## re-derived here without touching the stored blocks. The expensive part of the original run was
## 199 bootstrap REFITS per dataset for the residual map; that is NOT needed for Stute-Zhu, and for
## proj it is needed only through a rank-1 quadratic form:
##
##   ⭐ THE A MATRIX DEPENDS ONLY ON X, NEVER ON y.
##     A is built ONCE per dataset -- O(n^3 p), vectorised over the third index -- and every
##     bootstrap replicate then costs only O(n^2) for eps' A eps. That is what turns Liu et al.'s
##     reported ~2.5 h per setting into something that runs here.
##
## proj: T = eps-hat' A eps-hat / n^2, with
##   A_ij = sum_r ||X_i-X_r|| ||X_j-X_r|| (pi - arccos(a_ijr)),
##   a_ijr = (X_i-X_r)'(X_j-X_r)/(||X_i-X_r|| ||X_j-X_r||),   X standardised, intercept dropped.
## Calibrated by the SAME model-based bootstrap as every other arm: simulate y* from the fitted
## model, refit, recompute eps-hat*, exact rank p-value at B = 199.
##
## VALIDATION (LIU6, frozen): our proj must show high, MONOTONE power on Setting 2. If it does not,
## the implementation is wrong and no proj number may be reported.

ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def"
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
source(file.path(ROOT,"jasa_study","config.R")); source(file.path(ROOT,"R","scenarios.R"))
source(file.path(JASA$study,"R","scenarios_registry.R")); source(file.path(JASA$study,"R","methods_registry.R"))
source(file.path(JASA$study,"R","harness_core.R")); jasa_init()

B     <- as.integer(Sys.getenv("AU_B",    "199"))
NW    <- as.integer(Sys.getenv("AU_W",    "12"))
NOBS  <- as.integer(Sys.getenv("AU_N",    "100"))
SEED0 <- as.integer(Sys.getenv("AU_SEED", "101000"))
SET   <- as.integer(Sys.getenv("AU_SET",  "2"))
DOPROJ<- as.integer(Sys.getenv("AU_PROJ", "1"))
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)
OUT   <- file.path(PAP,"results","liu", sprintf("S%d_n%d_s%d", SET, NOBS, SEED0))
if (!dir.exists(OUT)) stop("no such block: ", OUT)
META  <- read.csv(file.path(OUT,"blk_meta.csv"))
NNULL <- META$nnull[1]; NRAND <- META$nrand[1]; NEVAL <- META$neval[1]

## ---- the four DGPs, identical to liu_grid.R ----------------------------------------------------
gen_liu <- function(n, setting, C) {
  if (setting == 1L) { x1 <- runif(n,-3,3); d <- data.frame(x1)
    eta <- -1.963 + 0.982*x1 + C*x1^2
  } else if (setting == 2L) { x1 <- runif(n,-3,3); x2 <- rnorm(n,0,1.5); x3 <- rt(n,4)
    d <- data.frame(x1,x2,x3); eta <- -2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2
  } else if (setting == 3L) { x1 <- runif(n,-3,3); x2 <- rbinom(n,1,.5); d <- data.frame(x1,x2)
    eta <- -1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2
  } else { x1 <- runif(n,-3,3); x2 <- runif(n,-3,3); d <- data.frame(x1,x2)
    eta <- 0.3*x1 + 0.3*x2 + C*x1*x2 }
  d$y <- rbinom(n, 1, plogis(eta)); d
}
gen_train_alt <- function(n, setting, fseed) {
  set.seed(fseed); d <- gen_liu(n, setting, 0)
  X <- as.matrix(d[, setdiff(names(d),"y"), drop=FALSE])
  p <- ncol(X); j <- sample(p,1); k <- sample(p,1)
  kind <- sample(c("quad","inter","bump","step","sine","abs"),1)
  xi <- X[,j]; xj <- X[,k]
  delta <- switch(kind, quad=xi^2, inter=xi*xj, bump=exp(-(xi-runif(1,-1.5,1.5))^2/runif(1,.15,.6)),
                  step=as.numeric(xi>runif(1,-1,1)), sine=sin(runif(1,.8,2.5)*xi+runif(1,0,pi)),
                  abs=abs(xi))
  sdd <- stats::sd(delta); amp <- runif(1,.4,2.0)*sqrt(100/n)
  if (sdd > 1e-8) delta <- delta*(amp/sdd)
  eta <- switch(as.character(setting),
    "1" = -1.963 + 0.982*X[,1], "2" = -2 + 0.3*X[,1] + 0.3*X[,2] + 0.3*X[,3],
    "3" = -1.972 + 0.135*X[,1] + 1.117*X[,2], "4" = 0.3*X[,1] + 0.3*X[,2])
  d$y <- rbinom(n, 1, plogis(eta + delta)); d
}
FRM <- function(d) stats::as.formula(paste("y ~", paste(setdiff(names(d),"y"), collapse=" + ")))
fit_pol <- function(frm, d) {
  g <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial())), error=function(e) NULL)
  if (is.null(g)) return(NULL); b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1])>8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(frm,data=d,family=binomial(),method="brglmFit",type="AS_mean")),
                  error=function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(f)
  }
  g
}

## ---- Escanciano's A matrix: built ONCE from X ---------------------------------------------------
proj_A <- function(X) {
  X <- scale(X); n <- nrow(X)
  A <- matrix(0, n, n)
  for (r in seq_len(n)) {
    D  <- sweep(X, 2, X[r, ], "-")                 # n x p
    nm <- sqrt(rowSums(D*D))
    G  <- tcrossprod(D)                            # n x n inner products
    NN <- outer(nm, nm)
    cs <- ifelse(NN > 1e-12, G/NN, 0)
    A  <- A + NN * (pi - acos(pmin(pmax(cs, -1), 1)))
  }
  A
}
proj_T <- function(e, A) as.numeric(crossprod(e, A %*% e)) / (length(e)^2)

one <- function(kind, i, C=NA_real_, ci=0L) {
  set.seed(SEED0 + SET*500000L + kind*50000L + ci*5000L + i)
  d <- if (kind == 1L) gen_train_alt(NOBS, SET, SEED0*29L + i)
       else gen_liu(NOBS, SET, if (kind==0L) 0 else C)
  if (length(unique(d$y)) < 2L) return(NULL)
  frm <- FRM(d); fit <- fit_pol(frm, d); if (is.null(fit)) return(NULL)
  sz <- tryCatch(suppressWarnings(
          run.all.gof(fit, tests="Stute-Zhu", include_slow=TRUE, install="no")$p_value[1]),
        error=function(e) NA_real_)
  pj <- NA_real_
  if (DOPROJ == 1L) {
    pj <- tryCatch({
      X <- stats::model.matrix(fit)[, -1, drop=FALSE]
      A <- proj_A(X)
      ph <- as.numeric(stats::fitted(fit))
      Tobs <- proj_T(d$y - ph, A)
      Tb <- numeric(B); d2 <- d
      for (b in seq_len(B)) {
        d2$y <- rbinom(NOBS, 1, ph)
        fb <- fit_pol(frm, d2)
        Tb[b] <- if (is.null(fb)) -Inf else proj_T(d2$y - as.numeric(stats::fitted(fb)), A)
      }
      (1 + sum(Tb >= Tobs)) / (B + 1)
    }, error = function(e) NA_real_)
  }
  c(kind=kind, rep=i, C=C, SZ=sz, proj=pj)
}

SRC <- c(file.path(ROOT,"jasa_study","config.R"), file.path(ROOT,"R","scenarios.R"))
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("SRC","B","NOBS","SEED0","SET","DOPROJ","gen_liu","gen_train_alt",
                              "FRM","fit_pol","proj_A","proj_T","one"), envir=environment())
invisible(parallel::clusterEvalQ(cl, {
  source(SRC[1]); source(file.path(JASA$study,"R","methods_registry.R"))
  source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
  source(SRC[2]); source(file.path(JASA$study,"R","scenarios_registry.R"))
  suppressMessages(library(brglm2)); NULL }))

jobs <- c(lapply(seq_len(NNULL), function(i) list(k=0L,i=i)),
          lapply(seq_len(NRAND), function(i) list(k=1L,i=i)))
for (ci in seq_along(CS)) if (CS[ci] > 0)
  jobs <- c(jobs, lapply(seq_len(NEVAL), function(i) list(k=2L,i=i,C=CS[ci],ci=ci)))
cat(sprintf("AUGMENT S%d n=%d seed=%d | Stute-Zhu%s | %d datasets\n",
            SET, NOBS, SEED0, if (DOPROJ==1L) " + proj" else "", length(jobs))); flush.console()
t0 <- Sys.time()
res <- Filter(Negate(is.null), parallel::parLapplyLB(cl, jobs, function(j)
  one(j$k, j$i, if (is.null(j$C)) NA_real_ else j$C, if (is.null(j$ci)) 0L else j$ci)))
parallel::stopCluster(cl)

R <- as.data.frame(do.call(rbind, res))
write.csv(R, file.path(OUT,"augment_pvalues.csv"), row.names=FALSE)
cat(sprintf("AUGMENT COMPLETE: %d rows, %.1f min | SZ NA %.3f | proj NA %.3f -> %s\n",
            nrow(R), as.numeric(difftime(Sys.time(),t0,units="mins")),
            mean(is.na(R$SZ)), mean(is.na(R$proj)), OUT))
