suppressMessages({ library(brglm2); library(ebrahim.gof) })
## liu_augment.R -- adds the Stute-Zhu test to an existing Liu block.
##
## Every dataset of a block is generated from a deterministic seed
## (set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i)), so the same data are re-derived here
## without touching the stored block. Stute-Zhu is ebrahim.gof's implementation (run.all.gof row
## "Stute-Zhu"; it calibrates by its own bootstrap and so draws random numbers after the fit).
##
## The deposited augment_pvalues.csv also carries `proj_superseded`: the projection statistic of
## the first run, whose weight had an extra length factor |X_i - X_l||X_j - X_l| on standardised
## covariates. It is kept only as the pairing key for the corrected projection test of
## theory/addon_arms/proj_addon.R (Liu et al.'s angle-only weight on the raw covariates), which
## replays the random-number state Stute-Zhu leaves and so uses the same bootstrap outcomes.
##
## Environment: AU_SET, AU_N, AU_SEED (the block), AU_W (2 workers), AU_DIR (default
## results/liu/S<set>_n<n>_s<seed> under DEEPGOF_ARCHIVE). Writes <block>/augment_sz.csv.

ARC   <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
NW    <- as.integer(Sys.getenv("AU_W",    "2"))
NOBS  <- as.integer(Sys.getenv("AU_N",    "100"))
SEED0 <- as.integer(Sys.getenv("AU_SEED", "101000"))
SET   <- as.integer(Sys.getenv("AU_SET",  "2"))
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)
OUT   <- Sys.getenv("AU_DIR", file.path(ARC, "results", "liu", sprintf("S%d_n%d_s%d", SET, NOBS, SEED0)))
if (!dir.exists(OUT)) stop("no such block: ", OUT)
META  <- read.csv(file.path(OUT, "blk_meta.csv"))
NNULL <- as.integer(Sys.getenv("AU_NULL", META$nnull[1]))
NRAND <- as.integer(Sys.getenv("AU_RAND", META$nrand[1]))
NEVAL <- as.integer(Sys.getenv("AU_EVAL", META$neval[1]))

## ---- the four DGPs and the fit policy, identical to liu_grid.R -----------------------------------
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

one <- function(kind, i, C=NA_real_, ci=0L) {
  set.seed(SEED0 + SET*500000L + kind*50000L + ci*5000L + i)
  d <- if (kind == 1L) gen_train_alt(NOBS, SET, SEED0*29L + i)
       else gen_liu(NOBS, SET, if (kind==0L) 0 else C)
  if (length(unique(d$y)) < 2L) return(NULL)
  frm <- FRM(d); fit <- fit_pol(frm, d); if (is.null(fit)) return(NULL)
  sz <- tryCatch(suppressWarnings(
          run.all.gof(fit, tests="Stute-Zhu", include_slow=TRUE, install="no")$p_value[1]),
        error=function(e) NA_real_)
  c(kind=kind, rep=i, C=C, SZ=sz)
}

cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("NOBS","SEED0","SET","gen_liu","gen_train_alt","FRM","fit_pol","one"),
                        envir=environment())
invisible(parallel::clusterEvalQ(cl, { suppressMessages({ library(brglm2); library(ebrahim.gof) }); NULL }))
jobs <- c(lapply(seq_len(NNULL), function(i) list(k=0L,i=i)),
          lapply(seq_len(NRAND), function(i) list(k=1L,i=i)))
for (ci in seq_along(CS)) if (CS[ci] > 0)
  jobs <- c(jobs, lapply(seq_len(NEVAL), function(i) list(k=2L,i=i,C=CS[ci],ci=ci)))
cat(sprintf("AUGMENT S%d n=%d seed=%d | Stute-Zhu | %d datasets\n", SET, NOBS, SEED0, length(jobs)))
t0 <- Sys.time()
res <- Filter(Negate(is.null), parallel::parLapplyLB(cl, jobs, function(j)
  one(j$k, j$i, if (is.null(j$C)) NA_real_ else j$C, if (is.null(j$ci)) 0L else j$ci)))
parallel::stopCluster(cl)
R <- as.data.frame(do.call(rbind, res))
write.csv(R, file.path(OUT, "augment_sz.csv"), row.names=FALSE)
cat(sprintf("AUGMENT COMPLETE: %d rows, %.1f min | SZ NA %.3f -> %s\n",
            nrow(R), as.numeric(difftime(Sys.time(),t0,units="mins")), mean(is.na(R$SZ)), OUT))
