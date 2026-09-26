suppressMessages({ library(brglm2); library(ebrahim.gof) })
## liu_grid.R -- Simulation II: the four settings of Liu et al. (2024, Statistics and Computing 34:175),
## one (setting, n, corpus seed) block per invocation.
##
## Their four DGPs, verbatim from the paper's Section 4 (C = departure parameter; the fitted model
## omits the C-term):
##   S1: -1.963 + 0.982 x1 + C x1^2                      x1 ~ U[-3,3]
##   S2: -2 + 0.3 x1 + 0.3 x2 + 0.3 x3 + C x1^2          x1 ~ U[-3,3], x2 ~ N(0,1.5^2), x3 ~ t4
##   S3: -1.972 + 0.135 x1 + 1.117 x2 + C x1 x2          x1 ~ U[-3,3], x2 ~ Bin(.5)
##   S4: 0.3 x1 + 0.3 x2 + C x1 x2   (NO intercept)      x1, x2 ~ U[-3,3]
##
## Per dataset this writes the p-values of the classical arms computed here (Hosmer-Lemeshow, HL
## equal-width, Pigeon-Heyse, Osius-Rojek: ebrahim.gof run.all.gof; Stukel: the marginal SstBoth
## form of Liu et al., ebrahim.gof:::gof_stukel(form = "marginal")) and the observed residual map
## with its B bootstrap maps (blk_maps.bin), which liu_score.py turns into the DeepGOF-1 p-value.
## Stute-Zhu comes from liu_augment.R; mHL, GiViTI and the projection test from
## theory/addon_arms/ (liu_addon.R, proj_addon.R); BAGofT from bagoft/confirm/.
##
## Datasets are generated from a deterministic seed per replicate,
##   set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i),
## kind 0 = null, 1 = random training departure (used only to train the block's network), 2 = the
## Liu departure at strength C. None of the tests draws random numbers, so the bootstrap maps that
## follow use the same random stream as in the deposited run.
##
## p = 1 MAP: Setting 1 has one covariate, so it uses a one-dimensional partition of rank(x1) into
## K^2 = 36 quantile cells, reshaped row-major to 6 x 6.
##
## Environment: LU_SET (2), LU_N (100), LU_SEED (101000), LU_B (199), LU_NULL (2000), LU_RAND (1000),
## LU_EVAL (1000), LU_W (2 worker processes), LU_OUT (default results/liu/S<set>_n<n>_s<seed>
## under DEEPGOF_ARCHIVE, the archive root; run from the root or set it).
## Smoke test (seconds):  LU_SET=1 LU_N=50 LU_NULL=3 LU_RAND=1 LU_EVAL=1 LU_B=9 LU_W=1 LU_OUT=smoke_liu Rscript benchmark/liu/liu_grid.R

ARC   <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
PART  <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Osius-Rojek")
ALLT  <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "Stukel", "Osius-Rojek")
K     <- 6L
B     <- as.integer(Sys.getenv("LU_B",    "199"))
NNULL <- as.integer(Sys.getenv("LU_NULL", "2000"))
NRAND <- as.integer(Sys.getenv("LU_RAND", "1000"))
NEVAL <- as.integer(Sys.getenv("LU_EVAL", "1000"))
NW    <- as.integer(Sys.getenv("LU_W",    "2"))
NOBS  <- as.integer(Sys.getenv("LU_N",    "100"))
SEED0 <- as.integer(Sys.getenv("LU_SEED", "101000"))
SET   <- as.integer(Sys.getenv("LU_SET",  "2"))
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)
OUT   <- Sys.getenv("LU_OUT", file.path(ARC, "results", "liu", sprintf("S%d_n%d_s%d", SET, NOBS, SEED0)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- Liu's four DGPs, verbatim -----------------------------------------------------------------
gen_liu <- function(n, setting, C) {
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
    eta <- 0.3*x1 + 0.3*x2 + C*x1*x2                      # no intercept in the truth
  }
  d$y <- rbinom(n, 1, plogis(eta)); d
}
FRM <- function(d) stats::as.formula(paste("y ~", paste(setdiff(names(d),"y"), collapse=" + ")))

## the fit policy: glm, with the bias-reduced fit (brglm2 AS_mean) when a slope exceeds 8 in
## absolute value, a coefficient is not finite, or the iteration does not converge
fit_pol <- function(frm, d) {
  g <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial())), error=function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1])>8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial(),
                                       method="brglmFit", type="AS_mean")), error=function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(list(fit=f, firth=1L))
  }
  list(fit=g, firth=0L)
}
pick_axes <- function(fit, d) {
  vars <- setdiff(names(d), "y")
  if (length(vars) == 1L) return(vars)                   # p = 1 -> the 1-D construction
  mm <- stats::model.matrix(fit); b <- stats::coef(fit)
  sc <- vapply(vars, function(v) {
    cols <- which(startsWith(colnames(mm), v)); if (!length(cols)) return(0)
    nm <- colnames(mm)[cols]
    val <- abs(b[nm]) * apply(mm[, cols, drop=FALSE], 2, stats::sd)
    val <- val[is.finite(val)]; if (!length(val)) 0 else max(val)
  }, numeric(1))
  sc[!is.finite(sc)] <- 0
  sel <- vars[order(sc, decreasing=TRUE)][1:2]
  vars[sort(match(sel, vars))]
}
## ties among covariate values go by row order here; the simulated rows are in random order with
## respect to y, so this equals the random tie rule of R/deepgof1.R in law (theory/tie_fix/)
resid_map <- function(fit, d, ax) {
  ph <- as.numeric(stats::fitted(fit)); r <- d$y - ph; n <- length(r)
  if (length(ax) == 1L) {                                # p = 1: 36 quantile cells, row-major
    i <- pmin(K*K, 1L + floor(K*K*(rank(d[[ax]], ties.method="first")-1)/n))
    idx <- factor(i, levels = 1:(K*K))
  } else {
    bn <- function(v) pmin(K, 1L + floor(K*(rank(as.numeric(v), ties.method="first")-1)/n))
    idx <- factor((bn(d[[ax[1]]])-1L)*K + bn(d[[ax[2]]]), levels = 1:(K*K))
  }
  s <- tapply(r, idx, sum); vv <- tapply(ph*(1-ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s)/sqrt(pmax(as.numeric(vv), 1e-8))
}
## the classical arms of this block (no test draws random numbers)
eval_liu <- function(fit) {
  out <- setNames(rep(NA_real_, length(ALLT)), ALLT)
  r <- tryCatch(suppressWarnings(run.all.gof(fit, tests=PART, include_slow=FALSE, install="no")),
                error=function(e) NULL)
  if (!is.null(r)) { pv <- setNames(r$p_value, r$Test)
    k <- intersect(PART, names(pv)); out[k] <- pv[k] }
  out["Stukel"] <- tryCatch({
    ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
    v <- ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))$p_value
    if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_
  }, error = function(e) NA_real_)
  out
}
## training alternative class: random departures on THIS setting's covariate design
gen_train_alt <- function(n, setting, fseed) {
  set.seed(fseed)
  d <- gen_liu(n, setting, 0)                            # the design, null outcome placeholder
  X <- as.matrix(d[, setdiff(names(d),"y"), drop=FALSE])
  p <- ncol(X); j <- sample(p,1); k <- sample(p,1)
  kind <- sample(c("quad","inter","bump","step","sine","abs"),1)
  xi <- X[,j]; xj <- X[,k]
  delta <- switch(kind, quad=xi^2, inter=xi*xj, bump=exp(-(xi-runif(1,-1.5,1.5))^2/runif(1,.15,.6)),
                  step=as.numeric(xi>runif(1,-1,1)), sine=sin(runif(1,.8,2.5)*xi+runif(1,0,pi)),
                  abs=abs(xi))
  sdd <- stats::sd(delta); amp <- runif(1, .4, 2.0)*sqrt(100/n)
  if (sdd > 1e-8) delta <- delta*(amp/sdd)
  eta <- switch(as.character(setting),
    "1" = -1.963 + 0.982*X[,1], "2" = -2 + 0.3*X[,1] + 0.3*X[,2] + 0.3*X[,3],
    "3" = -1.972 + 0.135*X[,1] + 1.117*X[,2], "4" = 0.3*X[,1] + 0.3*X[,2])
  d$y <- rbinom(n, 1, plogis(eta + delta)); d
}

one <- function(kind, i, C=NA_real_, ci=0L) {
  set.seed(SEED0 + SET*500000L + kind*50000L + ci*5000L + i)
  d <- if (kind == 1L) gen_train_alt(NOBS, SET, SEED0*29L + i) else gen_liu(NOBS, SET, if (kind==0L) 0 else C)
  if (length(unique(d$y)) < 2L) return(NULL)
  frm <- FRM(d); fp <- fit_pol(frm, d); if (is.null(fp)) return(NULL)
  pv <- eval_liu(fp$fit); if (all(is.na(pv))) return(NULL)
  ax <- tryCatch(pick_axes(fp$fit, d), error=function(e) NULL); if (is.null(ax)) return(NULL)
  M <- matrix(NA_real_, B+1L, K*K); M[1,] <- resid_map(fp$fit, d, ax)
  ph <- as.numeric(stats::fitted(fp$fit)); d2 <- d
  for (b in seq_len(B)) {
    d2$y <- rbinom(NOBS, 1, ph); fb <- fit_pol(frm, d2)
    M[b+1L,] <- if (is.null(fb)) rep(0, K*K) else resid_map(fb$fit, d2, ax)
  }
  list(kind=kind, rep=i, C=C, sep=as.numeric(fp$firth), p=pv, maps=M)
}

cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("K","B","NOBS","SEED0","SET","PART","ALLT","gen_liu","FRM","fit_pol",
                              "pick_axes","resid_map","eval_liu","gen_train_alt","one"), envir=environment())
invisible(parallel::clusterEvalQ(cl, { suppressMessages({ library(brglm2); library(ebrahim.gof) }); NULL }))

jobs <- c(lapply(seq_len(NNULL), function(i) list(k=0L,i=i)),
          lapply(seq_len(NRAND), function(i) list(k=1L,i=i)))
for (ci in seq_along(CS)) if (CS[ci] > 0)
  jobs <- c(jobs, lapply(seq_len(NEVAL), function(i) list(k=2L,i=i,C=CS[ci],ci=ci)))
cat(sprintf("LIU Setting %d | n=%d seed=%d | C in {%s} | %d jobs | K=%d B=%d\n",
            SET, NOBS, SEED0, paste(CS,collapse=","), length(jobs), K, B)); flush.console()
t0 <- Sys.time()
res <- Filter(Negate(is.null), parallel::parLapplyLB(cl, jobs, function(j)
  one(j$k, j$i, if (is.null(j$C)) NA_real_ else j$C, if (is.null(j$ci)) 0L else j$ci)))
parallel::stopCluster(cl)

idx <- data.frame(kind=vapply(res,function(z) z$kind,integer(1)),
                  rep =vapply(res,function(z) z$rep, integer(1)),
                  C   =vapply(res,function(z) z$C,   numeric(1)),
                  sep =vapply(res,function(z) z$sep, numeric(1)))
P <- do.call(rbind, lapply(res, function(z) z$p)); colnames(P) <- ALLT
write.csv(cbind(idx,P), file.path(OUT,"blk_pvalues.csv"), row.names=FALSE)
writeBin(as.numeric(t(do.call(rbind, lapply(res,function(z) z$maps)))),
         file.path(OUT,"blk_maps.bin"), size=4)
write.csv(data.frame(B=B,K=K,n=NOBS,seed=SEED0,setting=SET,ndatasets=length(res),
                     nnull=NNULL,nrand=NRAND,neval=NEVAL,ncell=K*K,Cs=paste(CS,collapse=",")),
          file.path(OUT,"blk_meta.csv"), row.names=FALSE)
cat(sprintf("LIU block COMPLETE: %d/%d datasets, %.1f min -> %s\n", length(res), length(jobs),
            as.numeric(difftime(Sys.time(),t0,units="mins")), OUT))
