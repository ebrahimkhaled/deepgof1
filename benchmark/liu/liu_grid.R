suppressMessages(library(brglm2))
## liu_grid.R -- GATE LIU: replicate Liu et al. (2024, Stat & Comput 34:175) settings exactly.
##
## Pre-declared in speckit/26_LIU_BENCHMARK.md. Read it before changing anything here.
##
## Their four DGPs, verbatim from the paper's section 4 (C = departure parameter, fitted model omits
## the C-term):
##   S1: -1.963 + 0.982 x1 + C x1^2                      x1 ~ U[-3,3]
##   S2: -2 + 0.3 x1 + 0.3 x2 + 0.3 x3 + C x1^2          x1 ~ U[-3,3], x2 ~ N(0,1.5^2), x3 ~ t4
##   S3: -1.972 + 0.135 x1 + 1.117 x2 + C x1 x2          x1 ~ U[-3,3], x2 ~ Bin(.5)
##   S4: 0.3 x1 + 0.3 x2 + C x1 x2   (NO intercept)      x1, x2 ~ U[-3,3]
##
## ⭐ S2 IS THE DECISIVE CELL. Liu et al: all nine classical tests fail there -- "their power
## functions are not even increasing with C" -- because C*x1^2 is nearly orthogonal to the fitted
## index (x1 enters with coefficient 0.3). Only their projection test gains power. That is exactly
## the covariate-space corner the residual map was built for, so LIU2 is this project's most direct
## external test of its own central mechanism claim.
##
## COMPARATORS: eight EXTERNAL published tests (the five partition ones + Stukel + Osius-Rojek +
## le-Cessie). The author's own EF/DEF/EDGE family is computed and deposited but is NEVER the
## comparator (26 section 2). `proj` and Stute-Zhu are ABSENT and that is a stated limitation.
##
## p=1 MAP (26 section 3, declared before use): Setting 1 has one covariate, so a two-axis grid
## degenerates. It uses a ONE-dimensional partition of rank(x1) into K^2 = 36 quantile cells,
## reshaped row-major to 6x6 -- a finer 1-D partition, reported as such.

ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def"
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
source(file.path(ROOT,"jasa_study","config.R")); source(file.path(ROOT,"R","scenarios.R"))
source(file.path(JASA$study,"R","scenarios_registry.R")); source(file.path(JASA$study,"R","methods_registry.R"))
source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
G <- 10L
source(file.path(PAP,"code","families_shared.R"))

EXT  <- c("HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","Stukel","Osius-Rojek","le-Cessie")
OWN  <- c("EF","DEF.poly2","DEF.poly3","DEF.stukel","A.poly","A.spline")
ALLT <- c(EXT, OWN)
BATT <- c("HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","EF",
          "DEF.poly2","DEF.poly3","DEF.stukel","Stukel","Osius-Rojek")

K     <- 6L
B     <- as.integer(Sys.getenv("LU_B",    "199"))
NNULL <- as.integer(Sys.getenv("LU_NULL", "2000"))
NRAND <- as.integer(Sys.getenv("LU_RAND", "1000"))
NEVAL <- as.integer(Sys.getenv("LU_EVAL", "1000"))
NW    <- as.integer(Sys.getenv("LU_W",    "12"))
NOBS  <- as.integer(Sys.getenv("LU_N",    "100"))
SEED0 <- as.integer(Sys.getenv("LU_SEED", "101000"))
SET   <- as.integer(Sys.getenv("LU_SET",  "2"))
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)
OUT   <- file.path(PAP,"results","liu", sprintf("S%d_n%d_s%d", SET, NOBS, SEED0))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

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
## the eight external tests + the author's family, from one battery call + le-Cessie + cdef
eval_liu <- function(fit) {
  out <- setNames(rep(NA_real_, length(ALLT)), ALLT)
  r <- tryCatch(suppressWarnings(run.all.gof(fit, tests=BATT, include_slow=FALSE, install="no")),
                error=function(e) NULL)
  if (!is.null(r)) { pv <- setNames(r$p_value, r$Test)
    k <- intersect(BATT, names(pv)); out[k] <- pv[k] }
  ## lc_fast(y, ph, X, covs) -- NOT a fit. Build its arguments from the fitted model.
  out["le-Cessie"] <- tryCatch({
    y <- as.numeric(fit$y); ph <- as.numeric(stats::fitted(fit))
    X <- stats::model.matrix(fit)
    cv <- lapply(seq_len(ncol(X))[-1], function(j) X[, j])      # covariates, intercept dropped
    unname(lc_fast(y, ph, X, cv)["p"])
  }, error = function(e) NA_real_)
  out["A.poly"]    <- tryCatch(suppressWarnings(cdef.gof(fit, basis="poly")$p_value),   error=function(e) NA_real_)
  out["A.spline"]  <- tryCatch(suppressWarnings(cdef.gof(fit, basis="spline")$p_value), error=function(e) NA_real_)
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

SRC <- c(file.path(ROOT,"jasa_study","config.R"), file.path(ROOT,"R","scenarios.R"),
         file.path(PAP,"code","families_shared.R"), file.path(PAP,"code","lecessie_algebraic.R"))
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("SRC","G","K","B","NOBS","SEED0","SET","EXT","OWN","ALLT","BATT",
                              "gen_liu","FRM","fit_pol","pick_axes","resid_map","eval_liu",
                              "gen_train_alt","one"), envir=environment())
invisible(parallel::clusterEvalQ(cl, {
  source(SRC[1]); source(file.path(JASA$study,"R","methods_registry.R"))
  source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
  source(SRC[2]); source(file.path(JASA$study,"R","scenarios_registry.R"))
  G <- 10L; suppressMessages(library(brglm2))
  source(SRC[3]); source(SRC[4]); NULL }))

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
