suppressMessages(library(brglm2))
## scn_grid.R -- GATE SCN: one block per (n, corpus) carrying ALL scenarios x severities.
##
## Pre-declared in speckit/25_SCENARIO_STUDY.md. Read it before changing anything here.
##
## S3 compared against a MIXTURE of departures and returned one number per cell. This gate resolves
## the mixture: nine single-kind scenarios, six severity levels each, so the reader sees WHICH
## departures the learned test wins on and AT WHAT SEVERITY -- the power-curve comparison a GOF
## paper is actually judged on.
##
## SEVERITY = sd(delta) on the LOGIT scale. The generators already parameterise exactly this
## (`delta <- delta*(amp/sd(delta))`), so severity is scenario-independent and interpretable, and
## s = 0 IS the null, giving every arm's own size as the left end of its power curve.
## Severity is NOT rescaled by n -- power rising with n down a column is the point.
##
## COMPARATORS: nine tests. `A.poly`/`A.spline` are EXCLUDED (25 section 2) -- they are the author's
## own unpublished EDGE/EDGES family, so benchmarking against them compares preprint to preprint.
## They are still COMPUTED and deposited (they cost nothing extra, eval_fast returns all eleven) so
## the exclusion is auditable, but no criterion or figure may use them.
##
## Config inherited from closed gates: K=6 (S4), deployment-matched mix=1.0, B=199, Firth policy,
## design D1 (p=2 textbook).

ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def"
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
source(file.path(ROOT,"jasa_study","config.R")); source(file.path(ROOT,"R","scenarios.R"))
source(file.path(JASA$study,"R","scenarios_registry.R")); source(file.path(JASA$study,"R","methods_registry.R"))
source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
G <- 10L
FEATS11 <- c("HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","EF",
             "DEF.poly2","DEF.poly3","DEF.stukel","A.poly","A.spline")
source(file.path(PAP,"code","families_shared.R"))
source(file.path(PAP,"code","fast_battery.R"))

K     <- 6L
B     <- as.integer(Sys.getenv("SC_B",    "199"))
NNULL <- as.integer(Sys.getenv("SC_NULL", "2000"))   # half trains, half measures size
NRAND <- as.integer(Sys.getenv("SC_RAND", "1000"))   # training alternative class (mixture)
NEVAL <- as.integer(Sys.getenv("SC_EVAL", "1000"))   # per scenario x severity
NW    <- as.integer(Sys.getenv("SC_W",    "12"))
NOBS  <- as.integer(Sys.getenv("SC_N",    "100"))
SEED0 <- as.integer(Sys.getenv("SC_SEED", "101000"))

## BLOCK A = localised in covariate space | BLOCK B = function of the index / smooth  (25 section 3)
SCN <- c("bump","step","inter","abs",  "indexq","indexs","quad","cubic","sine")
SEV <- c(0, 0.3, 0.6, 1.0, 1.5, 2.2)                  # sd(delta), logit scale
OUT <- file.path(PAP,"results","scn", sprintf("n%d_s%d", NOBS, SEED0))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

fit_pol <- function(d) {
  g <- tryCatch(suppressWarnings(glm(y~x1+x2, data=d, family=binomial())), error=function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1])>8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(y~x1+x2, data=d, family=binomial(),
                                       method="brglmFit", type="AS_mean")), error=function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(list(fit=f, firth=1L))
  }
  list(fit=g, firth=0L)
}
resid_map <- function(fit, d) {
  ph <- as.numeric(stats::fitted(fit)); r <- d$y - ph; n <- length(r)
  i1 <- pmin(K, 1L + floor(K*(rank(d$x1, ties.method="first")-1)/n))
  i2 <- pmin(K, 1L + floor(K*(rank(d$x2, ties.method="first")-1)/n))
  idx <- factor((i1-1L)*K + i2, levels=1:(K*K))
  s <- tapply(r, idx, sum); vv <- tapply(ph*(1-ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s)/sqrt(pmax(as.numeric(vv), 1e-8))
}
covs <- function(n) data.frame(x1=runif(n,-3,3), x2=rnorm(n))

## one scenario term, at a declared severity. sev = sd(delta) on the logit scale; sev 0 = the null.
gen_scn <- function(n, kind, sev, fseed) {
  d <- covs(n); eta0 <- 0.3 + 0.8*d$x1 - 0.5*d$x2
  if (sev <= 0) { d$y <- rbinom(n,1,plogis(eta0)); return(d) }
  set.seed(fseed)                              # the term's shape parameters, severity held separate
  xi <- if (runif(1) < .5) d$x1 else d$x2
  xj <- if (identical(xi, d$x1)) d$x2 else d$x1
  delta <- switch(kind,
    bump   = exp(-(xi - runif(1,-1.5,1.5))^2 / runif(1,.15,.6)),
    step   = as.numeric(xi > runif(1,-1,1)),
    inter  = xi*xj,
    abs    = abs(xi),
    indexq = (eta0-mean(eta0))^2/2,
    indexs = sin(runif(1,.8,2.0)*eta0),
    quad   = xi^2,
    cubic  = xi^3/3,
    sine   = sin(runif(1,.8,2.5)*xi + runif(1,0,pi)))
  sdd <- stats::sd(delta)
  if (sdd > 1e-8) delta <- delta*(sev/sdd)     # sd(delta) == sev exactly
  d$y <- rbinom(n, 1, plogis(eta0 + delta)); d
}
## the TRAINING alternative class stays the random MIXTURE (the shipped prior), never single kinds
gen_mix <- function(n, fam) {
  d <- covs(n); eta0 <- 0.3+0.8*d$x1-0.5*d$x2; X <- cbind(d$x1,d$x2); delta <- rep(0,n)
  for (tm in fam$terms) { xi <- X[,((tm$i-1L)%%2L)+1L]; xj <- X[,((tm$j-1L)%%2L)+1L]
    delta <- delta + tm$c*switch(tm$kind, quad=xi^2, cubic=xi^3/3, inter=xi*xj,
      sine=sin(tm$w*xi+tm$ph), bump=exp(-(xi-tm$mu)^2/tm$s), step=as.numeric(xi>tm$t),
      abs=abs(xi), indexq=(eta0-mean(eta0))^2/2, indexs=sin(tm$w*eta0)) }
  sdd <- sd(delta); if (sdd>1e-8) delta <- delta*(2*sqrt(100/n)*fam$amp/sdd)
  d$y <- rbinom(n,1,plogis(eta0+delta)); d }
tanf <- function(p) tan((0.5 - pmin(pmax(p,1e-12),1-1e-12))*pi)

## kind: 0 = null, 1 = training mixture, 2 = evaluation (carries scn + sev)
one <- function(kind, i, scn=NA_character_, sev=NA_real_, si=0L, vi=0L) {
  set.seed(SEED0 + kind*400000L + si*40000L + vi*4000L + i)
  d <- if (kind==0L) gen_scn(NOBS,"none",0,0L)
       else if (kind==1L) gen_mix(NOBS, make_family_p(SEED0*17L + i, 2L))
       else gen_scn(NOBS, scn, sev, SEED0*23L + si*7919L + i)
  if (length(unique(d$y)) < 2L) return(NULL)
  fp <- fit_pol(d); if (is.null(fp)) return(NULL)
  pv <- tryCatch(eval_fast(fp$fit), error=function(e) NULL); if (is.null(pv)) return(NULL)
  p11 <- pmin(pmax(as.numeric(pv[FEATS11]),1e-12),1-1e-12)
  M <- matrix(NA_real_, B+1L, K*K); M[1,] <- resid_map(fp$fit, d)
  ph <- as.numeric(stats::fitted(fp$fit)); d2 <- d
  for (b in seq_len(B)) {
    d2$y <- rbinom(NOBS,1,ph); fb <- fit_pol(d2)
    M[b+1L,] <- if (is.null(fb)) rep(0,K*K) else resid_map(fb$fit, d2)
  }
  list(kind=kind, rep=i, scn=scn, sev=sev, p11=p11, sep=as.numeric(fp$firth), maps=M)
}

SRC <- c(file.path(ROOT,"jasa_study","config.R"), file.path(ROOT,"R","scenarios.R"),
         file.path(PAP,"code","families_shared.R"), file.path(PAP,"code","fast_battery.R"),
         file.path(PAP,"code","universal_design.R"))
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("SRC","G","K","B","NOBS","SEED0","fit_pol","resid_map","covs",
                              "gen_scn","gen_mix","tanf","one","FEATS11"), envir=environment())
invisible(parallel::clusterEvalQ(cl, {
  source(SRC[1]); source(file.path(JASA$study,"R","methods_registry.R"))
  source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
  source(SRC[2]); source(file.path(JASA$study,"R","scenarios_registry.R"))
  G <- 10L
  FEATS11 <- c("HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","EF",
               "DEF.poly2","DEF.poly3","DEF.stukel","A.poly","A.spline")
  suppressMessages(library(brglm2))
  source(SRC[3]); source(SRC[4]); source(SRC[5]); NULL }))

jobs <- c(lapply(seq_len(NNULL), function(i) list(k=0L,i=i)),
          lapply(seq_len(NRAND), function(i) list(k=1L,i=i)))
for (si in seq_along(SCN)) for (vi in seq_along(SEV)) {
  if (SEV[vi] <= 0) next                               # severity 0 IS the null block, not duplicated
  jobs <- c(jobs, lapply(seq_len(NEVAL), function(i)
    list(k=2L,i=i,scn=SCN[si],sev=SEV[vi],si=si,vi=vi)))
}
cat(sprintf("SCN n=%d seed=%d | %d scenarios x %d severities | %d jobs | K=%d B=%d\n",
            NOBS, SEED0, length(SCN), sum(SEV>0), length(jobs), K, B)); flush.console()
t0 <- Sys.time()
res <- Filter(Negate(is.null), parallel::parLapplyLB(cl, jobs, function(j)
  one(j$k, j$i, if (is.null(j$scn)) NA_character_ else j$scn,
      if (is.null(j$sev)) NA_real_ else j$sev,
      if (is.null(j$si)) 0L else j$si, if (is.null(j$vi)) 0L else j$vi)))
parallel::stopCluster(cl)

idx <- data.frame(kind=vapply(res,function(z) z$kind,integer(1)),
                  rep =vapply(res,function(z) z$rep, integer(1)),
                  scn =vapply(res,function(z) z$scn, character(1)),
                  sev =vapply(res,function(z) z$sev, numeric(1)),
                  sep =vapply(res,function(z) z$sep, numeric(1)))
P <- do.call(rbind, lapply(res,function(z) z$p11)); colnames(P) <- FEATS11
write.csv(cbind(idx,P), file.path(OUT,"blk_pvalues.csv"), row.names=FALSE)
writeBin(as.numeric(t(do.call(rbind, lapply(res,function(z) z$maps)))),
         file.path(OUT,"blk_maps.bin"), size=4)
write.csv(data.frame(B=B,K=K,n=NOBS,seed=SEED0,ndatasets=length(res),nnull=NNULL,nrand=NRAND,
                     neval=NEVAL,ncell=K*K,scn=paste(SCN,collapse=","),sev=paste(SEV,collapse=",")),
          file.path(OUT,"blk_meta.csv"), row.names=FALSE)
cat(sprintf("SCN block COMPLETE: %d/%d datasets, %.1f min -> %s\n", length(res), length(jobs),
            as.numeric(difftime(Sys.time(),t0,units="mins")), OUT))
