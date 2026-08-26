suppressMessages(library(brglm2))
## s3_grid.R -- GATE S3: the confirmatory grid, one (design, n, seed) block per invocation.
##
## Pre-declared in speckit/21_S3_PREDECLARATION.md. Read it before changing anything here.
##
## Configuration inherited from CLOSED gates and NOT re-opened: K=6 (gate S4), deployment-matched
## training with mix pinned at 1.0, B=199 exact rank p-value, Firth on separation, amplitude
## 2*sqrt(100/n) (gate S4b's regime).
##
## THE p>2 PROJECTION (first-time construction, declared in 21 section 2): grid on the two
## predictors with the largest |beta_j| * sd(x_j) in the OBSERVED fit, by rank. The chosen pair is
## computed ONCE per dataset and reused for every bootstrap replicate of that dataset -- letting it
## move between the observed row and its replicates would break the exchangeability Proposition 1
## rests on. A factor's score is the largest |beta*sd| among its dummies; if selected, its level
## index is rank-binned like any other axis, which populates exactly nlev distinct bands.
## When p = 2 this reduces EXACTLY to the construction used in every earlier gate.

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
source(file.path(PAP,"code","universal_design.R"))

K      <- 6L                                            # gate S4's selection
B      <- as.integer(Sys.getenv("S3_B",    "199"))
NNULL  <- as.integer(Sys.getenv("S3_NULL", "2000"))
NRAND  <- as.integer(Sys.getenv("S3_RAND", "1000"))
NEVAL  <- as.integer(Sys.getenv("S3_EVAL", "2000"))
NW     <- as.integer(Sys.getenv("S3_W",    "12"))
NOBS   <- as.integer(Sys.getenv("S3_N",    "100"))
SEED0  <- as.integer(Sys.getenv("S3_SEED", "101000"))
DES    <- Sys.getenv("S3_DESIGN", "D1")
AMPMULT <- 2*sqrt(100/NOBS)                             # gate S4b's regime, not tunable per cell

## the four designs of 21 section 4
CFG <- switch(DES,
  D1 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.50),   # textbook
  D2 = list(p=3L, law="indep", rho=0.5, cat=3L, rate=0.50),   # one 3-level factor
  D3 = list(p=3L, law="ar1",   rho=0.7, cat=0L, rate=0.50),   # correlated
  D4 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.15))   # 15% events
OUT <- file.path(PAP,"results","s3", sprintf("%s_n%d_s%d", DES, NOBS, SEED0))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

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

## the two grid axes: largest |beta| * sd, scored per PREDICTOR (factor = max over its dummies)
pick_axes <- function(fit, d) {
  vars <- setdiff(names(d), "y")
  mm <- stats::model.matrix(fit); b <- stats::coef(fit)
  sc <- vapply(vars, function(v) {
    cols <- which(startsWith(colnames(mm), v))
    if (!length(cols)) return(0)
    nm <- colnames(mm)[cols]
    vals <- abs(b[nm]) * apply(mm[, cols, drop=FALSE], 2, stats::sd)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) 0 else max(vals)
  }, numeric(1))
  sc[!is.finite(sc)] <- 0
  sel <- vars[order(sc, decreasing=TRUE)][1:2]
  ## return in DATA-COLUMN order, not score order. Score order would transpose the grid whenever
  ## the second predictor outscored the first, and every closed gate used a fixed (x1,x2) order --
  ## so score order would silently make S3 incomparable with them. Data order keeps the p=2 case
  ## EXACTLY the earlier construction, which is what 21 section 2 claims.
  vars[sort(match(sel, vars))]
}
axis_num <- function(x) if (is.factor(x)) as.numeric(x) else as.numeric(x)

resid_map <- function(fit, d, ax) {
  ph <- as.numeric(stats::fitted(fit)); r <- d$y - ph; n <- length(r)
  bin <- function(v) pmin(K, 1L + floor(K*(rank(axis_num(v), ties.method="first")-1)/n))
  i1 <- bin(d[[ax[1]]]); i2 <- bin(d[[ax[2]]])
  idx <- factor((i1-1L)*K + i2, levels=1:(K*K))
  s <- tapply(r, idx, sum); vv <- tapply(ph*(1-ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  as.numeric(s)/sqrt(pmax(as.numeric(vv), 1e-8))
}
tanf <- function(p) tan((0.5 - pmin(pmax(p,1e-12),1-1e-12))*pi)

one <- function(kind, i) {
  set.seed(SEED0 + kind*300000L + i)
  g <- tryCatch({
    if (kind == 0L) gen_design_null(NOBS, CFG$p, CFG$law, CFG$rho, CFG$cat, CFG$rate)
    else {
      fs <- if (kind == 1L) SEED0*17L + i else SEED0*17L + 700000L + i   # DISJOINT eval stream
      fam <- make_family_p(fs, CFG$p); fam$amp <- fam$amp * AMPMULT
      gen_design_alt(fam, NOBS, CFG$p, CFG$law, CFG$rho, CFG$cat, CFG$rate)
    }}, error=function(e) NULL)
  if (is.null(g) || length(unique(g$data$y)) < 2L) return(NULL)
  d <- g$data; frm <- stats::formula(g$fit)
  fp <- fit_pol(frm, d); if (is.null(fp)) return(NULL)
  fit <- fp$fit
  pv <- tryCatch(eval_fast(fit), error=function(e) NULL); if (is.null(pv)) return(NULL)
  p11 <- pmin(pmax(as.numeric(pv[FEATS11]), 1e-12), 1-1e-12)
  ax <- tryCatch(pick_axes(fit, d), error=function(e) NULL); if (is.null(ax)) return(NULL)

  M <- matrix(NA_real_, B+1L, K*K)
  M[1,] <- resid_map(fit, d, ax)                 # axes FIXED for the whole block
  ph <- as.numeric(stats::fitted(fit)); d2 <- d
  for (b in seq_len(B)) {
    d2$y <- stats::rbinom(NOBS, 1, ph)
    fb <- fit_pol(frm, d2)
    M[b+1L,] <- if (is.null(fb)) rep(0, K*K) else resid_map(fb$fit, d2, ax)
  }
  list(kind=kind, rep=i, p11=p11, pcct=0.5-atan(mean(tanf(p11)))/pi, sep=as.numeric(fp$firth),
       ax=paste(ax, collapse="|"), maps=M)
}

SRC <- c(file.path(ROOT,"jasa_study","config.R"), file.path(ROOT,"R","scenarios.R"),
         file.path(PAP,"code","families_shared.R"), file.path(PAP,"code","fast_battery.R"),
         file.path(PAP,"code","universal_design.R"))
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("SRC","G","K","B","NOBS","SEED0","CFG","AMPMULT","fit_pol",
                              "pick_axes","axis_num","resid_map","tanf","one","FEATS11"),
                        envir=environment())
invisible(parallel::clusterEvalQ(cl, {
  source(SRC[1]); source(file.path(JASA$study,"R","methods_registry.R"))
  source(file.path(JASA$study,"R","harness_core.R")); jasa_init()
  source(SRC[2]); source(file.path(JASA$study,"R","scenarios_registry.R"))
  G <- 10L
  FEATS11 <- c("HL","HL-equalwidth","Pigeon-Heyse","Tsiatis","Xie","EF",
               "DEF.poly2","DEF.poly3","DEF.stukel","A.poly","A.spline")
  suppressMessages(library(brglm2))
  source(SRC[3]); source(SRC[4]); source(SRC[5]); NULL }))

cat(sprintf("S3 %s n=%d seed=%d | p=%d law=%s rho=%.1f cat=%d rate=%.2f | amp x%.2f K=%d B=%d\n",
            DES, NOBS, SEED0, CFG$p, CFG$law, CFG$rho, CFG$cat, CFG$rate, AMPMULT, K, B))
flush.console(); t0 <- Sys.time()
jobs <- c(lapply(seq_len(NNULL), function(i) list(k=0L,i=i)),
          lapply(seq_len(NRAND), function(i) list(k=1L,i=i)),
          lapply(seq_len(NEVAL), function(i) list(k=2L,i=i)))
res <- Filter(Negate(is.null), parallel::parLapplyLB(cl, jobs, function(j) one(j$k, j$i)))
parallel::stopCluster(cl)

idx <- data.frame(kind=vapply(res,function(z) z$kind,integer(1)),
                  rep =vapply(res,function(z) z$rep, integer(1)),
                  sep =vapply(res,function(z) z$sep, numeric(1)),
                  cct =vapply(res,function(z) z$pcct,numeric(1)),
                  axes=vapply(res,function(z) z$ax,  character(1)))
P <- do.call(rbind, lapply(res,function(z) z$p11)); colnames(P) <- FEATS11
write.csv(cbind(idx,P), file.path(OUT,"blk_pvalues.csv"), row.names=FALSE)
writeBin(as.numeric(t(do.call(rbind, lapply(res,function(z) z$maps)))),
         file.path(OUT,"blk_maps.bin"), size=4)
write.csv(data.frame(B=B,K=K,n=NOBS,seed=SEED0,design=DES,ndatasets=length(res),
                     nnull=NNULL,nrand=NRAND,neval=NEVAL,ampmult=AMPMULT,ncell=K*K),
          file.path(OUT,"blk_meta.csv"), row.names=FALSE)
cat(sprintf("S3 block COMPLETE: %d/%d datasets, %.1f min -> %s\n",
            length(res), length(jobs), as.numeric(difftime(Sys.time(),t0,units="mins")), OUT))
