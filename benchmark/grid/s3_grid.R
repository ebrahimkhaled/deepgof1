suppressMessages({ library(brglm2); library(ebrahim.gof) })
## s3_grid.R -- Simulation I: the 60-cell design grid, one (design, n, corpus seed) block per invocation.
##
## Fixed design choices: K = 6, B = 199 exact rank p-value, the bias-reduced fit on separation,
## departure amplitude 2 sqrt(100/n) times a random family amplitude.
##
## THE p > 2 PROJECTION: the grid is built on the two predictors with the largest |beta_j| sd(x_j) in
## the OBSERVED fit, by rank. The pair is chosen once per dataset and reused for every bootstrap
## replicate of that dataset. A factor's score is the largest |beta sd| among its dummies; if
## selected, its level index is rank-binned like any other axis. When p = 2 this reduces exactly to
## the two-covariate construction.
##
## Per dataset: the p-values of the three published partition tests (Hosmer-Lemeshow, HL equal-width,
## Pigeon-Heyse; ebrahim.gof run.all.gof, clipped to [1e-12, 1 - 1e-12] as in the deposited files) and
## the observed map with its B bootstrap maps (blk_maps.bin), which s3_score.py scores. None of the
## tests draws random numbers, so the bootstrap maps use the same random stream as the deposited run.
##
## Environment: S3_DESIGN (D1..D4), S3_N, S3_SEED, S3_B (199), S3_NULL (2000), S3_RAND (1000),
## S3_EVAL (2000), S3_W (2 workers), S3_OUT (default results/s3/<design>_n<n>_s<seed> under
## DEEPGOF_ARCHIVE). Smoke test: S3_DESIGN=D2 S3_N=50 S3_NULL=3 S3_RAND=1 S3_EVAL=2 S3_B=9 S3_W=1 S3_OUT=smoke_D2
##   Rscript benchmark/grid/s3_grid.R

ARC    <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
HERE_G <- local({ fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
                  if (length(fa)) dirname(normalizePath(fa[1])) else file.path(ARC, "benchmark", "grid") })
source(file.path(HERE_G, "s3_designs.R"))
PART   <- c("HL", "HL-equalwidth", "Pigeon-Heyse")
K      <- 6L
B      <- as.integer(Sys.getenv("S3_B",    "199"))
NNULL  <- as.integer(Sys.getenv("S3_NULL", "2000"))
NRAND  <- as.integer(Sys.getenv("S3_RAND", "1000"))
NEVAL  <- as.integer(Sys.getenv("S3_EVAL", "2000"))
NW     <- as.integer(Sys.getenv("S3_W",    "2"))
NOBS   <- as.integer(Sys.getenv("S3_N",    "100"))
SEED0  <- as.integer(Sys.getenv("S3_SEED", "101000"))
DES    <- Sys.getenv("S3_DESIGN", "D1")
AMPMULT <- 2*sqrt(100/NOBS)

## the four designs
CFG <- switch(DES,
  D1 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.50),   # textbook
  D2 = list(p=3L, law="indep", rho=0.5, cat=3L, rate=0.50),   # one 3-level factor
  D3 = list(p=3L, law="ar1",   rho=0.7, cat=0L, rate=0.50),   # correlated
  D4 = list(p=2L, law="indep", rho=0.5, cat=0L, rate=0.15))   # 15% events
OUT <- Sys.getenv("S3_OUT", file.path(ARC, "results", "s3", sprintf("%s_n%d_s%d", DES, NOBS, SEED0)))
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

## the two grid axes: largest |beta| * sd, scored per PREDICTOR (factor = max over its dummies),
## returned in data-column order
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
eval_part <- function(fit) {
  out <- setNames(rep(NA_real_, length(PART)), PART)
  r <- tryCatch(suppressWarnings(run.all.gof(fit, tests = PART, include_slow = FALSE, install = "no")),
                error = function(e) NULL)
  if (!is.null(r)) { pv <- setNames(r$p_value, r$Test); k <- intersect(PART, names(pv)); out[k] <- pv[k] }
  out
}

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
  pv <- tryCatch(eval_part(fit), error=function(e) NULL); if (is.null(pv)) return(NULL)
  p3 <- pmin(pmax(as.numeric(pv[PART]), 1e-12), 1-1e-12)
  ax <- tryCatch(pick_axes(fit, d), error=function(e) NULL); if (is.null(ax)) return(NULL)

  M <- matrix(NA_real_, B+1L, K*K)
  M[1,] <- resid_map(fit, d, ax)                 # axes FIXED for the whole block
  ph <- as.numeric(stats::fitted(fit)); d2 <- d
  for (b in seq_len(B)) {
    d2$y <- stats::rbinom(NOBS, 1, ph)
    fb <- fit_pol(frm, d2)
    M[b+1L,] <- if (is.null(fb)) rep(0, K*K) else resid_map(fb$fit, d2, ax)
  }
  list(kind=kind, rep=i, p3=p3, sep=as.numeric(fp$firth), ax=paste(ax, collapse="|"), maps=M)
}

cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("HERE_G","PART","K","B","NOBS","SEED0","CFG","AMPMULT","fit_pol",
                              "pick_axes","axis_num","resid_map","eval_part","one"), envir=environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages({ library(brglm2); library(ebrahim.gof) })
  source(file.path(HERE_G, "s3_designs.R")); NULL }))

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
                  axes=vapply(res,function(z) z$ax,  character(1)))
P <- do.call(rbind, lapply(res,function(z) z$p3)); colnames(P) <- PART
write.csv(cbind(idx,P), file.path(OUT,"blk_pvalues.csv"), row.names=FALSE)
writeBin(as.numeric(t(do.call(rbind, lapply(res,function(z) z$maps)))),
         file.path(OUT,"blk_maps.bin"), size=4)
write.csv(data.frame(B=B,K=K,n=NOBS,seed=SEED0,design=DES,ndatasets=length(res),
                     nnull=NNULL,nrand=NRAND,neval=NEVAL,ampmult=AMPMULT,ncell=K*K),
          file.path(OUT,"blk_meta.csv"), row.names=FALSE)
cat(sprintf("S3 block COMPLETE: %d/%d datasets, %.1f min -> %s\n",
            length(res), length(jobs), as.numeric(difftime(Sys.time(),t0,units="mins")), OUT))
