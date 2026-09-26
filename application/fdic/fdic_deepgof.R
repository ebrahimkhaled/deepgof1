## fdic_deepgof.R -- DeepGOF-1 on the FDIC bank-failure cells. Two products, both deposited:
##
##   1. ADEQUACY: deepgof1() on each state cell's fitted model (the column "DG-1 observed p" of the
##      paper's FDIC table), next to the classical tests' observed p-values in observed_pvalues.csv.
##   2. SIZE ON THE REAL DESIGN: simulate REPS null outcomes from the fitted model on the real
##      covariate matrix, run the deployed test on each, and read the rejection rate at .05.
##
## Uses the archived deepgof1() (R/deepgof1.R), i.e. the shipped weights. Seeds are deterministic
## per (state, replicate).
##
## TIES. The adequacy calls use the fixed rule of ebrahim.gof 2.8.0 and of R/deepgof1.R by default
## (FDG_TIES=random): observed p-values .350, .175, .785, .655, .330, .905 for GA, FL, IL, CA, WA, MN,
## the values of the paper's FDIC table. FDG_TIES=roworder gives the values computed before the
## tie-breaking fix (.320, .185, .805, .710, .350, .905); the difference is bootstrap Monte Carlo
## movement only (no cell rejects under either rule; theory/tie_fix/ reports 30 random row orders
## per cell). fdic_tie_check.R checks both against the deposit. The covariates eq_ta and roa have
## no ties; npa_ta, bro_dep and cre_ta are tied only at zero.
##
## The deposited size column (deepgof_size.csv, deepgof_pvalues_<S>.csv; the "DG-1 size" column of
## the paper) was computed with ties broken by row order, so the size part keeps that rule by default
## (FDG_SIZE_TIES=roworder) and reproduces the deposit exactly. FDG_SIZE_TIES=random reruns it under
## the fixed rule.
##
## Run from the archive root (or set DEEPGOF_ARCHIVE). Writes results/fdic/deepgof_observed.csv,
## results/fdic/deepgof_size.csv and results/fdic/deepgof_pvalues_<S>.csv. FDG_REPS=0 runs the
## adequacy calls only (seconds).

ARC  <- normalizePath(Sys.getenv("DEEPGOF_ARCHIVE", "."))
REL  <- file.path(ARC, "R", "deepgof1.R")
REPS <- as.integer(Sys.getenv("FDG_REPS", "2000"))
NW   <- as.integer(Sys.getenv("FDG_W",    "2"))
TIES <- Sys.getenv("FDG_TIES", "random")
SIZE_TIES <- Sys.getenv("FDG_SIZE_TIES", "roworder")
SEED0<- 771000L
OUTD <- file.path(ARC, "results", "fdic"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)

source(REL)

d  <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_2008_banks.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_state_counts.csv"), stringsAsFactors = FALSE)
CELLS <- c("GA", setdiff(st$state[st$failures >= 15], "GA"))
FRM <- failed ~ eq_ta + npa_ta + roa + bro_dep + cre_ta

## ---- 1. the adequacy calls (seconds) -----------------------------------------------------------
obs <- list()
for (S in CELLS) {
  ds  <- d[d$state == S, ]
  fit <- suppressWarnings(glm(FRM, data = ds, family = binomial()))
  set.seed(SEED0 + match(S, CELLS))
  r <- deepgof1(fit, ties = TIES)
  obs[[S]] <- data.frame(state = S, n = nrow(ds), events = sum(ds$failed),
                         S_obs = r$statistic, p_deepgof = r$p.value,
                         axis1 = r$axes[1], axis2 = r$axes[2])
  cat(sprintf("%s: n=%3d events=%3d  DeepGOF-1 p = %.4f  (axes %s, %s; ties %s)\n",
              S, nrow(ds), sum(ds$failed), r$p.value, r$axes[1], r$axes[2], TIES)); flush.console()
}
OB <- do.call(rbind, obs)
write.csv(OB, file.path(OUTD, "deepgof_observed.csv"), row.names = FALSE)
if (REPS < 1L) quit(save = "no")

## ---- 2. the size on the real design (the long part; wall times of the deposited run in DEEPGOF_RUN.log)
cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("REL", "FRM", "SEED0", "SIZE_TIES", "ARC"), envir = environment())
invisible(parallel::clusterEvalQ(cl, { source(REL); NULL }))

SZ <- list()
for (S in CELLS) {
  ds  <- d[d$state == S, ]
  fit <- suppressWarnings(glm(FRM, data = ds, family = binomial()))
  ph  <- as.numeric(fitted(fit))
  parallel::clusterExport(cl, c("ds", "ph", "S"), envir = environment())
  t0 <- Sys.time()
  pv <- unlist(parallel::parLapplyLB(cl, seq_len(REPS), function(i) {
    set.seed(SEED0 + 1000L * which(LETTERS %in% substr(S, 1, 1)) + i)
    d2 <- ds; d2$failed <- rbinom(nrow(d2), 1L, ph)
    f2 <- tryCatch(suppressWarnings(glm(FRM, data = d2, family = binomial())),
                   error = function(e) NULL)
    if (is.null(f2)) return(NA_real_)
    tryCatch(deepgof1(f2, ties = SIZE_TIES)$p.value, error = function(e) NA_real_)
  }))
  write.csv(data.frame(rep = seq_len(REPS), p = pv),
            file.path(OUTD, paste0("deepgof_pvalues_", S, ".csv")), row.names = FALSE)
  sz <- mean(ifelse(is.na(pv), 1, pv) <= 0.05)      # a failed test cannot reject
  SZ[[S]] <- data.frame(state = S, n = nrow(ds), events = sum(ds$failed),
                        size = sz, na_rate = mean(is.na(pv)),
                        mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(sprintf("%s: DeepGOF-1 size on the real design = %.4f  (NA %.3f, %.1f min; ties %s)\n",
              S, sz, mean(is.na(pv)), SZ[[S]]$mins, SIZE_TIES)); flush.console()
}
parallel::stopCluster(cl)
write.csv(do.call(rbind, SZ), file.path(OUTD, "deepgof_size.csv"), row.names = FALSE)
cat("FDIC DEEPGOF COMPLETE\n")
