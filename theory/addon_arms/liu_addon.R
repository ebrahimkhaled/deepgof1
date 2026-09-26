# =============================================================================
# liu_addon.R -- mHL and GiViTI on the Liu et al. (2024) grid, on EXACTLY the datasets of the
# stored blocks in release/deepgof1/results/liu/ (24 blocks: S1-S4 x n 50/100/500 x corpus
# 101000/102000).
#
# Same data. Every dataset of a block comes from its own seed (liu_grid.R, one()):
#   set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i)
# and is generated and fitted here by gen_liu(), gen_train_alt() and fit_pol() copied VERBATIM
# from release/deepgof1/benchmark/liu/liu_grid.R -- the same mechanism liu_augment.R used to add
# Stute-Zhu and proj. Nothing after the fit uses the RNG, so the arms cannot shift the data.
#
# Rows. One row per dataset, keyed (setting, n, seed0, kind, rep, C) exactly as blk_pvalues.csv
# and augment_pvalues.csv, plus ci and the dataset seed. HL_check is the HL p-value recomputed
# here: it must equal the stored HL column (verify_pairing.R checks it).
#   kind 0 = null (C = 0), rep 1-2000; kind 2 = departure C in {.25,.5,.75,1}, rep 1-1000.
#   kind 1 (the 1,000 random departures used only to TRAIN the map network) is skipped unless
#   AD_KIND1 = 1: no reported size or power uses it.
# A dataset liu_grid.R dropped (degenerate y, failed fit) is kept here with status set, so a join
# on the key simply leaves it out.
#
# Resumable. Work is cut into chunks of AD_CHUNK datasets; each chunk is written to
# liu_rows/<block>/chunk_<kind>_<ci>_<first rep>.csv (.part, then rename) as soon as it is done.
# When every chunk of a block is on disk they are bound into liu_rows/<block>.csv. A block whose
# final CSV exists is skipped; a chunk whose file exists is skipped.
#
# Run:   Rscript liu_addon.R
# Env:   AD_W       workers (4)
#        AD_BLOCKS  comma list of block names, e.g. S2_n100_s101000 (default: all 24)
#        AD_GIVITI  1 = compute GiViTI (default 1); 0 = HL/Stukel/mHL only
#        AD_KIND1   0 (default) / 1
#        AD_CHUNK   datasets per chunk (50)
#        AD_MAXREP  cap on rep per kind, for smoke tests (default: none)
#        AD_OUT     output folder (default liu_rows; smoke tests use liu_rows_smoke)
# =============================================================================
HERE  <- "/path/to/project/paper_deepgof/theory/addon_arms"
RES   <- "/path/to/project/paper_deepgof/release/deepgof1/results/liu"
NW    <- as.integer(Sys.getenv("AD_W", "4"))
DOGV  <- Sys.getenv("AD_GIVITI", "1") == "1"
KIND1 <- Sys.getenv("AD_KIND1", "0") == "1"
CHUNK <- as.integer(Sys.getenv("AD_CHUNK", "50"))
MAXR  <- as.integer(Sys.getenv("AD_MAXREP", "100000"))
OUTD  <- file.path(HERE, Sys.getenv("AD_OUT", "liu_rows")); dir.create(OUTD, showWarnings = FALSE)
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)
source(file.path(HERE, "addon_common.R"))            # write_atomic() for the bind step

BLOCKS <- if (nzchar(Sys.getenv("AD_BLOCKS"))) strsplit(Sys.getenv("AD_BLOCKS"), ",")[[1]] else
  sort(basename(list.dirs(RES, recursive = FALSE)))

liu_defs <- quote({
  ## ---- VERBATIM from benchmark/liu/liu_grid.R (lines 55-84) ----
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
  ## ---- VERBATIM from liu_grid.R (lines 131-147) ----
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
  ## ---- one dataset: the data of liu_grid.R::one(), then the add-on arms ----
  liu_one <- function(SET, NOBS, SEED0, kind, i, C, ci, dogv) {
    s <- SEED0 + SET*500000L + kind*50000L + ci*5000L + i
    set.seed(s)
    d <- if (kind == 1L) gen_train_alt(NOBS, SET, SEED0*29L + i) else gen_liu(NOBS, SET, if (kind==0L) 0 else C)
    key <- data.frame(setting = SET, n = NOBS, seed0 = SEED0, kind = kind, rep = i, C = C, ci = ci,
                      dseed = s, events = sum(d$y))
    na <- c(HL_check = NA, Stukel_joint = NA, Stukel_marginal = NA, mHL = NA, GiViTI = NA,
            giviti_status = NA, giviti_secs = NA)
    if (length(unique(d$y)) < 2L) return(cbind(key, sep = NA, t(na), status = "dropped_degenerate_y"))
    fp <- fit_pol(FRM(d), d)
    if (is.null(fp)) return(cbind(key, sep = NA, t(na), status = "dropped_no_fit"))
    a <- arms_for_fit(fp$fit, do_giviti = dogv)
    cbind(key, sep = fp$firth, t(a), status = "ok")
  }
})

## ---- the job list: chunks of datasets, blocks with a final CSV skipped ----
block_parts <- function(b) {
  m <- regmatches(b, regexec("^S([0-9])_n([0-9]+)_s([0-9]+)$", b))[[1]]
  meta <- utils::read.csv(file.path(RES, b, "blk_meta.csv"))
  list(SET = as.integer(m[2]), NOBS = as.integer(m[3]), SEED0 = as.integer(m[4]),
       nnull = meta$nnull[1], nrand = meta$nrand[1], neval = meta$neval[1])
}
chunks <- list(); todo_blocks <- character(0)
for (b in BLOCKS) {
  if (file.exists(file.path(OUTD, paste0(b, ".csv")))) next
  bp <- block_parts(b); bd <- file.path(OUTD, b); dir.create(bd, showWarnings = FALSE)
  todo_blocks <- c(todo_blocks, b)
  specs <- list(list(kind = 0L, ci = 0L, C = NA_real_, R = bp$nnull))
  if (KIND1) specs <- c(specs, list(list(kind = 1L, ci = 0L, C = NA_real_, R = bp$nrand)))
  for (ci in seq_along(CS)) if (CS[ci] > 0)
    specs <- c(specs, list(list(kind = 2L, ci = ci, C = CS[ci], R = bp$neval)))
  for (sp in specs) {
    R <- min(sp$R, MAXR)
    for (s in seq(1L, R, by = CHUNK)) {
      f <- file.path(bd, sprintf("chunk_k%d_c%d_r%04d.csv", sp$kind, sp$ci, s))
      if (!file.exists(f))
        chunks[[length(chunks) + 1L]] <- c(bp[c("SET", "NOBS", "SEED0")], sp[c("kind", "ci", "C")],
                                           list(reps = s:min(R, s + CHUNK - 1L), f = f))
    }
  }
}
cat(sprintf("liu_addon | %d blocks to finish | %d chunks to run | GiViTI %s | kind1 %s | %d workers | %s\n",
            length(todo_blocks), length(chunks), DOGV, KIND1, NW, format(Sys.time())))
flush.console()

if (length(chunks)) {
  t0 <- Sys.time()
  cl <- parallel::makeCluster(min(NW, length(chunks)))
  parallel::clusterExport(cl, c("HERE", "liu_defs", "DOGV"))
  invisible(parallel::clusterEvalQ(cl, { source(file.path(HERE, "addon_common.R")); eval(liu_defs); NULL }))
  done <- 0L
  ## chunks go out in batches so progress is logged; each chunk writes its own file
  for (bt in split(chunks, ceiling(seq_along(chunks) / (4L * NW)))) {
    parallel::parLapplyLB(cl, bt, function(ch) {
      R <- do.call(rbind, lapply(ch$reps, function(i)
        liu_one(ch$SET, ch$NOBS, ch$SEED0, ch$kind, i, ch$C, ch$ci, DOGV)))
      write_atomic(R, ch$f); NULL
    })
    done <- done + length(bt)
    el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    cat(sprintf("%s  %d / %d chunks | %.2f h | ETA %s\n", format(Sys.time(), "%Y-%m-%d %H:%M"), done,
                length(chunks), el, format(Sys.time() + (length(chunks) - done) * el / done * 3600, "%H:%M")))
    flush.console()
  }
  invisible(parallel::clusterEvalQ(cl, giviti_close()))
  parallel::stopCluster(cl)
}

## ---- bind finished blocks (only when every expected chunk is on disk) ----
for (b in todo_blocks) {
  bd <- file.path(OUTD, b); bp <- block_parts(b)
  need <- min(bp$nnull, MAXR) + (if (KIND1) min(bp$nrand, MAXR) else 0) + 4L * min(bp$neval, MAXR)
  fs <- list.files(bd, pattern = "^chunk_.*\\.csv$", full.names = TRUE)
  R  <- do.call(rbind, lapply(fs, utils::read.csv))
  if (is.null(R) || nrow(R) != need) { cat(sprintf("  %s: %d / %d rows, not bound yet\n", b, NROW(R), need)); next }
  R <- R[order(R$kind, R$ci, R$rep), ]
  write_atomic(R, file.path(OUTD, paste0(b, ".csv")))
  cat(sprintf("  %s: bound %d rows | GiViTI NA %.4f | mHL NA %.4f\n", b, nrow(R),
              mean(is.na(R$GiViTI[R$status == "ok"])), mean(is.na(R$mHL[R$status == "ok"]))))
}
cat("liu_addon DONE", format(Sys.time()), "\n")
