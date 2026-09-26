# =============================================================================
# run_trial.R -- STEP 2: the EXPLORATORY BAGofT small-n trial.
#
# *** EXPLORATORY, NOT CONFIRMATORY. *** The numbers this script produces may
# never be published as a head-to-head. If the result is promising, a separate
# PRE-DECLARED confirmatory run on FRESH seeds must follow.
#
# WHAT IS REUSED, NOT REIMPLEMENTED
#   * the BAGofT call  = jasa_study/R/methods_registry.R :: bagoft_p(), at the
#     PUBLISHED DEFAULTS (nsplits=100, nsim=100, parRF(ntree=60),
#     ne=floor(5*sqrt(n))). Nothing is lowered.
#   * the datasets     = paper_deepgof/code/liu_grid.R :: gen_liu() and its exact
#     per-dataset seeding, so these are BIT-IDENTICAL to the datasets the
#     existing DeepGOF-1 / comparator numbers in results/liu/VERDICT.csv were
#     computed on (seed 101000).
#
# CHECKPOINTING: every finished dataset writes its OWN one-row csv into rows/
# from inside the worker, immediately. The master CSV is appended after each
# chunk. Re-running skips every (setting, n, C, rep) already on disk.
#
# Run: Rscript run_trial.R        Env: EB_WORKERS EB_REPS EB_BUDGET_MIN EB_CELLS EB_SCOPE
#   EB_SCOPE=full  -> Phase 3 (declared 2026-09-16): the whole Liu grid, one rep-major phase.
#                     Setting 1 uses parRF_dropfix.R (BAGofT 1.0.0 crashes on one covariate).
# =============================================================================
source("/path/to/project/paper_deepgof/explore_bagoft/common.R")
source("/path/to/project/paper_deepgof/explore_bagoft/parRF_dropfix.R")

SEED0   <- 101000L
NW      <- as.integer(Sys.getenv("EB_WORKERS", "8"))       # machine is busy: <= 8
REPS    <- as.integer(Sys.getenv("EB_REPS", "60"))         # target reps per cell
BUDGET  <- as.numeric(Sys.getenv("EB_BUDGET_MIN", "90")) * 60
ROWDIR  <- file.path(EXPL, "rows"); dir.create(ROWDIR, showWarnings = FALSE)
MASTER  <- file.path(EXPL, "bagoft_pvalues.csv")

## ---- the pre-declared cell list --------------------------------------------
## Setting 2 = Liu's decisive cell (all nine of their tests fail there).
## Setting 3 = the binary-covariate cell where DeepGOF-1 is weakest.
## C = 0 is the SIZE cell and is never dropped.
CELLS <- expand.grid(SET = c(2L, 3L), n = c(50L, 100L), C = c(0, 0.5, 1.0),
                     KEEP.OUT.ATTRS = FALSE)
## order the cells so that, whatever the budget cuts off, coverage is balanced
## and the size cells are never the ones left short: C=0 first, then C=1, C=0.5.
CELLS <- CELLS[order(match(CELLS$C, c(0, 1.0, 0.5)), CELLS$SET, CELLS$n), ]

## Phase 3: all 60 cells of Liu's grid. Ordered by n so that each parallel chunk holds datasets of
## similar cost -- a chunk waits for its slowest job, and n = 500 costs ~1.5x n = 50.
SCOPE <- Sys.getenv("EB_SCOPE", "phase12")
if (SCOPE == "full") {
  CELLS <- expand.grid(SET = 1:4, n = c(50L, 100L, 500L), C = c(0, 0.25, 0.5, 0.75, 1.0),
                       KEEP.OUT.ATTRS = FALSE)
  CELLS$SET <- as.integer(CELLS$SET)
  CELLS <- CELLS[order(CELLS$n, CELLS$SET, CELLS$C), ]
}
if (nzchar(Sys.getenv("EB_CELLS", ""))) {                  # e.g. "2,3" settings filter
  keep <- as.integer(strsplit(Sys.getenv("EB_CELLS"), ",")[[1]])
  CELLS <- CELLS[CELLS$SET %in% keep, ]
}

rowfile <- function(SET, n, C, i)
  file.path(ROWDIR, sprintf("S%d_n%d_C%s_r%04d.csv", SET, n, sub("\\.", "p", format(C)), i))

## ---- one dataset ------------------------------------------------------------
one_bag <- function(j) {
  f <- rowfile(j$SET, j$n, j$C, j$i)
  if (file.exists(f)) return(NULL)
  dd <- liu_dataset(SEED0, j$SET, j$n, j$C, j$i)
  if (is.null(dd$data)) {
    row <- data.frame(setting = j$SET, n = j$n, C = j$C, rep = j$i, seed = dd$seed,
                      p = NA_real_, p_med = NA_real_, bag_pmin_stat = NA_real_,
                      time_s = 0, status = "dropped_degenerate_y")
  } else {
    frm <- FRM(dd$data)
    set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)   # BAGofT's own RNG stream
    ## Setting 1 has a single covariate, which crashes BAGofT 1.0.0's parRF. It uses the one-line
    ## drop = FALSE fix in parRF_dropfix.R, bit-identical to the original wherever that one runs.
    r <- if (j$SET == 1L) bagoft_p_dropfix(dd$data, frm, nsplits = 100L, nsim = 100L, ntree = 60L)
         else bagoft_p(dd$data, frm, nsplits = 100L, nsim = 100L, ntree = 60L)
    row <- data.frame(setting = j$SET, n = j$n, C = j$C, rep = j$i, seed = dd$seed,
                      p = r$p, p_med = r$p_med, bag_pmin_stat = r$stat_min,
                      time_s = round(r$time_s, 2),
                      status = if (is.na(r$p)) "error_or_na" else "ok")
  }
  row$nsplits <- 100L; row$nsim <- 100L; row$ntree <- 60L
  row$label <- if (j$SET == 1L) "EXPLORATORY_NOT_CONFIRMATORY_DROPFIX" else "EXPLORATORY_NOT_CONFIRMATORY"
  utils::write.csv(row, f, row.names = FALSE)     # CHECKPOINT: immediate, per dataset
  row
}

## ---- job list: PHASED, rep-major round robin inside each phase --------------
## Phase 1 = the 8 priority cells, C in {0, 1.0}: the SIZE cell (mandatory) and
##           the maximum-departure power cell, where a difference is most likely
##           to be visible at the tiny rep counts this budget allows.
## Phase 2 = the 4 C = 0.5 cells, run only if the budget survives phase 1.
## Within a phase the loop is rep-major, so every cell in that phase grows by
## one rep before any cell grows by two -- the budget can be cut at ANY point
## and the coverage stays balanced.
PH <- if (SCOPE == "full") list(CELLS) else
        list(CELLS[CELLS$C %in% c(0, 1.0), ], CELLS[CELLS$C == 0.5, ])
jobs <- list()
for (ph in PH) {
  if (!nrow(ph)) next
  for (i in seq_len(REPS))
    for (k in seq_len(nrow(ph)))
      jobs[[length(jobs) + 1L]] <- list(SET = ph$SET[k], n = ph$n[k],
                                        C = ph$C[k], i = i)
}
todo <- Filter(function(j) !file.exists(rowfile(j$SET, j$n, j$C, j$i)), jobs)

cat(sprintf("EXPLORATORY BAGofT trial | %d cells x %d reps = %d datasets | %d already on disk | %d to do | %d workers | budget %.0f min\n",
            nrow(CELLS), REPS, length(jobs), length(jobs) - length(todo), length(todo), NW, BUDGET/60))
print(CELLS, row.names = FALSE); flush.console()
if (!length(todo)) { cat("nothing to do\n"); quit(save = "no") }

cl <- parallel::makePSOCKcluster(min(NW, length(todo)))
parallel::clusterExport(cl, c("SEED0", "ROWDIR", "rowfile"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  source("/path/to/project/paper_deepgof/explore_bagoft/common.R")
  source("/path/to/project/paper_deepgof/explore_bagoft/parRF_dropfix.R")
  NULL }))
parallel::clusterExport(cl, "one_bag", envir = environment())

t0 <- Sys.time(); ndone <- 0L
for (chunk in split(todo, ceiling(seq_along(todo) / NW))) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (el > BUDGET) { cat(sprintf("BUDGET REACHED at %.1f min -- stopping cleanly.\n", el/60)); break }
  rows <- Filter(Negate(is.null), parallel::parLapplyLB(cl, chunk, one_bag))
  if (length(rows)) {
    df <- do.call(rbind, rows)
    utils::write.table(df, MASTER, sep = ",", row.names = FALSE,
                       append = file.exists(MASTER), col.names = !file.exists(MASTER))
    ndone <- ndone + nrow(df)
  }
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %4d datasets | %.1f min elapsed | %.0f s/dataset wall | eta-to-budget %.0f more\n",
              ndone, el/60, el/max(1,ndone), (BUDGET-el)/max(1e-9, el/max(1,ndone))))
  flush.console()
}
parallel::stopCluster(cl)
cat(sprintf("DONE: %d datasets this run, %.1f min. Master: %s\n", ndone,
            as.numeric(difftime(Sys.time(), t0, units="mins")), MASTER))
