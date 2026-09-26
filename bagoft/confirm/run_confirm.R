# =============================================================================
# run_confirm.R -- the BAGofT confirmatory arm (PREDECLARATION.md is the contract).
#
# 48 alternative cells x replicates 13-112 and 12 null blocks x replicates 13-212 of the
# benchmark's corpus 101000; BAGofT by bagoft_fast() at the published defaults. One row file per
# dataset, written to .part and renamed, so an interrupted run resumes where it stopped and a
# half-written file can never be mistaken for a finished one.
#
# This script computes NO rejection rate and NO comparison. It reports progress and errors only;
# analyze_confirm.py reads the rows once every dataset is on disk.
#
# Run:  Rscript run_confirm.R            Env: BC_WORKERS (20), BC_BATCH (200), BC_SMOKE (0)
#       BC_SMOKE=1 runs 4 datasets at nsplits = nsim = 2 into rows_smoke/ -- plumbing only.
# =============================================================================
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
EXPL <- file.path(PAP, "explore_bagoft")
HERE <- file.path(PAP, "confirm_bagoft")
source(file.path(EXPL, "common.R"))
source(file.path(EXPL, "bagoft_fast.R"))

SEED0   <- 101000L
NW      <- as.integer(Sys.getenv("BC_WORKERS", "20"))
BATCH   <- as.integer(Sys.getenv("BC_BATCH", "200"))
SMOKE   <- Sys.getenv("BC_SMOKE", "0") == "1"
NSPLITS <- if (SMOKE) 2L else 100L
NSIM    <- if (SMOKE) 2L else 100L
ROWDIR  <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows")
dir.create(ROWDIR, showWarnings = FALSE)

ALT_REPS  <- 13:112                     # replicates 1-12 were seen by the exploratory trial
NULL_REPS <- 13:212
CS        <- c(0.25, 0.5, 0.75, 1.0)

rowfile <- function(SET, n, C, i)
  file.path(ROWDIR, sprintf("S%d_n%d_C%s_r%04d.csv", SET, n, sub("\\.", "p", format(C)), i))

## ---- the job list: rep-major round robin, so coverage stays balanced whenever it is cut ---------
jobs <- list()
for (i in sort(union(ALT_REPS, NULL_REPS)))
  for (n in c(50L, 100L, 500L)) for (SET in 1:4) {
    if (i %in% NULL_REPS) jobs[[length(jobs) + 1L]] <- list(SET = SET, n = n, C = 0, i = i)
    if (i %in% ALT_REPS) for (C in CS)
      jobs[[length(jobs) + 1L]] <- list(SET = SET, n = n, C = C, i = i)
  }
if (SMOKE) jobs <- jobs[c(1, 2, 13, 50)]

## ---- one dataset ------------------------------------------------------------------------------
one_job <- function(j) {
  f <- rowfile(j$SET, j$n, j$C, j$i)
  if (file.exists(f)) return(NULL)
  dd  <- liu_dataset(SEED0, j$SET, j$n, j$C, j$i)
  row <- data.frame(setting = j$SET, n = j$n, C = j$C, rep = j$i, seed = dd$seed,
                    p = NA_real_, p_med = NA_real_, p_min = NA_real_, bag_pmin_stat = NA_real_,
                    time_s = 0, status = "dropped_degenerate_y", err = "")
  if (!is.null(dd$data)) {
    dat <- as.data.frame(lapply(dd$data, as.vector))       # what bagoft_p() does first
    frm <- FRM(dat); environment(frm) <- globalenv()
    set.seed(dd$seed + 500000000L)                           # the trial's BAGofT RNG stream
    t0 <- proc.time()[["elapsed"]]
    r  <- tryCatch(suppressWarnings(bagoft_fast(dat, frm, nsplits = NSPLITS, nsim = NSIM,
                                                ntree = 60L)),
                   error = function(e) e)
    row$time_s <- round(proc.time()[["elapsed"]] - t0, 2)
    if (inherits(r, "error")) {
      row$status <- "error"; row$err <- substr(conditionMessage(r), 1, 200)
    } else {
      row$p <- r$p.value; row$p_med <- r$p.value2; row$p_min <- r$p.value3
      row$bag_pmin_stat <- r$pmin; row$status <- "ok"
    }
  }
  row$nsplits <- NSPLITS; row$nsim <- NSIM; row$ntree <- 60L
  row$impl  <- "bagoft_fast (identical to BAGofT 1.0.0, VERIFY_FAST.log)"
  row$label <- if (SMOKE) "SMOKE_TEST_NOT_DATA" else "CONFIRMATORY"
  tmp <- paste0(f, ".part")
  utils::write.csv(row, tmp, row.names = FALSE)
  file.rename(tmp, f)                                         # atomic on the same volume
  row
}

todo <- Filter(function(j) !file.exists(rowfile(j$SET, j$n, j$C, j$i)), jobs)
cat(sprintf("BAGofT CONFIRMATORY%s | %d datasets | %d on disk | %d to do | %d workers | %s\n",
            if (SMOKE) " -- SMOKE TEST" else "", length(jobs), length(jobs) - length(todo),
            length(todo), NW, format(Sys.time())))
flush.console()
if (!length(todo)) { cat("nothing to do\n"); quit(save = "no") }

cl <- parallel::makeCluster(min(NW, length(todo)))
parallel::clusterExport(cl, c("PAP", "EXPL", "HERE", "SEED0", "NSPLITS", "NSIM", "ROWDIR",
                              "SMOKE", "rowfile", "one_job"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(EXPL, "common.R")); source(file.path(EXPL, "bagoft_fast.R")); NULL }))

t0 <- Sys.time(); done <- 0L; nerr <- 0L
for (chunk in split(todo, ceiling(seq_along(todo) / BATCH))) {
  rows <- Filter(Negate(is.null), parallel::parLapplyLB(cl, chunk, one_job))
  done <- done + length(rows)
  nerr <- nerr + sum(vapply(rows, function(r) r$status == "error", TRUE))
  el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  cat(sprintf("%s  %5d / %d written | errors %d | %.2f h elapsed | ETA %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M"), done, length(todo), nerr, el,
              format(Sys.time() + (length(todo) - done) * el / max(done, 1) * 3600, "%H:%M")))
  flush.console()
}
parallel::stopCluster(cl)
cat(sprintf("DONE %s: %d datasets this run, %d errors, %.2f h\n", format(Sys.time()), done, nerr,
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
