## 13_run.R -- the study. D1, n = 100, B = 199, gamma = .01.
##
## Arms: kind 0 (null, gen_design_null) and kind 2 (alternative, the evaluation stream of
## code/s3_grid.R with amplitude 2*sqrt(100/n)). Seeds exactly as s3_grid.R. Per dataset the
## plug-in test and the exact test see the SAME data and the SAME common random numbers U.
##
## CHECKPOINTED. Each dataset's row is written to rows/k<kind>_r<rep>.csv the moment it finishes
## and appended to per_dataset.csv under a lock, so a cut-off loses nothing and a restart with the
## same command resumes. Jobs are interleaved (null 1, alt 1, null 2, alt 2, ...) so that any
## prefix of the queue is a balanced pair of arms.
##
## env: EX_W (workers, default 6), EX_REPS (datasets per arm, default 300),
##      EX_DEADLINE_MIN (stop queuing new blocks after this many minutes, default 300)

EX <- "/path/to/project/paper_deepgof/theory/exact_mmc"
source(file.path(EX, "10_core.R"))

NW    <- as.integer(Sys.getenv("EX_W", "6"))
NREPS <- as.integer(Sys.getenv("EX_REPS", "300"))
DLMIN <- as.numeric(Sys.getenv("EX_DEADLINE_MIN", "300"))
ROWS  <- file.path(EX, "rows")
PDCSV <- file.path(EX, "per_dataset.csv")
PROG  <- file.path(EX, "progress.log")
dir.create(ROWS, recursive = TRUE, showWarnings = FALSE)

with_lock <- function(dir, expr) {
  lk <- file.path(dir, ".pdlock"); got <- FALSE
  for (k in 1:6000) {
    if (isTRUE(suppressWarnings(dir.create(lk, showWarnings = FALSE)))) { got <- TRUE; break }
    Sys.sleep(0.02)
  }
  on.exit(if (got) unlink(lk, recursive = TRUE), add = TRUE)
  force(expr)
}
hb <- function(msg) cat(sprintf("%s pid%-6d %s\n", format(Sys.time(), "%H:%M:%S"), Sys.getpid(), msg),
                        file = PROG, append = TRUE)
save_row <- function(row) {
  df <- as.data.frame(as.list(row), check.names = FALSE)
  names(df) <- sub("\\..*$", "", names(df))
  rf <- file.path(ROWS, sprintf("k%d_r%d.csv", as.integer(row[["kind"]]), as.integer(row[["rep"]])))
  tf <- sprintf("%s.tmp%d", rf, Sys.getpid())
  utils::write.csv(df, tf, row.names = FALSE)
  if (file.exists(rf)) unlink(rf)
  file.rename(tf, rf)
  with_lock(EX, {
    if (!file.exists(PDCSV)) utils::write.csv(df, PDCSV, row.names = FALSE)
    else utils::write.table(df, PDCSV, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
  })
  invisible(TRUE)
}
job <- function(j) {
  hb(sprintf("START kind=%d rep=%d", j$k, j$i))
  r <- tryCatch(run_dataset(j$k, j$i), error = function(e) { hb(sprintf("ERR %s", conditionMessage(e))); NULL })
  if (is.null(r)) { hb(sprintf("SKIP  kind=%d rep=%d", j$k, j$i)); return(NULL) }
  save_row(r)
  hb(sprintf("DONE  kind=%d rep=%d  p_plug %.4f p_exact %.4f  n_eval %d  %.0f s",
             j$k, j$i, r[["p_plugin"]], r[["p_exact"]], r[["n_eval"]], r[["secs"]]))
  r
}

done <- character(0)
fs <- list.files(ROWS, pattern = "^k\\d+_r\\d+\\.csv$")
if (length(fs)) done <- sub("\\.csv$", "", fs)
JOBS <- Filter(function(j) !(sprintf("k%d_r%d", j$k, j$i) %in% done),
               unlist(lapply(seq_len(NREPS), function(i) list(list(k = 0L, i = i), list(k = 2L, i = i))),
                      recursive = FALSE))

cat(sprintf("exact_mmc | D1 n=%d B=%d gamma=%.2f alpha=%.2f | %d reps/arm | %d workers\n",
            NOBS, BREP, GAMMA, ALPHA, NREPS, NW))
cat(sprintf("grid %d^3 = %d, NM starts %d x maxit %d | RCAP %g Wald SE\n",
            NGRID, NGRID^3, NMSTART, NMMAXIT, RCAP))
cat(sprintf("resume: %d rows on disk, %d jobs to run, deadline %.0f min\n",
            length(done), length(JOBS), DLMIN))
flush.console()

if (length(JOBS)) {
  cl <- parallel::makeCluster(NW)
  parallel::clusterExport(cl, c("EX", "ROWS", "PDCSV", "PROG", "with_lock", "hb", "save_row", "job"),
                          envir = environment())
  invisible(parallel::clusterEvalQ(cl, { source(file.path(EX, "10_core.R")); NULL }))
  t0 <- Sys.time(); BLK <- NW * 3L; ndone <- 0L
  for (st in seq(1L, length(JOBS), by = BLK)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (el > DLMIN) { cat(sprintf("deadline reached at %.1f min; stopping with %d done\n", el, ndone)); break }
    id <- st:min(st + BLK - 1L, length(JOBS))
    invisible(parallel::parLapplyLB(cl, JOBS[id], job))
    ndone <- ndone + length(id)
    cat(sprintf("%s  block %d-%d done (%d/%d jobs, %.1f min elapsed)\n",
                format(Sys.time(), "%H:%M:%S"), min(id), max(id), ndone, length(JOBS),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }
  try(parallel::stopCluster(cl), silent = TRUE)
}

## rebuild per_dataset.csv from the row files
fs <- list.files(ROWS, pattern = "^k\\d+_r\\d+\\.csv$", full.names = TRUE)
A <- if (length(fs)) do.call(rbind, lapply(fs, function(f) tryCatch(utils::read.csv(f), error = function(e) NULL))) else NULL
if (!is.null(A) && nrow(A)) {
  A <- A[!duplicated(A[, c("kind", "rep")]), , drop = FALSE]
  A <- A[order(A$kind, A$rep), , drop = FALSE]
  with_lock(EX, utils::write.csv(A, PDCSV, row.names = FALSE))
  cat(sprintf("\n%d rows on disk: null %d, alt %d\n", nrow(A), sum(A$kind == 0), sum(A$kind == 2)))
  cat(sprintf("median seconds/dataset %.0f, median n_eval %.0f\n", median(A$secs), median(A$n_eval)))
}
cat("RUN INVOCATION COMPLETE\n")
