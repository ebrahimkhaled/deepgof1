## run_rate.R -- level-rate study for DeepGOF-1 (2026-09-24)
## 3 null designs (D1, D4, D3) x 5 sample sizes x REPS null replicates, B = 199, K = 6,
## released deepgof1() with the analyst's plain-glm fit. Every 100-replicate chunk is written to
## RATE_RAW/<design>_n<n>_c<chunk>.csv as soon as it finishes; chunks already on disk are skipped,
## so the run can be killed and restarted (with a different core count) at any time.
##
##   Rscript run_rate.R [cores] [reps_per_cell]
##
## Chunk seeds are deterministic: seed = 20260924 + 1e6*design_index + 1e3*n_index + chunk.

args  <- commandArgs(trailingOnly = TRUE)
CORES <- if (length(args) >= 1) as.integer(args[1]) else 4L
REPS  <- if (length(args) >= 2) as.integer(args[2]) else 7600L
CHUNK <- as.integer(Sys.getenv("RATE_CHUNK", "100"))
HERE  <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
RAW   <- Sys.getenv("RATE_RAW", file.path(HERE, "RATE_RAW"))
dir.create(RAW, showWarnings = FALSE)

DES  <- c("D1", "D4", "D3")
NS   <- c(50L, 100L, 200L, 400L, 800L)
NCH  <- ceiling(REPS / CHUNK)

## task list, round-robin: chunk 1 of every cell, then chunk 2 of every cell, ...
tasks <- list()
for (ch in seq_len(NCH)) for (di in seq_along(DES)) for (ni in seq_along(NS))
  tasks[[length(tasks) + 1L]] <- list(di = di, ni = ni, ch = ch)
fname <- function(t) file.path(RAW, sprintf("%s_n%d_c%03d.csv", DES[t$di], NS[t$ni], t$ch))
todo  <- Filter(function(t) !file.exists(fname(t)), tasks)
cat(sprintf("%d chunks total, %d already done, %d to run on %d cores\n",
            length(tasks), length(tasks) - length(todo), length(todo), CORES))
if (length(todo) == 0L) quit(status = 0)

do_chunk <- function(t) {
  seed <- 20260924L + 1000000L * t$di + 1000L * t$ni + t$ch
  set.seed(seed)
  out <- matrix(NA_real_, CHUNK, 4L, dimnames = list(NULL, c("p", "stat", "sep", "ninf")))
  for (i in seq_len(CHUNK)) out[i, ] <- one_rep(DES[t$di], NS[t$ni])
  df <- data.frame(design = DES[t$di], n = NS[t$ni], chunk = t$ch, rep = seq_len(CHUNK),
                   seed = seed, out)
  tmp <- paste0(fname(t), ".tmp")
  write.csv(df, tmp, row.names = FALSE)
  file.rename(tmp, fname(t))          # atomic: a chunk file is either complete or absent
  sprintf("%s n=%d chunk %d done", DES[t$di], NS[t$ni], t$ch)
}

cl <- parallel::makeCluster(CORES)
parallel::clusterExport(cl, c("HERE", "CHUNK", "DES", "NS", "RAW", "fname", "do_chunk"))
invisible(parallel::clusterEvalQ(cl, source(file.path(HERE, "rate_common.R"))))
t0 <- Sys.time()
res <- parallel::parLapplyLB(cl, todo, do_chunk)
parallel::stopCluster(cl)
cat(sprintf("finished %d chunks in %.1f h\n", length(res), as.numeric(difftime(Sys.time(), t0, units = "hours"))))
