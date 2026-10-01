# run_rivals_n1000.R -- the projection test and BAGofT on the SAME stored corrupted data sets at n = 1,000 that
# run_dg_on_edge_data.R gave DeepGOF-1 (logit truth; clean and x multiplied by 4 in k = 1, 2, 5, 10 records).
# Data sets are rebuilt from their stored seeds and checked against the stored events and fitted b.x.
#   projection.gof(fit, B = 199), ebrahim.gof 2.9.0: the first 200 replicates of each of the five cells
#   BAGofT 1.0.0 at its package defaults (parRF), as in the paper: the first 100 replicates of clean, k = 5, k = 10
# One row per replicate and test in rows_rivals/ (resumable; the slow BAGofT jobs are queued first).
BAT  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/simulations/battery"
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/corruption"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
OUT  <- file.path(HERE, "rows_rivals"); dir.create(OUT, showWarnings = FALSE)
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"

CELLS <- data.frame(file = c("logit_clean_n1000", sprintf("logit_C1_r%s_n1000", c("001", "002", "005", "010"))),
                    k = c(0L, 1L, 2L, 5L, 10L), n_proj = 200L, n_bag = c(100L, 0L, 0L, 0L, 100L))
jobs <- list()
for (q in seq_len(nrow(CELLS))) {
  S <- utils::read.csv(gzfile(file.path(BAT, "9", paste0(CELLS$file[q], "_pvalues.csv.gz"))))
  S <- S[order(S$rep), c("rep", "seed", "events", "b.x")]
  for (tst in c("BAGofT", "proj")) {
    m <- if (SMOKE) 1L else if (tst == "proj") CELLS$n_proj[q] else CELLS$n_bag[q]
    if (SMOKE && tst == "BAGofT") m <- 0L
    for (i in seq_len(m)) {
      f <- file.path(OUT, sprintf("%s_%s_r%04d.csv", tst, CELLS$file[q], S$rep[i]))
      if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], s = S[i, ], test = tst, f = f)
    }
  }
}

jobs <- jobs[order(vapply(jobs, function(j) j$test != "BAGofT", TRUE))]   # the slow BAGofT jobs first

one_job <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  RNGkind("L'Ecuyer-CMRG"); set.seed(j$s$seed)
  n <- 1000L
  x <- stats::runif(n, -3, 3); d <- stats::rbinom(n, 1, 0.5)
  y <- stats::rbinom(n, 1, stats::plogis(0.6 * x + 0.5 * d))
  if (j$cc$k > 0) { idx <- sample.int(n, j$cc$k); x[idx] <- 4 * x[idx] }
  dat <- data.frame(y = y, x = x, d = d)
  fit <- suppressWarnings(stats::glm(y ~ x + d, data = dat, family = stats::binomial()))
  ok <- sum(y) == j$s$events && abs(stats::coef(fit)[["x"]] - j$s$b.x) < 1e-8
  p <- NA_real_; secs <- NA_real_
  if (ok) {
    RNGkind("Mersenne-Twister"); set.seed(j$s$seed + 13L)
    t0 <- proc.time()[["elapsed"]]
    p <- tryCatch(if (j$test == "proj") projection.gof(fit, B = 199)$p.value else
           suppressWarnings(suppressMessages(BAGofT::BAGofT(
             testModel = BAGofT::testGlmBi(formula = y ~ x + d, link = "logit"),
             data = dat, parFun = BAGofT::parRF())))$p.value, error = function(e) NA_real_)
    secs <- proc.time()[["elapsed"]] - t0
  }
  row <- data.frame(test = j$test, file = j$cc$file, k = j$cc$k, rep = j$s$rep, seed = j$s$seed, match = ok,
                    p = p, secs = secs)
  tmp <- paste0(j$f, ".part"); utils::write.csv(row, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("rivals-n1000%s | %d jobs | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("RV_WORKERS", if (SMOKE) "6" else "20")))
  parallel::clusterExport(cl, c("LIB", "one_job"))
  invisible(parallel::parLapplyLB(cl, jobs, one_job, chunk.size = 1L))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
