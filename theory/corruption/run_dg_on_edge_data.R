# run_dg_on_edge_data.R -- DeepGOF-1 on the SAME corrupted data sets as the stored corrupted-record studies,
# so its false-alarm rate is compared with the stored p-values of the other tests, replicate by replicate.
# Each data set is rebuilt from its stored seed with the stored generator (L'Ecuyer-CMRG, draw order x, d, y,
# then the corruption), and checked against a stored fingerprint before DeepGOF-1 is run:
#   block 9  (logit truth, n = 1,000, 1,000 reps): clean; C1 = x multiplied by 4 in k = 1, 2, 5, 10 records;
#            C3 = x divided by 10 in k = 1, 2, 5, 10 records.   Fingerprint: events and the fitted b.x.
#   block 9b (n = 500, 100 reps, the projection test stored): clean, x4, x8 in k = 3 records.   Fingerprint: events.
# DeepGOF-1: one set of B = 199 refits, axes = all-pairs = combined here (two covariates, one pair); ties in the
# binary d broken by a random order as deepgof1() does. Output rows/<cell>_r<first>.csv (resumable).
BAT  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/simulations/battery"
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/corruption"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, Sys.getenv("CO_OUT", if (SMOKE) "rows_edge_smoke" else "rows_edge")); dir.create(OUT, showWarnings = FALSE)

CELLS <- rbind(
  data.frame(block = "9", file = c("logit_clean_n1000", sprintf("logit_C1_r%s_n1000", c("001", "002", "005", "010")),
                                   sprintf("logit_C3_r%s_n1000", c("001", "002", "005", "010"))),
             n = 1000L, corr = c("clean", rep("C1", 4), rep("C3", 4)), k = c(0L, 1L, 2L, 5L, 10L, 1L, 2L, 5L, 10L)),
  data.frame(block = "9b", file = c("clean", "x4", "x8"), n = 500L, corr = c("clean", "x4", "x8"), k = c(0L, 3L, 3L)))
jobs <- list()
if (nzchar(Sys.getenv("CO_CELL"))) CELLS <- CELLS[CELLS$file == Sys.getenv("CO_CELL"), ]
for (q in seq_len(nrow(CELLS))) {
  S <- utils::read.csv(gzfile(file.path(BAT, CELLS$block[q], paste0(CELLS$file[q], "_pvalues.csv.gz"))))
  keep <- intersect(c("rep", "seed", "events", "b.x"), names(S)); S <- S[order(S$rep), keep]
  if (SMOKE) S <- S[1:2, ]
  for (ch in split(seq_len(nrow(S)), ceiling(seq_len(nrow(S)) / 25))) {
    f <- file.path(OUT, sprintf("%s_%s_r%04d.csv", CELLS$block[q], CELLS$file[q], S$rep[ch[1]]))
    if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[q, ], S = S[ch, , drop = FALSE], f = f)
  }
}

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(z) ebrahim.gof:::.dg_score((z - M$mu) / M$sd, M)
  RNGkind("L'Ecuyer-CMRG")
  R <- do.call(rbind, lapply(seq_len(nrow(j$S)), function(i) {
    s <- j$S[i, ]; n <- j$cc$n
    RNGkind("L'Ecuyer-CMRG"); set.seed(s$seed)
    x <- stats::runif(n, -3, 3); d <- stats::rbinom(n, 1, 0.5)
    y <- stats::rbinom(n, 1, stats::plogis(0.6 * x + 0.5 * d))
    if (j$cc$block == "9") {
      if (j$cc$corr == "C1") { idx <- sample.int(n, j$cc$k); x[idx] <- 4 * x[idx] }
      if (j$cc$corr == "C3") { idx <- sample.int(n, j$cc$k); x[idx] <- x[idx] / 10 }
    } else if (j$cc$corr != "clean") {
      idx <- sample.int(n, j$cc$k); x[idx] <- (if (j$cc$corr == "x4") 4 else 8) * x[idx]
    }
    dat <- data.frame(y = y, x = x, d = d)
    fit <- suppressWarnings(stats::glm(y ~ x + d, data = dat, family = stats::binomial()))
    ok <- sum(y) == s$events && (is.null(s$b.x) || abs(stats::coef(fit)[["x"]] - s$b.x) < 1e-8)
    p <- NA_real_
    if (ok) {
      mm <- stats::model.matrix(fit); X <- ebrahim.gof:::.dg_covariates(fit)
      tcols <- ebrahim.gof:::.dg_term_cols(fit, mm, colnames(X))
      if (Sys.getenv("CO_MT") == "1") RNGkind("Mersenne-Twister")
      set.seed(s$seed + 7L)
      tb <- if (any(apply(X, 2, anyDuplicated) > 0L)) sample.int(nrow(mm)) else NULL
      E <- ebrahim.gof:::.dg_engine(fit, 199L, 6L, X, mm, tb, score, need = c("axes", "allpairs"), tcols = tcols)
      p <- ebrahim.gof:::.dg_pvalues(E$obs$v, E$boot)[["combined"]]
    }
    data.frame(block = j$cc$block, file = j$cc$file, rep = s$rep, seed = s$seed, match = ok, DG = p)
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("dg-on-edge%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(as.integer(Sys.getenv("CO_WORKERS", if (SMOKE) "4" else "10")))
  parallel::clusterExport(cl, c("LIB", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
