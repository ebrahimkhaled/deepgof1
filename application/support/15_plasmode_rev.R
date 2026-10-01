# plasmode2_SUPPORT.R -- referee B5: level and power on the SUPPORT design itself, with ebrahim.gof 2.9.0.
# Each replicate draws n patients (without replacement) from the 8,873 complete cases, keeps their real
# covariates, and simulates in-hospital death from a TRUE model fitted on all 8,873:
#   null : M0, linear in every covariate          alt : M1, a natural spline (3 df) in mean arterial pressure
# The analyst fits M0 and runs DeepGOF-1 in its three readings (axis rule, all-pairs, combined; one set of
# B = 199 refits), the projection test (B = 199; n <= 1,000), Stukel's test (joint and SstBoth), the
# Hosmer-Lemeshow test and GiViTI. Null replicates 300 at n = 250 and 500, 200 at 1,000; alternatives 200 and 150; and
# 100 null replicates at the full n = 8,873 (DeepGOF-1 only), the level on the application's own design.
# One CSV row per replicate in plasmode2_rows/ (resumable). Usage: Rscript plasmode2_SUPPORT.R [smoke]
LIB <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
suppressMessages({ library(ebrahim.gof, lib.loc = LIB); library(splines) })
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data_prep.R")
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
ROWS <- file.path(HERE, if (SMOKE) "plasmode2_rows_smoke" else "plasmode2_rows")
dir.create(ROWS, showWarnings = FALSE)
p <- prep(2L); D <- p$data; F0 <- p$formula
truth <- list(null = readRDS(file.path(HERE, "screen_SUPPORT", "M0.rds")),
              alt = readRDS(file.path(HERE, "screen_SUPPORT", "M1_b.rds"))$fit)
PT <- lapply(truth, function(m) as.numeric(stats::predict(m, newdata = D, type = "response")))

nrep <- function(k) if (SMOKE) 1L else k
jobs <- rbind(expand.grid(rep = seq_len(nrep(300L)), n = c(250L, 500L), truth = "null", stringsAsFactors = FALSE),
              expand.grid(rep = seq_len(nrep(200L)), n = 1000L, truth = "null", stringsAsFactors = FALSE),
              expand.grid(rep = seq_len(nrep(200L)), n = c(250L, 500L), truth = "alt", stringsAsFactors = FALSE),
              expand.grid(rep = seq_len(nrep(150L)), n = 1000L, truth = "alt", stringsAsFactors = FALSE),
              expand.grid(rep = seq_len(nrep(100L)), n = nrow(D), truth = "null", stringsAsFactors = FALSE))
jobs$seed <- 30261000L + 100000L * (jobs$truth == "alt") + 1000L * match(jobs$n, c(250L, 500L, 1000L, nrow(D))) + jobs$rep
jobs$file <- file.path(ROWS, sprintf("%s_n%04d_r%04d.csv", jobs$truth, jobs$n, jobs$rep))
todo <- jobs[!file.exists(jobs$file), ]
todo <- todo[order(-todo$n == nrow(D)), ]           # the slow full-n nulls first, so they do not trail
todo <- todo[rev(seq_len(nrow(todo))), ]           # second runner: work from the other end of the queue

one_job <- function(jb) {
  suppressMessages({ library(ebrahim.gof, lib.loc = LIB); library(splines) })
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
  set.seed(jb$seed)
  i <- if (jb$n == nrow(D)) seq_len(nrow(D)) else sample(nrow(D), jb$n)
  d <- D[i, ]
  d$y <- stats::rbinom(jb$n, 1L, PT[[jb$truth]][i])
  fit <- tryCatch(suppressWarnings(stats::glm(F0, data = d, family = stats::binomial())), error = function(e) NULL)
  out <- c(DG_axes = NA, DG_all = NA, DG_comb = NA, proj = NA, StukelJ = NA, StukelS = NA, HL = NA, GiViTI = NA)
  secs <- NA
  if (!is.null(fit) && fit$converged) {
    t0 <- Sys.time()
    f <- function(expr) tryCatch(expr, error = function(e) NA_real_)
    pv <- function(z) { v <- z$p_value; if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
    mm <- stats::model.matrix(fit); X <- ebrahim.gof:::.dg_covariates(fit)
    tcols <- ebrahim.gof:::.dg_term_cols(fit, mm, colnames(X))
    set.seed(jb$seed + 7L)
    tb <- if (any(apply(X, 2, anyDuplicated) > 0L)) sample.int(nrow(mm)) else NULL
    E <- ebrahim.gof:::.dg_engine(fit, 199L, 6L, X, mm, tb, score, need = c("axes", "allpairs"), tcols = tcols)
    P <- ebrahim.gof:::.dg_pvalues(E$obs$v, E$boot)
    out[c("DG_axes", "DG_all", "DG_comb")] <- P[c("axes", "allpairs", "combined")]
    if (jb$n <= 1000L) {
      ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
      out[["HL"]] <- f(pv(ebrahim.gof:::gof_hl(ctx)))
      out[["StukelJ"]] <- f(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "joint"))))
      out[["StukelS"]] <- f(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))))
      out[["GiViTI"]] <- f({ v <- givitiR::givitiCalibrationTest(o = as.numeric(ctx$y), e = as.numeric(ctx$ph),
                                                                 devel = "internal")$p.value
                            if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ })
      set.seed(jb$seed + 11L)
      out[["proj"]] <- f(projection.gof(fit, B = 199)$p.value)
    }
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  row <- data.frame(truth = jb$truth, n = jb$n, rep = jb$rep, seed = jb$seed, events = sum(d$y),
                    converged = !is.null(fit) && fit$converged, secs = secs, t(out))
  tmp <- paste0(jb$file, ".part"); utils::write.csv(row, tmp, row.names = FALSE); file.rename(tmp, jb$file)
  jb$file
}

cat(sprintf("plasmode2%s | %d jobs of %d | %s\n", if (SMOKE) " SMOKE" else "", nrow(todo), nrow(jobs), format(Sys.time())))
if (nrow(todo)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("PL_WORKERS", "5")))
  parallel::clusterExport(cl, c("D", "PT", "F0", "LIB", "one_job"))
  invisible(parallel::parLapplyLB(cl, split(todo, seq_len(nrow(todo))), one_job))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
