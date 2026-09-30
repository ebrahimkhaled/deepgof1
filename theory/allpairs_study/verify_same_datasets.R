# =============================================================================
# verify_same_datasets.R -- prove run_phd.R regenerates the benchmark's OWN datasets.
#
# WHY. The all-pairs run of this same afternoon indexed the departure strengths 1..4 when
# liu_grid.R indexes them 2..5 inside c(0,.25,.5,.75,1). The datasets it produced came from the
# right distribution with the wrong draws, so every number would have looked reasonable and the arm
# would have sat in the table beside columns computed on other data. Only this check caught it.
#
# HOW. Self-consistency is not evidence -- run_phd.R copies its generator from liu_grid.R, so
# comparing it against itself proves nothing. Instead this regenerates the dataset and recomputes a
# CLASSICAL test that liu_grid.R already deposited (blk_pvalues.csv), through the same package call.
# If the dataset is the same one, the p-value matches to machine precision. If it is not, it will
# not, whatever the reason.
#
# Run:  Rscript verify_same_datasets.R
# =============================================================================
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
HERE <- file.path(PAP, "phd")
suppressMessages(library(brglm2))
suppressMessages(library(ebrahim.gof))
source(file.path(HERE, "run_phd_defs.R"))

CS_ALL <- c(0, 0.25, 0.5, 0.75, 1.0)
PER <- as.integer(Sys.getenv("VS_PER", "3"))
TEST <- "HL"

rows <- list()
for (SET in c(2L, 3L, 4L)) for (NOBS in c(50L, 100L, 500L)) for (SEED0 in c(101000L, 102000L)) {
  dep <- utils::read.csv(file.path(PAP, "results", "liu",
                                   sprintf("S%d_n%d_s%d", SET, NOBS, SEED0), "blk_pvalues.csv"))
  dep$C[is.na(dep$C)] <- 0
  for (ci in seq_along(CS_ALL)) {
    kind <- if (ci == 1L) 0L else 2L
    for (i in seq_len(PER)) {
      set.seed(SEED0 + SET * 500000L + kind * 50000L + (if (ci == 1L) 0L else ci) * 5000L + i)
      d <- gen_liu(NOBS, SET, if (kind == 0L) 0 else CS_ALL[ci])
      if (length(unique(d$y)) < 2L) next
      xn <- setdiff(names(d), "y")
      frm <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
      fp <- fit_pol(frm, d); if (is.null(fp)) next
      got <- tryCatch({
        r <- suppressWarnings(suppressMessages(
          run.all.gof(fp$fit, tests = TEST, include_slow = FALSE, install = "no")))
        setNames(r$p_value, r$Test)[[TEST]]
      }, error = function(e) NA_real_)
      k <- which(dep$kind == kind & dep$C == (if (kind == 0L) 0 else CS_ALL[ci]) & dep$rep == i)
      want <- if (length(k) == 1L) dep[[TEST]][k] else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        setting = SET, n = NOBS, seed = SEED0, C = CS_ALL[ci], rep = i,
        recomputed = got, deposited = want, absdiff = abs(got - want))
    }
  }
}
R <- do.call(rbind, rows)
R$match <- is.finite(R$absdiff) & R$absdiff < 1e-10
utils::write.csv(R, file.path(HERE, "VERIFY_SAME_DATASETS.csv"), row.names = FALSE)
cat(sprintf("\n%d of %d datasets reproduce the deposited %s p-value (worst |diff| = %.3g)\n",
            sum(R$match), nrow(R), TEST, max(R$absdiff, na.rm = TRUE)))
if (any(!R$match)) { cat("\nMISMATCHES:\n"); print(R[!R$match, ], row.names = FALSE) } else
  cat("PASS -- run_phd.R runs on the benchmark's own datasets.\n")
