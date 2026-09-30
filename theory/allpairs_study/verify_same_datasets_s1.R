# Gate for run_s1.R: prove it regenerates the benchmark's OWN Setting 1 datasets.
# Same design as the Settings 2-4 gate: recompute a classical test that liu_grid.R already
# deposited, through the same package call. Self-consistency proves nothing, so the generator is
# sourced from run_phd_defs.R and the comparison is against an independently deposited quantity.
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
HERE <- file.path(PAP, "phd")
suppressMessages(library(brglm2)); suppressMessages(library(ebrahim.gof))
source(file.path(HERE, "run_phd_defs.R"))

CS_ALL <- c(0, 0.25, 0.5, 0.75, 1.0)
PER <- as.integer(Sys.getenv("VS_PER", "3")); TEST <- "HL"; SET <- 1L
rows <- list()
for (NOBS in c(50L, 100L, 500L)) for (SEED0 in c(101000L, 102000L)) {
  dep <- utils::read.csv(file.path(PAP, "results", "liu",
                                   sprintf("S%d_n%d_s%d", SET, NOBS, SEED0), "blk_pvalues.csv"))
  dep$C[is.na(dep$C)] <- 0
  for (ci in seq_along(CS_ALL)) {
    kind <- if (ci == 1L) 0L else 2L
    for (i in seq_len(PER)) {
      set.seed(SEED0 + SET * 500000L + kind * 50000L + (if (ci == 1L) 0L else ci) * 5000L + i)
      d <- gen_liu(NOBS, SET, if (kind == 0L) 0 else CS_ALL[ci])
      if (length(unique(d$y)) < 2L) next
      fp <- fit_pol(y ~ x1, d); if (is.null(fp)) next
      got <- tryCatch({
        r <- suppressWarnings(suppressMessages(
          run.all.gof(fp$fit, tests = TEST, include_slow = FALSE, install = "no")))
        setNames(r$p_value, r$Test)[[TEST]]
      }, error = function(e) NA_real_)
      k <- which(dep$kind == kind & dep$C == (if (kind == 0L) 0 else CS_ALL[ci]) & dep$rep == i)
      rows[[length(rows) + 1L]] <- data.frame(
        n = NOBS, seed = SEED0, C = CS_ALL[ci], rep = i, recomputed = got,
        deposited = if (length(k) == 1L) dep[[TEST]][k] else NA_real_)
    }
  }
}
R <- do.call(rbind, rows); R$absdiff <- abs(R$recomputed - R$deposited)
R$match <- is.finite(R$absdiff) & R$absdiff < 1e-10
utils::write.csv(R, file.path(HERE, "VERIFY_SAME_DATASETS_S1.csv"), row.names = FALSE)
cat(sprintf("\n%d of %d Setting 1 datasets reproduce the deposited %s p-value (worst %.3g)\n",
            sum(R$match), nrow(R), TEST, max(R$absdiff, na.rm = TRUE)))
if (any(!R$match)) { cat("\nMISMATCHES:\n"); print(R[!R$match, ], row.names = FALSE) } else
  cat("PASS -- run_s1.R runs on the benchmark's own datasets.\n")
