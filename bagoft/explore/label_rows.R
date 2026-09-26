# label_rows.R -- one-off: stamp the EXPLORATORY label onto every CSV already
# written by run_trial.R (rows written before the label column existed).
# *** EXPLORATORY, NOT CONFIRMATORY. *** Values are not touched; one column is added.
EXPL <- "/path/to/project/paper_deepgof/explore_bagoft"
LAB  <- "EXPLORATORY_NOT_CONFIRMATORY"
stamp <- function(f) {
  r <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (!is.null(r$label)) return(0L)
  before <- r
  r$label <- LAB
  utils::write.csv(r, f, row.names = FALSE)
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  stopifnot(isTRUE(all.equal(back[, names(before)], before)))   # values unchanged
  1L
}
fs <- c(list.files(file.path(EXPL, "rows"), pattern = "\\.csv$", full.names = TRUE),
        file.path(EXPL, "bagoft_pvalues.csv"))
n <- sum(vapply(fs, stamp, integer(1)))
cat(sprintf("stamped %d of %d csv files\n", n, length(fs)))
