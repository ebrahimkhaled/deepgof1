# scan_candidates.R -- METADATA ONLY (dimensions and variable types). No model is fitted, no test is run.
# Lists every data set in the installed packages named in the plan, and for each data frame reports
# n, number of columns, the two-level (binary) variables, and the count of numeric variables with
# more than 10 distinct values ("continuous").
pk <- c("aplore3", "MASS", "mlbench", "carData", "survival", "rms")
out <- list()
for (p in pk) {
  if (!requireNamespace(p, quietly = TRUE)) { cat("not installed:", p, "\n"); next }
  ds <- utils::data(package = p)$results[, "Item"]
  ds <- sub(" .*", "", ds)
  for (nm in ds) {
    e <- new.env()
    ok <- tryCatch({ utils::data(list = nm, package = p, envir = e); TRUE }, error = function(err) FALSE, warning = function(w) TRUE)
    if (!ok || !exists(nm, envir = e)) next
    d <- get(nm, envir = e)
    if (!is.data.frame(d)) { out[[length(out) + 1]] <- data.frame(pkg = p, data = nm, class = class(d)[1], n = NA, ncol = NA,
                             n_cc = NA, binary = "", n_cont = NA); next }
    bin <- names(d)[vapply(d, function(v) length(unique(stats::na.omit(v))) == 2, TRUE)]
    cont <- sum(vapply(d, function(v) is.numeric(v) && length(unique(stats::na.omit(v))) > 10, TRUE))
    out[[length(out) + 1]] <- data.frame(pkg = p, data = nm, class = "data.frame", n = nrow(d), ncol = ncol(d),
                                         n_cc = sum(stats::complete.cases(d)), binary = paste(bin, collapse = ","), n_cont = cont)
  }
}
R <- do.call(rbind, out)
utils::write.csv(R, "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche/candidate_scan.csv", row.names = FALSE)
options(width = 250)
print(R, row.names = FALSE)
