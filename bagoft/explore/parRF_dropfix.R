# =============================================================================
# parRF_dropfix.R -- BAGofT's partition function with its single-covariate crash fixed.
#
# *** EXPLORATORY, NOT CONFIRMATORY. ***
#
# THE BUG (BAGofT 1.0.0, parRF, line 20 of the deparsed function):
#     datRf <- Train.data[, -which(names(Train.data) == Rsp)]
# With ONE covariate the result is a single column, which R silently drops to a plain
# vector. The next line (datRf$res <- res) coerces that vector to a list, and dcPre()
# then fails with "argument is of length zero". So BAGofT cannot run on any model with a
# single covariate -- Liu et al.'s Setting 1.
#
# THE FIX: the same line with drop = FALSE. Nothing else changes. With two or more
# covariates the subset is already a data frame, so drop = FALSE is a no-op and the
# patched function returns bit-identical results; verify_dropfix.R proves that on
# Settings 2-4 before Setting 1 is run. A duplicated covariate column would NOT be a
# neutral workaround: parRF sets maxnodes = 5 x (number of covariates), so a second
# column would double the tree size.
# =============================================================================
suppressMessages(library(BAGofT))

make_parRF_dropfix <- function() {
  src <- deparse(BAGofT::parRF)
  bad <- "datRf <- Train.data[, -which(names(Train.data) == Rsp)]"
  i <- grep(bad, src, fixed = TRUE)
  if (length(i) != 1L) stop("parRF source changed: expected exactly one match for the buggy line")
  src[i] <- sub("== Rsp)]", "== Rsp), drop = FALSE]", src[i], fixed = TRUE)
  # evaluate inside the BAGofT namespace so dcPre() and the other internals resolve as usual
  eval(parse(text = src), envir = asNamespace("BAGofT"))
}
parRF_dropfix <- make_parRF_dropfix()

# bagoft_p() from jasa_study/R/methods_registry.R, identical except for parFun.
# Used ONLY for Setting 1; Settings 2-4 keep the registry wrapper unchanged.
bagoft_p_dropfix <- function(dat, formula, nsplits = 100L, nsim = 100L, ntree = 60L) {
  environment(formula) <- globalenv()
  dat <- as.data.frame(lapply(dat, as.vector))
  t0 <- Sys.time()
  b <- tryCatch({
    invisible(capture.output(
      res <- suppressMessages(suppressWarnings(
        BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = formula, link = "logit"),
                       parFun = parRF_dropfix(ntree = ntree),
                       data = dat, nsplits = nsplits,
                       ne = floor(5 * sqrt(nrow(dat))), nsim = nsim)))))
    res
  }, error = function(e) NULL)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(b)) return(list(p = NA_real_, p_med = NA_real_, p_min = NA_real_,
                              stat_min = NA_real_, time_s = el))
  getnum <- function(x) if (!is.null(x) && is.numeric(x) && length(x) == 1) as.numeric(x) else NA_real_
  list(p = getnum(b$p.value), p_med = getnum(b$p.value2), p_min = getnum(b$p.value3),
       stat_min = getnum(b$pmin), time_s = el)
}
