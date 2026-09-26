# probe_setting1.R -- why does BAGofT fail on Liu setting 1 (a single covariate)?
# EXPLORATORY. Calls BAGofT directly, without the wrapper's tryCatch, so the real error shows.
source("/path/to/project/paper_deepgof/explore_bagoft/common.R")

dd <- liu_dataset(101000L, 1L, 50L, 1.0, 1L)
dat <- as.data.frame(lapply(dd$data, as.vector))
cat("columns:", paste(names(dat), collapse = ", "), "| n =", nrow(dat), "\n\n")

try_call <- function(label, expr) {
  cat("---", label, "---\n")
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat("ERROR:", conditionMessage(r), "\n")
    cat("call :", paste(deparse(conditionCall(r)), collapse = " "), "\n\n")
  } else cat("ok: p.value =", r$p.value, "\n\n")
}

set.seed(1)
try_call("BAGofT, parRF() default, one covariate",
  suppressMessages(BAGofT(testModel = testGlmBi(formula = y ~ x1, link = "logit"),
                          parFun = parRF(ntree = 60), data = dat,
                          nsplits = 4, ne = floor(5 * sqrt(nrow(dat))), nsim = 4)))

## where inside BAGofT does it break?
set.seed(1)
tr <- tryCatch(withCallingHandlers(
  suppressMessages(BAGofT(testModel = testGlmBi(formula = y ~ x1, link = "logit"),
                          parFun = parRF(ntree = 60), data = dat,
                          nsplits = 4, ne = floor(5 * sqrt(nrow(dat))), nsim = 4)),
  error = function(e) { cat("traceback frames:\n"); print(sys.calls()) }),
  error = function(e) NULL)

## does a second, pure-noise covariate make it run? (diagnostic only)
set.seed(1)
dat2 <- dat; dat2$z <- rnorm(nrow(dat2))
try_call("diagnostic: add one noise column to the data, model still y ~ x1",
  suppressMessages(BAGofT(testModel = testGlmBi(formula = y ~ x1, link = "logit"),
                          parFun = parRF(ntree = 60), data = dat2,
                          nsplits = 4, ne = floor(5 * sqrt(nrow(dat2))), nsim = 4)))

cat("parRF source (first lines):\n")
cat(head(deparse(BAGofT::parRF), 40), sep = "\n")
