# verify_dropfix.R -- prove the drop = FALSE patch is inert wherever unpatched BAGofT runs,
# then confirm it runs on the single-covariate Setting 1. EXPLORATORY; cheap settings.
source("/path/to/project/paper_deepgof/explore_bagoft/common.R")
source("/path/to/project/paper_deepgof/explore_bagoft/parRF_dropfix.R")

src <- deparse(parRF_dropfix)
cat("patched line:", trimws(grep("datRf <- Train.data", src, value = TRUE, fixed = TRUE)), "\n\n")

run <- function(parFun, dat, seed) {
  set.seed(seed)
  suppressWarnings(suppressMessages(BAGofT::BAGofT(
    testModel = BAGofT::testGlmBi(formula = FRM(dat), link = "logit"), parFun = parFun,
    data = as.data.frame(lapply(dat, as.vector)), nsplits = 6, ne = floor(5 * sqrt(nrow(dat))), nsim = 6)))
}

cat("=== 1. identity on the settings where unpatched BAGofT works ===\n")
ok_all <- TRUE
for (SET in 2:4) for (n in c(50L, 500L)) for (i in 1:2) {
  dd <- liu_dataset(101000L, SET, n, 1.0, i)
  if (is.null(dd$data)) next
  a <- run(BAGofT::parRF(ntree = 60), dd$data, dd$seed)
  b <- run(parRF_dropfix(ntree = 60), dd$data, dd$seed)
  same <- identical(c(a$p.value, a$p.value2, a$p.value3, a$pmean, a$pmedian, a$pmin),
                    c(b$p.value, b$p.value2, b$p.value3, b$pmean, b$pmedian, b$pmin))
  ok_all <- ok_all && same
  cat(sprintf("  S%d n=%-3d rep %d: original p=%.3f pmin=%.3g | patched p=%.3f pmin=%.3g | identical=%s\n",
              SET, n, i, a$p.value, a$pmin, b$p.value, b$pmin, same))
}
cat("ALL IDENTICAL:", ok_all, "\n\n")

cat("=== 2. Setting 1 (one covariate): original vs patched ===\n")
for (n in c(50L, 100L, 500L)) {
  dd <- liu_dataset(101000L, 1L, n, 1.0, 1L)
  a <- tryCatch(run(BAGofT::parRF(ntree = 60), dd$data, dd$seed), error = function(e) conditionMessage(e))
  b <- tryCatch(run(parRF_dropfix(ntree = 60), dd$data, dd$seed), error = function(e) conditionMessage(e))
  cat(sprintf("  S1 n=%-3d: original -> %s | patched -> %s\n", n,
              if (is.character(a)) paste("ERROR:", a) else sprintf("p=%.3f", a$p.value),
              if (is.character(b)) paste("ERROR:", b) else sprintf("p=%.3f pmin=%.3g", b$p.value, b$pmin)))
}

cat("\n=== 3. the wrapper used for Setting 1 returns the same fields as the registry wrapper ===\n")
dd <- liu_dataset(101000L, 1L, 50L, 1.0, 2L)
set.seed(dd$seed + JASA$METHOD_RNG_OFFSET)
w <- bagoft_p_dropfix(dd$data, FRM(dd$data), nsplits = 4L, nsim = 4L, ntree = 60L)
cat("  names:", paste(names(w), collapse = ", "), "| p =", w$p, "| error-free:", !is.na(w$p), "\n")
