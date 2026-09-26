## time_one.R -- time one replicate at each n and design (single core)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "rate_common.R"))
cat("ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "\n")
set.seed(1)
for (des in c("D1", "D4", "D3")) for (n in c(50, 100, 200, 400, 800)) {
  t0 <- proc.time()[["elapsed"]]
  r <- one_rep(des, n)
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("%s n=%4d  p=%.3f stat=%.3f sep=%d ninf=%d  %.2fs\n", des, n, r["p"], r["stat"], r["sep"], r["ninf"], t1 - t0))
}
