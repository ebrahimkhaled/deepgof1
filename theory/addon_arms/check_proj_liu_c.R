# check_proj_liu_c.R -- the compiled weight (proj_A_liu_fast) against the R definition (proj_A_liu):
# agreement to rounding on the benchmark designs, including duplicated rows, and timing.
# Output: check_proj_liu_c.log
HERE <- "/path/to/project/paper_deepgof/theory/addon_arms"
source(file.path(HERE, "proj_liu.R"))
sink(file.path(HERE, "check_proj_liu_c.log"), split = TRUE)
cat("check_proj_liu_c |", format(Sys.time()), "\n")
set.seed(7)
des <- list(
  "Liu S1 n=50"            = function() cbind(runif(50, -3, 3)),
  "Liu S2 n=100"           = function() cbind(runif(100, -3, 3), rnorm(100, 0, 1.5), rt(100, 4)),
  "Liu S3 n=50"            = function() cbind(runif(50, -3, 3), rbinom(50, 1, .5)),
  "Liu S4 n=500"           = function() cbind(runif(500, -3, 3), runif(500, -3, 3)),
  "Liu S2 n=500"           = function() cbind(runif(500, -3, 3), rnorm(500, 0, 1.5), rt(500, 4)),
  "niche n=200 p=10"       = function() matrix(rnorm(2000), 200, 10),
  "niche n=500 p=10"       = function() matrix(rnorm(5000), 500, 10),
  "duplicated rows n=60"   = function() { X <- cbind(rnorm(50), rnorm(50)); rbind(X, X[1:10, ]) })
for (nm in names(des)) {
  X <- des[[nm]]()
  t1 <- system.time(A1 <- proj_A_liu(X))[[3]]
  t2 <- system.time(A2 <- proj_A_liu_fast(X))[[3]]
  A3 <- proj_A_liu(X, std = FALSE); A4 <- proj_A_liu_fast(X, std = FALSE)
  e <- rnorm(nrow(X))
  cat(sprintf("%-22s max|A_c - A_R| std %.2e, raw %.2e | rel diff in T std %.2e, raw %.2e | R %.2f s, C %.2f s\n",
              nm, max(abs(A1 - A2)), max(abs(A3 - A4)), abs(proj_T(e, A1) - proj_T(e, A2)) / proj_T(e, A1),
              abs(proj_T(e, A3) - proj_T(e, A4)) / proj_T(e, A3), t1, t2))
}
sink()
