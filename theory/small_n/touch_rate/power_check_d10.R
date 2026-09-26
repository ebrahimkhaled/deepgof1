# =============================================================================
# power_check_d10.R -- W5(ii), the COST of a false fire: when the covariate-only rule deletes a
# clean record, Theorem S.17 says the level is untouched, so the only cost is power at n-1 in
# place of n. This measures that cost per fire, paired: on each dataset the base test is run on
# all n rows and on the n-1 rows left after deleting the record with the largest robust
# outlyingness (the record the rule WOULD delete, deleted unconditionally so every replicate
# pays the cost; multiply by the touch rate for the expected cost).
# d = 10 standard normal covariates, n = 100, eta = -0.5 + x1 + .6 x2 + .3 x3 - .3 x4 + C (x1^2 - 1),
# C = 1, fitted model y ~ x1 + ... + x10. DeepGOF-1 (ebrahim.gof::deepgof1, B = 199, K = 6; the ten
# covariates are continuous, so the tie rule does not arise). 400 paired replicates, rows written per job,
# resumable.  Run: Rscript power_check_d10.R      Env: PC_WORKERS (4), PC_REPS (400)
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/small_n/touch_rate"
PAP  <- "/path/to/project/paper_deepgof"
ROWD <- file.path(HERE, "rows_power"); dir.create(ROWD, showWarnings = FALSE)
NW   <- as.integer(Sys.getenv("PC_WORKERS", "4"))
REPS <- as.integer(Sys.getenv("PC_REPS", "400"))
CHUNK <- 10L; SEEDB <- 92000L; N <- 100L; D <- 10L; CC <- 1.0; B <- 199L

gen <- function(n) {
  X <- matrix(stats::rnorm(n * D), n, D, dimnames = list(NULL, paste0("x", 1:D)))
  eta <- -0.5 + X[, 1] + 0.6 * X[, 2] + 0.3 * X[, 3] - 0.3 * X[, 4] + CC * (X[, 1]^2 - 1)
  d <- as.data.frame(X); d$y <- stats::rbinom(n, 1, stats::plogis(eta)); d
}
omax_index <- function(d) {            # the record the rule would delete (covariates only)
  X <- as.matrix(d[, paste0("x", 1:D)])
  o <- apply(X, 2, function(v) abs(v - stats::median(v)) / stats::mad(v))
  om <- apply(o, 1, max)
  c(which.max(om), max(om))
}
two_p <- function(d) {
  fit <- tryCatch(suppressWarnings(stats::glm(y ~ ., data = d, family = stats::binomial())),
                  error = function(e) NULL)
  if (is.null(fit) || !fit$converged || any(!is.finite(stats::coef(fit)))) return(c(DG = NA))
  dg <- tryCatch(ebrahim.gof::deepgof1(fit, B = B, K = 6L)$p.value, error = function(e) NA_real_)
  c(DG = dg)
}
one_rep <- function(r) {
  set.seed(SEEDB + r)
  d <- gen(N)
  k <- omax_index(d)
  full <- two_p(d); trm <- two_p(d[-k[1], , drop = FALSE])
  data.frame(rep = r, deleted = k[1], omax = k[2], events = sum(d$y),
             p_DG = full[["DG"]], t_DG = trm[["DG"]])
}
jobs <- list()
for (s in seq(1L, REPS, by = CHUNK)) {
  f <- file.path(ROWD, sprintf("r%04d.csv", s))
  if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(reps = s:min(REPS, s + CHUNK - 1L), f = f)
}
cat(sprintf("power check d=%d n=%d C=%.1f | %d jobs | %d workers | %s\n", D, N, CC, length(jobs), NW, format(Sys.time())))
if (!length(jobs)) { cat("nothing to do\n"); quit(save = "no") }
run_job <- function(jb) {
  R <- do.call(rbind, lapply(jb$reps, one_rep))
  tmp <- paste0(jb$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, jb$f)
  sprintf("%s reps %d-%d", format(Sys.time()), min(jb$reps), max(jb$reps))
}
cl <- parallel::makeCluster(min(NW, length(jobs)))
parallel::clusterExport(cl, c("PAP", "SEEDB", "N", "D", "CC", "B", "gen", "omax_index", "two_p", "one_rep"))
invisible(parallel::clusterEvalQ(cl, { suppressMessages(library(ebrahim.gof)); NULL }))
t0 <- Sys.time()
msg <- parallel::parLapplyLB(cl, jobs, run_job)
parallel::stopCluster(cl)
cat(paste(unlist(msg), collapse = "\n"), "\n")
cat(sprintf("DONE %s: %.1f min\n", format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
