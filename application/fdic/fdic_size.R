## fdic_size.R -- the size of the published partition tests on the REAL covariate geometry.
##
## For each state cell (Georgia first, then every state with at least 15 failures; the rule was
## fixed before any test was run): fit the 5-covariate logit, then run a parametric bootstrap
## UNDER THE FITTED NULL at the real design matrix -- simulate y ~ Bernoulli(pi-hat), refit, and
## score the Hosmer-Lemeshow, HL equal-width and Pigeon-Heyse tests. The rejection rate at nominal
## .05 is the size of each test on that real geometry. Only the covariates and the fitted
## probabilities enter; the real outcome does not.
##
## Tests: ebrahim.gof (run.all.gof), the same implementations as the benchmark. p-values are clipped
## to [1e-12, 1 - 1e-12] as in the deposited files. The paper's table reads the size with a missing
## p-value counted as no rejection (tables/make_tables.py); the summary written here drops missing
## values, as the original run did. Seeds are 88000 + replicate, as in the deposited run.
##
## Run from the archive root (or set DEEPGOF_ARCHIVE). Writes results/fdic/boot_pvalues_<S>.csv,
## results/fdic/observed_pvalues.csv and results/fdic/real_design_size.csv.

suppressMessages(library(ebrahim.gof))
ARC  <- Sys.getenv("DEEPGOF_ARCHIVE", ".")
REPS <- as.integer(Sys.getenv("FD_REPS", "2000"))
NW   <- as.integer(Sys.getenv("FD_W",    "2"))
MINF <- as.integer(Sys.getenv("FD_MINFAIL", "15"))
OUTD <- file.path(ARC, "results", "fdic"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
TESTS <- c("HL", "HL-equalwidth", "Pigeon-Heyse")

d  <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_2008_banks.csv"), stringsAsFactors = FALSE)
st <- read.csv(file.path(ARC, "application", "fdic", "data", "fdic_state_counts.csv"), stringsAsFactors = FALSE)
CELLS <- c("GA", setdiff(st$state[st$failures >= MINF], "GA"))
FRM <- failed ~ eq_ta + npa_ta + roa + bro_dep + cre_ta

eval_part <- function(fit) {
  out <- setNames(rep(NA_real_, length(TESTS)), TESTS)
  r <- tryCatch(suppressWarnings(run.all.gof(fit, tests = TESTS, include_slow = FALSE, install = "no")),
                error = function(e) NULL)
  if (!is.null(r)) { pv <- setNames(r$p_value, r$Test); k <- intersect(TESTS, names(pv)); out[k] <- pv[k] }
  pmin(pmax(out, 1e-12), 1 - 1e-12)
}

cl <- parallel::makeCluster(NW)
parallel::clusterExport(cl, c("FRM", "TESTS", "eval_part"), envir = environment())
invisible(parallel::clusterEvalQ(cl, { suppressMessages(library(ebrahim.gof)); NULL }))

rows <- list(); obs <- list()
for (S in CELLS) {
  ds <- d[d$state == S, ]
  fit <- suppressWarnings(glm(FRM, data = ds, family = binomial()))
  ph  <- as.numeric(fitted(fit))
  obs[[S]] <- data.frame(state = S, n = nrow(ds), events = sum(ds$failed), rate = mean(ds$failed),
                         t(eval_part(fit)), check.names = FALSE)
  cat(sprintf("%s: n=%3d events=%3d  bootstrapping %d nulls on the real design ...\n",
              S, nrow(ds), sum(ds$failed), REPS)); flush.console()
  parallel::clusterExport(cl, c("ds", "ph"), envir = environment())
  P <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(REPS), function(i) {
    set.seed(88000L + i)
    d2 <- ds; d2$failed <- rbinom(length(ph), 1, ph)
    if (length(unique(d2$failed)) < 2L) return(rep(NA_real_, length(TESTS)))
    f2 <- tryCatch(suppressWarnings(glm(FRM, data = d2, family = binomial())), error = function(e) NULL)
    if (is.null(f2)) rep(NA_real_, length(TESTS)) else eval_part(f2)
  }))
  colnames(P) <- TESTS
  write.csv(P, file.path(OUTD, paste0("boot_pvalues_", S, ".csv")), row.names = FALSE)
  rows[[S]] <- data.frame(state = S, n = nrow(ds), events = sum(ds$failed), test = TESTS,
                          size = as.numeric(colMeans(P <= .05, na.rm = TRUE)),
                          na_rate = as.numeric(colMeans(!is.finite(P))), row.names = NULL)
}
parallel::stopCluster(cl)
write.csv(do.call(rbind, rows), file.path(OUTD, "real_design_size.csv"), row.names = FALSE)
write.csv(do.call(rbind, obs), file.path(OUTD, "observed_pvalues.csv"), row.names = FALSE)
print(do.call(rbind, rows), row.names = FALSE)
