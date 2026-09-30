# run_shipped_liu.R -- the SHIPPED DeepGOF-1 on the benchmark of Section 6 (Liu et al. 2024 settings), for the
# referee's M1: the benchmark scored DeepGOF-1 with a network retrained per block; does the deployed network
# do as well on the same datasets?
# Settings 2-4 only: Setting 1 has one covariate and deepgof1() needs two (the benchmark used a 1-D variant).
# Datasets regenerated from their seeds by the archived harness (benchmark/liu/liu_grid.R, gen_liu,
# set.seed(SEED0 + SET*500000 + kind*50000 + ci*5000 + i)), SEED0 = 101000: the block's held-out nulls
# (kind 0, i = 1001..2000) and the first 500 alternatives at each C in {.25,.5,.75,1} (kind 2).
# Plain glm fit and deepgof1(fit, B = 199), as an analyst runs it. Output: rows_liu/<S>_n<n>_<k>_<C>_c<chunk>.csv
suppressMessages(library(ebrahim.gof))
stopifnot(utils::packageDescription("ebrahim.gof")$Version == "2.8.0", exists("projection.gof"))
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
ROWS <- file.path(HERE, if (SMOKE) "rows_liu_smoke" else "rows_liu")
dir.create(ROWS, showWarnings = FALSE, recursive = TRUE)
SEED0 <- 101000L; CS <- c(0, 0.25, 0.5, 0.75, 1.0)

gen_liu <- function(n, setting, C) {                   # verbatim from benchmark/liu/liu_grid.R
  if (setting == 1L) {
    x1 <- runif(n,-3,3); d <- data.frame(x1)
    eta <- -1.963 + 0.982*x1 + C*x1^2
  } else if (setting == 2L) {
    x1 <- runif(n,-3,3); x2 <- rnorm(n,0,1.5); x3 <- rt(n,4)
    d <- data.frame(x1,x2,x3)
    eta <- -2 + 0.3*x1 + 0.3*x2 + 0.3*x3 + C*x1^2
  } else if (setting == 3L) {
    x1 <- runif(n,-3,3); x2 <- rbinom(n,1,.5); d <- data.frame(x1,x2)
    eta <- -1.972 + 0.135*x1 + 1.117*x2 + C*x1*x2
  } else {
    x1 <- runif(n,-3,3); x2 <- runif(n,-3,3); d <- data.frame(x1,x2)
    eta <- 0.3*x1 + 0.3*x2 + C*x1*x2
  }
  d$y <- rbinom(n, 1, plogis(eta)); d
}

jobs <- list()
for (S in 2:4) for (n in c(50L, 100L, 500L)) {
  sets <- list(list(k = 0L, ci = 0L, C = NA, ids = if (SMOKE) 1001:1002 else 1001:2000))
  for (ci in 2:5) sets[[length(sets) + 1]] <- list(k = 2L, ci = ci, C = CS[ci], ids = if (SMOKE) 1:2 else 1:500)
  for (st in sets) {
    ch <- split(st$ids, ceiling(seq_along(st$ids) / if (SMOKE) 2L else 50L))
    for (c in seq_along(ch))
      jobs[[length(jobs) + 1]] <- list(S = S, n = n, k = st$k, ci = st$ci, C = st$C, ids = ch[[c]],
        file = file.path(ROWS, sprintf("S%d_n%03d_k%d_ci%d_c%02d.csv", S, n, st$k, st$ci, c)))
  }
}
todo <- Filter(function(j) !file.exists(j$file), jobs)

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof))
  r <- do.call(rbind, lapply(j$ids, function(i) {
    set.seed(SEED0 + j$S * 500000L + j$k * 50000L + j$ci * 5000L + i)
    d <- gen_liu(j$n, j$S, if (j$k == 0L) 0 else j$C)
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste("y ~", paste(setdiff(names(d), "y"), collapse = " + "))),
                                       data = d, family = stats::binomial()))
    set.seed(880000L + j$S * 100000L + j$n * 10L + j$k * 7L + j$ci + i * 13L)
    dg <- tryCatch(deepgof1(fit, B = 199), error = function(e) NULL)
    data.frame(setting = j$S, n = j$n, kind = j$k, C = j$C, rep = i, HLcheck_events = sum(d$y),
               converged = fit$converged, p_shipped = if (is.null(dg)) NA else dg$p.value)
  }))
  tmp <- paste0(j$file, ".part"); utils::write.csv(r, tmp, row.names = FALSE); file.rename(tmp, j$file)
  j$file
}

cat(sprintf("shipped liu%s | %d of %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(todo), length(jobs), format(Sys.time())))
if (length(todo)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else 20L)
  parallel::clusterExport(cl, c("gen_liu", "SEED0", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, todo, one_chunk))
  parallel::stopCluster(cl)
}
R <- do.call(rbind, lapply(Filter(file.exists, vapply(jobs, `[[`, "", "file")), utils::read.csv))
cat("rows:", nrow(R), " failed:", sum(is.na(R$p_shipped)), "|", format(Sys.time()), "\n")
if (SMOKE) print(head(R, 20))
