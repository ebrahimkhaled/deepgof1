# run_ablation.R -- referee B1 (does the network add anything over simple statistics of the same map?) and
# referee M2 (a reading that combines the axis rule and the all-pairs reading). On every dataset, ONE set of
# B = 199 parametric-bootstrap refits gives five p-values (ebrahim.gof 2.9.0 internals, private library):
#   axes      the shipped network on the axis-rule map (the default test)
#   allpairs  the largest network score over every pair of covariates
#   combined  the smaller of the two p-values above, calibrated exactly over the B + 1 draws
#   ss        the sum of squared cells of the axis-rule map (a covariate-partition chi-square)
#   max       the largest |cell| of the axis-rule map
# Designs (fitted model linear in all covariates):
#   liu     Liu et al. (2024) Settings 2-4, n 50/100/500 (generator of the benchmark); nulls 500, alternatives
#           250 per departure strength C in {.25, .5, .75, 1}
#   niche   ten N(0,1) covariates, eta = -0.5 + x1 + 0.8 x2 + 0.15 sum x3..x10 + C f, f = x1 x2 or x1^2 - 1,
#           C in {.3, .5, .8}, n 200/500; nulls 600, alternatives 150 per cell
#   ushape  five N(0,1) covariates, eta = -0.5 + 0.8 x2 + 0.6 x3 + 0.4 x4 + 0.3 x5 + C (x1^2 - 1),
#           C in {.3, .6}, n 200/500; nulls 1000, alternatives 400 per cell
# Output rows/<design>_<cell>_r<first>.csv, 25 datasets per file (resumable). Usage: Rscript run_ablation.R [smoke]
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/ablation"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L

cells <- list()
for (S in 2:4) for (n in c(50L, 100L, 500L)) {
  cells[[length(cells) + 1]] <- list(design = "liu", S = S, fam = "null", C = 0, n = n, reps = 500L)
  for (C in c(.25, .5, .75, 1)) cells[[length(cells) + 1]] <- list(design = "liu", S = S, fam = "alt", C = C, n = n, reps = 250L)
}
for (n in c(200L, 500L)) {
  cells[[length(cells) + 1]] <- list(design = "niche", S = 0, fam = "null", C = 0, n = n, reps = 600L)
  for (f in c("inter", "quad")) for (C in c(.3, .5, .8))
    cells[[length(cells) + 1]] <- list(design = "niche", S = 0, fam = f, C = C, n = n, reps = 150L)
  cells[[length(cells) + 1]] <- list(design = "ushape", S = 0, fam = "null", C = 0, n = n, reps = 1000L)
  for (C in c(.3, .6)) cells[[length(cells) + 1]] <- list(design = "ushape", S = 0, fam = "quad", C = C, n = n, reps = 400L)
}
for (k in seq_along(cells)) cells[[k]]$id <- k
if (SMOKE) cells <- lapply(cells[c(1, 2, 37, 39, 45, 46)], function(z) { z$reps <- 2L; z })

jobs <- list()
for (cc in cells) for (ch in split(seq_len(cc$reps), ceiling(seq_len(cc$reps) / 25))) {
  f <- file.path(OUT, sprintf("%s_c%03d_r%04d.csv", cc$design, cc$id, ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = cc, reps = ch, f = f)
}
if (Sys.getenv("AB_REV") == "1") jobs <- rev(jobs)

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
  gen <- function(cc) {
    n <- cc$n
    if (cc$design == "liu") {
      C <- cc$C
      if (cc$S == 2L) { x1 <- stats::runif(n, -3, 3); x2 <- stats::rnorm(n, 0, 1.5); x3 <- stats::rt(n, 4)
        d <- data.frame(x1, x2, x3); eta <- -2 + 0.3 * x1 + 0.3 * x2 + 0.3 * x3 + C * x1^2
      } else if (cc$S == 3L) { x1 <- stats::runif(n, -3, 3); x2 <- stats::rbinom(n, 1, .5); d <- data.frame(x1, x2)
        eta <- -1.972 + 0.135 * x1 + 1.117 * x2 + C * x1 * x2
      } else { x1 <- stats::runif(n, -3, 3); x2 <- stats::runif(n, -3, 3); d <- data.frame(x1, x2)
        eta <- 0.3 * x1 + 0.3 * x2 + C * x1 * x2 }
    } else if (cc$design == "niche") {
      X <- matrix(stats::rnorm(n * 10), n, 10, dimnames = list(NULL, paste0("x", 1:10))); d <- as.data.frame(X)
      f <- switch(cc$fam, null = 0, inter = X[, 1] * X[, 2], quad = X[, 1]^2 - 1)
      eta <- -0.5 + X[, 1] + 0.8 * X[, 2] + 0.15 * rowSums(X[, 3:10]) + cc$C * f
    } else {
      X <- matrix(stats::rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5))); d <- as.data.frame(X)
      eta <- -0.5 + X %*% c(0, 0.8, 0.6, 0.4, 0.3) + cc$C * (X[, 1]^2 - 1)
    }
    d$y <- stats::rbinom(n, 1, stats::plogis(eta)); d
  }
  R <- do.call(rbind, lapply(j$reps, function(r) {
    s <- 8000000L + j$cc$id * 10000L + r
    set.seed(s)
    d <- gen(j$cc)
    base <- data.frame(design = j$cc$design, S = j$cc$S, fam = j$cc$fam, C = j$cc$C, n = j$cc$n, rep = r)
    na <- data.frame(p_axes = NA, p_allpairs = NA, p_combined = NA, p_ss = NA, p_max = NA, axes_hit = NA, top_hit = NA)
    if (length(unique(d$y)) < 2L) return(cbind(base, na))
    xn <- setdiff(names(d), "y")
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = d,
                                       family = stats::binomial()))
    mm <- stats::model.matrix(fit); X <- ebrahim.gof:::.dg_covariates(fit)
    tcols <- ebrahim.gof:::.dg_term_cols(fit, mm, colnames(X))
    set.seed(s + 1L)
    tb <- if (any(apply(X, 2, anyDuplicated) > 0L)) sample.int(nrow(mm)) else NULL
    E <- tryCatch(ebrahim.gof:::.dg_engine(fit, B, 6L, X, mm, tb, score, need = c("axes", "allpairs"), tcols = tcols),
                  error = function(e) NULL)
    if (is.null(E)) return(cbind(base, na))
    P <- ebrahim.gof:::.dg_pvalues(E$obs$v, E$boot)
    active <- if (j$cc$design == "liu") "x1" else "x1"            # the covariate that carries the misfit
    top <- if (!is.na(E$obs$top)) E$cells$pairs[E$obs$top, ] else NA
    cbind(base, data.frame(p_axes = P[["axes"]], p_allpairs = P[["allpairs"]], p_combined = P[["combined"]],
                           p_ss = P[["ss"]], p_max = P[["max"]],
                           axes_hit = as.integer(active %in% E$obs$axes), top_hit = as.integer(active %in% top)))
  }))
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("ablation%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("AB_WORKERS", "12")))
  parallel::clusterExport(cl, c("LIB", "B", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
