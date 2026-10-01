# check_clean_level.R -- DeepGOF-1 read .072 (1,000 reps, SE .007) on the CLEAN n = 1,000 design of the corrupted-
# record study (x ~ U(-3, 3), binary d). Is the excess size real, and does the binary covariate cause it?
#   A  the same design, fresh seeds                                     (replicates the .072 or not)
#   B  d replaced by a continuous N(0, 1) covariate, same coefficient   (is it the binary covariate?)
#   C  as A at n = 500                                                   (does it shrink with n?)
# 1,000 replicates each, B = 199. Output check_clean_level.csv
LIB <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290d"
one <- function(r) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  M <- ebrahim.gof:::deepgof1_weights
  score <- function(z) ebrahim.gof:::.dg_score((z - M$mu) / M$sd, M)
  res <- c(A = NA, B = NA, C = NA)
  for (v in names(res)) {
    n <- if (v == "C") 500L else 1000L
    set.seed(5100000L + match(v, names(res)) * 100000L + r)
    x <- stats::runif(n, -3, 3); d <- if (v == "B") stats::rnorm(n) else stats::rbinom(n, 1, .5)
    y <- stats::rbinom(n, 1, stats::plogis(.6 * x + .5 * d))
    fit <- suppressWarnings(stats::glm(y ~ x + d, family = stats::binomial()))
    mm <- stats::model.matrix(fit); X <- ebrahim.gof:::.dg_covariates(fit)
    tc <- ebrahim.gof:::.dg_term_cols(fit, mm, colnames(X))
    tb <- if (any(apply(X, 2, anyDuplicated) > 0L)) sample.int(nrow(mm)) else NULL
    E <- ebrahim.gof:::.dg_engine(fit, 199L, 6L, X, mm, tb, score, need = c("axes", "allpairs"), tcols = tc)
    res[[v]] <- ebrahim.gof:::.dg_pvalues(E$obs$v, E$boot)[["axes"]]
  }
  res
}
cl <- parallel::makeCluster(20L)
parallel::clusterExport(cl, "LIB")
P <- do.call(rbind, parallel::parLapplyLB(cl, 1:1000, one))
parallel::stopCluster(cl)
utils::write.csv(P, "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/corruption/check_clean_level.csv", row.names = FALSE)
print(round(colMeans(P <= .05), 3)); print(round(colMeans(P <= .10), 3))
