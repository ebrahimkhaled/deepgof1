## 11_bench.R -- everything that must pass before a single dataset of the study is run.
##
##  (1) gen_D1_null_theta() is RNG-identical to gen_design_null(n, 2, "indep", .5, 0, .5) of
##      code/universal_design.R, and its captured theta0 really is the data-generating parameter.
##  (2) the fast fit/map/score path agrees with the shipped .deepgof_map + deepgof_score path.
##  (3) pB(theta-hat) computed from the CRN matrix agrees with the same p-value computed through
##      the shipped path on the SAME replicate responses.
##  (4) the Bernstein set covers theta0 at least 1 - gamma of the time, and its box contains it.
##  (5) timings: one pB evaluation, one grid, one whole dataset.

EX <- "/path/to/project/paper_deepgof/theory/exact_mmc"
source(file.path(EX, "10_core.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a
.deepgof_map <- local({
  txt <- readLines(file.path(RD, "deepgof.R"))
  i0  <- grep("^\\.deepgof_map <- function", txt)
  i1  <- grep("^\\}", txt); i1 <- min(i1[i1 > i0])
  e   <- new.env(parent = globalenv())
  eval(parse(text = txt[i0:i1]), envir = e)
  get(".deepgof_map", envir = e)
})
fit_pol <- function(frm, d) {
  g <- tryCatch(suppressWarnings(glm(frm, data = d, family = binomial())), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  b <- coef(g)
  if (any(!is.finite(b)) || any(abs(b[-1]) > 8) || !g$converged) {
    f <- tryCatch(suppressWarnings(glm(frm, data = d, family = binomial(),
                                       method = "brglmFit", type = "AS_mean")), error = function(e) NULL)
    if (!is.null(f) && all(is.finite(coef(f)))) return(list(fit = f, firth = 1L))
  }
  list(fit = g, firth = 0L)
}
score_shipped <- function(m) deepgof_score((m - MU) / SDV, W, 6L)

cat("=== (0) shipped forward pass vs PyTorch =========================================\n")
Xt <- as.matrix(utils::read.csv(file.path(RD, "test_input.csv"), header = FALSE))
Sp <- as.numeric(readLines(file.path(RD, "test_output.csv")))
cat(sprintf("batched net vs test_output.csv : max|diff| = %.3e\n",
            max(abs(deepgof_score_batch(Xt, W) - Sp))))

cat("\n=== (1) generator RNG-identity and theta0 =======================================\n")
dmax <- 0; tmax <- 0
for (i in 1:25) {
  set.seed(4000 + i); a <- gen_design_null(NOBS, 2L, "indep", 0.5, 0L, 0.5)
  set.seed(4000 + i); b <- gen_D1_null_theta(NOBS, 0.5)
  dmax <- max(dmax, max(abs(as.matrix(a$data) - as.matrix(b$data))))
  ## theta0 must reproduce the fitted-probability vector the generator actually sampled from
  mmb <- model.matrix(~ x1 + x2, data = b$data)
  ## rebuild the generator's p0 independently
  set.seed(4000 + i)
  Z <- matrix(rnorm(NOBS * 2L), NOBS, 2L); Z[, 2] <- runif(NOBS, -3, 3)
  dd <- as.data.frame(Z); names(dd) <- c("x1", "x2")
  MM <- model.matrix(~ ., data = dd); bet <- c(0, runif(2, -0.8, 0.8))
  lin <- as.numeric(MM %*% bet); lin <- lin / max(sd(lin), 1e-8) * 1.2
  b0 <- uniroot(function(z) mean(plogis(z + lin)) - 0.5, c(-25, 25))$root
  tmax <- max(tmax, max(abs(plogis(b0 + lin) - plogis(mmb %*% b$theta0))))
}
cat(sprintf("25 seeds: max|data diff| = %.3e ; max|pi(theta0) - p0| = %.3e\n", dmax, tmax))

cat("\n=== (2) fast path vs shipped path ===============================================\n")
set.seed(7)
dif <- numeric(0); mism <- 0L; ntot <- 0L
for (r in 1:6) {
  g <- gen_D1_null_theta(NOBS); P <- prep_D1(g$data)
  ph0 <- plogis(as.numeric(P$mm %*% g$theta0))
  for (b in 1:40) {
    yb <- if (b == 1L) g$data$y else rbinom(NOBS, 1, ph0)
    if (length(unique(yb)) < 2L) next
    dd <- g$data; dd$y <- yb
    fp <- fit_pol(y ~ x1 + x2, dd); if (is.null(fp)) next
    Ss <- score_shipped(.deepgof_map(fp$fit, 6L)$map)
    f  <- fast_fit(as.numeric(yb), P$mm)
    Sf <- score_one(fast_map(as.numeric(yb), f$ph, P))
    dif <- c(dif, abs(Ss - Sf)); mism <- mism + (fp$firth != f$firth); ntot <- ntot + 1L
  }
}
cat(sprintf("%d comparisons: max|score diff| = %.3e, Firth-flag mismatches = %d\n",
            ntot, max(dif), mism))

cat("\n=== (3) pB() vs the shipped path on the same replicate responses ================\n")
set.seed(13)
g <- gen_D1_null_theta(NOBS); P <- prep_D1(g$data); y <- as.numeric(g$data$y)
fo <- fast_fit(y, P$mm); Tobs <- score_one(fast_map(y, fo$ph, P))
U <- matrix(runif(NOBS * BREP), NOBS, BREP)
r <- pB(as.numeric(fo$coef), Tobs, U, P)
pin <- plogis(as.numeric(P$mm %*% fo$coef))
cnt <- 0L
for (b in 1:BREP) {
  yb <- as.numeric(U[, b] <= pin)
  dd <- g$data; dd$y <- yb
  fp <- fit_pol(y ~ x1 + x2, dd)
  Sb <- if (is.null(fp)) Inf else score_shipped(.deepgof_map(fp$fit, 6L)$map)
  cnt <- cnt + (Sb >= Tobs)
}
cat(sprintf("pB = %.6f (count %d) ; shipped-path replication = %.6f (count %d)\n",
            r$p, r$cnt, (1 + cnt) / (BREP + 1), cnt))

cat("\n=== (3b) early abandonment does not change the count ============================\n")
r2 <- pB(as.numeric(fo$coef), Tobs, U, P, floor_count = -1L, chunk = 17L)
r3 <- pB(as.numeric(fo$coef), Tobs, U, P, floor_count = BREP + 5L)   # must abandon immediately
cat(sprintf("chunk 17: count %d (same = %s) ; floor B+5: nrep = %d, full = %s\n",
            r2$cnt, identical(r2$cnt, r$cnt), r3$nrep, r3$full))

cat("\n=== (4) Bernstein set: coverage of theta0, and of the box =======================\n")
nin <- 0L; nbox <- 0L; NCOV <- 400L
relh <- matrix(NA_real_, NCOV, 3)
for (i in 1:NCOV) {
  gg <- gen_dataset(0L, 100000L + i)
  PP <- prep_D1(gg$data); yy <- as.numeric(gg$data$y)
  ff <- fast_fit(yy, PP$mm); if (is.null(ff) || is.null(ff$ph)) next
  nin <- nin + (in_set_slack(matrix(gg$theta0, 1L), yy, PP$mm) <= 1)
  if (i <= 60L) {
    bx <- set_box(yy, PP$mm, as.numeric(ff$coef), ff$ph, 5000L + i, nmc = 8000L)
    nbox <- nbox + all(gg$theta0 >= bx$lo & gg$theta0 <= bx$hi)
    relh[i, ] <- bx$half / bx$se
  }
}
cat(sprintf("theta0 in C_gamma : %d / %d = %.4f   (nominal >= %.2f)\n", nin, NCOV, nin / NCOV, 1 - GAMMA))
cat(sprintf("theta0 in the box : %d / 60\n", nbox))
cat(sprintf("median box half-width / Wald SE : %.2f %.2f %.2f\n",
            median(relh[, 1], na.rm = TRUE), median(relh[, 2], na.rm = TRUE),
            median(relh[, 3], na.rm = TRUE)))

cat("\n=== (5) timing ==================================================================\n")
t0 <- Sys.time(); invisible(pB(as.numeric(fo$coef), Tobs, U, P)); e1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("one full pB (B = %d) : %.3f s\n", BREP, e1))
t0 <- Sys.time(); bx <- set_box(y, P$mm, as.numeric(fo$coef), fo$ph, 99L); e2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("set_box              : %.3f s ; half/SE = %.2f %.2f %.2f ; F = %.2f ; in = %d/20000\n",
            e2, bx$half[1]/bx$se[1], bx$half[2]/bx$se[2], bx$half[3]/bx$se[3], bx$F, bx$nin))
ax <- lapply(1:3, function(k) seq(bx$lo[k], bx$hi[k], length.out = NGRID))
G <- as.matrix(expand.grid(ax[[1]], ax[[2]], ax[[3]]))
sl <- in_set_slack(G, y, P$mm)
cat(sprintf("grid: %d of %d points inside C_gamma\n", sum(sl <= 1), nrow(G)))

cat("\n--- three whole datasets, end to end ---\n")
for (kk in c(0L, 0L, 2L)) {
  ii <- if (kk == 0L) sample.int(1000, 1) else 1L
  rr <- run_dataset(kk, ii)
  cat(sprintf("kind %d rep %d : p_plugin %.4f  p_exact %.4f  n_eval %d  nrep %d  at_hat %d  nm+ %d  %.1f s\n",
              kk, ii, rr[["p_plugin"]], rr[["p_exact"]], rr[["n_eval"]], rr[["nrep_tot"]],
              rr[["at_hat"]], rr[["nm_improved"]], rr[["secs"]]))
}
cat("\nBENCH COMPLETE\n")
