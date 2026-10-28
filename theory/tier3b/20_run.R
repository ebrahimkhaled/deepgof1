# 20_run.R -- Tier 3b runner: non-polynomial departures (threshold, bump, oscillation, corner).
#
# Designs (fitted model linear in all covariates; truth has C * (h - E h) added on the logit scale):
#   D2   x ~ N(0, I_2),  eta = -0.5 + x1 + 0.6 x2 + C h(x1, x2)                      (paper's Table 2 design)
#   D10  x ~ N(0, I_10), eta = -0.5 + x1 + 0.8 x2 + 0.15 sum_{j>=3} x_j + C h(x1, x2)  (ten-covariate design)
# DeepGOF-1 and the arms that share its bootstrap (ONE parametric bootstrap, B = 199, per data set):
#   D2   theory/tier3/20_run_AB.R combo_liu(): deepgof1() of ebrahim.gof 2.8.0 replicated step for step
#        (plain glm, same RNG use; the smoke run asserts p_dg == deepgof1(fit)$p.value on every smoke row).
#   D10  theory/tier3/20_run_AB.R combo_niche(): theory/niche/run_niche.R's loop (fit_pol refits, Firth when
#        not converged, the batched forward pass of the distributed network), selected pair.
#   p_dg (DeepGOF-1, default selected-axes reading), p_ss (chi-square of the 36 map cells), p_score (score test
#   for xa^2, xb^2, xa*xb of the selected pair, re-selected per refit), p_score_asym (chi-square(df) p-value).
# Other arms on the same fit: p_hl (ebrahim.gof::run.all.gof HL, G = 10), p_stukel_fb (Stukel score test,
#   theory/tier3/30_run_C.R stukel(), one-sided fallback), p_grp (GRPtests::GRPtest(X, y, fam = "binomial"),
#   package defaults; NA when no finite p-value; seeded per data set, run inside keep_rng(), RNG asserted).
# Fingerprints per data set: events and xfp = sum of all covariates (printed to 17 digits); 30_verify.R
#   regenerates a sample of rows from their recorded seeds and asserts both.
# Modes: smoke (nulls and alternatives, 2 reps per cell, pilot-grid C, + deepgof1() check on D2)
#        pilot (n = 316, C grid, 100 data sets per point, no GRP)
#        main  (n = 200, 500; C from tables/C_CHOSEN.csv; 500 alternatives per cell, 2,000 nulls per design x n)
# One CSV per chunk of 25 data sets, each row appended as it finishes; reruns skip rows already on disk.
# Usage: Rscript 20_run.R <smoke|pilot|main> [workers]
args <- commandArgs(TRUE)
MODE <- if (length(args) >= 1) args[1] else stop("mode?")
NW   <- if (length(args) >= 2) as.integer(args[2]) else 22L
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b/b_common.R")
ROWD <- file.path(T3B, paste0("rows_", MODE)); dir.create(ROWD, showWarnings = FALSE)
FAILF <- file.path(T3B, "logs", "FINGERPRINT_FAIL.txt")
COLS <- c("mode", "design", "shape", "n", "ci", "C", "rep", "dseed", "bseed", "events", "xfp", "firth",
          "p_dg", "p_dg_pkg", "p_ss", "p_score", "p_score_asym", "score_df", "n_fail", "sc_na_boot", "axes",
          "p_hl", "p_stukel_fb", "stukel_fallback", "p_grp", "grp_msg", "t_combo", "t_grp", "t_other")
CH <- 25L
HM <- h_means()

## ---------------- pilot grid ----------------
HT <- utils::read.csv(file.path(TABB, "H_MEANS.csv"))
GRID_F <- 2^((-4:6) / 2)                                   # 11 points, factor sqrt(2) apart (0.25x to 8x)
pilot_grid <- function(shape) { c0 <- 0.6 / HT$sd_resid_linear[HT$shape == shape]; signif(c0 * GRID_F, 3) }

## ---------------- jobs ----------------
cells <- list()
add <- function(design, shape, n, ci, C, reps) cells[[length(cells) + 1]] <<- list(design = design, shape = shape, n = n, ci = ci, C = C, reps = reps)
if (MODE == "pilot") {
  for (d in DESIGNS) for (s in SHAPES) { g <- pilot_grid(s); for (k in seq_along(g)) add(d, s, N_PILOT, 10L + k, g[k], 1:100) }
  ## added after the first pilot: with C > 0 the corner departure saturates (DeepGOF-1 raw power <= .08 up to
  ## C = 25.6, because fitted probabilities in that corner are already high); the corner is re-piloted with C < 0
  for (d in DESIGNS) { g <- -pilot_grid("corner"); for (k in seq_along(g)) add(d, "corner", N_PILOT, 21L + k, g[k], 1:100) }
} else if (MODE == "smoke") {
  for (d in DESIGNS) for (n in NS) {
    add(d, "null", n, 0L, 0, 1:2)
    for (s in SHAPES) add(d, s, n, 1L, pilot_grid(s)[7], 1:2)          # a large C: alternatives should reject
  }
} else if (MODE == "main") {
  CC <- utils::read.csv(file.path(TABB, "C_CHOSEN.csv"))
  for (d in DESIGNS) for (n in NS) {
    add(d, "null", n, 0L, 0, 1:2000)
    for (s in SHAPES) for (ci in 1:3) add(d, s, n, ci, CC$C[CC$design == d & CC$shape == s & CC$ci == ci], 1:500)
  }
} else stop("unknown mode")
jobs <- list()
for (cl in cells) for (ch in split(cl$reps, ceiling(seq_along(cl$reps) / CH)))
  jobs[[length(jobs) + 1]] <- c(cl[c("design", "shape", "n", "ci", "C")], list(ids = ch,
    file = file.path(ROWD, sprintf("%s_%s_n%03d_c%02d_r%04d.csv", cl$design, cl$shape, cl$n, cl$ci, ch[1]))))
done_reps <- function(f) {
  if (!file.exists(f)) return(integer(0))
  L <- readLines(f, warn = FALSE); if (length(L) <= 1) return(integer(0))
  ok <- c(TRUE, lengths(regmatches(L[-1], gregexpr(",", L[-1]))) == length(COLS) - 1L)
  if (!all(ok)) writeLines(L[ok], f)                     # drop a line cut short by a crash
  if (sum(ok) <= 1) return(integer(0))
  utils::read.csv(f)$rep
}
todo <- Filter(function(j) length(setdiff(j$ids, done_reps(j$file))) > 0, jobs)
todo <- todo[order(-vapply(todo, function(j) j$n * (if (j$design == "D10") 1.3 else 1), 1))]   # slowest first
cat(sprintf("Tier3b %s | %d data sets in %d chunks | %d chunks to do | %d workers | %s\n", MODE,
            sum(lengths(lapply(jobs, `[[`, "ids"))), length(jobs), length(todo), NW, format(Sys.time())))
flush.console()

## ---------------- worker ----------------
worker_init <- function() {
  source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3b/b_common.R")
  suppressMessages({ library(ebrahim.gof); library(GRPtests); library(brglm2); library(statmod) })
  stopifnot(utils::packageDescription("ebrahim.gof")$Version == "2.8.0")
  ex <- parse(file.path(PAP, "theory/niche/run_niche.R"))
  HERE <<- file.path(PAP, "theory/niche")
  for (e in ex) { s <- paste(deparse(e), collapse = " "); if (grepl("^defs <-", s)) eval(e, globalenv()) }
  eval(defs, globalenv())                                # fit_pol, load_score_batch, ... verbatim
  SB <<- load_score_batch()                              # the distributed network, batched (D10)
  HM <<- h_means()
  NULL
}

run_chunk <- function(j) {
  B <- 199L; K <- 6L
  M <- ebrahim.gof:::deepgof1_weights
  score1 <- function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M)
  sc3 <- function(mm, f, ax) {
    ph <- as.numeric(stats::fitted(f)); xa <- mm[, ax[1]]; xb <- mm[, ax[2]]
    score_stat(mm, ph, f$y - ph, cbind(xa^2, xb^2, xa * xb))
  }
  pv <- function(x) { o <- x[1]; if (is.na(o)) return(NA_real_); b <- x[-1]; b[is.na(b)] <- Inf; (1 + sum(b >= o)) / length(x) }
  ## ---- D2: deepgof1() 2.8.0 step for step (theory/tier3/20_run_AB.R combo_liu, verbatim) ----
  combo_liu <- function(fit) {
    mm <- stats::model.matrix(fit)
    tied <- any(apply(mm[, -1, drop = FALSE], 2, anyDuplicated) > 0L)
    tb <- if (tied) sample.int(nrow(mm)) else NULL
    one <- function(f) { o <- ebrahim.gof:::.dg_map(f, K, tb); st <- sc3(mm, f, o$axes)
      list(v = c(score1(o$map), sum(o$map^2), st[["T"]]), df = st[["df"]], map = o$map, axes = o$axes) }
    o <- one(fit); ph <- as.numeric(stats::fitted(fit))
    dat <- stats::model.frame(fit); yn <- names(dat)[1]; frm <- stats::formula(fit)
    V <- matrix(NA_real_, B + 1L, 3); V[1, ] <- o$v; nf <- 0L
    for (b in seq_len(B)) {
      dat[[yn]] <- stats::rbinom(length(ph), 1L, ph)
      fb <- tryCatch(suppressWarnings(stats::glm(frm, data = dat, family = stats::binomial())), error = function(e) NULL)
      if (is.null(fb)) { V[b + 1L, ] <- Inf; nf <- nf + 1L } else { s <- one(fb); V[b + 1L, ] <- s$v }
    }
    pdg <- { s <- V[, 1]; (1 + sum(s[-1] >= s[1])) / (B + 1) }
    list(p = c(p_dg = pdg, p_ss = pv(V[, 2]), p_score = pv(V[, 3]),
               p_score_asym = if (is.na(o$v[3])) NA else stats::pchisq(o$v[3], o$df, lower.tail = FALSE)),
         df = o$df, nf = nf, scna = sum(is.na(V[-1, 3])), axes = paste(o$axes, collapse = "+"))
  }
  ## ---- D10: run_niche.R dg_both(), selected pair only (theory/tier3/20_run_AB.R combo_niche, verbatim) ----
  combo_niche <- function(fit, dat) {
    mm <- stats::model.matrix(fit); vars <- colnames(mm)[-1]; n <- nrow(mm)
    bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
    BI <- sapply(vars, function(v) bn(mm[, v]))
    PR <- utils::combn(length(vars), 2L); np <- ncol(PR)
    Zq <- lapply(seq_len(np), function(q) { Z <- matrix(0, n, K * K)
      Z[cbind(seq_len(n), (BI[, PR[1, q]] - 1L) * K + BI[, PR[2, q]])] <- 1; Z })
    one <- function(f) {
      ph <- as.numeric(stats::fitted(f)); r <- f$y - ph
      scv <- abs(stats::coef(f)[vars]) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
      ax <- sort(order(scv, decreasing = TRUE)[1:2]); q <- which(PR[1, ] == ax[1] & PR[2, ] == ax[2])
      S <- crossprod(Zq[[q]], cbind(r, ph * (1 - ph))); m <- S[, 1] / sqrt(pmax(S[, 2], 1e-8))
      st <- sc3(mm, f, vars[ax])
      list(v = c(sum(m^2), st[["T"]]), df = st[["df"]], map = m, axes = vars[ax])
    }
    o <- one(fit); ph0 <- as.numeric(stats::fitted(fit)); d2 <- dat
    V <- matrix(NA_real_, B + 1L, 2); MP <- matrix(NA_real_, B + 1L, K * K); V[1, ] <- o$v; MP[1, ] <- o$map; nf <- 0L
    for (b in seq_len(B)) {
      d2$y <- stats::rbinom(n, 1, ph0); fb <- fit_pol(d2)
      if (is.null(fb)) { V[b + 1L, ] <- Inf; nf <- nf + 1L } else { s <- one(fb); V[b + 1L, ] <- s$v; MP[b + 1L, ] <- s$map }
    }
    sdg <- SB(MP); sdg[is.na(sdg)] <- Inf
    list(p = c(p_dg = (1 + sum(sdg[-1] >= sdg[1])) / (B + 1), p_ss = pv(V[, 1]), p_score = pv(V[, 2]),
               p_score_asym = if (is.na(o$v[2])) NA else stats::pchisq(o$v[2], o$df, lower.tail = FALSE)),
         df = o$df, nf = nf, scna = sum(is.na(V[-1, 2])), axes = paste(o$axes, collapse = "+"))
  }
  stukel <- function(fit) {                              # theory/tier3/30_run_C.R, verbatim
    ph <- as.numeric(stats::fitted(fit)); g <- stats::qlogis(ph)
    Z <- cbind(0.5 * g^2 * (ph >= 0.5), -0.5 * g^2 * (ph < 0.5))
    z <- suppressWarnings(as.numeric(statmod::glm.scoretest(fit, Z)))
    ok <- is.finite(z) & apply(Z, 2, function(v) any(v != 0))
    z[!ok] <- NA
    both <- if (all(ok)) stats::pchisq(sum(z^2), 2, lower.tail = FALSE) else NA
    fb <- if (all(ok)) 0L else if (any(ok)) 1L else 2L
    pf <- if (fb == 0L) both else if (fb == 1L) stats::pchisq(z[ok]^2, 1, lower.tail = FALSE) else NA
    c(p = pf, fallback = fb)
  }
  hl <- function(fit) {
    before <- .Random.seed
    r <- keep_rng(tryCatch({ rr <- suppressWarnings(suppressMessages(ebrahim.gof::run.all.gof(fit, tests = "HL", include_slow = FALSE, install = "no")))
                             rr$p_value[rr$Test == "HL"][1] }, error = function(e) NA_real_))
    if (!identical(before, .Random.seed)) stop("run.all.gof moved the RNG stream")
    if (length(r) != 1L || !is.finite(r)) NA_real_ else r
  }
  grp <- function(X, y, seed) {
    before <- .Random.seed
    r <- keep_rng({ seed_at(seed)
      tryCatch(list(p = as.numeric(suppressWarnings(suppressMessages(GRPtests::GRPtest(X, y, fam = "binomial")))), msg = ""),
               error = function(e) list(p = NA_real_, msg = conditionMessage(e))) })
    if (!identical(before, .Random.seed)) stop("GRPtest moved the RNG stream despite keep_rng")
    if (length(r$p) != 1L || !is.finite(r$p)) r <- list(p = NA_real_, msg = if (nzchar(r$msg)) r$msg else "no finite p-value")
    r
  }
  fail <- function(msg) { cat(msg, "\n", file = FAILF, append = TRUE); stop(msg) }
  done <- done_reps(j$file)
  if (!file.exists(j$file) || !length(done)) writeLines(paste(COLS, collapse = ","), j$file)
  for (i in j$ids) {
    if (i %in% done) next
    dseed <- data_seed(j$design, j$shape, j$n, j$ci, i); bseed <- boot_seed(dseed)
    seed_at(dseed); dat <- gen_b(j$design, j$shape, j$n, j$C, HM)
    xn <- setdiff(names(dat), "y"); ev <- sum(dat$y); xfp <- sum(as.matrix(dat[, xn]))
    t0 <- proc.time()[["elapsed"]]
    if (j$design == "D2") {
      fit <- suppressWarnings(stats::glm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = dat, family = stats::binomial()))
      firth <- 0
      seed_at(bseed); cb <- combo_liu(fit)
      pkg <- NA_real_
      if (MODE == "smoke") { seed_at(bseed); pkg <- ebrahim.gof::deepgof1(fit, B = 199L)$p.value
        if (pkg != cb$p[["p_dg"]]) fail(sprintf("%s %s n%d rep%d: replicated p_dg %.4f != deepgof1() %.4f", j$design, j$shape, j$n, i, cb$p[["p_dg"]], pkg)) }
    } else {
      fit <- fit_pol(dat); firth <- as.numeric(inherits(fit, "brglmFit") || identical(fit$method, "brglmFit"))
      seed_at(bseed); cb <- combo_niche(fit, dat); pkg <- NA_real_
    }
    t_combo <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    ph_hl <- hl(fit)
    sk <- tryCatch(stukel(fit), error = function(e) c(p = NA, fallback = NA))
    t_other <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    g <- if (MODE == "pilot") list(p = NA_real_, msg = "not run in pilot") else grp(dat[, xn, drop = FALSE], dat$y, grp_seed(dseed))
    t_grp <- proc.time()[["elapsed"]] - t0
    if (sum(dat$y) != ev) fail(sprintf("%s %s n%d rep%d: events changed after third-party calls", j$design, j$shape, j$n, i))
    msg <- substr(gsub("[,\"\r\n]", " ", g$msg), 1, 90)
    row <- c(MODE, j$design, j$shape, j$n, j$ci, sprintf("%.6g", j$C), i, dseed, bseed, ev, sprintf("%.17g", xfp), firth,
             sprintf("%.17g", c(cb$p[["p_dg"]], pkg, cb$p[["p_ss"]], cb$p[["p_score"]], cb$p[["p_score_asym"]])),
             cb$df, cb$nf, cb$scna, cb$axes, sprintf("%.17g", c(ph_hl, sk[["p"]])), sk[["fallback"]],
             sprintf("%.17g", g$p), msg, sprintf("%.3f", c(t_combo, t_grp, t_other)))
    stopifnot(length(row) == length(COLS))
    cat(paste(row, collapse = ","), "\n", sep = "", file = j$file, append = TRUE)
  }
  j$file
}

if (length(todo)) {
  t0 <- Sys.time()
  cl <- parallel::makeCluster(min(NW, length(todo)))
  parallel::clusterExport(cl, c("COLS", "FAILF", "MODE", "done_reps"))
  invisible(parallel::clusterCall(cl, worker_init))
  res <- tryCatch(parallel::parLapplyLB(cl, todo, run_chunk), error = function(e) { cat("RUN STOPPED:", conditionMessage(e), "\n"); NULL })
  parallel::stopCluster(cl)
  cat(sprintf("wall %.1f min | %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time())))
  if (is.null(res)) quit(save = "no", status = 1)
}
R <- do.call(rbind, lapply(Filter(file.exists, vapply(jobs, `[[`, "", "file")), utils::read.csv))
cat(sprintf("rows on disk: %d of %d | GRP NA: %d | HL NA: %d | Stukel NA: %d | refit failures: %d | %s\n", nrow(R),
            sum(lengths(lapply(jobs, `[[`, "ids"))), sum(is.na(R$p_grp)), sum(is.na(R$p_hl)), sum(is.na(R$p_stukel_fb)),
            sum(R$n_fail), format(Sys.time())))
if (MODE == "smoke") {
  cat(sprintf("D2: replicated p_dg == deepgof1() on %d of %d rows\n", sum(R$p_dg == R$p_dg_pkg, na.rm = TRUE), sum(R$design == "D2")))
  print(R[, c("design", "shape", "n", "C", "rep", "events", "p_dg", "p_ss", "p_score", "p_score_asym", "p_hl", "p_stukel_fb", "p_grp", "t_combo", "t_grp")], digits = 3)
}
