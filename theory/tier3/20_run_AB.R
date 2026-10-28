# 20_run_AB.R -- Tier 3 items A and B, on the data sets of the published results.
#
# Data sets (regenerated from their archived seeds; a fingerprint is asserted for every one):
#   liu   Settings 2-4, n in {50,100,500}, seed block 101000: nulls 1001-2000 and the first 500 alternatives
#         at each C in {.25,.5,.75,1} -- the 27,000 data sets of theory/shipped_level/rank_shipped_liu.R.
#         set.seed(101000 + S*500000 + kind*50000 + ci*5000 + i); gen_liu() verbatim; plain glm.
#         Fingerprint: events == stored HLcheck_events AND p_dg == stored p_shipped (rows_liu).
#   niche the ten-covariate design (theory/niche/run_niche.R), all 14 scenarios, 2,000 nulls per n and
#         1,000 alternatives per cell (16,000). set.seed(123000 + id*100000 + r); gen() and fit_pol() verbatim.
#         Fingerprint: events == stored events AND p_dg == stored DG (theory/niche/rows).
# Arms, all from ONE parametric bootstrap (B = 199) per data set:
#   p_dg          DeepGOF-1 default (selected-axes) reading. liu: deepgof1() of ebrahim.gof 2.8.0 replicated
#                 step for step (same RNG use), as 02_pilot_A_C.R; niche: run_niche.R's dg_both() loop
#                 (fit_pol refits, the batched forward pass), restricted to the selected pair.
#   p_ss          chi-square of the 36 map cells (sum of squares), calibrated by the same refits
#   p_score       score test for xa^2, xb^2, xa*xb of the pair DeepGOF-1 selects (re-selected in every
#                 refit, as DeepGOF-1 is), calibrated by the same refits; aliased columns dropped
#                 (Setting 3: xb binary, df 2). p_score_asym: chi-square(df) p-value of the observed statistic.
#   p_<net>       the 12 Tier 3 networks (nets/, item B), scoring the SAME 200 maps.
#   p_grp         GRPtests::GRPtest(X, y, fam = "binomial"), package defaults; NA when it returns no p-value.
#                 Seeded per data set (data seed + 2e8), run inside keep_rng(); RNG state asserted unchanged.
# Failed refits score +Inf (deepgof1's rule) for every statistic; a bootstrap score statistic that is NA
# (all columns aliased) also counts as +Inf (counted in sc_na_boot).
# Output: one row per data set, appended as it finishes, to rows/<chunk>.csv; reruns skip done rows.
# Usage: Rscript 20_run_AB.R [smoke] [workers]
args  <- commandArgs(TRUE)
SMOKE <- length(args) >= 1 && args[1] == "smoke"
NW    <- if (length(args) >= 2) as.integer(args[2]) else 22L
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
ROWD  <- file.path(T3, if (SMOKE) "rows_AB_smoke" else "rows_AB"); dir.create(ROWD, showWarnings = FALSE)
FAILF <- file.path(T3, "logs", "FINGERPRINT_FAIL.txt")
B <- 199L; K <- 6L; CS <- c(0, .25, .5, .75, 1)
COLS <- c("design", "S", "n", "kind", "ci", "C", "rep", "scen", "events", "firth", "p_dg", "p_dg_stored",
          "p_ss", "p_score", "p_score_asym", "score_df", "n_fail", "sc_na_boot", "axes", "p_grp", "grp_msg",
          paste0("p_", NETS), "p_grp_pilotseed", "t_combo", "t_nets", "t_grp")

## ---------------- stored fingerprints ----------------
SL <- do.call(rbind, lapply(list.files(file.path(PAP, "theory/shipped_level/rows_liu"), "\\.csv$", full.names = TRUE), utils::read.csv))
SL$key <- sprintf("%d_%d_%d_%.2f_%d", SL$setting, SL$n, SL$kind, ifelse(is.na(SL$C), 0, SL$C), SL$rep)
SN <- do.call(rbind, lapply(list.files(file.path(PAP, "theory/niche/rows"), "\\.csv$", full.names = TRUE), function(f)
  utils::read.csv(f, stringsAsFactors = FALSE, na.strings = "NA")))
SN$key <- paste(SN$scen, SN$rep)
stopifnot(nrow(SL) == 27000, !anyDuplicated(SL$key), nrow(SN) == 16000, !anyDuplicated(SN$key))

## ---------------- jobs ----------------
ex <- parse(file.path(PAP, "theory/niche/run_niche.R"))
ne <- new.env(); ne$SMOKE <- FALSE                       # the STUDY's scenario table, never its smoke one
for (e in ex) { s <- paste(deparse(e), collapse = " ")
  if (grepl("^(SC|SEEDB|NNULL|NALT) <-", s) || grepl("^SC\\$", s)) eval(e, ne) }
SC <- ne$SC
stopifnot(ne$SEEDB == 123000L, nrow(SC) == 14, sum(SC$reps) == 16000)
CH <- 25L
jobs <- list()
for (S in 2:4) for (n in c(50L, 100L, 500L)) {
  sets <- list(list(k = 0L, ci = 0L, ids = 1001:2000))
  for (ci in 2:5) sets[[length(sets) + 1]] <- list(k = 2L, ci = ci, ids = 1:500)
  for (st in sets) {
    ids <- if (SMOKE) st$ids[1:2] else st$ids
    for (ch in split(ids, ceiling(seq_along(ids) / CH))) {
      keys <- sprintf("%d_%d_%d_%.2f_%d", S, n, st$k, CS[st$ci + (st$ci == 0)], ch)
      if (st$ci == 0) keys <- sprintf("%d_%d_%d_%.2f_%d", S, n, st$k, 0, ch)
      fp <- SL[match(keys, SL$key), ]; stopifnot(!anyNA(fp$p_shipped))
      jobs[[length(jobs) + 1]] <- list(design = "liu", S = S, n = n, k = st$k, ci = st$ci, C = if (st$ci) CS[st$ci] else 0,
        ids = ch, ev = fp$HLcheck_events, pst = fp$p_shipped,
        file = file.path(ROWD, sprintf("liu_S%d_n%03d_k%d_ci%d_r%04d.csv", S, n, st$k, st$ci, ch[1])))
    }
  }
}
for (i in seq_len(nrow(SC))) {
  sc <- SC[i, ]; ids <- seq_len(sc$reps); if (SMOKE) ids <- 1:2
  for (ch in split(ids, ceiling(seq_along(ids) / CH))) {
    fp <- SN[match(paste(sc$scen, ch), SN$key), ]; stopifnot(!anyNA(fp$DG), nrow(fp) == length(ch))
    jobs[[length(jobs) + 1]] <- list(design = "niche", sc = sc, S = 0L, n = sc$n, k = as.integer(sc$fam != "null"), ci = sc$id,
      C = sc$C, ids = ch, ev = fp$events, pst = fp$DG,
      file = file.path(ROWD, sprintf("niche_%s_r%04d.csv", sc$scen, ch[1])))
  }
}
done_reps <- function(f) {
  if (!file.exists(f)) return(integer(0))
  L <- readLines(f, warn = FALSE); if (length(L) <= 1) return(integer(0))
  ok <- c(TRUE, lengths(regmatches(L[-1], gregexpr(",", L[-1]))) == length(COLS) - 1L & grepl("\n?$", L[-1]))
  if (!all(ok)) writeLines(L[ok], f)                    # drop a line cut short by a crash
  if (sum(ok) <= 1) return(integer(0))
  utils::read.csv(f)$rep
}
todo <- Filter(function(j) length(setdiff(j$ids, done_reps(j$file))) > 0, jobs)
todo <- todo[order(-vapply(todo, function(j) j$n, 1))]          # largest n first
cat(sprintf("Tier3 A+B%s | %d data sets in %d chunks | %d chunks to do | %d workers | %s\n", if (SMOKE) " SMOKE" else "",
            sum(lengths(lapply(jobs, `[[`, "ids"))), length(jobs), length(todo), NW, format(Sys.time())))
flush.console()

## ---------------- worker ----------------
worker_init <- function() {
  source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
  suppressMessages({ library(ebrahim.gof); library(GRPtests); library(brglm2) })
  stopifnot(utils::packageDescription("ebrahim.gof")$Version == "2.8.0")
  ex <- parse(file.path(PAP, "theory/niche/run_niche.R"))
  HERE <<- file.path(PAP, "theory/niche")
  for (e in ex) { s <- paste(deparse(e), collapse = " "); if (grepl("^defs <-", s)) eval(e, globalenv()) }
  eval(defs, globalenv())                               # gen, frm, fit_pol, load_score_batch, ... verbatim
  SB    <<- load_score_batch()                          # the distributed network, batched (niche DG)
  XNETS <<- load_extra_nets()
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
  ## ---- liu: deepgof1() 2.8.0 step for step ----
  combo_liu <- function(fit) {
    mm <- stats::model.matrix(fit)
    tied <- any(apply(mm[, -1, drop = FALSE], 2, anyDuplicated) > 0L)
    tb <- if (tied) sample.int(nrow(mm)) else NULL
    one <- function(f) { o <- ebrahim.gof:::.dg_map(f, K, tb); st <- sc3(mm, f, o$axes)
      list(v = c(score1(o$map), sum(o$map^2), st[["T"]]), df = st[["df"]], map = o$map, axes = o$axes) }
    o <- one(fit); ph <- as.numeric(stats::fitted(fit))
    dat <- stats::model.frame(fit); yn <- names(dat)[1]; frm <- stats::formula(fit)
    V <- matrix(NA_real_, B + 1L, 3); MP <- matrix(NA_real_, B + 1L, K * K); V[1, ] <- o$v; MP[1, ] <- o$map; nf <- 0L
    for (b in seq_len(B)) {
      dat[[yn]] <- stats::rbinom(length(ph), 1L, ph)
      fb <- tryCatch(suppressWarnings(stats::glm(frm, data = dat, family = stats::binomial())), error = function(e) NULL)
      if (is.null(fb)) { V[b + 1L, ] <- Inf; nf <- nf + 1L } else { s <- one(fb); V[b + 1L, ] <- s$v; MP[b + 1L, ] <- s$map }
    }
    pdg <- { s <- V[, 1]; (1 + sum(s[-1] >= s[1])) / (B + 1) }       # exactly deepgof1's p (failures are Inf)
    list(p = c(p_dg = pdg, p_ss = pv(V[, 2]), p_score = pv(V[, 3]),
               p_score_asym = if (is.na(o$v[3])) NA else stats::pchisq(o$v[3], o$df, lower.tail = FALSE)),
         df = o$df, nf = nf, scna = sum(is.na(V[-1, 3])), maps = MP, axes = paste(o$axes, collapse = "+"))
  }
  ## ---- niche: run_niche.R dg_both(), selected pair only ----
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
         df = o$df, nf = nf, scna = sum(is.na(V[-1, 2])), maps = MP, axes = paste(o$axes, collapse = "+"))
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
  for (ii in seq_along(j$ids)) {
    i <- j$ids[ii]; if (i %in% done) next
    if (j$design == "liu") {
      dseed <- liu_seed(101000L, j$S, j$k, j$ci, i)
      seed_at(dseed); d <- gen_liu(j$n, j$S, j$C); xn <- setdiff(names(d), "y")
      if (sum(d$y) != j$ev[ii]) fail(sprintf("liu S%d n%d k%d ci%d rep%d: events %d != stored %d", j$S, j$n, j$k, j$ci, i, sum(d$y), j$ev[ii]))
      fit <- suppressWarnings(stats::glm(stats::as.formula(paste("y ~", paste(xn, collapse = " + "))), data = d, family = stats::binomial()))
      t0 <- proc.time()[["elapsed"]]; seed_at(liu_dgseed(j$S, j$n, j$k, j$ci, i)); cb <- combo_liu(fit)
      t_combo <- proc.time()[["elapsed"]] - t0; firth <- 0; scen <- sprintf("S%d_n%d", j$S, j$n)
      X <- d[, xn, drop = FALSE]; y <- d$y
    } else {
      dseed <- 123000L + j$sc$id * 100000L + i
      seed_at(dseed); dat <- gen(j$sc)
      if (sum(dat$y) != j$ev[ii]) fail(sprintf("niche %s rep%d: events %d != stored %d", j$sc$scen, i, sum(dat$y), j$ev[ii]))
      fit <- fit_pol(dat); firth <- as.numeric(inherits(fit, "brglmFit") || identical(fit$method, "brglmFit"))
      t0 <- proc.time()[["elapsed"]]; cb <- combo_niche(fit, dat); t_combo <- proc.time()[["elapsed"]] - t0
      scen <- j$sc$scen; X <- dat[, setdiff(names(dat), "y")]; y <- dat$y
    }
    if (!isTRUE(all.equal(cb$p[["p_dg"]], j$pst[ii], tolerance = 0)) || cb$p[["p_dg"]] != j$pst[ii])
      fail(sprintf("%s %s rep%d: p_dg %.4f != stored %.4f", j$design, scen, i, cb$p[["p_dg"]], j$pst[ii]))
    t0 <- proc.time()[["elapsed"]]
    pn <- vapply(XNETS, function(f) boot_p(f(cb$maps)), 0)
    t_nets <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]; g <- grp(X, y, dseed + 200000000L); t_grp <- proc.time()[["elapsed"]] - t0
    gp <- if (SMOKE) grp(X, y, 77L + i)$p else NA
    if (sum(y) != j$ev[ii]) fail(sprintf("%s %s rep%d: events changed after GRPtest", j$design, scen, i))
    msg <- substr(gsub("[,\"\r\n]", " ", g$msg), 1, 90)
    row <- c(j$design, j$S, j$n, j$k, j$ci, sprintf("%.2f", j$C), i, scen, sum(y), firth,
             sprintf("%.17g", c(cb$p[["p_dg"]], j$pst[ii], cb$p[["p_ss"]], cb$p[["p_score"]], cb$p[["p_score_asym"]])),
             cb$df, cb$nf, cb$scna, cb$axes, sprintf("%.17g", g$p), msg, sprintf("%.17g", pn), sprintf("%.17g", gp),
             sprintf("%.3f", c(t_combo, t_nets, t_grp)))
    stopifnot(length(row) == length(COLS))
    cat(paste(row, collapse = ","), "\n", sep = "", file = j$file, append = TRUE)
  }
  j$file
}

if (length(todo)) {
  t0 <- Sys.time()
  cl <- parallel::makeCluster(min(NW, length(todo)))
  parallel::clusterExport(cl, c("COLS", "FAILF", "SMOKE", "done_reps"))
  invisible(parallel::clusterCall(cl, worker_init))
  res <- tryCatch(parallel::parLapplyLB(cl, todo, run_chunk), error = function(e) { cat("RUN STOPPED:", conditionMessage(e), "\n"); NULL })
  parallel::stopCluster(cl)
  cat(sprintf("wall %.1f min | %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time())))
  if (is.null(res)) quit(save = "no", status = 1)
}
R <- do.call(rbind, lapply(Filter(file.exists, vapply(jobs, `[[`, "", "file")), utils::read.csv))
cat(sprintf("rows on disk: %d of %d | p_dg == stored on all: %s | GRP NA: %d | %s\n", nrow(R),
            sum(lengths(lapply(jobs, `[[`, "ids"))), all(R$p_dg == R$p_dg_stored), sum(is.na(R$p_grp)), format(Sys.time())))
