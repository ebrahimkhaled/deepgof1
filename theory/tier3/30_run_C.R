# 30_run_C.R -- Tier 3 item C: Stukel (Liu et al.'s SstBoth form, with a one-sided fallback) and
# Osius-Rojek on the whole Liu benchmark: Settings 1-4, n in {50,100,500}, seed blocks 101000 and 102000,
# nulls 1-2000 and alternatives 1-1000 at each C in {.25,.5,.75,1} (6,000 per block, 144,000 in all).
#
# Data sets regenerated exactly as code/liu_grid.R did (set.seed(SEED0 + S*500000 + kind*50000 + ci*5000 + i),
# gen_liu() verbatim, fit_pol() verbatim: plain glm, Firth (brglm2 AS_mean) when not converged, a
# non-finite coefficient or |b| > 8). Fingerprint for EVERY data set: the stored le Cessie p-value
# (results/liu/<block>/blk_pvalues.csv, computed by lc_fast on the same fit) must be reproduced to 1e-9
# (relative), and the stored Firth flag must agree; data sets the harness dropped are recorded, not scored.
#
# Stukel: z_ge, z_lt = statmod::glm.scoretest() for the two Stukel columns, each added alone (Liu et al.'s
#   and LogisticDx's SstBoth); stat = z_ge^2 + z_lt^2 on 2 df when both exist. When one column is identically
#   zero (no fitted probability on that side of .5) its score is undefined, and p uses the other term alone
#   on 1 df (fallback = 1). fallback = 2: neither term defined.
# Osius-Rojek: hand-coded (Hosmer & Lemeshow 3rd ed., 5.2.2; 02_pilot_A_C.R) and LogisticDx::gof()
#   (independent code, also gives SstBoth); both compared with the stored value.
# Output: rows_C/<block>_c<chunk>.csv, one row per data set appended as it finishes; reruns skip done rows.
# Usage: Rscript 30_run_C.R [smoke] [workers]
args  <- commandArgs(TRUE)
SMOKE <- length(args) >= 1 && args[1] == "smoke"
NW    <- if (length(args) >= 2) as.integer(args[2]) else 22L
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
ROWD  <- file.path(T3, if (SMOKE) "rows_C_smoke" else "rows_C"); dir.create(ROWD, showWarnings = FALSE)
FAILF <- file.path(T3, "logs", "FINGERPRINT_FAIL_C.txt")
CS <- c(0, .25, .5, .75, 1)
COLS <- c("block", "S", "n", "seed0", "kind", "ci", "C", "rep", "status", "events", "firth", "lc", "lc_stored",
          "z_ge", "z_lt", "p_sst_both", "p_stukel_fb", "fallback", "p_stukel_stored", "p_sst_ldx",
          "p_or_hand", "p_or_ldx", "p_or_stored", "min_pq")
CH <- 250L
jobs <- list()
for (S in 1:4) for (n in c(50L, 100L, 500L)) for (seed0 in c(101000L, 102000L)) {
  blk <- sprintf("S%d_n%d_s%d", S, n, seed0)
  ids <- rbind(data.frame(k = 0L, ci = 0L, i = 1:2000), do.call(rbind, lapply(2:5, function(ci) data.frame(k = 2L, ci = ci, i = 1:1000))))
  if (SMOKE) ids <- ids[c(1, 2, 1001, 1002, 2001, 2002, 3001, 4001, 5001, 6000), ]
  for (c in split(seq_len(nrow(ids)), ceiling(seq_len(nrow(ids)) / CH)))
    jobs[[length(jobs) + 1]] <- list(blk = blk, S = S, n = n, seed0 = seed0, ids = ids[c, ],
      file = file.path(ROWD, sprintf("%s_c%02d.csv", blk, c[1] %/% CH + 1L)))
}
done_keys <- function(f) {
  if (!file.exists(f)) return(character(0))
  L <- readLines(f, warn = FALSE); if (length(L) <= 1) return(character(0))
  ok <- c(TRUE, lengths(regmatches(L[-1], gregexpr(",", L[-1]))) == length(COLS) - 1L)
  if (!all(ok)) writeLines(L[ok], f)
  if (sum(ok) <= 1) return(character(0))
  z <- utils::read.csv(f); paste(z$kind, z$ci, z$rep)
}
todo <- Filter(function(j) length(setdiff(paste(j$ids$k, j$ids$ci, j$ids$i), done_keys(j$file))) > 0, jobs)
todo <- todo[order(-vapply(todo, function(j) j$n, 1))]
cat(sprintf("Tier3 C%s | %d data sets in %d chunks | %d to do | %d workers | %s\n", if (SMOKE) " SMOKE" else "",
            sum(vapply(jobs, function(j) nrow(j$ids), 1)), length(jobs), length(todo), NW, format(Sys.time()))); flush.console()

worker_init <- function() {
  source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
  suppressMessages({ library(brglm2); library(statmod); library(LogisticDx) })
  ## lc_fast only (the file's self-check section below it re-seeds and is not run)
  src <- readLines(file.path(PAP, "code/lecessie_algebraic.R"))
  a <- grep("^lc_fast <- function", src); b <- a - 1 + grep("^}", src[a:length(src)])[1]
  eval(parse(text = src[a:b]), envir = globalenv())
  NULL
}
run_chunk <- function(j) {
  fit_pol <- function(frm, d) {                          # code/liu_grid.R, verbatim
    g <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial())), error=function(e) NULL)
    if (is.null(g)) return(NULL)
    b <- coef(g)
    if (any(!is.finite(b)) || any(abs(b[-1])>8) || !g$converged) {
      f <- tryCatch(suppressWarnings(glm(frm, data=d, family=binomial(),
                                         method="brglmFit", type="AS_mean")), error=function(e) NULL)
      if (!is.null(f) && all(is.finite(coef(f)))) return(list(fit=f, firth=1L))
    }
    list(fit=g, firth=0L)
  }
  or_hand <- function(fit) {                             # 02_pilot_A_C.R hand()
    X <- stats::model.matrix(fit); y <- fit$y; ph <- as.numeric(stats::fitted(fit))
    key <- apply(X, 1, paste, collapse = "|")
    m <- as.numeric(tapply(y, key, length)); yy <- as.numeric(tapply(y, key, sum))
    pj <- as.numeric(tapply(ph, key, mean)); Xj <- X[!duplicated(key), , drop = FALSE][order(unique(key)), , drop = FALSE]
    vj <- m * pj * (1 - pj); X2 <- sum((yy - m * pj)^2 / vj); J <- length(m)
    cj <- (1 - 2 * pj) / vj
    wf <- stats::lm.wfit(Xj, cj, vj); RSS <- sum(vj * wf$residuals^2)
    A <- 2 * (J - sum(1 / m)); z <- (X2 - (J - ncol(X))) / sqrt(A + RSS)
    2 * stats::pnorm(-abs(z))
  }
  stukel <- function(fit) {
    ph <- as.numeric(stats::fitted(fit)); g <- stats::qlogis(ph)
    Z <- cbind(0.5 * g^2 * (ph >= 0.5), -0.5 * g^2 * (ph < 0.5))
    z <- suppressWarnings(as.numeric(statmod::glm.scoretest(fit, Z)))
    ok <- is.finite(z) & apply(Z, 2, function(v) any(v != 0))
    z[!ok] <- NA
    both <- if (all(ok)) stats::pchisq(sum(z^2), 2, lower.tail = FALSE) else NA
    fb <- if (all(ok)) 0L else if (any(ok)) 1L else 2L
    pf <- if (fb == 0L) both else if (fb == 1L) stats::pchisq(z[ok]^2, 1, lower.tail = FALSE) else NA
    c(z_ge = z[1], z_lt = z[2], p_sst_both = both, p_stukel_fb = pf, fallback = fb)
  }
  fail <- function(msg) { cat(msg, "\n", file = FAILF, append = TRUE); stop(msg) }
  st <- utils::read.csv(file.path(PAP, "results/liu", j$blk, "blk_pvalues.csv"), check.names = FALSE)
  st <- st[st$kind != 1, ]; st$key <- paste(st$kind, ifelse(is.na(st$C), 0, st$C), st$rep)
  done <- done_keys(j$file)
  if (!file.exists(j$file) || !length(done)) writeLines(paste(COLS, collapse = ","), j$file)
  for (r in seq_len(nrow(j$ids))) {
    k <- j$ids$k[r]; ci <- j$ids$ci[r]; i <- j$ids$i[r]; C <- if (ci) c(0, .25, .5, .75, 1)[ci] else 0
    if (paste(k, ci, i) %in% done) next
    seed_at(liu_seed(j$seed0, j$S, k, ci, i)); d <- gen_liu(j$n, j$S, C)
    s <- st[match(paste(k, C, i), st$key), ]
    out <- setNames(rep(NA_real_, 14), c("lc", "lc_stored", "z_ge", "z_lt", "p_sst_both", "p_stukel_fb", "fallback",
                                         "p_stukel_stored", "p_sst_ldx", "p_or_hand", "p_or_ldx", "p_or_stored", "firth", "min_pq"))
    status <- "ok"
    if (is.na(s$rep)) {
      status <- if (length(unique(d$y)) < 2L) "dropped_one_class" else "not_stored"
    } else {
      frm <- stats::as.formula(paste("y ~", paste(setdiff(names(d), "y"), collapse = " + ")))
      assign("d", d, envir = globalenv())
      fp <- fit_pol(frm, d); fit <- fp$fit
      X <- stats::model.matrix(fit); ph <- as.numeric(stats::fitted(fit))
      lc <- tryCatch(unname(lc_fast(as.numeric(fit$y), ph, X, lapply(seq_len(ncol(X))[-1], function(q) X[, q]))["p"]), error = function(e) NA)
      if (fp$firth != s$sep || !isTRUE(abs(lc - s[["le-Cessie"]]) <= 1e-9 * max(1e-12, abs(s[["le-Cessie"]]))))
        fail(sprintf("%s k%d ci%d rep%d: firth %d/%d, le Cessie %.12g vs stored %.12g", j$blk, k, ci, i, fp$firth, s$sep, lc, s[["le-Cessie"]]))
      sk <- tryCatch(stukel(fit), error = function(e) c(z_ge = NA, z_lt = NA, p_sst_both = NA, p_stukel_fb = NA, fallback = NA))
      oh <- tryCatch(or_hand(fit), error = function(e) NA)
      lx <- keep_rng(tryCatch({ G <- as.data.frame(suppressWarnings(LogisticDx::gof(fit, plotROC = FALSE))$gof)
                                nm <- as.character(G[[1]]); pc <- grep("^p", names(G), ignore.case = TRUE)[1]
                                c(G[match("OsRo", nm), pc], G[match("SstBoth", nm), pc]) }, error = function(e) c(NA, NA)))
      out[c("lc", "lc_stored", names(sk), "p_stukel_stored", "p_sst_ldx", "p_or_hand", "p_or_ldx", "p_or_stored", "firth", "min_pq")] <-
        c(lc, s[["le-Cessie"]], sk, s$Stukel, lx[2], oh, lx[1], s[["Osius-Rojek"]], fp$firth, min(ph, 1 - ph))
    }
    row <- c(j$blk, j$S, j$n, j$seed0, k, ci, sprintf("%.2f", C), i, status, sum(d$y), out[["firth"]],
             sprintf("%.17g", out[c("lc", "lc_stored", "z_ge", "z_lt", "p_sst_both", "p_stukel_fb")]), out[["fallback"]],
             sprintf("%.17g", out[c("p_stukel_stored", "p_sst_ldx", "p_or_hand", "p_or_ldx", "p_or_stored", "min_pq")]))
    stopifnot(length(row) == length(COLS))
    cat(paste(row, collapse = ","), "\n", sep = "", file = j$file, append = TRUE)
  }
  j$file
}
if (length(todo)) {
  t0 <- Sys.time()
  cl <- parallel::makeCluster(min(NW, length(todo)))
  parallel::clusterExport(cl, c("COLS", "FAILF", "done_keys"))
  invisible(parallel::clusterCall(cl, worker_init))
  res <- tryCatch(parallel::parLapplyLB(cl, todo, run_chunk), error = function(e) { cat("RUN STOPPED:", conditionMessage(e), "\n"); NULL })
  parallel::stopCluster(cl)
  cat(sprintf("wall %.1f min | %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time())))
  if (is.null(res)) quit(save = "no", status = 1)
}
R <- do.call(rbind, lapply(Filter(file.exists, vapply(jobs, `[[`, "", "file")), utils::read.csv))
cat(sprintf("rows on disk: %d | status: %s | %s\n", nrow(R), paste(names(table(R$status)), table(R$status), collapse = " "), format(Sys.time())))
