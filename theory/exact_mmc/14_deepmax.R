## 14_deepmax.R -- how much does the protocol maximiser (11^3 box grid + 3 Nelder-Mead runs)
## under-estimate the true supremum?
##
## A failed search breaks exactness: p_exact is only valid if it really is the sup over C_gamma.
## This audit re-maximises on a subset of datasets with a search ~20x heavier --
##   * the six coordinate extremes of C_gamma found by set_box(),
##   * NDEEP points drawn uniformly from C_gamma by rejection (not from its bounding box),
##   * eight Nelder-Mead restarts from the best eight,
## using the IDENTICAL data, U, T_obs and theta-hat as 13_run.R (everything is seeded), so the two
## numbers are directly comparable dataset by dataset.
##
## env: EXD_W (workers), EXD_REPS (datasets per arm), EXD_NPTS, EXD_NMS

EX <- "/path/to/project/paper_deepgof/theory/exact_mmc"
source(file.path(EX, "10_core.R"))

NW    <- as.integer(Sys.getenv("EXD_W", "6"))
NREPS <- as.integer(Sys.getenv("EXD_REPS", "12"))
NDEEP <- as.integer(Sys.getenv("EXD_NPTS", "300"))
NMS   <- as.integer(Sys.getenv("EXD_NMS", "8"))
DBUDGET <- as.numeric(Sys.getenv("EXD_BUDGET", "1200"))
ROWS  <- file.path(EX, "deep_rows")
DCSV  <- file.path(EX, "deepmax.csv")
PROG  <- file.path(EX, "deep_progress.log")
dir.create(ROWS, recursive = TRUE, showWarnings = FALSE)

hbd <- function(msg) cat(sprintf("%s pid%-6d %s\n", format(Sys.time(), "%H:%M:%S"), Sys.getpid(), msg),
                         file = PROG, append = TRUE)

## uniform-by-rejection sample of C_gamma inside the truncation box
sample_set <- function(y, mm, bx, want, seed, chunk = 20000L, maxdraw = 600000L) {
  L <- bern_L(); Mj <- apply(abs(mm), 2, max)
  set.seed(seed)
  acc <- NULL; drawn <- 0L
  while (drawn < maxdraw && (is.null(acc) || nrow(acc) < want)) {
    Th <- cbind(runif(chunk, bx$lo[1], bx$hi[1]), runif(chunk, bx$lo[2], bx$hi[2]),
                runif(chunk, bx$lo[3], bx$hi[3]))
    ok <- which(in_set_slack(Th, y, mm, L, Mj) <= 1)
    if (length(ok)) acc <- rbind(acc, Th[ok, , drop = FALSE])
    drawn <- drawn + chunk
  }
  if (is.null(acc)) return(NULL)
  if (nrow(acc) > want) acc <- acc[sample.int(nrow(acc), want), , drop = FALSE]
  list(pts = acc, drawn = drawn)
}

## the six coordinate extremes of C_gamma (feasible points, on the faces of the bounding box)
extremes <- function(y, mm, thetahat, se, bx) {
  L <- bern_L(); Mj <- apply(abs(mm), 2, max)
  clo <- thetahat - RCAP * se; chi <- thetahat + RCAP * se
  pen <- function(th, k, s) {
    ex <- max(c(0, (clo - th) / se, (th - chi) / se))
    if (ex > 0) return(1e6 + 1e3 * ex)
    sl <- in_set_slack(matrix(th, 1L), y, mm, L, Mj)
    if (sl > 1) return(1e3 + 1e3 * (sl - 1))
    -s * th[k]
  }
  out <- list()
  st0 <- if (!is.null(bx$mc_in)) bx$mc_in else matrix(thetahat, 1L)
  for (k in 1:3) for (s in c(-1, 1)) {
    st <- st0[which.max(s * st0[, k]), ]
    o <- tryCatch(optim(st, pen, k = k, s = s, method = "Nelder-Mead",
                        control = list(maxit = 400, reltol = 1e-10, parscale = pmax(se, 1e-8))),
                  error = function(e) NULL)
    if (!is.null(o) && o$value < 0) out[[length(out) + 1L]] <- as.numeric(o$par)
  }
  out
}

deep_one <- function(j) {
  kind <- j$k; i <- j$i
  t0 <- Sys.time(); hbd(sprintf("START kind=%d rep=%d", kind, i))
  g <- tryCatch(gen_dataset(kind, i), error = function(e) NULL)
  if (is.null(g) || length(unique(g$data$y)) < 2L) return(NULL)
  P <- prep_D1(g$data); y <- as.numeric(g$data$y)
  fo <- fast_fit(y, P$mm); if (is.null(fo) || is.null(fo$ph)) return(NULL)
  th <- as.numeric(fo$coef)
  Tobs <- score_one(fast_map(y, fo$ph, P))
  set.seed(SEED_U + kind * 300000L + i); U <- matrix(runif(P$n * BREP), P$n, BREP)
  bx <- set_box(y, P$mm, th, fo$ph, SEED_BOX + kind * 300000L + i)

  L <- bern_L(); Mj <- apply(abs(P$mm), 2, max)
  st <- new.env(parent = emptyenv()); st$bc <- -1L; st$best <- th; st$ne <- 0L; st$nr <- 0L
  st$trunc <- 0L; tstop <- Sys.time() + DBUDGET
  ev <- function(x) {
    if (st$trunc || Sys.time() > tstop) { st$trunc <- 1L; return((1 + max(st$bc, 0L)) / (BREP + 1)) }
    r <- pB(x, Tobs, U, P, floor_count = st$bc)
    st$ne <- st$ne + 1L; st$nr <- st$nr + r$nrep
    if (r$full && r$cnt > st$bc) { st$bc <- r$cnt; st$best <- x }
    r$p
  }
  r0 <- pB(th, Tobs, U, P, floor_count = -1L)
  st$bc <- r0$cnt; st$best <- th; st$ne <- 1L; st$nr <- r0$nrep
  p_hat <- r0$p
  ex <- extremes(y, P$mm, th, bx$se, bx)
  for (e in ex) ev(e)
  cnt_ext <- st$bc
  S <- sample_set(y, P$mm, bx, NDEEP, 991000L + kind * 300000L + i)
  pp <- numeric(0)
  if (!is.null(S)) { pp <- numeric(nrow(S$pts)); for (q in seq_len(nrow(S$pts))) pp[q] <- ev(S$pts[q, ]) }
  cnt_rej <- st$bc
  negp <- function(x) {
    s <- in_set_slack(matrix(x, 1L), y, P$mm, L, Mj)
    if (s > 1) return(1 + (s - 1))
    -ev(x)
  }
  psc <- pmax(bx$half / 20, 1e-8)
  if (length(pp)) {
    top <- order(pp, decreasing = TRUE)[seq_len(min(NMS, length(pp)))]
    for (q in top) tryCatch(optim(S$pts[q, ], negp, method = "Nelder-Mead",
                                  control = list(maxit = 120, reltol = 1e-12, parscale = psc)),
                            error = function(e) NULL)
  }
  row <- c(kind = kind, rep = i, Tobs = Tobs,
           p_hat = p_hat, p_ext = (1 + cnt_ext) / (BREP + 1),
           p_rej = (1 + cnt_rej) / (BREP + 1), p_deep = (1 + st$bc) / (BREP + 1),
           n_eval = st$ne, nrep_tot = st$nr, n_pts = if (is.null(S)) 0L else nrow(S$pts),
           n_drawn = if (is.null(S)) 0L else S$drawn,
           accept = if (is.null(S)) 0 else nrow(S$pts) / S$drawn,
           sup_0 = st$best[1], sup_1 = st$best[2], sup_2 = st$best[3],
           box_unbdd = bx$unbdd, trunc = st$trunc,
           secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  df <- as.data.frame(as.list(row), check.names = FALSE); names(df) <- sub("\\..*$", "", names(df))
  rf <- file.path(ROWS, sprintf("k%d_r%d.csv", kind, i))
  tf <- sprintf("%s.tmp%d", rf, Sys.getpid())
  utils::write.csv(df, tf, row.names = FALSE); if (file.exists(rf)) unlink(rf); file.rename(tf, rf)
  hbd(sprintf("DONE  kind=%d rep=%d  p_hat %.3f p_ext %.3f p_rej %.3f p_deep %.3f  ne %d  %.0f s",
              kind, i, p_hat, (1 + cnt_ext)/(BREP+1), (1 + cnt_rej)/(BREP+1),
              (1 + st$bc)/(BREP+1), st$ne, row[["secs"]]))
  row
}

done <- sub("\\.csv$", "", list.files(ROWS, pattern = "^k\\d+_r\\d+\\.csv$"))
JOBS <- Filter(function(j) !(sprintf("k%d_r%d", j$k, j$i) %in% done),
               unlist(lapply(seq_len(NREPS), function(i) list(list(k = 0L, i = i), list(k = 2L, i = i))),
                      recursive = FALSE))
cat(sprintf("deep-max audit | %d reps/arm | %d rejection points + 6 extremes + %d NM | %d workers\n",
            NREPS, NDEEP, NMS, NW))
cat(sprintf("resume: %d done, %d to run\n", length(done), length(JOBS))); flush.console()
if (length(JOBS)) {
  cl <- parallel::makeCluster(NW)
  parallel::clusterExport(cl, c("EX", "ROWS", "PROG", "NDEEP", "NMS", "hbd", "sample_set",
                                "extremes", "deep_one", "DBUDGET"), envir = environment())
  invisible(parallel::clusterEvalQ(cl, { source(file.path(EX, "10_core.R")); NULL }))
  invisible(parallel::parLapplyLB(cl, JOBS, deep_one))
  try(parallel::stopCluster(cl), silent = TRUE)
}
fs <- list.files(ROWS, pattern = "^k\\d+_r\\d+\\.csv$", full.names = TRUE)
A <- if (length(fs)) do.call(rbind, lapply(fs, function(f) tryCatch(utils::read.csv(f), error = function(e) NULL))) else NULL
if (!is.null(A)) { A <- A[order(A$kind, A$rep), ]; utils::write.csv(A, DCSV, row.names = FALSE)
  cat(sprintf("\n%d deep rows; median p_deep %.3f, median n_eval %.0f, median secs %.0f\n",
              nrow(A), median(A$p_deep), median(A$n_eval), median(A$secs))) }
cat("DEEPMAX COMPLETE\n")
