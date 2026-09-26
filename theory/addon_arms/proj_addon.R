# =============================================================================
# proj_addon.R -- the projection test with the weight of Liu et al. (2024) (proj_liu.R), on EXACTLY
# the datasets AND the bootstrap draws of the stored proj arm, for the Liu grid and the niche.
#
# Same data. Datasets are regenerated from their own seeds with the generator code of the add-on
# scripts (liu_defs parsed out of liu_addon.R, niche_defs and SC parsed out of niche_addon.R; both
# verbatim copies of the originals, verified in verify_pairing.R).
#
# Same bootstrap draws. In the stored runs the proj bootstrap started from the RNG state left by
# the arm computed just before it on the same dataset:
#   Liu grid (code/liu_augment.R, one()):  set.seed -> data -> fit -> Stute-Zhu (run.all.gof,
#            200 replicates, each rbinom(n, 1, p): 200 n uniforms) -> proj bootstrap.
#   niche    (theory/niche/run_niche.R, one_rep()):  set.seed -> data -> fit -> the closed-form
#            arms (no random numbers) -> DeepGOF-1 bootstrap (199 replicates of rbinom(n, 1, p):
#            199 n uniforms, only when the fit exists) -> proj bootstrap.
# rbinom(n, 1, p) with 0 < p < 1 uses exactly one uniform per draw, whatever p is, so the runner
# replays those draws with runif() (AD_REPLAY = "draws", default) instead of recomputing the arm.
# AD_REPLAY = "arm" recomputes Stute-Zhu through run.all.gof instead (Liu grid only; slower), to
# show that both give the same state.
#
# Per dataset the bootstrap is run ONCE: the same y* and the same refits feed the quadratic forms,
#   proj_liu_raw  the weight of the definition ((pi - arccos) / 2 pi, zero-vector rule) on the raw
#                 covariates: Liu et al.'s statistic exactly -- the PRIMARY proj arm
#   proj_liu_std  the same weight on standardised covariates (the deployed arm's scaling) -- secondary
#   proj_old  the deployed weight (NN * (pi - arccos)), recomputed -> must equal the stored proj.
#             Computed on every dataset with n <= 100 and on reps <= AD_OLDREP (25) otherwise
#             (in R it is O(n^3) and costs about 15 s at n = 500); NA elsewhere.
# So the two arms are paired on data and on bootstrap draws, and the rows carrying proj_old are the
# pairing check (analyze_proj.R joins proj_old to the stored proj column).
#
# Resumable, every dataset on disk as soon as it is done: a chunk's rows are appended to
# <chunk>.part one dataset at a time; a finished chunk is renamed to <chunk>.csv. A rerun skips
# finished chunks and, inside a .part, the datasets already written.
#
# Run:   Rscript proj_addon.R
# Env:   AD_BENCH   liu (default) / niche
#        AD_W       workers (4)
#        AD_BLOCKS  (liu) comma list of blocks, e.g. S2_n100_s101000; default all 24
#        AD_SCEN    (niche) comma list of scenarios; default all 14
#        AD_KIND1   (liu) 0 (default) / 1: also the 1,000 kind-1 training departures
#        AD_CHUNK   datasets per chunk (liu 25, niche 20 = the stored row blocks)
#        AD_MAXREP  cap on rep per kind / scenario (smoke and timing runs)
#        AD_REPLAY  draws (default) / arm
#        AD_OLDREP  reps (n > 100) that also recompute the deployed weight (25)
#        AD_OUT     output folder (default proj_rows_liu / proj_rows_niche)
# =============================================================================
HERE  <- "/path/to/project/paper_deepgof/theory/addon_arms"
RES   <- "/path/to/project/paper_deepgof/release/deepgof1/results/liu"
BENCH <- Sys.getenv("AD_BENCH", "liu"); stopifnot(BENCH %in% c("liu", "niche"))
NW    <- as.integer(Sys.getenv("AD_W", "4"))
KIND1 <- Sys.getenv("AD_KIND1", "0") == "1"
CHUNK <- as.integer(Sys.getenv("AD_CHUNK", if (BENCH == "liu") "25" else "20"))
MAXR  <- as.integer(Sys.getenv("AD_MAXREP", "100000"))
REPLAY<- Sys.getenv("AD_REPLAY", "draws"); stopifnot(REPLAY %in% c("draws", "arm"))
OUTD  <- file.path(HERE, Sys.getenv("AD_OUT", paste0("proj_rows_", BENCH))); dir.create(OUTD, showWarnings = FALSE)
OLDREP<- as.integer(Sys.getenv("AD_OLDREP", "25"))
B     <- 199L                                   # the stored proj arm: B = 199 in both benchmarks
CS    <- c(0, 0.25, 0.5, 0.75, 1.0)

grab <- function(file, names) {                 # evaluate only the named top-level assignments
  ex <- parse(file)
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                    ((is.name(e[[2]]) && as.character(e[[2]]) %in% names) ||
                     (is.call(e[[2]]) && identical(e[[2]][[2]], as.name("SC")) && "SC" %in% names)))
    eval(e, globalenv())
}
if (BENCH == "liu") grab(file.path(HERE, "liu_addon.R"), "liu_defs") else {
  grab(file.path(HERE, "niche_addon.R"), c("CH", "SEEDB", "NNULL", "NALT", "SC", "niche_defs"))
  if (nzchar(Sys.getenv("AD_SCEN"))) SC <- SC[SC$scen %in% strsplit(Sys.getenv("AD_SCEN"), ",")[[1]], ]
}

proj_defs <- quote({
  suppressMessages(library(brglm2))
  source(file.path(HERE, "proj_liu.R"))
  na_row <- c(proj_old = NA, proj_liu_raw = NA, proj_liu_std = NA, T_old = NA, T_liu_raw = NA, T_liu_std = NA,
              nfail = NA, secs_A = NA, secs = NA)
  ## with_old: also recompute the deployed weight (pairing check). It is O(n^3) in R, so by default
  ## it is done on every dataset with n <= 100 and on reps <= OLDREP otherwise.
  run_proj <- function(d, frm, fit, with_old) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch({
      X  <- stats::model.matrix(fit)[, -1, drop = FALSE]
      AL <- list(raw = proj_A_liu_fast(X, std = FALSE), std = proj_A_liu_fast(X, std = TRUE))
      if (with_old) AL$old <- proj_A_nn(X)
      ta <- proc.time()[["elapsed"]] - t0
      o  <- proj_p_multi(d, frm, fit, B, AL)
      c(proj_old = if (with_old) o$p[["old"]] else NA, proj_liu_raw = o$p[["raw"]], proj_liu_std = o$p[["std"]],
        T_old = if (with_old) o$T[["old"]] else NA, T_liu_raw = o$T[["raw"]], T_liu_std = o$T[["std"]],
        nfail = o$nfail, secs_A = ta)
    }, error = function(e) NULL)
    if (is.null(r)) { r <- na_row[-length(na_row)] }
    c(r, secs = proc.time()[["elapsed"]] - t0)
  }
  ## ---- Liu grid: the data and the RNG path of code/liu_augment.R::one() ----
  liu_proj_one <- function(SET, NOBS, SEED0, kind, i, C, ci) {
    s <- SEED0 + SET*500000L + kind*50000L + ci*5000L + i
    set.seed(s)
    d <- if (kind == 1L) gen_train_alt(NOBS, SET, SEED0*29L + i) else gen_liu(NOBS, SET, if (kind==0L) 0 else C)
    key <- data.frame(setting = SET, n = NOBS, seed0 = SEED0, kind = kind, rep = i, C = C, ci = ci, dseed = s)
    if (length(unique(d$y)) < 2L) return(cbind(key, t(na_row), status = "dropped_degenerate_y"))
    frm <- FRM(d); fit <- proj_fit_pol(frm, d)
    if (is.null(fit)) return(cbind(key, t(na_row), status = "dropped_no_fit"))
    if (REPLAY == "draws") invisible(stats::runif(200L * NOBS)) else
      invisible(tryCatch(suppressWarnings(ebrahim.gof::run.all.gof(fit, tests = "Stute-Zhu", include_slow = TRUE,
                                                                   install = "no")), error = function(e) NULL))
    cbind(key, t(run_proj(d, frm, fit, NOBS <= 100L || i <= OLDREP)), status = "ok")
  }
  ## ---- niche: the data and the RNG path of theory/niche/run_niche.R::one_rep() ----
  niche_proj_one <- function(sc, r) {
    s <- SEEDB + sc$id * 100000L + r
    set.seed(s); dat <- gen(sc); fit <- fit_pol(dat)
    key <- data.frame(scen = sc$scen, fam = sc$fam, C = sc$C, n = sc$n, rep = r, dseed = s)
    if (!is.null(fit)) invisible(stats::runif(199L * sc$n))       # the DeepGOF-1 bootstrap draws
    f2 <- frm(dat); fp <- proj_fit_pol(f2, dat)
    if (is.null(fp)) return(cbind(key, t(na_row), status = "no_proj_fit"))
    cbind(key, t(run_proj(dat, f2, fp, sc$n <= 100L || r <= OLDREP)), status = "ok")
  }
  ## append one dataset's row to the chunk's .part; rename when the chunk is complete
  run_chunk <- function(ch) {
    part <- paste0(ch$f, ".part"); have <- integer(0)
    if (file.exists(part)) have <- tryCatch(utils::read.csv(part)$rep, error = function(e) integer(0))
    for (i in setdiff(ch$reps, have)) {
      row <- if (BENCH == "liu") liu_proj_one(ch$SET, ch$NOBS, ch$SEED0, ch$kind, i, ch$C, ch$ci) else
        niche_proj_one(ch$sc, i)
      utils::write.table(row, part, sep = ",", row.names = FALSE, col.names = !file.exists(part),
                         append = file.exists(part), qmethod = "double")
    }
    R <- utils::read.csv(part)
    if (all(ch$reps %in% R$rep)) {
      R <- R[!duplicated(R$rep), ]; R <- R[order(R$rep), ]
      utils::write.csv(R, part, row.names = FALSE)
      if (file.exists(ch$f)) file.remove(ch$f); file.rename(part, ch$f)
    }
    NULL
  }
})

## ---- the job list ----
chunks <- list()
if (BENCH == "liu") {
  BLOCKS <- if (nzchar(Sys.getenv("AD_BLOCKS"))) strsplit(Sys.getenv("AD_BLOCKS"), ",")[[1]] else
    sort(basename(list.dirs(RES, recursive = FALSE)))
  for (b in BLOCKS) {
    m <- regmatches(b, regexec("^S([0-9])_n([0-9]+)_s([0-9]+)$", b))[[1]]
    meta <- utils::read.csv(file.path(RES, b, "blk_meta.csv"))
    bp <- list(SET = as.integer(m[2]), NOBS = as.integer(m[3]), SEED0 = as.integer(m[4]))
    bd <- file.path(OUTD, b); dir.create(bd, showWarnings = FALSE)
    specs <- list(list(kind = 0L, ci = 0L, C = NA_real_, R = meta$nnull[1]))
    if (KIND1) specs <- c(specs, list(list(kind = 1L, ci = 0L, C = NA_real_, R = meta$nrand[1])))
    for (ci in seq_along(CS)) if (CS[ci] > 0)
      specs <- c(specs, list(list(kind = 2L, ci = ci, C = CS[ci], R = meta$neval[1])))
    for (sp in specs) {
      R <- min(sp$R, MAXR)
      for (s in seq(1L, R, by = CHUNK)) {
        f <- file.path(bd, sprintf("chunk_k%d_c%d_r%04d.csv", sp$kind, sp$ci, s))
        if (!file.exists(f))
          chunks[[length(chunks) + 1L]] <- c(bp, sp[c("kind", "ci", "C")], list(reps = s:min(R, s + CHUNK - 1L), f = f,
                                             cost = bp$NOBS))
      }
    }
  }
} else {
  for (i in seq_len(nrow(SC))) {
    sc <- SC[i, ]; R <- min(sc$reps, MAXR)
    for (s in seq(1L, R, by = CHUNK)) {
      f <- file.path(OUTD, sprintf("%s_r%04d.csv", sc$scen, s))
      if (!file.exists(f)) chunks[[length(chunks) + 1L]] <- list(sc = sc, reps = s:min(R, s + CHUNK - 1L), f = f,
                                                                 cost = sc$n)
    }
  }
}
## largest n first, so the long chunks do not straggle at the end
chunks <- chunks[order(-vapply(chunks, function(ch) ch$cost, 1))]
cat(sprintf("proj_addon %s | %d chunks to run | replay %s | %d workers | %s\n", BENCH, length(chunks), REPLAY, NW,
            format(Sys.time())))
flush.console()

if (length(chunks)) {
  t0 <- Sys.time()
  cl <- parallel::makeCluster(min(NW, length(chunks)))
  exp <- c("HERE", "proj_defs", "B", "BENCH", "REPLAY", "OLDREP")
  exp <- c(exp, if (BENCH == "liu") "liu_defs" else c("niche_defs", "SEEDB"))
  parallel::clusterExport(cl, exp)
  invisible(parallel::clusterEvalQ(cl, {
    if (BENCH == "liu") eval(liu_defs) else eval(niche_defs)
    eval(proj_defs)
    if (REPLAY == "arm") suppressMessages(library(ebrahim.gof))
    NULL }))
  done <- 0L
  for (bt in split(chunks, ceiling(seq_along(chunks) / (2L * NW)))) {
    parallel::parLapplyLB(cl, bt, function(ch) run_chunk(ch))
    done <- done + length(bt)
    el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    cat(sprintf("%s  %d / %d chunks | %.2f h | ETA %s\n", format(Sys.time(), "%Y-%m-%d %H:%M"), done, length(chunks), el,
                format(Sys.time() + (length(chunks) - done) * el / done * 3600, "%Y-%m-%d %H:%M")))
    flush.console()
  }
  parallel::stopCluster(cl)
}
cat("proj_addon", BENCH, "DONE", format(Sys.time()), "\n")
