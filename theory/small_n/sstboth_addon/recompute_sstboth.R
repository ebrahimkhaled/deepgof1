## NOTE (archive): this script sources sim_defs.R of the two corrupted-record studies (gen_clin, corrupt,
## fit_pol, p_classical). That file is covered by the studies' declaration hashes (Online Resource 1,
## Section 19) and is available from the author together with the declarations; it is not in this archive.
# =============================================================================
# recompute_sstboth.R -- the Stukel arm of the two corrupted-record studies in the SstBoth form
# (Liu et al. 2024: the sum of the squares of the two marginal Stukel score z-statistics, referred to
# chi-square(2); ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))), on EXACTLY the stored datasets.
#
# Datasets are regenerated from the seeds of run_small_n.R (Study 1, base 81000) and run_trim.R
# (Study 2, base 91000): set.seed, gen_clin, sample.int, corrupt -- the draws that precede every
# bootstrap, so DeepGOF-1 need not be rerun. Pairing is verified per replicate by recomputing the
# stored Hosmer-Lemeshow and joint-form Stukel p-values (full data, and for Study 1 the LOO maximum;
# for Study 2 the trimmed index and the trimmed p-values).
#
# SstBoth is undefined (NaN) when one Stukel column is identically zero, i.e. no fitted risk on one
# side of .5; as in the package tables it counts as non-rejection (p = 1), and the event is recorded.
# A failed fit gives p = 1, as in the original scripts.
#
# One job = (study, cell, block of 25 replicates), written atomically to rows_s1/ or rows_s2/; rerun
# skips finished blocks. Run:  Rscript recompute_sstboth.R [workers=4]
# =============================================================================
HERE <- "/path/to/project/paper_deepgof/theory/small_n"
OUT  <- file.path(HERE, "sstboth_addon")
args <- commandArgs(trailingOnly = TRUE)
NW <- if (length(args)) as.integer(args[1]) else 4L
CHUNK <- 25L

CELLS <- rbind(
  expand.grid(n = c(30L, 50L), C = 0, type = c("none", "decimal", "sign", "outcome"),
              reps = 1000L, stringsAsFactors = FALSE),
  expand.grid(n = c(30L, 50L), C = c(0.5, 1.0, 1.5), type = c("none", "decimal"),
              reps = 300L, stringsAsFactors = FALSE))
CELLS$cell <- seq_len(nrow(CELLS))

jobs <- list()
for (study in 1:2) {
  rd <- file.path(OUT, sprintf("rows_s%d", study)); dir.create(rd, showWarnings = FALSE)
  for (k in seq_len(nrow(CELLS))) {
    cc <- CELLS[k, ]
    for (s in seq(1L, cc$reps, by = CHUNK)) {
      f <- file.path(rd, sprintf("c%02d_r%04d.csv", cc$cell, s))
      if (!file.exists(f))
        jobs[[length(jobs) + 1L]] <- list(study = study, cc = cc, reps = s:min(cc$reps, s + CHUNK - 1L), f = f)
    }
  }
}
if (length(args) > 1) {                         # smoke: the first job of each study only
  pick <- c(which(vapply(jobs, function(j) j$study == 1, TRUE))[1], which(vapply(jobs, function(j) j$study == 2, TRUE))[1])
  jobs <- jobs[pick[!is.na(pick)]]
}
cat(sprintf("SstBoth add-on | %d jobs to run | %d workers | %s\n", length(jobs), NW, format(Sys.time())))
flush.console()
if (!length(jobs)) { cat("nothing to do\n"); quit(save = "no") }

worker_init <- function(HERE) {
  source(file.path(HERE, "sim_defs.R"))          # gen_clin, corrupt, fit_pol, p_classical (unchanged)
  CUT <<- 4.5
  trim_index <<- function(d) {                   # verbatim from run_trim.R
    o <- pmax(abs(d$x1 - stats::median(d$x1)) / stats::mad(d$x1),
              abs(d$x2 - stats::median(d$x2)) / stats::mad(d$x2))
    if (max(o) > CUT) which.max(o) else 0L
  }
  ## both Stukel forms from the package, on the context run.all.gof builds; NaN/NA recorded as undefined
  st2 <<- function(fit) {
    ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
    pv <- function(z) { v <- z$p_value; if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
    c(j = tryCatch(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "joint"))), error = function(e) NA_real_),
      m = tryCatch(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))), error = function(e) NA_real_),
      side = as.numeric(any(ctx$ph >= 0.5) && any(ctx$ph < 0.5)))
  }
  NULL
}

run_job <- function(jb) {
  cc <- jb$cc
  R <- vector("list", length(jb$reps))
  for (ii in seq_along(jb$reps)) {
    r <- jb$reps[ii]
    set.seed((if (jb$study == 1) 81000L else 91000L) + cc$cell * 100000L + r)
    d <- gen_clin(cc$n, cc$C)
    j <- sample.int(cc$n, 1L)
    if (cc$type != "none") d <- corrupt(d, cc$type, j)
    fit <- fit_pol(d)
    if (is.null(fit)) { pc <- c(HL = 1, Stukel = 1); s <- c(j = 1, m = 1, side = NA) }
    else { pc <- p_classical(fit, d); s <- st2(fit) }
    row <- data.frame(cell = cc$cell, rep = r, corrupted = if (cc$type == "none") NA_integer_ else j,
                      events = sum(d$y), chk_p_HL = pc[["HL"]], chk_p_Stj = pc[["Stukel"]],
                      p_Stj = s[["j"]], p_Stm = s[["m"]], full_undef = as.integer(is.na(s[["m"]])),
                      full_twoside = s[["side"]])
    if (jb$study == 1) {
      lj <- lm_ <- numeric(cc$n); nund <- 0L
      for (i in seq_len(cc$n)) {
        fi <- fit_pol(d[-i, , drop = FALSE])
        if (is.null(fi)) { lj[i] <- 1; lm_[i] <- 1; next }
        si <- st2(fi)
        lj[i] <- si[["j"]]; lm_[i] <- si[["m"]]
        if (is.na(si[["m"]])) nund <- nund + 1L
      }
      lj[is.na(lj)] <- 1; lm_[is.na(lm_)] <- 1               # as all_p: loo[is.na(loo)] <- 1
      row$chk_loo_Stj <- max(lj)
      row$loo_Stm <- max(lm_); row$nres_Stm <- sum(lm_ > 0.05); row$amax_Stm <- which.max(lm_)
      row$loo_undef <- nund
    } else {
      k <- trim_index(d); row$trimmed <- k
      if (k == 0L) { row$chk_t_Stj <- row$p_Stj; row$t_Stm <- row$p_Stm; row$t_undef <- row$full_undef
                     row$chk_t_HL <- row$chk_p_HL }
      else {
        dk <- d[-k, , drop = FALSE]; fk <- fit_pol(dk)
        if (is.null(fk)) { row$chk_t_Stj <- 1; row$t_Stm <- 1; row$t_undef <- 0L; row$chk_t_HL <- 1 }
        else { sk <- st2(fk); pk <- p_classical(fk, dk)
               row$chk_t_Stj <- pk[["Stukel"]]; row$t_Stm <- sk[["m"]]; row$t_undef <- as.integer(is.na(sk[["m"]]))
               row$chk_t_HL <- pk[["HL"]] }
      }
    }
    R[[ii]] <- row
  }
  R <- do.call(rbind, R)
  tmp <- paste0(jb$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, jb$f)
  sprintf("%s  study %d cell %d reps %d-%d", format(Sys.time()), jb$study, cc$cell, min(jb$reps), max(jb$reps))
}

if (NW <= 1L) {
  worker_init(HERE)
  for (jb in jobs) cat(run_job(jb), "\n")
} else {
  cl <- parallel::makeCluster(min(NW, length(jobs)))
  parallel::clusterCall(cl, worker_init, HERE)
  t0 <- Sys.time()
  msg <- parallel::parLapplyLB(cl, jobs, run_job)
  parallel::stopCluster(cl)
  cat(paste(unlist(msg), collapse = "\n"), "\n")
  cat(sprintf("DONE %s: %.2f h\n", format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}
