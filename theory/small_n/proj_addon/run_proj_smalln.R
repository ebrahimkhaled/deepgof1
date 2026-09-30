# run_proj_smalln.R -- the projection test of Liu et al. (2024) added to the two corrupted-record studies of
# Section 8 (referee M9: the robustness section omitted the recommended default test).
# The datasets are regenerated exactly from their seeds with the studies' own definitions (sim_defs.R):
#   study "small_n" (run_small_n.R, SEEDB 81000): the size cells (C = 0), four corruptions, n = 30, 50;
#   study "trim"    (run_trim.R,    SEEDB 91000): the size cells of "none" and "decimal", with and without
#                                                 the covariate-only deletion rule (cut 4.5).
# Each regenerated dataset is checked against the deposited row (events and the corrupted index).
# projection.gof(fit, B = 199) on the study's own fit (Firth when the plain fit separates). Its bootstrap
# has a seed of its own, drawn after the data, so the data stream is untouched.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/small_n"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT  <- file.path(HERE, "proj_addon", if (SMOKE) "rows_smoke" else "rows")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BP <- 199L; CUT <- 4.5

CELLS <- rbind(
  expand.grid(n = c(30L, 50L), C = 0, type = c("none", "decimal", "sign", "outcome"),
              reps = 1000L, stringsAsFactors = FALSE),
  expand.grid(n = c(30L, 50L), C = c(0.5, 1.0, 1.5), type = c("none", "decimal"),
              reps = 300L, stringsAsFactors = FALSE))
CELLS$cell <- seq_len(nrow(CELLS))
SIZE <- CELLS[CELLS$C == 0, ]
jobs <- list()
for (st in c("small_n", "trim")) for (k in seq_len(nrow(SIZE))) {
  cc <- SIZE[k, ]
  if (st == "trim" && !(cc$type %in% c("none", "decimal"))) next
  reps <- if (SMOKE) 1:3 else seq_len(cc$reps)
  for (ch in split(reps, ceiling(seq_along(reps) / 100))) {
    f <- file.path(OUT, sprintf("%s_c%02d_r%04d.csv", st, cc$cell, ch[1]))
    if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(st = st, cc = cc, reps = ch, f = f)
  }
}

one_chunk <- function(j) {
  source(file.path(HERE, "sim_defs.R"))                  # gen_clin, corrupt, fit_pol
  seedb <- if (j$st == "small_n") 81000L else 91000L
  rows <- lapply(j$reps, function(r) {
    set.seed(seedb + j$cc$cell * 100000L + r)
    d <- gen_clin(j$cc$n, j$cc$C)
    jj <- sample.int(j$cc$n, 1L)
    if (j$cc$type != "none") d <- corrupt(d, j$cc$type, jj)
    ## the deposited row for this replicate: files hold 25 replicates each
    fr <- file.path(HERE, if (j$st == "small_n") "rows" else "rows_trim",
                    sprintf("c%02d_r%04d.csv", j$cc$cell, 25L * ((r - 1L) %/% 25L) + 1L))
    dr <- utils::read.csv(fr); dr <- dr[dr$rep == r, ]
    match_ok <- nrow(dr) == 1L && dr$events == sum(d$y) &&
      (j$cc$type == "none" || dr$corrupted == jj)
    pp <- function(dd, s) {
      fit <- fit_pol(dd); if (is.null(fit)) return(1)
      set.seed(s)
      tryCatch(projection.gof(fit, B = BP)$p.value, error = function(e) NA_real_)
    }
    s0 <- 7000000L + (if (j$st == "trim") 500000L else 0L) + j$cc$cell * 10000L + r
    p_full <- pp(d, s0)
    k <- 0L; p_trim <- NA_real_
    if (j$st == "trim") {
      o <- pmax(abs(d$x1 - stats::median(d$x1)) / stats::mad(d$x1),
                abs(d$x2 - stats::median(d$x2)) / stats::mad(d$x2))
      k <- if (max(o) > CUT) which.max(o) else 0L
      p_trim <- if (k == 0L) p_full else pp(d[-k, , drop = FALSE], s0 + 5000L)
      match_ok <- match_ok && dr$trimmed == k
    }
    data.frame(study = j$st, cell = j$cc$cell, n = j$cc$n, type = j$cc$type, rep = r,
               match = match_ok, trimmed = k, p_proj = p_full, t_proj = p_trim)
  })
  R <- do.call(rbind, rows)
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("proj add-on%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else 10L)
  parallel::clusterExport(cl, c("HERE", "BP", "CUT", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
R <- do.call(rbind, lapply(list.files(OUT, "\\.csv$", full.names = TRUE), utils::read.csv))
cat("rows:", nrow(R), " regenerated exactly:", sum(R$match), " failed p:", sum(is.na(R$p_proj)), "\n")
A <- aggregate(cbind(full = p_proj <= .05, trimmed = t_proj <= .05) ~ study + n + type, R, mean, na.action = na.pass)
print(A, digits = 3)
utils::write.csv(A, file.path(HERE, "proj_addon", if (SMOKE) "SMOKE.csv" else "PROJ_SIZE.csv"), row.names = FALSE)
