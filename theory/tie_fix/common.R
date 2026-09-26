# common.R -- the two versions of deepgof1() compared in the tie-breaking study.
#   dg_new(): ebrahim.gof 2.8.0 (worktree feature-proj-bagoft), ties broken at random, once per call.
#   dg_old(): the 2.7.0 / release statistic, ties broken by row order (rank(ties.method = "first")).
# Both use the same shipped weights and the same bootstrap loop; only the tie-breaking differs.
WT <- "/path/to/ebrahim.gof"
PAP <- "/path/to/project/paper_deepgof"
OUTD <- file.path(PAP, "theory", "tie_fix")
suppressMessages(pkgload::load_all(WT, quiet = TRUE, export_all = TRUE, helpers = FALSE))

dg_new <- function(fit, B = 199L) deepgof1(fit, B = B)

dg_old <- function(fit, B = 199L, K = 6L) {
  M <- deepgof1_weights
  score <- function(m) .dg_score((m - M$mu) / M$sd, M)
  obs <- .dg_map(fit, K, NULL); Sobs <- score(obs$map)
  ph <- as.numeric(stats::fitted(fit)); dat <- stats::model.frame(fit); yn <- names(dat)[1]
  frm <- stats::formula(fit); Sb <- numeric(B)
  for (b in seq_len(B)) {
    dat[[yn]] <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(frm, data = dat, family = stats::binomial())),
                   error = function(e) NULL)
    Sb[b] <- if (is.null(fb)) Inf else score(.dg_map(fb, K, NULL)$map)
  }
  list(statistic = Sobs, p.value = (1 + sum(Sb >= Sobs)) / (B + 1), axes = obs$axes)
}

# resumable append: one row per finished unit
done_keys <- function(path) if (file.exists(path)) read.csv(path)$key else character(0)
append_row <- function(path, row) {
  write.table(row, path, sep = ",", row.names = FALSE, col.names = !file.exists(path),
              append = file.exists(path))
}
