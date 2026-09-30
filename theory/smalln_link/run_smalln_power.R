# run_smalln_power.R -- power at very small n against misfit in the covariates, the companion of
# run_smalln_link.R (same design, same tests, and its null datasets give the matched-level thresholds).
# Two covariates x1, x2 ~ N(0, 1); eta = 0.25 + 1.5 x1 + 1.0 x2 + C f, with f = x1^2 - 1 (omitted quadratic)
# or f = x1 x2 (omitted interaction), C in {0.5, 1.0}; y ~ Bern(plogis(eta)); n in {20, 30, 50}; 1,000 per cell.
# The fitted model is y ~ x1 + x2 (plain glm). Tests: deepgof1() (two covariates: both readings coincide),
# Hosmer-Lemeshow, Stukel (joint and SstBoth), GiViTI (givitiR) and the projection test, B = 199.
# ebrahim.gof 2.9.0 from the private library. Output: rows_power/<cell>_r<first>.csv, 50 datasets per file.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/smalln_link"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290c"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_power_smoke" else "rows_power"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L

CELLS <- expand.grid(fam = c("quad", "inter"), C = c(0.5, 1.0), n = c(20L, 30L, 50L), reps = 1000L,
                     stringsAsFactors = FALSE)
CELLS$id <- 100L + seq_len(nrow(CELLS))            # ids apart from the link study's 1..20
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (k in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[k]), ceiling(seq_len(CELLS$reps[k]) / 50))) {
  f <- file.path(OUT, sprintf("c%03d_%s_C%.1f_n%03d_r%04d.csv", CELLS$id[k], CELLS$fam[k], CELLS$C[k], CELLS$n[k], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[k, ], reps = ch, f = f)
}
if (Sys.getenv("SP_REV") == "1") jobs <- rev(jobs)

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  gv_env <- new.env()
  giviti_p <- function(o, e) {
    if (is.null(gv_env$rs) || !identical(gv_env$rs$get_state(), "idle")) {
      if (!is.null(gv_env$rs)) try(gv_env$rs$kill(), silent = TRUE)
      gv_env$rs <- tryCatch({ r <- callr::r_session$new(wait = TRUE)
        r$run(function() { suppressMessages(requireNamespace("givitiR")); NULL }); r }, error = function(e) NULL)
      if (is.null(gv_env$rs)) return(NA_real_)
    }
    rs <- gv_env$rs
    rs$call(function(o, e) givitiR::givitiCalibrationTest(o = o, e = e, devel = "internal")$p.value,
            args = list(o = o, e = e))
    if (!identical(rs$poll_process(120000), "ready")) { try(rs$kill(), silent = TRUE); gv_env$rs <- NULL; return(NA_real_) }
    r <- tryCatch(rs$read(), error = function(e) NULL)
    if (is.null(r) || !is.null(r$error) || length(r$result) != 1 || !is.finite(r$result)) return(NA_real_)
    as.numeric(r$result)
  }
  pv <- function(z) { v <- tryCatch(z$p_value, error = function(e) NA); if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
  R <- do.call(rbind, lapply(j$reps, function(r) {
    s <- 7000000L + j$cc$id * 10000L + r
    set.seed(s)
    n <- j$cc$n; x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    f <- if (j$cc$fam == "quad") x1^2 - 1 else x1 * x2
    eta <- 0.25 + 1.5 * x1 + 1.0 * x2 + j$cc$C * f
    d <- data.frame(x1, x2, y = stats::rbinom(n, 1, stats::plogis(eta)))
    base <- data.frame(id = j$cc$id, fam = j$cc$fam, C = j$cc$C, n = n, rep = r, events = sum(d$y))
    if (length(unique(d$y)) < 2L)
      return(cbind(base, sep = NA, p_DG = NA, p_HL = NA, p_StukelJ = NA, p_StukelS = NA, p_GiViTI = NA, p_proj = NA))
    fit <- suppressWarnings(stats::glm(y ~ x1 + x2, data = d, family = stats::binomial()))
    sep <- as.integer(!fit$converged || any(abs(stats::coef(fit)[-1]) > 8))
    set.seed(s + 1L); pdg <- tryCatch(deepgof1(fit, B = B)$p.value, error = function(e) NA_real_)
    ctx <- tryCatch(ebrahim.gof:::.gof_context(fit, G = 10), error = function(e) NULL)
    phl <- psj <- pss <- pgv <- NA_real_
    if (!is.null(ctx)) {
      phl <- pv(tryCatch(ebrahim.gof:::gof_hl(ctx), error = function(e) NULL))
      psj <- pv(tryCatch(ebrahim.gof:::gof_stukel(ctx, list(form = "joint")), error = function(e) NULL))
      pss <- pv(tryCatch(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal")), error = function(e) NULL))
      pgv <- giviti_p(as.numeric(ctx$y), as.numeric(ctx$ph))
    }
    set.seed(s + 2L); ppj <- tryCatch(projection.gof(fit, B = B)$p.value, error = function(e) NA_real_)
    cbind(base, sep = sep, p_DG = pdg, p_HL = phl, p_StukelJ = psj, p_StukelS = pss, p_GiViTI = pgv, p_proj = ppj)
  }))
  if (!is.null(gv_env$rs)) try(gv_env$rs$close(), silent = TRUE)
  tmp <- paste0(j$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, j$f)
  j$f
}

cat(sprintf("smalln_power%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("SP_WORKERS", "12")))
  parallel::clusterExport(cl, c("LIB", "B", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
cat("done", format(Sys.time()), "\n")
