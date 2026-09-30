# run_smalln_link.R -- DeepGOF-1 at very small n and under link misfit.
# Two covariates x1, x2 ~ N(0, 1); linear predictor eta = 0.25 + 1.5 x1 + 1.0 x2 (about half events).
#   nulls: y ~ Bern(plogis(eta)), n in {20, 30, 50, 100}; 2,000 per n at 20 and 30, 1,000 at 50 and 100
#   links: y ~ Bern(F(eta)) with F = probit, cloglog, Stukel heavier tails (a1 = a2 = -1) and Stukel
#          asymmetric (a1 = 0.5, a2 = -1); 1,000 per link and n
# The fitted model is the logistic regression y ~ x1 + x2 with a plain glm(), as an analyst fits it.
# Tests on every dataset: deepgof1() (with two covariates the axis rule and the all-pairs reading are the
# same map), Hosmer-Lemeshow, Stukel (joint score test and the SstBoth form), GiViTI (givitiR, as in the
# paper) and the projection test of Liu et al. (2024), B = 199. Recorded per dataset: events, whether the
# fit separated (not converged or some |b| > 8), and every p-value (NA when a test returns none).
# ebrahim.gof 2.9.0 from the private library. Output: rows/<cell>_r<first>.csv, one file per 50 datasets.
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/smalln_link"
LIB  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/shipped_level/lib290c"
SMOKE <- length(commandArgs(TRUE)) > 0 && commandArgs(TRUE)[1] == "smoke"
OUT <- file.path(HERE, if (SMOKE) "rows_smoke" else "rows"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
B <- 199L

NS <- c(20L, 30L, 50L, 100L)
LINKS <- c("probit", "cloglog", "stukel_heavy", "stukel_asym")
CELLS <- rbind(data.frame(link = "logit", n = NS, reps = c(2000L, 2000L, 1000L, 1000L)),
               expand.grid(link = LINKS, n = NS, reps = 1000L, stringsAsFactors = FALSE))
CELLS$id <- seq_len(nrow(CELLS))
if (SMOKE) CELLS$reps <- 2L
jobs <- list()
for (k in seq_len(nrow(CELLS))) for (ch in split(seq_len(CELLS$reps[k]), ceiling(seq_len(CELLS$reps[k]) / 50))) {
  f <- file.path(OUT, sprintf("c%02d_%s_n%03d_r%04d.csv", CELLS$id[k], CELLS$link[k], CELLS$n[k], ch[1]))
  if (!file.exists(f)) jobs[[length(jobs) + 1]] <- list(cc = CELLS[k, ], reps = ch, f = f)
}

one_chunk <- function(j) {
  suppressMessages(library(ebrahim.gof, lib.loc = LIB))
  stopifnot(utils::packageVersion("ebrahim.gof", lib.loc = LIB) == "2.9.0")
  ## Stukel's (1988) generalized logistic: p = plogis(h(eta)); a < 0 gives heavier tails, a > 0 lighter
  h1 <- function(e, a) if (a > 0) (exp(a * e) - 1) / a else if (a < 0) -log(1 - a * e) / a else e
  stukel_p <- function(eta, a1, a2) {
    h <- ifelse(eta >= 0, h1(pmax(eta, 0), a1), -h1(pmax(-eta, 0), a2))
    stats::plogis(h)
  }
  prob <- function(link, eta) switch(link,
    logit = stats::plogis(eta), probit = stats::pnorm(eta), cloglog = 1 - exp(-exp(eta)),
    stukel_heavy = stukel_p(eta, -1, -1), stukel_asym = stukel_p(eta, 0.5, -1))
  ## GiViTI in a long-lived callr session, as in the paper's add-on arms: a crash in givitiR's compiled
  ## code kills the session, not the worker; 120 s timeout; NA with a status otherwise
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
    s <- 6000000L + j$cc$id * 100000L + r
    set.seed(s)
    n <- j$cc$n; x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    eta <- 0.25 + 1.5 * x1 + 1.0 * x2
    d <- data.frame(x1, x2, y = stats::rbinom(n, 1, prob(j$cc$link, eta)))
    base <- data.frame(id = j$cc$id, link = j$cc$link, n = n, rep = r, events = sum(d$y))
    if (length(unique(d$y)) < 2L)                   # no model can be fitted: every test is undefined
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

## a second runner can take the chunks from the other end (SL_REV=1, SL_WORKERS=n); a chunk both reach is
## computed twice with the same seeds and written atomically, so the result is the same file
if (Sys.getenv("SL_REV") == "1") jobs <- rev(jobs)
cat(sprintf("smalln_link%s | %d chunks | %s\n", if (SMOKE) " SMOKE" else "", length(jobs), format(Sys.time())))
if (length(jobs)) {
  cl <- parallel::makeCluster(if (SMOKE) 4L else as.integer(Sys.getenv("SL_WORKERS", "16")))
  parallel::clusterExport(cl, c("LIB", "B", "one_chunk"))
  invisible(parallel::parLapplyLB(cl, jobs, one_chunk))
  parallel::stopCluster(cl)
}
R <- do.call(rbind, lapply(list.files(OUT, "\\.csv$", full.names = TRUE), utils::read.csv))
cat("rows:", nrow(R), "\n"); print(colSums(is.na(R[, grep("^p_", names(R))])))
if (SMOKE) print(R)
