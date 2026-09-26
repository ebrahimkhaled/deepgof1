# =============================================================================
# addon_common.R -- the three arm computations shared by liu_addon.R and niche_addon.R.
#
# Every arm is the INSTALLED ebrahim.gof 2.8.0 code, called on the same context object that
# run.all.gof() builds (fitted risks clamped to [1e-6, 1 - 1e-6], G = 10):
#   HL      ebrahim.gof:::gof_hl     -- recomputed only as the pairing check
#   Stukel  ebrahim.gof:::gof_stukel, form "joint" (the 2.8.0 default) and form "marginal"
#           (the 2.7.0-and-earlier LogisticDx SstBoth sum) -- both, for the Stukel-version check
#   mHL     ebrahim.gof:::gof_ftest  -- lm(deviance residual ~ decile-of-risk group), ANOVA F
#                                       (Liu et al. 2024, section 2.2; run.all.gof row "F-test")
#   GiViTI  givitiR::givitiCalibrationTest(o = y, e = clamped risk, devel = "internal")$p.value
#           -- exactly the call inside ebrahim.gof:::gof_giviti. gof_giviti starts a fresh callr
#           process for every dataset (about a second of start-up on Windows). Here each worker
#           keeps ONE callr::r_session alive and sends it every call: the same isolation (a crash
#           in givitiR's compiled code kills the session, not the worker), the same 120 s timeout,
#           without the start-up per dataset. A dead or timed-out session is replaced and that
#           dataset is recorded as NA with its status.
# The package is only loaded, never edited.
# =============================================================================
suppressMessages({ library(brglm2); library(ebrahim.gof) })
stopifnot(as.character(utils::packageVersion("ebrahim.gof")) == "2.8.0")

.GV <- new.env()
giviti_session <- function() {
  if (is.null(.GV$rs) || !identical(.GV$rs$get_state(), "idle")) {
    if (!is.null(.GV$rs)) try(.GV$rs$kill(), silent = TRUE)
    .GV$rs <- callr::r_session$new(wait = TRUE)
    .GV$rs$run(function() { suppressMessages(requireNamespace("givitiR")); NULL })
  }
  .GV$rs
}
giviti_close <- function() { if (!is.null(.GV$rs)) try(.GV$rs$close(), silent = TRUE); .GV$rs <- NULL }

## returns c(p, status) with status 0 ok, 1 givitiR error, 2 timeout, 3 session died
giviti_p <- function(o, e, devel = "internal", timeout_s = 120) {
  rs <- tryCatch(giviti_session(), error = function(err) NULL)
  if (is.null(rs)) return(c(p = NA_real_, status = 3))
  rs$call(function(o, e, devel) givitiR::givitiCalibrationTest(o = o, e = e, devel = devel)$p.value,
          args = list(o = o, e = e, devel = devel))
  st <- rs$poll_process(timeout_s * 1000)
  if (!identical(st, "ready")) { try(rs$kill(), silent = TRUE); .GV$rs <- NULL; return(c(p = NA_real_, status = 2)) }
  r <- tryCatch(rs$read(), error = function(err) NULL)
  if (is.null(r) || !identical(rs$get_state(), "idle")) {
    try(rs$kill(), silent = TRUE); .GV$rs <- NULL; return(c(p = NA_real_, status = 3))
  }
  if (!is.null(r$error)) return(c(p = NA_real_, status = 1))
  pv <- r$result
  if (length(pv) != 1 || !is.finite(pv)) return(c(p = NA_real_, status = 1))
  c(p = as.numeric(pv), status = 0)
}

## HL (check), mHL and, when asked, GiViTI for one fitted glm
arms_for_fit <- function(fit, do_giviti = TRUE) {
  ctx <- ebrahim.gof:::.gof_context(fit, G = 10)
  pv  <- function(z) { v <- z$p_value; if (length(v) == 1 && is.finite(v)) as.numeric(v) else NA_real_ }
  hl  <- tryCatch(pv(ebrahim.gof:::gof_hl(ctx)),    error = function(e) NA_real_)
  mhl <- tryCatch(pv(ebrahim.gof:::gof_ftest(ctx)), error = function(e) NA_real_)
  ## both Stukel forms: milliseconds, and they settle which form any stored column holds
  stj <- tryCatch(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "joint"))),    error = function(e) NA_real_)
  stm <- tryCatch(pv(ebrahim.gof:::gof_stukel(ctx, list(form = "marginal"))), error = function(e) NA_real_)
  gv  <- c(p = NA_real_, status = NA_real_); tg <- NA_real_
  if (do_giviti) {
    t0 <- proc.time()[["elapsed"]]
    gv <- giviti_p(as.numeric(ctx$y), as.numeric(ctx$ph))
    tg <- proc.time()[["elapsed"]] - t0
  }
  c(HL_check = hl, Stukel_joint = stj, Stukel_marginal = stm, mHL = mhl, GiViTI = unname(gv["p"]), giviti_status = unname(gv["status"]),
    giviti_secs = tg)
}

## write a data.frame atomically: .part, then rename
write_atomic <- function(df, f) {
  tmp <- paste0(f, ".part"); utils::write.csv(df, tmp, row.names = FALSE)
  if (file.exists(f)) file.remove(f)
  file.rename(tmp, f)
}
