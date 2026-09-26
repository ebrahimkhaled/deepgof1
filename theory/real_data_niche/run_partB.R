# =============================================================================
# run_partB.R -- Part B of PREDECLARATION.md (the confirmatory real-design simulation).
# Reads scenarios.csv (written by calibrate.R). Replicate r of scenario id: set.seed(31000000 + 100000 id + r);
# y ~ Bernoulli(plogis(eta0 + C g)) on the fixed real design X; fitted model = main effects on all covariates
# (fit_pol). Arms, in this RNG order: HL, Stukel marginal/joint, mHL, Osius-Rojek, GiViTI (no draws),
# Stute-Zhu (B = 200), DeepGOF-1 (dg_fast, B = 199), projection (Liu weight, raw X, B = 199).
# One CSV per (scenario, block of 20 reps), written .part then renamed: resumable.
# Env: PB_W workers (default 3); PB_SMOKE=1 -> rows_smoke/, seeds 999000000+, 2 reps, B = 19.
# =============================================================================
RD <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/real_data_niche"
SMOKE <- Sys.getenv("PB_SMOKE", "0") == "1"
NW <- as.integer(Sys.getenv("PB_W", "3"))
ROWD <- file.path(RD, if (SMOKE) "rows_smoke" else "rows"); dir.create(ROWD, showWarnings = FALSE)
SEEDB <- if (SMOKE) 999000000L else 31000000L
BB <- if (SMOKE) 19L else 199L
CH <- if (SMOKE) 2L else 20L
SC <- utils::read.csv(file.path(RD, "scenarios.csv"), stringsAsFactors = FALSE)
if (SMOKE) { SC$reps <- 2L; SC <- SC[SC$fam == "null" | SC$level == 3, ] }

jobs <- list()
for (i in seq_len(nrow(SC))) { sc <- SC[i, ]
  for (s in seq(1L, sc$reps, by = CH)) {
    f <- file.path(ROWD, sprintf("%s_r%04d.csv", sc$scen, s))
    if (!file.exists(f)) jobs[[length(jobs) + 1L]] <- list(sc = sc, reps = s:min(sc$reps, s + CH - 1L), f = f)
  }
}
## rep-block-major order, so a partial run is balanced across scenarios
jobs <- jobs[order(vapply(jobs, function(j) j$reps[1], 1), vapply(jobs, function(j) j$sc$id, 1))]
cat(sprintf("partB%s | %d scenarios | %d jobs | B = %d | %d workers | %s\n", if (SMOKE) " SMOKE" else "", nrow(SC), length(jobs), BB, NW, format(Sys.time())))
flush.console()
if (!length(jobs)) quit(save = "no")

cl <- parallel::makeCluster(min(NW, length(jobs)))
parallel::clusterExport(cl, c("SEEDB", "BB", "RD"))
invisible(parallel::clusterEvalQ(cl, {
  source(file.path(RD, "common.R"))
  .CACHE <- new.env(); .SB <- load_score_batch()
  get_design <- function(k) {
    key <- paste0("d", k)
    if (is.null(.CACHE[[key]])) { D <- design_setup(k); D$P <- dg_prep(D$X); D$AL <- list(raw = proj_A_liu_fast(as.matrix(D$X), std = FALSE))
      D$frm <- stats::as.formula(paste("y ~", paste(D$vars, collapse = " + "))); .CACHE[[key]] <- D }
    .CACHE[[key]]
  }
  one_rep <- function(sc, r) {
    D <- get_design(sc$k)
    set.seed(SEEDB + 100000L * sc$id + r)
    dat <- data.frame(y = gen_y(D, sc$design, sc$fam, sc$C), D$X)
    t0 <- proc.time()[[3]]
    fit <- fit_pol(dat)
    out <- c(HL = NA, Stukel_marginal = NA, Stukel_joint = NA, mHL = NA, OsiusRojek = NA, GiViTI = NA, StuteZhu = NA,
             DG = NA, proj = NA, sel_ab = NA, firth = NA, secs = NA)
    if (!is.null(fit)) {
      out[["firth"]] <- as.numeric(inherits(fit, "brglmFit") || identical(fit$method, "brglmFit"))
      a <- arms_cheap(fit); out[names(a)] <- a
      dg <- tryCatch(dg_fast(fit, dat, .SB, BB, D$P), error = function(e) NULL)
      if (!is.null(dg)) { out[["DG"]] <- dg$p; out[["sel_ab"]] <- as.numeric(setequal(dg$ax, c(D$a, D$b))) }
      fp <- proj_fit_pol(D$frm, dat)
      if (!is.null(fp)) out[["proj"]] <- tryCatch(unname(proj_p_multi(dat, D$frm, fp, B = BB, AL = D$AL)$p["raw"]), error = function(e) NA)
    }
    out[["secs"]] <- proc.time()[[3]] - t0
    data.frame(scen = sc$scen, id = sc$id, data = sc$data, design = sc$design, fam = sc$fam, level = sc$level, C = sc$C,
               rep = r, events = sum(dat$y), t(out))
  }
  NULL
}))
res <- parallel::parLapplyLB(cl, jobs, function(jb) {
  R <- do.call(rbind, lapply(jb$reps, function(r) one_rep(jb$sc, r)))
  tmp <- paste0(jb$f, ".part"); utils::write.csv(R, tmp, row.names = FALSE); file.rename(tmp, jb$f); jb$f
})
parallel::stopCluster(cl)
cat("DONE", format(Sys.time()), "\n")
