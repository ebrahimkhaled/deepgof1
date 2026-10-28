# 00_niche_repro.R -- can the stored niche DeepGOF-1 p-values (theory/niche/rows, column DG) be reproduced?
# The niche run (theory/niche/run_niche.R) did NOT call deepgof1(): it used its own loop (fit_pol with a
# Firth fallback, batched forward pass of code/rdeploy, the bootstrap continuing the data-generation
# stream after run.all.gof + lc_fast). This script evaluates run_niche.R's `defs` verbatim and checks
# (1) whether run.all.gof / lc_fast touch the RNG, (2) whether dg_both() reproduces the stored DG and
# DG_all on a few data sets per scenario, (3) the time of one dg_both call.
PAP  <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof"
NSRC <- file.path(PAP, "theory/niche/run_niche.R")
ex <- parse(NSRC)
for (e in ex) { s <- paste(deparse(e), collapse = " ")
  if (grepl("^(HERE|defs|SC|SEEDB|B|SMOKE|NNULL|NALT) <-", s) || grepl("^SC\\$", s)) eval(e, globalenv()) }
eval(defs)
sb <- load_score_batch()
st <- do.call(rbind, lapply(list.files(file.path(PAP, "theory/niche/rows"), full.names = TRUE), read.csv))
cat("stored rows:", nrow(st), "\n")
res <- list()
for (i in seq_len(nrow(SC))) {
  sc <- SC[i, ]
  for (r in c(1L, 2L)) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(SEEDB + sc$id * 100000L + r)
    dat <- gen(sc); fit <- fit_pol(dat)
    s0 <- .Random.seed
    rr <- suppressWarnings(suppressMessages(run.all.gof(fit, tests = c("HL", "Stukel"), include_slow = FALSE, install = "no")))
    s1 <- .Random.seed
    ph <- as.numeric(stats::fitted(fit)); X <- stats::model.matrix(fit)
    l1 <- lc_fast(dat$y, ph, X, lapply(seq_len(ncol(X))[-1], function(j) X[, j]))
    s2 <- .Random.seed
    t0 <- proc.time()[["elapsed"]]
    dg <- dg_both(fit, dat, sb, 199L)
    tt <- proc.time()[["elapsed"]] - t0
    z <- st[st$scen == sc$scen & st$rep == r, ]
    res[[length(res) + 1]] <- data.frame(scen = sc$scen, rep = r, ev = sum(dat$y), ev_st = z$events,
      rng_touched_by_gof = !identical(s0, s1), rng_touched_by_lc = !identical(s1, s2),
      DG = dg[["DG"]], DG_st = z$DG, DGall = dg[["DG_all"]], DGall_st = z$DG_all, secs = tt)
  }
}
R <- do.call(rbind, res); print(R)
cat("events match:", all(R$ev == R$ev_st), "| DG match:", sum(R$DG == R$DG_st), "of", nrow(R),
    "| DG_all match:", sum(R$DGall == R$DGall_st), "\n")
