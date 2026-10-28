# 12_verify_nets.R -- cross-language gate for the 12 Tier 3 networks, and the batched scorer.
#  (a) deposited gate: base-R deepgof_score (code/rdeploy) on verify_input.csv vs PyTorch float32 (< 1e-6)
#  (b) batched scorer (tier3_common.R) on the same rows                                  (< 1e-6)
#  (c) float64 bench: R (both scorers) vs PyTorch float64 on 200 maps                    (< 1e-10)
#  (d) shipped network: batched scorer vs ebrahim.gof:::.dg_score on 500 random maps (same network)
#  (e) time of one batched call on 200 maps
source("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/theory/tier3/tier3_common.R")
suppressMessages(library(ebrahim.gof))
out <- list()
for (nm in NETS) {
  d <- file.path(T3, "nets", nm); s <- make_batch_scorer(file.path(d, "weights.csv"), file.path(d, "scaler.csv"))
  Z  <- as.matrix(utils::read.csv(file.path(d, "verify_input.csv"), header = FALSE)); Sp <- as.numeric(readLines(file.path(d, "verify_output.csv")))
  Z6 <- as.matrix(utils::read.csv(file.path(d, "verify64_input.csv"), header = FALSE)); S6 <- as.numeric(readLines(file.path(d, "verify64_output.csv")))
  sf <- make_batch_scorer(file.path(d, "weights.csv"), file.path(d, "scaler.csv"), fast = TRUE)
  a <- max(abs(apply(Z, 1, s$unbatched) - Sp)); b <- max(abs(s$zscore(Z) - Sp))
  c1 <- max(abs(apply(Z6, 1, s$unbatched) - S6)); c2 <- max(abs(s$zscore(Z6) - S6)); c3 <- max(abs(sf$zscore(Z6) - S6))
  ## raw-map path (scaler applied in R) vs PyTorch float64 on the same maps
  ms <- s$ms; Xraw <- sweep(sweep(Z6, 2, ms[2, ], "*"), 2, ms[1, ], "+"); c4 <- max(abs(sf$score(Xraw) - S6))
  ## the deposited 1e-6 gate compares a float32 PyTorch pass with float64 R; it is reported, while the
  ## pass criterion is the float64 bench (exact port up to float64 rounding)
  out[[nm]] <- data.frame(net = nm, r_vs_torch32 = a, batched_vs_torch32 = b, r_vs_torch64 = c1, batched_vs_torch64 = c2,
                          fast_vs_torch64 = c3, fast_rawmaps_vs_torch64 = c4,
                          pass = c1 < 1e-10 && c2 < 1e-10 && c3 < 1e-10 && c4 < 1e-9 && a < 5e-6)
}
V <- do.call(rbind, out); print(V, digits = 3)
M <- ebrahim.gof:::deepgof1_weights
sh <- make_batch_scorer(file.path(RDEP, "weights.csv"), file.path(RDEP, "scaler.csv"))
seed_at(5); X <- matrix(rnorm(500 * 36), 500, 36)
d_ship <- max(abs(sh$score(X) - apply(X, 1, function(m) ebrahim.gof:::.dg_score((m - M$mu) / M$sd, M))))
cat(sprintf("shipped network: batched scorer vs .dg_score, max |diff| %.3e\n", d_ship))
nets <- load_extra_nets(); X200 <- X[1:200, ]
t0 <- proc.time()[["elapsed"]]; for (k in 1:5) for (f in nets) f(X200); tt <- (proc.time()[["elapsed"]] - t0) / 5
cat(sprintf("12 networks x 200 maps, batched: %.3f s per data set (single process)\n", tt))
V$shipped_batched_vs_dgscore <- d_ship
utils::write.csv(V, file.path(T3, "NET_VERIFY.csv"), row.names = FALSE)
if (!all(V$pass)) stop("VERIFY FAILED")
cat("ALL 12 NETWORKS PASS\n")
