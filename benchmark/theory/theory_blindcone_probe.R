## Run with the working directory theory/_support/rdeploy of this archive (the script reads
## deepgof_forward.R, weights.csv and scaler.csv from there).
## theory_blindcone_probe.R -- characterize the exceptional cone {v : Sinf(v) <= 0} of the shipped
## weights: multi-start descent for a better global minimum, and the worst direction's 6x6 pattern.
source("deepgof_forward.R")
W  <- load_weights("weights.csv"); ms <- as.matrix(utils::read.csv("scaler.csv", header=FALSE))
sdv <- ms[2, ]; K <- 6L
W0 <- W; for (nm in names(W0)) if (grepl("bias", nm)) W0[[nm]] <- W0[[nm]] * 0
Sinf <- function(u) deepgof_score(u / sdv, W0, K)

set.seed(7)
N <- 50000L
w <- matrix(rnorm(N*K*K), N); w <- w / sqrt(rowSums(w^2))
vals <- apply(w, 1, Sinf)
starts <- w[order(vals)[1:20], ]
best <- Inf; bv <- NULL
for (s in 1:20) {
  v <- starts[s, ]; cur <- Sinf(v)
  for (r in c(.5,.2,.08,.03,.01,.004)) for (it in 1:1500) {
    cand <- v + r*rnorm(K*K); cand <- cand/sqrt(sum(cand^2))
    sc <- Sinf(cand); if (sc < cur) { cur <- sc; v <- cand }
  }
  if (cur < best) { best <- cur; bv <- v }
}
cat(sprintf("multi-start global-ish minimum of Sinf on the sphere: %.4f\n", best))
cat(sprintf("sampled measure of {Sinf<=0}: %d of %d draws in THIS run
", sum(vals <= 0), N))
cat("pooled with theory_sphere_check.R 200k run (0 hits): one-sided 95% CP bound
")
cat(sprintf("on the spherical measure = 1 - 0.05^(1/250000) = %.3e
", 1 - 0.05^(1/250000)))
cat("\nworst direction, as the 6x6 map pattern (rows = x1 bins, cols = x2 bins):\n")
M <- matrix(bv, K, K, byrow = TRUE)
for (i in 1:K) cat(sprintf("% .2f ", M[i, ]), "\n")
cat(sprintf("\npattern diagnostics: mean %.3f | sd of row means %.3f | sd of col means %.3f\n",
            mean(M), sd(rowMeans(M)), sd(colMeans(M))))
cat(sprintf("lag-1 autocorr along rows: %.3f  (natural departures are SMOOTH: near +1)\n",
    cor(as.numeric(M[, 1:(K-1)]), as.numeric(M[, 2:K]))))
write.csv(data.frame(v = bv, Sinf = best), "../../results/theory_blind_direction.csv", row.names = FALSE)
## margins of the natural departure directions, for the theorem's applicability table
vec_rows <- function(val) as.numeric(sapply(1:K, function(i) rep(val[i], K)))
nat <- list("curvature stripe" = vec_rows(c(1.2,-.4,-1.1,-1.1,-.4,1.2)),
            "step/threshold"   = vec_rows(c(-1,-1,-1,1,1,1)),
            "interaction saddle"= as.numeric(outer(c(-1,-1,-1,1,1,1), c(-1,-1,-1,1,1,1))),
            "diffuse drift"    = rep(1, K*K),
            "local bump"       = { z <- rep(0,K*K); z[c(15,16,21,22)] <- 1; z })
cat("\nSinf at the natural departure directions (the theorem's applicability margins):\n")
for (nm in names(nat)) { v <- nat[[nm]]/sqrt(sum(nat[[nm]]^2))
  cat(sprintf("  %-19s %6.3f   (mirror %6.3f)\n", nm, Sinf(v), Sinf(-v))) }
