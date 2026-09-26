## Run with the working directory theory/_support/rdeploy of this archive (the script reads
## deepgof_forward.R, weights.csv and scaler.csv from there).
## theory_verify_numbers.R -- verify the two corrected numbers before they enter the supplement:
## (1) the recession closed form vs slopes measured FAR past the last kink (t = 1e6 .. 2e6);
## (2) the 50k probe draw's actual count of nonpositive directions (the script had a literal 0).
source("deepgof_forward.R")
W  <- load_weights("weights.csv"); ms <- as.matrix(utils::read.csv("scaler.csv", header=FALSE))
mu <- ms[1, ]; sdv <- ms[2, ]; K <- 6L
W0 <- W; for (nm in names(W0)) if (grepl("bias", nm)) W0[[nm]] <- W0[[nm]] * 0
Sinf  <- function(u) deepgof_score(u / sdv, W0, K)
Sfull <- function(m) deepgof_score((m - mu) / sdv, W, K)
vec_rows <- function(val) as.numeric(sapply(1:K, function(i) rep(val[i], K)))
stripe <- c(1.2,-.4,-1.1,-1.1,-.4,1.2); step <- c(-1,-1,-1,1,1,1)
dirs <- list(stripe_x1=vec_rows(stripe), step_x1=vec_rows(step),
             saddle=as.numeric(outer(step,step)), diffuse=rep(1,K*K),
             random={set.seed(11); rnorm(K*K)})
cat("(1) closed form vs slope on t = [1e6, 2e6] (past every kink):\n")
mx <- 0
for (nm in names(dirs)) {
  v <- dirs[[nm]]/sqrt(sum(dirs[[nm]]^2))
  sl <- (Sfull(2e6*v) - Sfull(1e6*v)) / 1e6
  d <- abs(Sinf(v) - sl); mx <- max(mx, d)
  cat(sprintf("   %-10s |Sinf - slope| = %.3e\n", nm, d))
}
cat(sprintf("   MAX over five directions: %.3e\n", mx))
## and the honest number on the ORIGINAL t=2000..4000 interval
mx2 <- max(sapply(dirs, function(d) { v <- d/sqrt(sum(d^2))
  abs(Sinf(v) - (Sfull(4000*v)-Sfull(2000*v))/2000) }))
cat(sprintf("   on the shipped t=[2000,4000] interval: max = %.3e\n\n", mx2))

cat("(2) the 50k probe draw (seed 7), counted for real:\n")
set.seed(7); N <- 50000L
w <- matrix(rnorm(N*K*K), N); w <- w / sqrt(rowSums(w^2))
vals <- apply(w, 1, Sinf)
cat(sprintf("   min = %.4f | count(Sinf <= 0) = %d of %d\n", min(vals), sum(vals <= 0), N))
cat(sprintf("   pooled with the 200k run (also 0 hits): one-sided 95%% CP bound on the\n"))
cat(sprintf("   spherical measure of the blind cone = 1 - 0.05^(1/250000) = %.3e\n",
            1 - 0.05^(1/250000)))
